#!/usr/bin/env python
from InitialStateAnalysis.Limits.Limits import Limits
import os
import logging
import sys
import argparse
import datetime
import subprocess
import numpy as np
from InitialStateAnalysis.Plotters.plotUtils import _3L_MASSES, _4L_MASSES, ZMASS
from InitialStateAnalysis.Plotters.Plotter import Plotter
from multiprocessing import Pool
from mklimits import *
import ROOT
import pickle
from array import array

def wrapper(args):
    analysis = args[0]
    region = args[1]
    period = args[2]
    optVar = args[3]
    binRange = args[4]
    mass = args[5]
    bp = args[6]
    bgMode = args[7]
    scaleFactor = args[8]
    unblind = args[9]
    # iterative optimization
    numBins = 10
    prevLimit = -1
    currLimit = 0
    limits = [-1000,100,1000]
    optLevel = 0
    while abs(limits[0]-limits[1])/limits[1]>0.01 or abs(limits[2]-limits[1])/limits[1]>0.01:
        optLevel += 1
        binning = [numBins,binRange[0],binRange[1]]
        # process datacards and get limitdata
        #for mass in masses:
        #    optimizer(analysis,region,period,optVar,mass,bp,bgMode,scaleFactor,unblind,binning)
        optimizerWrapper([analysis,region,period,optVar,mass,bp,bgMode,scaleFactor,unblind,binning])
        # maximize expected limits
        limits, binRange = minimizer(analysis,region,period,optVar,bp,bgMode,binning,mass)
        if optLevel > 3: break # to prevent endless recursion
        if limits[1]<0.00001: limits[1]=0.00001
    with open('optimization.log','a') as f:
        f.write('%s %s %i %s - Optimal bin: %f, Expected Limit: %f\n' %(str(datetime.datetime.now()),optVar,mass,bp,sum(binRange)/2,limits[1]))
    return (optVar,mass,bp,sum(binRange)/2,limits[1])


def optimizerWrapper(args):
    analysis = args[0]
    region = args[1]
    period = args[2]
    optVar = args[3]
    mass = args[4]
    bp = args[5]
    bgMode = args[6]
    scaleFactor = args[7]
    unblind = args[8]
    binning = args[9]
    optimizer(analysis,region,period,optVar,mass,bp,bgMode,scaleFactor,unblind,binning)

def minimizer(analysis,region,period,optVar,bp,bgMode,binning,mass):
    # get limitdata
    bins = [(binning[1]+x*(binning[2]-binning[1])/binning[0]) for x in range(binning[0])]
    expectedLimits = [0] * len(bins)
    for c in range(len(bins)):
        quartiles = np.empty((6, 1), dtype=float)
        destDir = 'datacards_optimization/%s/%s_%itev_%s/%s/%i' % (optVar,analysis,period,region,bp,mass)
        rootName = '%s_%s%i.root' %(bp,optVar,c)
        fname = os.path.join(destDir, rootName)
        file = ROOT.TFile(fname,"READ")
        tree = file.Get("limit")
        if not tree: continue
        for i, row in enumerate(tree):
            quartiles[i,0] = row.limit
        n = 1
        expectedLimits[c] = quartiles[2][0]
    for b,l in zip(bins,expectedLimits):
        with open('optimization.log','a') as f:
            f.write('%s %s %i %s - %f: %f\n' %(str(datetime.datetime.now()),optVar,mass,bp,b,l))
    minExpectedLimit = min(expectedLimits)
    minBin = expectedLimits.index(minExpectedLimit)
    if minBin>0 and minBin<(len(expectedLimits)-1): limits = [expectedLimits[minBin-1], expectedLimits[minBin], expectedLimits[minBin+1]]
    elif minBin>0: limits = [expectedLimits[minBin-1], expectedLimits[minBin], expectedLimits[minBin]]
    else: limits = [expectedLimits[minBin], expectedLimits[minBin], expectedLimits[minBin+1]]
    lowRangeBin = bins[minBin-1] if minBin else bins[0]
    highRangeBin = bins[minBin+1] if minBin<(len(bins)-1) else bins[-1]
    return (limits, [lowRangeBin,highRangeBin])

def optimizer(analysis,region,period,optVar,mass,bp,bgMode,scaleFactor,unblind,binning):
    # first, generate datacards
    bins = [binning[1]+x*(binning[2]-binning[1])/binning[0] for x in range(binning[0])]
    cutMap = {
        'pre'  : '1',
        'st'   : 'finalstate.sT>1.1*%f+60.' %mass,
        'zmass': 'fabs(z1.mass-%f)>80.' %ZMASS,
        'dphi' : 'h1.dPhi<%f/600.+1.95' %mass,
        'hmass': 'h1.mass>0.9*%f & h1.mass<1.1*%f' %(mass,mass),
    }
    cut = 'select.PassTight'
    for c in range(len(bins)):
        denomCut = cut
        if True: # add cut on top of rest or do it alone
            if 'h1' in optVar or 'h2' in optVar:
                denomCut = cut + ' & '.join([cutMap[x] for x in cutMap if x!='hmass'])
            elif 'z1' in optVar or 'z2' in optVar:
                denomCut = cut + ' & '.join([cutMap[x] for x in cutMap if x!='zmass'])
            elif 'dP' in optVar or 'dR' in optVar:
                denomCut = cut + ' & '.join([cutMap[x] for x in cutMap if x!='dphi'])
            elif 'st' == optVar:
                denomCut = cut + ' & '.join([cutMap[x] for x in cutMap if x!='st'])
            else:
                denomCut = cut + ' & '.join([cutMap[x] for x in cutMap])
        if optVar=='h1mass'     : passCut = 'fabs(h1.mass-%f)<%f&&%s' % (mass, bins[c], denomCut)
        if optVar=='h1massUnder': passCut = '%f-h1.mass<%f&&h1.mass<%f&&%s' % (mass, bins[c], mass, denomCut)
        if optVar=='h1massOver' : passCut = 'h1.mass-%f<%f&&h1.mass>%f&&%s' % (mass, bins[c], mass, denomCut)
        if optVar=='h2mass'     : passCut = 'fabs(h2.mass-%f)<%f&&%s' % (mass, bins[c], denomCut)
        if optVar=='h2massUnder': passCut = '%f-h2.mass<%f&&h2.mass<%f&&%s' % (mass, bins[c], mass, denomCut)
        if optVar=='h2massOver' : passCut = 'h2.mass-%f<%f&&h2.mass>%f&&%s' % (mass, bins[c], mass, denomCut)
        if optVar=='z1mass'     : passCut = 'fabs(z1.mass-%f)>%f&&%s' % (ZMASS, bins[c], denomCut)
        if optVar=='z2mass'     : passCut = 'fabs(z2.mass-%f)>%f&&%s' % (ZMASS, bins[c], denomCut)
        if optVar=='dPhi1'      : passCut = 'fabs(h1.dPhi)<%f&&%s' % (bins[c], denomCut)
        if optVar=='dPhi2'      : passCut = 'fabs(h2.dPhi)<%f&&%s' % (bins[c], denomCut)
        if optVar=='dR1'        : passCut = 'h1.dR<%f&&%s' % (bins[c], denomCut)
        if optVar=='dR2'        : passCut = 'h2.dR<%f&&%s' % (bins[c], denomCut)
        if optVar=='met'        : passCut = 'finalstate.met>%f&&%s' % (bins[c], denomCut)
        if optVar=='st'         : passCut = 'finalstate.sT>%f&&%s' % (bins[c], denomCut)
        datacardDir = './datacards_optimization/%s' % optVar
        # make datacards
        BP(analysis,region,period,mass,bp,mode=bgMode,scalefactor=scaleFactor,unblind=unblind,datacardDir=datacardDir,fullCut=passCut,optVar=optVar,optCut=c)
        # merge the datacards
        pipe = subprocess.PIPE
        dn = '> /dev/null 2>&1' # send to devnull
        command = ''
        combineDir = '/cms/dntaylor/HIGGSCOMBINE/CMSSW_6_1_1/src/'
        command += 'pushd %s; eval `scramv1 runtime -sh`; popd;' % (combineDir)
        cardDir = '%s/%s_%itev_%s/%s/%s' % (datacardDir, analysis, period, region, bp, mass)
        command += 'pushd %s;' % (cardDir)
        cardname = '%s_%s%i' %(bp,optVar,c)
        command += 'combineCards.py %s_[em][em][em].txt > %s_comb.txt;'  %(cardname,cardname)
        #os.system(command)
        out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]
        # run the limit
        destDir = 'datacards_optimization/%s/%i' % (bp,mass)
        command = ''
        command += 'pushd %s; mkdir -p %s/%s;' % (cardDir, combineDir, destDir)
        command += 'cp %s_comb.txt %s/%s/%s.txt;' %(cardname,combineDir,destDir,cardname)
        command += 'pushd %s/%s; eval `scramv1 runtime -sh`;' % (combineDir,destDir)
        command += 'combine -m %i -M Asymptotic %s.txt;' %(mass,cardname)
        command += 'popd; cp %s/%s/higgsCombineTest.Asymptotic.mH%i.root %s.root;' % (combineDir,destDir,mass,cardname)
        #os.system(command)
        out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]

class MultiOptimizer( ROOT.TPyMultiGenFunction ):
    '''Class to run optimization on multiple variables using datacards'''
    def __init__(self,args):
        self.analysis = args[0]
        self.region = args[1]
        self.period = args[2]
        self.mass = args[3]
        self.bp = args[4]
        self.bgMode = args[5]
        self.scaleFactor = args[6]
        self.unblind = args[7]
        ROOT.TPyMultiGenFunction.__init__(self,self)

    def NDim(self): return 3

    def DoEval(self,args):
        '''The optimizer to create datacard, run the limit, return the expected limit'''
        # get variables
        st = args[0]
        z1mass = args[1]
        #dPhi1 = args[2]
        dR1 = args[2]
        # set selection
        cutMap = {
            'st'   : 'finalstate.sT>%f' %st,
            'zmass': 'fabs(z1.mass-%f)>%f' %(ZMASS,z1mass),
            #'dphi' : 'h1.dPhi<%f' %dPhi1,
            'dR' : 'h1.dR<%f' %dR1,
            'hmass': 'h1.mass>0.9*%f & h1.mass<1.1*%f' %(self.mass,self.mass),
        }
        passCut = 'select.PassTight & ' + ' & '.join([cutMap[x] for x in cutMap])
        datacardDir = './datacards_optimization_multi'
        # make datacards
        BP(self.analysis,self.region,self.period,self.mass,self.bp,mode=self.bgMode,scalefactor=self.scaleFactor,unblind=self.unblind,datacardDir=datacardDir,fullCut=passCut)
        # merge the datacards
        pipe = subprocess.PIPE
        dn = '> /dev/null 2>&1' # send to devnull
        command = ''
        combineDir = '/cms/dntaylor/HIGGSCOMBINE/CMSSW_6_1_1/src/'
        command += 'pushd %s; eval `scramv1 runtime -sh`; popd;' % (combineDir)
        cardDir = '%s/%s_%itev_%s/%s/%s' % (datacardDir, self.analysis, self.period, self.region, self.bp, self.mass)
        command += 'pushd %s;' % (cardDir)
        cardname = '%s_%s%i' %(self.bp,'',0)
        command += 'combineCards.py %s_[em][em][em].txt > %s_comb.txt;'  %(cardname,cardname)
        #os.system(command)
        out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]
        # run the limit
        destDir = 'datacards_optimization_multi/%s/%i' % (self.bp,self.mass)
        command = ''
        command += 'pushd %s; mkdir -p %s/%s;' % (cardDir, combineDir, destDir)
        command += 'cp %s_comb.txt %s/%s/%s.txt;' %(cardname,combineDir,destDir,cardname)
        command += 'pushd %s/%s; eval `scramv1 runtime -sh`;' % (combineDir,destDir)
        command += 'combine -m %i -M Asymptotic %s.txt;' %(self.mass,cardname)
        command += 'popd; cp %s/%s/higgsCombineTest.Asymptotic.mH%i.root %s.root;' % (combineDir,destDir,self.mass,cardname)
        #os.system(command)
        out = subprocess.Popen(command, shell=True,stdout=pipe,stderr=subprocess.STDOUT).communicate()[0]
        # get limitdata
        quartiles = np.empty((6, 1), dtype=float)
        rootName = '%s.root' % cardname
        fname = os.path.join(cardDir, rootName)
        file = ROOT.TFile(fname,"READ")
        tree = file.Get("limit")
        for i, row in enumerate(tree):
            quartiles[i,0] = row.limit
        expectedLimit = quartiles[2][0]
        #print 'Optimization values: st: %f, z1mass: %f, dPhi1: %f' % (st, z1mass, dPhi1)
        print 'Optimization values: st: %f, z1mass: %f, dR1: %f' % (st, z1mass, dR1)
        print 'Expected limit: %f' % expectedLimit
        with open('optimization_multi.log','a') as f:
            #f.write('%s %i %s - st: %f, z1mass: %f, dPhi1: %f, Limit: %f\n' %(str(datetime.datetime.now()),self.mass,self.bp,st,z1mass,dPhi1,expectedLimit))
            f.write('%s %i %s - st: %f, z1mass: %f, dR1: %f, Limit: %f\n' %(str(datetime.datetime.now()),self.mass,self.bp,st,z1mass,dR1,expectedLimit))
        return expectedLimit

def multiWrapper(args):
    print 'Setup minimizer'
    mass = args[3]
    bp = args[4]
    # setup minimizer
    #mini = ROOT.Math.Factory.CreateMinimizer("GSLMultiMin", "BFGS") # wont go past first iteration, via derivatives, problematic
    mini = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Simplex")
    #mini = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Migrad")
    mini.SetMaxFunctionCalls(1000000)
    mini.SetMaxIterations(100000)
    mini.SetTolerance(0.01)
    mini.SetPrintLevel(1)

    miniFunct = MultiOptimizer(args)
    mini.SetFunction(miniFunct)

    print 'Define variables'
    # setup variables
    currCuts = [
        ['st', 1.1*mass+60., 50., 100., 1500.],
        ['z1mass', 80., 5., 1., 120.],
        #['dPhi1', mass/600.+1.95, .1, 1, 3.14],
        ['dR1', 3., .1, 1., 5.],
    ]

    for v,params in enumerate(currCuts):
        mini.SetLimitedVariable(v,params[0],params[1],params[2],params[3],params[4])

    print 'Minimize'
    # minimize
    mini.Minimize()

    # return values
    result = {}
    errors = mini.Errors()
    minimum = mini.X()
    for params,minVals,minErr in zip(currCuts,minimum,errors):
        result[params[0]] = [minVals,minErr]
    with open('optimization_multi.log','a') as f:
        #f.write('%s %i %s - Optimal: st: %f, z1mass: %f, dPhi1: %f\n' %(str(datetime.datetime.now()),mass,bp,result['st'][0],result['z1mass'][0],result['dPhi1'][0]))
        f.write('%s %i %s - Optimal: st: %f, z1mass: %f, dR1: %f\n' %(str(datetime.datetime.now()),mass,bp,result['st'][0],result['z1mass'][0],result['dR1'][0]))
    return [mass, bp, result]


def BP(analysis,region,period,mass,bp,**kwargs):
    optVar = kwargs.pop('optVar','')
    optCut = kwargs.pop('optCut',0)
    if bp == 'ee100':
        s = Scales(1., 0., 0., 0., 0., 0.)
    elif bp == 'em100':
        s = Scales(0., 1., 0., 0., 0., 0.)
    elif bp == 'mm100':
        s = Scales(0., 0., 0., 1., 0., 0.)
    else:
        print 'Unknown branching point: %s' %bp
    sf = getattr(s,'scale_%s'%analysis)
    print 'Processing %s mass %i optimization %s cut %i' %(bp, mass, optVar, optCut)
    chanMap = {
        'Hpp3l': {
             'names': ['ee','em','mm'],
             'ee'   : ['eee','eem'],
             'em'   : ['eme','emm','mee','mem'],
             'mm'   : ['mme','mmm'],
             'eegen': ['eee','eem','eet'],
             'emgen': ['eme','emm','emt'],
             'mmgen': ['mme','mmm','mmt'],
        },
        'Hpp4l': {
             'names': ['eeee','eeem','eemm','emem','emmm','mmmm'],
             'eeee' : ['eeee'],
             'eeem' : ['eeem','eeme','emee','meee'],
             'eemm' : ['eemm','mmee'],
             'emem' : ['emem','emme','meem','meme'],
             'emmm' : ['emmm','memm','mmem','mmme'],
             'mmmm' : ['mmmm'],
        },
    }

    chanCuts = []
    chanScales = []
    for c in chanMap[analysis]['names']:
        thisScale = sf(c[:2],c[2:])
        if thisScale==0: continue
        chanCut = '('+' | '.join(['channel=="%s"'%x for x in chanMap[analysis][c]])+')'
        chanCut += ' & (' + ' | '.join(['genChannel=="%s"'%x for x in chanMap[analysis][c+'gen']+['aaa']]) + ')'
        #limit(analysis,region,period,mass,bp=bp,name='%s_%s%i%s'%(bp,optVar,optCut,c),directory=bp,channels=[chanCut],scale=[thisScale],**kwargs)
        # now do individual channels
        for theChan in chanMap[analysis][c]:
            thisChanCut = 'channel=="%s"' % theChan
            thisChanCut += ' & (' + ' | '.join(['genChannel=="%s"'%x for x in chanMap[analysis][c+'gen']+['aaa']]) + ')'
            limit(analysis,region,period,mass,bp=bp,name='%s_%s%i_%s'%(bp,optVar,optCut,theChan),directory=bp,channels=[thisChanCut],scale=[thisScale],doIndividualChannel=True,**kwargs)
        chanCuts += [chanCut]
        chanScales += [thisScale]
    #limit(analysis,region,period,mass,name='%s_%s%i'%(bp,optVar,optCut),bp=bp,directory=bp,channels=chanCuts,scale=chanScales,**kwargs)

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Produce datacards")

    parser.add_argument('analysis', type=str, choices=['Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('region', type=str, choices=['Hpp3l','Hpp4l','WZ'], help='Analysis to run')
    parser.add_argument('period', type=int, choices=[8, 13], help='Energy (TeV)')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500,help='Mass for signal')
    parser.add_argument('-am','--allMasses',action='store_true',help='Run over all masses for signal')
    parser.add_argument('-ub','--unblind',action='store_true',help='unblind')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='ee100',default='ee100',choices=['ee100','em100','mm100'],help='Choose branching point for H++')
    parser.add_argument('-ov','--optimizationVariable',nargs='?',type=str,const='z1mass',default='z1mass',choices=['z1mass','dPhi1','dR1','st'],help='Choose variable to optimize')
    parser.add_argument('-ab','--allBranchingPoints',action='store_true',help='Run over all branching points for H++')
    parser.add_argument('-av','--allOptimizationVariables',action='store_true',help='Run over all optimization variables')
    parser.add_argument('-bg','--bgMode',nargs='?',type=str,const='sideband',default='sideband',choices=['mc','sideband'],help='Choose BG estimation')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for MC.')
    parser.add_argument('-dm','--doMultidimensional',action='store_true',help='Run the multidimensional fit')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    # plan
    # for each variable
    #   iterate over cuts
    #   for each mass
    #     for each 100% BP
    #       create datacard with this selection (the rest at current default)
    # process all datacards
    # run limits script
    # two optimization variables
    #   1) S/sqrt(B)
    #   2) expected limit

    branchingPoints = ['ee100','em100','mm100']
    masses = _3L_MASSES if args.analysis=='Hpp3l' else _4L_MASSES
    plotTypes = {
        #'h1mass'     : [0,500],
        #'h1massUnder': [0,500],
        #'h1massOver' : [0,500],
        #'h2mass'     : [0,500],
        #'h2massUnder': [0,500],
        #'h2massOver' : [0,500],
        'z1mass'     : [5,105],
        #'z2mass'     : [0,100],
        'dPhi1'      : [1.5,3.],
        #'dPhi2'      : [0,3.2],
        'dR1'        : [1.,4.],
        #'dR2'        : [0,6.],
        #'met'        : [0,100],
        'st'         : [100,1500],
    }
    optimizationVariables = plotTypes.keys()

    if not args.allMasses: masses = [args.mass]
    if not args.allBranchingPoints: branchingPoints = [args.branchingPoint]
    if not args.allOptimizationVariables: optimizationVariables = [args.optimizationVariable]

    if args.doMultidimensional: # run the multi fit with root minimizer
        with open('optimization_multi.log','w') as f:
            f.write('')

        poolArgs = [[m,b] for b in branchingPoints for m in masses]

        if len(poolArgs)==1:
            job = poolArgs[0]
            result = multiWrapper([args.analysis,args.region,args.period,job[0],job[1],args.bgMode,args.scaleFactor,args.unblind])
            pickle.dump(result,open('optimization_multi.pkl','wb'))
        else:
            p = Pool(8)
            try:
                result = p.map_async(multiWrapper, [(args.analysis,args.region,args.period,job[0],job[1],args.bgMode,args.scaleFactor,args.unblind) for job in poolArgs]).get(999999)
                pickle.dump(result,open('optimization_multi.pkl','wb'))
            except KeyboardInterrupt:
                p.terminate()
                print 'limits cancelled'
                sys.exit(1)

    else: # single var minimizer
        with open('optimization.log','w') as f:
            f.write('')

        poolArgs = [[v,b,m] for v in optimizationVariables for b in branchingPoints for m in masses]

        if len(poolArgs)==1:
            job = poolArgs[0]
            result = wrapper([args.analysis,args.region,args.period,job[0],plotTypes[job[0]],job[2],job[1],args.bgMode,args.scaleFactor,args.unblind])
            pickle.dump(result,open('optimization.pkl','wb'))
        else:
            p = Pool(8)
            try:
                result = p.map_async(wrapper, [(args.analysis,args.region,args.period,job[0],plotTypes[job[0]],job[2],job[1],args.bgMode,args.scaleFactor,args.unblind) for job in poolArgs]).get(999999)
                pickle.dump(result,open('optimization.pkl','wb'))
            except KeyboardInterrupt:
                p.terminate()
                print 'limits cancelled'
                sys.exit(1)
    
    return 0


if __name__ == "__main__":
    main()
