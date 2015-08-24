'''
Optimizer.py

Optimize selections.

Author: Devin N. Taylor, UW-Madison
'''
import math
import ROOT
import pickle
import logging
import sys
from multiprocessing import Pool

class FunctionOptimizer( ROOT.TPyMultiGenFunction ):
    '''Multi-variable optimizer'''
    def __init__(self,args):
        self.plotter = args[0]
        self.sigSelection = args[1]
        self.bgSelection = args[2]
        self.cuts = args[3]
        self.formula = args[4]
        self.maximize = args[5]
        ROOT.TPyMultiGenFunction.__init__(self,self)

    def NDim(self): return len(self.cuts)

    def DoEval(self,args):
        sigSel = self.sigSelection
        bgSel = self.bgSelection
        for cut,arg in zip(self.cuts,args):
            sigSel += ' & %s %f' % (cut,arg)
            bgSel += ' & %s %f' % (cut,arg)
        sig = self.plotter.getSignalEntries(sigSel)
        bg = self.plotter.getBackgroundEntries(bgSel)
        #result = self.formula.EvalPar(sig,bg) if bg else sig
        result = sig/math.sqrt(bg) if bg else 0.
        print 'Current args:', [args[x] for x in range(len(self.cuts))], 'result:', result
        return -result if self.maximize else result

class Optimizer(object):
    '''Class to optimize a series of cuts'''
    def __init__(self,plotter):
        self.plotter = plotter
        self.cuts = []
        self.formula = ROOT.TFormula('f','[0]/sqrt([1])')

    def setSelection(self,selection,**kwargs):
        '''Set preselection for both background and signal.'''
        signalSelection = kwargs.pop('signalSelection','1')
        backgroundSelection = kwargs.pop('backgroundSelection','1')
        self.selection = selection
        self.signalSelection = signalSelection
        self.backgroundSelection = backgroundSelection

    def addCut(self,name,func,min,max,step):
        '''Add a cut to optimize'''
        self.cuts += [{'name': name, 'func': func, 'min': min, 'max': max, 'step': step}]

    def setSignificance(self,formula):
        '''Set the significance formula, variables are SIG and BG.'''
        formula.replace('SIG','[0]')
        formula.replace('BG','[1]')
        self.formula = ROOT.TFormula('f',formula)
        if self.formula.Compile()!=0: print 'Formula failed to compile: %s' % formula

    def optimize(self,**kwargs):
        '''Optimize selections'''
        masses = kwargs.pop('masses',[])
        sigSel = self.selection + ' & ' + self.signalSelection
        bgSel = self.selection + ' & ' + self.backgroundSelection
        # do the simple N-D optimization with minuit simplex method
        if False: # broken
            mini = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Simplex")
            mini.SetMaxFunctionCalls(1000)
            mini.SetMaxIterations(1000)
            mini.SetTolerance(0.0001)
            mini.SetPrintLevel(1)
            args = [self.plotter,sigSel,bgSel,[cut['func'] for cut in self.cuts],self.formula,True]
            miniFunct = FunctionOptimizer(args)
            mini.SetFunction(miniFunct)
            for v,cut in enumerate(self.cuts):
                mini.SetLimitedVariable(v,cut['name'],float(cut['max']-cut['min'])/2.,10.*cut['step'],cut['min'],cut['max'])
            mini.Minimize()
        # get efficiencies
        if False:
            for cut in self.cuts:
                #print cut
                cutRange = [cut['min']+x*cut['step'] for x in range(int((cut['max']-cut['min'])/cut['step']))]
                effs = []
                denom = self.plotter.getSignalEntries(sigSel)
                for cutVal in cutRange:
                    num = self.plotter.getSignalEntries('%s & %s %f' %(sigSel, cut['func'], cutVal))
                    eff = num/denom
                    effs += [eff]
                    #print '  %f: %f' %(cutVal, eff)
                #effs = [self.plotter.getSignalEntries('%s & %s %f' %(sigSel, cut['func'], cutVal))/self.plotter.getSignalEntries(sigSel) for cutVal in cutRange]
                #print effs
                effCutoffs = [0.7,0.8,0.9,0.95,0.99,0.999]
                if effs[0]>effs[-1]: # reversed list
                    for c in effCutoffs:
                        highEffs = [i for i,e in enumerate(effs) if e<c]
                        if highEffs: 
                            val = cutRange[highEffs[0]]
                        else:
                            val = -1
                        print '%s: Eff=%f: %f' % (cut['name'],c,val)
                else:
                    for c in effCutoffs:
                        highEffs = [i for i,e in enumerate(effs) if e>c]
                        if highEffs: 
                            val = cutRange[highEffs[0]]
                        else:
                            val = -1
                        print '%s: Eff=%f: %f' % (cut['name'],c,val)

        optimizationVals = {}
        for cut in self.cuts:
            optimizationVals[cut['name']] = {}
            for mass in masses:
                optimizationVals[cut['name']][mass] = self.getIndividualMassCut(cut,mass,sigSel,bgSel)
        with open('optimize.pkl','wb') as file:
            pickle.dump(optimizationVals,file)


    def getIndividualMassCut(self,cut,mass,sigSel,bgSel):
        optimizationVals = {}
        cutname = cut['name']
        print '%s:%i' % (cutname,mass)
        cutRange = [cut['min']+x*cut['step'] for x in range(int((cut['max']-cut['min'])/cut['step']))]
        thisFunc = cut['func'].replace('MASS',str(mass))
        sigSample = 'HPlusPlusHMinusHTo3L_M-%i_8TeV-calchep-pythia6' % mass
        print '%s:%s Passing signal' % (cutname,mass)
        sigpass = [self.plotter.getSignalEntries('%s & %s %f' %(sigSel, thisFunc, cutVal),signal=sigSample,doError=True) for cutVal in cutRange]
        print '%s:%s All signal' % (cutname,mass)
        sigall = self.plotter.getSignalEntries(sigSel,signal=sigSample,doError=True)
        print '%s:%s Passing background' % (cutname,mass)
        bgpass = [self.plotter.getBackgroundEntries('%s & %s %f' %(bgSel, thisFunc, cutVal),doError=True) for cutVal in cutRange]
        print '%s:%s All background' % (cutname,mass)
        bgall = self.plotter.getBackgroundEntries(bgSel,doError=True)
        sigEff = []
        bgEff = []
        significance1 = []
        significance2 = []
        significance3 = []
        for s,b in zip(sigpass,bgpass):
            sigEffVal = s[0]/sigall[0] if sigall[0] else -1.
            sigEffErr = sigEffVal * (s[1]**2/s[0]**2 + sigall[1]**2/sigall[0]**2)**0.5 if sigall[0] and s[0] else -1.
            sigEff += [[sigEffVal, sigEffErr]]
            bgEffVal = b[0]/bgall[0] if bgall[0] else -1.
            bgEffErr = bgEffVal * (b[1]**2/b[0]**2 + bgall[1]**2/bgall[0]**2)**0.5 if b[0] and bgall[0] else -1.
            bgEff += [[bgEffVal, bgEffErr]]
            significance1Val = s[0]/b[0]**0.5 if b[0] else -1.
            significance1Err = significance1Val * (s[1]**2/s[0]**2 + 0.5**2 * b[1]**2/b[0]**2)**0.5 if s[0] and b[0] else -1.
            significance1 += [[significance1Val, significance1Err]]
            significance2Val = s[0]/(s[0] + b[0])**0.5 if (s[0]+b[0]) else -1.
            significance2Err = significance2Val * (s[1]**2/s[0]**2 + 0.5**2 * (s[1]**2 + b[1]**2)/(s[0] + b[0])**2)**0.5 if s[0] and (s[0]+b[0]) else -1.
            significance2 += [[significance2Val, significance2Err]]
            significance3Val = (s[0]+b[0])**0.5 - b[0]**0.5
            significance3Err = 0.5 * ((s[1]**2+b[1]**2)/(s[0]+b[0]) + b[1]**2/b[0])**0.5 if (s[0]+b[0]) and b[0] else -1.
            significance3 += [[significance3Val, significance3Err]]
        optimizationVals['cuts'] = cutRange
        optimizationVals['sigPass'] = sigpass
        optimizationVals['bgPass'] = bgpass
        optimizationVals['sigAll'] = sigall
        optimizationVals['bgAll'] = bgall
        optimizationVals['sigEff'] = sigEff
        optimizationVals['bgEff'] = bgEff
        optimizationVals['significance1'] = significance1
        optimizationVals['significance2'] = significance2
        optimizationVals['significance3'] = significance3
        return optimizationVals
