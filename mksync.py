#!/usr/bin/env python

from plotters.Plotter import Plotter
from plotters.plotUtils import *
import argparse
import itertools
import sys

ZMASS = 91.1876

def sync(analysis,channel,period,**kwargs):
    '''Print sync information to file.'''
    runTau = kwargs.pop('runTau',False)
    blind = kwargs.pop('blind',True)
    doBjetVeto = kwargs.pop('doBjetVeto',False)
    doEfficiency = kwargs.pop('doEfficiency',False)
    doCutflow = kwargs.pop('doCutflow',False)
    cut = kwargs.pop('cut','1')

    # WZ only for now
    if not (analysis == channel == 'WZ'): return

    fs = ['eee','eem','mme','mmm']

    print ''
    print '%s:%s:%iTeV' % (analysis, channel, period)
    print 'Selection to be used:'
    print cut
    if doBjetVeto:
        print 'Bjet veto to be applied: finalstate.bjetVeto30==0'
    print ''

    # sync on WZ sample
    # sync on channels: eee, eem, emm, mmm
    ntuples = 'ntuples%s_%stev_%s' % (analysis,period,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,period)
    mergeDict = getMergeDict(period)
    nl = 3 if analysis == 'WZ' or analysis == 'Hpp3l' else 4
    sigMap = getSigMap(nl,0)
    channelBackground = {
        'WZ' : ['T','TT', 'TTV', 'Z', 'ZZ','WZ'],
        'Hpp3l' : ['T', 'TT', 'TTV','Z','DB'],
        'Hpp4l' : ['T', 'TT', 'Z', 'TTV','DB']
    }
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in channelBackground[channel]])
    intLumi = 19700
    plotter.setIntLumi(intLumi)
    cutflow = {
        'pre' : cut,
        'bjet' : 'finalstate.bjetVeto30Medium==0',
        'zpt' : '(z1.Pt1>20.&z1.Pt2>10.)',
        'zmass' : 'fabs(z1.mass-%f)<20.' % ZMASS,
        'wpt' : 'w1.Pt1>20.',
        'wdr' : 'w1.dR1_z1_1>0.1 & w1.dR1_z1_2>0.1',
        'met' : 'w1.met>30.',
        'm3l' : 'finalstate.Mass>100.'
    }
    cutflowMap = defineCutFlowMap('WZ',0,0)
    cutflows = ['pre','zpt','zmass','wdr','wpt','met']
    if doBjetVeto: cutflows = ['pre','zpt','zmass','bjet','wdr','wpt','met']

    print 'WZ event counts'
    tempCut = cut
    if doBjetVeto: tempCut += ' & finalstate.bjetVeto30Medium==0'
    for chan in fs:
        num = plotter.getNumEntries('%s&channel=="%s"' %(tempCut,chan), sigMap[period]['WZ'], doUnweighted=True)
        print '%s: %i' % (chan, num)
    print ''

    # yields in each channel
    print 'Yields (scaled to %i pb-1)' % intLumi
    theCut = '&'.join([cutflow[x] for x in cutflows])
    print 'Selection applied: %s' % theCut
    for chan in fs:
        print '%8s |    Channel |      Yield' % chan
        for b in channelBackground[channel]:
            val = plotter.getNumEntries('%s&channel=="%s"' %(theCut,chan), sigMap[period][b])
            print '         | %10s | %10.2f' %(b,val)
        print ''

    # cut flow
    if doCutflow:
        print 'Cutflows'
        for chan in fs:
            print '%8s |         WZ |         BG |        S/B' % chan
            for c in range(len(cutflows)):
                wz = 0
                bg = 0
                theCut = '&'.join([cutflow[x] for x in cutflows[0:c+1]])
                for b in channelBackground[channel]:
                    val = plotter.getNumEntries('%s&channel=="%s"' %(theCut,chan), sigMap[period][b])
                    if b=='WZ': wz += val
                    else: bg += val
                print '%8s | %10.2f | %10.2f | %10.2f' %(cutflows[c],wz,bg,wz/bg)
            print ''

    # efficiency of cuts
    if doEfficiency:
        print 'Cut efficiencies'
        wzPreCuts = {}
        bgPreCuts = {}
        wzFullCuts = {}
        bgFullCuts = {}
        for chan in fs:
            wzPre = 0
            bgPre = 0
            wzFull = 0
            bgFull = 0
            theFullCut = '&'.join([cutflow[x] for x in cutflows])
            for b in channelBackground[channel]:
                valPre = plotter.getNumEntries('%s&channel=="%s"' %(cut,chan), sigMap[period][b])
                valFull = plotter.getNumEntries('%s&channel=="%s"' %(theFullCut,chan), sigMap[period][b])
                if b=='WZ':
                    wzPre += valPre
                    wzFull += valFull
                else:
                    bgPre += valPre
                    bgFull += valFull
            wzPreCuts[chan] = wzPre
            bgPreCuts[chan] = bgPre
            wzFullCuts[chan] = wzFull
            bgFullCuts[chan] = bgFull

        for c in cutflows[1:]:
            print '%8s |  WZ Pre Eff |  BG Pre Eff | WZ Post Eff | BG Post Eff' % c
            for chan in fs:
                wzAllbut = 0
                bgAllbut = 0
                wzOnly = 0
                bgOnly = 0
                theCut = '&'.join([cutflow[x] for x in cutflows if x != c])
                theOnlyCut = '%s&%s' % (cut, cutflow[c])
                for b in channelBackground[channel]:
                    valAllbut = plotter.getNumEntries('%s&channel=="%s"' %(theCut,chan), sigMap[period][b])
                    valOnly = plotter.getNumEntries('%s&channel=="%s"' %(theOnlyCut,chan), sigMap[period][b])
                    if b=='WZ':
                        wzAllbut += valAllbut
                        wzOnly += valOnly
                    else:
                        bgAllbut += valAllbut
                        bgOnly += valOnly
                wzEffPre = wzOnly/wzPreCuts[chan]
                bgEffPre = bgOnly/bgPreCuts[chan]
                wzEffPost = wzFullCuts[chan]/wzAllbut
                bgEffPost = bgFullCuts[chan]/bgAllbut
                print '%8s | %11.4f | %11.4f | %11.4f | %11.4f' % (chan, wzEffPre, bgEffPre, wzEffPost, bgEffPost)
            print ''



def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot a given channel and period")

    parser.add_argument('analysis', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to plot')
    parser.add_argument('channel', type=str, choices=['WZ','Hpp3l','Hpp4l','FakeRate'], help='Channel in analysis')
    parser.add_argument('period', type=int, choices=[7,8,13], help='Energy (TeV)')
    parser.add_argument('-rt','--runTau',action='store_true',help='Run Tau finalStates (not implemented)')
    parser.add_argument('-ub','--unblind',action='store_false',help='Unblind signal channel')
    parser.add_argument('-db','--doBjetVeto',action='store_true',help='Add a bjet veto')
    parser.add_argument('-c','--cut',type=str,default='select.passTight',help='Cut to be applied to plots (default = "select.passTight").')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    if args.period == 7:
        print "7 TeV not implemented"

    # loose no iso
    # e: loose CBID
    # m: isLooseMuon

    # loose w/ iso
    # e: loose CBID and iso<0.2
    # m: isLooseMuon and iso<0.2

    # tight
    # e: medium CBID and iso<0.15
    # m: isTightMuon and iso<0.12

    sync(args.analysis,args.channel,args.period,runTau=args.runTau,blind=args.unblind,cut=args.cut,doBjetVeto=args.doBjetVeto)

    return 0


if __name__ == "__main__":
    main()
