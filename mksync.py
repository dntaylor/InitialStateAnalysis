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
    doEfficiency = kwargs.pop('doEfficiency',False)
    doCutflow = kwargs.pop('doCutflow',False)
    doYields = kwargs.pop('doYields',False)
    doCounts = kwargs.pop('doCounts',False)
    doCorrelation = kwargs.pop('doCorrelation',False)
    cut = kwargs.pop('cut','1')
    mass = kwargs.pop('mass',500)

    print ''
    print '%s:%s:%iTeV' % (analysis, channel, period)
    print 'Selection to be used:'
    print cut
    if analysis in ['Hpp3l','Hpp4l']:
        print 'Mass %i' % mass
    print ''

    # sync on WZ sample
    # sync on channels: eee, eem, emm, mmm
    ntuples = 'ntuples/%s_%sTeV_%s' % (analysis,period,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,period)
    mergeDict = getMergeDict(period)
    nl = 3 if analysis == 'WZ' or analysis == 'Hpp3l' else 4
    finalStates, leptons = getChannels(nl)
    sigMap = getSigMap(nl,mass)
    channelBackground =  getChannelBackgrounds(period)
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in channelBackground[channel]])
    if analysis in ['Hpp3l', 'Hpp4l']: plotter.initializeSignalSamples([sigMap[period]['Sig']])
    intLumi = getIntLumiMap()[period]
    plotter.setIntLumi(intLumi)
    cutflowMap = defineCutFlowMap(analysis,finalStates,mass)
    cutflow = {}
    for c,l in zip(cutflowMap['cuts'],cutflowMap['labels_simple']):
        cutflow[l] = c
    cutflows = cutflowMap['labels_simple']
    allMC = channelBackground[channel]
    s = 'Sig' if analysis in ['Hpp3l', 'Hpp4l'] else 'WZ'
    if analysis in ['Hpp3l', 'Hpp4l']: allMC += ['Sig']

    if doCounts:
        print '{0} event counts'.format(analysis)
        tempCut = cut
        for chan in finalStates:
            num = plotter.getNumEntries('%s&channel=="%s"' %(tempCut,chan), sigMap[period][s], doUnweighted=True)
            print '%s: %i' % (chan, num)
        print ''

    # yields in each channel
    if doYields:
        print 'Yields (scaled to %i pb-1)' % intLumi
        #theCut = '&'.join([cutflow[x] for x in cutflows])
        print 'Selection applied: %s' % cut
        for chan in finalStates:
            print '%8s |    Channel |      Yield' % chan
            for b in allMC:
                val = plotter.getNumEntries('%s&channel=="%s"' %(cut,chan), sigMap[period][b])
                print '         | %10s | %10.2f' %(b,val)
            print ''

    # cut flow
    if doCutflow:
        print 'Cutflows'
        for chan in finalStates:
            print '%15s |        Sig |         BG |        S/B' % chan
            for c in range(len(cutflows)):
                sig = 0
                bg = 0
                theCut = '&'.join([cutflow[x] for x in cutflows[0:c+1]])
                for b in allMC:
                    val = plotter.getNumEntries('%s&channel=="%s"' %(theCut,chan), sigMap[period][b])
                    if b==s: sig += val
                    else: bg += val
                print '%15s | %10.2f | %10.2f | %10.2f' %(cutflows[c],sig,bg,sig/bg)
            print ''

    # efficiency of cuts
    if doEfficiency:
        print 'Cut efficiencies'
        sigPreCuts = {}
        bgPreCuts = {}
        sigFullCuts = {}
        bgFullCuts = {}
        for chan in finalStates:
            sigPre = 0
            bgPre = 0
            sigFull = 0
            bgFull = 0
            theFullCut = '&'.join([cutflow[x] for x in cutflows])
            for b in allMC:
                valPre = plotter.getNumEntries('%s&channel=="%s"' %(cut,chan), sigMap[period][b])
                valFull = plotter.getNumEntries('%s&channel=="%s"' %(theFullCut,chan), sigMap[period][b])
                if b==s:
                    sigPre += valPre
                    sigFull += valFull
                else:
                    bgPre += valPre
                    bgFull += valFull
            sigPreCuts[chan] = sigPre
            bgPreCuts[chan] = bgPre
            sigFullCuts[chan] = sigFull
            bgFullCuts[chan] = bgFull

        for c in cutflows[1:]:
            print '%15s |  Sg Pre Eff |  BG Pre Eff | Sg Post Eff | BG Post Eff' % c
            for chan in finalStates:
                sigAllbut = 0
                bgAllbut = 0
                sigOnly = 0
                bgOnly = 0
                theCut = '&'.join([cutflow[x] for x in cutflows if x != c])
                theOnlyCut = '%s&%s' % (cut, cutflow[c])
                for b in allMC:
                    valAllbut = plotter.getNumEntries('%s&channel=="%s"' %(theCut,chan), sigMap[period][b])
                    valOnly = plotter.getNumEntries('%s&channel=="%s"' %(theOnlyCut,chan), sigMap[period][b])
                    if b==s:
                        sigAllbut += valAllbut
                        sigOnly += valOnly
                    else:
                        bgAllbut += valAllbut
                        bgOnly += valOnly
                sigEffPre = sigOnly/sigPreCuts[chan] if sigPreCuts[chan] else -1
                bgEffPre = bgOnly/bgPreCuts[chan] if bgPreCuts[chan] else -1
                sigEffPost = sigFullCuts[chan]/sigAllbut if sigAllbut else -1
                bgEffPost = bgFullCuts[chan]/bgAllbut if bgAllbut else -1
                print '%15s | %11.4f | %11.4f | %11.4f | %11.4f' % (chan, sigEffPre, bgEffPre, sigEffPost, bgEffPost)
            print ''

    if doCorrelation:
        print 'Correlation matrix'
        for chan in finalStates:
            print ' | '.join(['%15s'%(chan+' denom.')] + ['%15s' % x for x in cutflows[1:]])
            for dc in cutflows[1:]:
                row = []
                for nc in cutflows[1:]:
                    numCut = '%s & %s & %s' % (cut, cutflow[dc], cutflow[nc])
                    denomCut = '%s & %s' % (cut, cutflow[dc])
                    numSig = plotter.getNumEntries('%s&channel=="%s"' %(numCut,chan), sigMap[period][s])
                    denomSig = plotter.getNumEntries('%s&channel=="%s"' %(denomCut,chan), sigMap[period][s])
                    eff = numSig/denomSig if denomSig else -1
                    row += [eff]
                print ' | '.join(['%15s'%dc] + ['%15.4f' % x for x in row])
            print ''


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot a given channel and period")

    parser.add_argument('analysis', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to plot')
    parser.add_argument('channel', type=str, choices=['WZ','Hpp3l','Hpp4l','FakeRate'], help='Channel in analysis')
    parser.add_argument('period', type=int, choices=[7,8,13], help='Energy (TeV)')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500)
    parser.add_argument('-rt','--runTau',action='store_true',help='Run Tau finalStates (not implemented)')
    parser.add_argument('-ub','--unblind',action='store_false',help='Unblind signal channel')
    parser.add_argument('-dc','--doCounts',action='store_true',help='run sig counts')
    parser.add_argument('-dcf','--doCutflow',action='store_true',help='run cutflow')
    parser.add_argument('-de','--doEfficiency',action='store_true',help='run efficiencies')
    parser.add_argument('-dy','--doYields',action='store_true',help='run yields')
    parser.add_argument('-dco','--doCorrelation',action='store_true',help='run correlation')
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

    sync(args.analysis,args.channel,args.period,
         runTau=args.runTau,
         blind=args.unblind,
         cut=args.cut,
         mass=args.mass,
         doEfficiency=args.doEfficiency,
         doCutflow=args.doCutflow,
         doYields=args.doYields,
         doCounts=args.doCounts,
         doCorrelation=args.doCorrelation
    )

    return 0


if __name__ == "__main__":
    main()
