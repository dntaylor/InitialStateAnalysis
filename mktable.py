#!/usr/bin/env python

from plotters.Plotter import Plotter
from plotters.plotUtils import *
from plotters.plotUtils import ZMASS
from plotters.xsec import xsecs
import argparse
import itertools
import sys

# Want to crate a table of the form
#
# | Sample Group | Sample         | Number of events MC | xsec    | Scaled to 19.7 fb-1 |
# |              |                | Loose ID | Tight ID |         | Loose ID | Tight ID |
# |ZZ            |                | xxxxxxx  | xxxxx    |         | xxxxx    | xxxx     |
# |              | ZZTo4mu        | xxxxx    | xxxxx    | xxxxx   | xxxxx    | xxxx     |
# ....

def table(analysis,channel,period,**kwargs):
    myCut = kwargs.pop('myCut','1')
    mass = kwargs.pop('mass',500)
    scaleFactor = kwargs.pop('scaleFactor','event.pu_weight*event.lep_scale*event.trig_scale')

    numleps = {
        'Hpp2l': 2,
        'Z'    : 2,
        'TT'   : 2,
        'WZ'   : 3,
        'Hpp3l': 3,
        'Hpp4l': 4,
    }
    nl = numleps[analysis]
    ntuples = 'ntuples%s_%stev_%s' % (analysis,period,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,period)
    intLumiMap = getIntLumiMap()
    sigMap = getSigMap(nl,mass)
    channelBackground = {
        'Hpp2l' : ['T', 'TT', 'TTV', 'W', 'Z', 'ZG', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Z'     : ['T', 'TT', 'TTV', 'W', 'Z', 'ZG', 'VVV', 'ZZ', 'WW', 'WZ'],
        'WZ'    : ['T', 'TT', 'TTV', 'Z', 'ZG', 'VVV', 'ZZ', 'WW', 'WZ'],
        'TT'    : ['T', 'TT', 'TTV', 'W', 'Z', 'ZG', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Hpp3l' : ['T', 'TT', 'TTV', 'Z', 'ZG', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Hpp4l' : ['TT', 'Z', 'DB']
    }
    if period==13:
        channelBackground = {
            'WZ'    : ['T', 'TT', 'TTV', 'Z', 'ZZ', 'WZ'],
            'TT'    : ['T', 'TT', 'TTV', 'Z', 'ZZ', 'WZ'],
            'Hpp3l' : ['T', 'TT', 'TTV', 'Z', 'DB'],
            'Hpp4l' : ['T', 'TT', 'Z', 'TTV', 'DB']
        }

    finalStates, leptons = getChannels(nl)
    mergeDict = getMergeDict(period)

    titleString = '| {0:12} | {1:60} | {2:21} | {3:11} | {4:25} |'.format('Sample Group','Sample','Number of events: MC','xsec (pb)','Scaled to {0:4.1f} fb-1'.format(intLumiMap[period]/1000.))
    subtitleString = '| {0:12} | {0:60} | {1:9} | {2:9} | {0:11} | {1:11} | {2:11} |'.format('','Loose ID','Tight ID')
    entryString = '| {group:12} | {sample:60} | {loosemc:9} | {tightmc:9} | {xsec:11} | {loose:11.4f} | {tight:11.4f} |'
    subentryString = '| {group:12} | {sample:60} | {loosemc:9} | {tightmc:9} | {xsec:11.4f} | {loose:11.4f} | {tight:11.4f} |'

    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,scaleFactor=scaleFactor)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in channelBackground[channel]])
    plotter.setIntLumi(intLumiMap[period])

    # print everything out
    print 'Cut to be applied: {0}'.format(myCut)
    print '-'*145
    print titleString
    print subtitleString
    print '-'*145
    for c in channelBackground[channel]:
        numloose = plotter.getNumEntries('{0}'.format(myCut), sigMap[period][c])
        numtight = plotter.getNumEntries('{0} & select.PassTight'.format(myCut), sigMap[period][c])
        print entryString.format(group=sigMap[period][c],loosemc='',tightmc='',loose=numloose,tight=numtight,sample='',xsec='')
        for s in sorted(mergeDict[sigMap[period][c]]):
            sCut = mergeDict[sigMap[period][c]][s]
            nummcloose = plotter.getNumEntries('{0} & {1}'.format(myCut,sCut), s, doUnweighted=True)
            nummctight = plotter.getNumEntries('{0} & {1} & select.PassTight'.format(myCut,sCut), s, doUnweighted=True)
            numloose = plotter.getNumEntries('{0} & {1}'.format(myCut,sCut), s)
            numtight = plotter.getNumEntries('{0} & {1} & select.PassTight'.format(myCut,sCut), s)
            xsec = xsecs[period][s]
            print subentryString.format(sample=s,loosemc=nummcloose,tightmc=nummctight,loose=numloose,tight=numtight,xsec=xsec,group='')
        print '-'*145
        

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Make a table of MC events")

    parser.add_argument('analysis', type=str, choices=['WZ','Hpp3l','Hpp4l','Hpp2l'], help='Analysis')
    parser.add_argument('channel', type=str, choices=['WZ','Hpp3l','Hpp4l','ZZ','TT','Z','Hpp2l'], help='Channel in analysis')
    parser.add_argument('period', type=int, choices=[7,8,13], help='Energy (TeV)')
    parser.add_argument('-c','--cut',type=str,default='1',help='Cut to be applied to events.')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    if args.period == 7:
        print "7 TeV not implemented"

    table(args.analysis,args.channel,args.period,myCut=args.cut)

    return 0


if __name__ == "__main__":
    main()
              
