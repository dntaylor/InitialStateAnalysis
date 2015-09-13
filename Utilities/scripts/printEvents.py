#!/usr/bin/env python

import os
import sys
import glob
import argparse
import logging

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

from InitialStateAnalysis.Analyzers.ntuples import buildNtuple

def printEvent(row,treeTuple,**kwargs):
    analysis = kwargs.pop('analysis','Hpp3l')
    eventList = kwargs.pop('eventList',[])
    mode = kwargs.pop('mode','isa')
    channel = kwargs.pop('channel','')
    isMC = kwargs.pop('isMC',False)
    branches = kwargs.pop('branches',{})

    printThisEvent = True
    if eventList: # check to see if event is in eventlist
        if mode=='fsa':
            thisEvent = [row.run, row.lumi, row.evt]
        elif mode=='isa':
            eventStruct = branches['event']
            thisEvent = [eventStruct.run, eventStruct.lumi, eventStruct.evt]
        else:
            thisEvent = []
        printThisEvent = True if thisEvent in eventList else False

    if not printThisEvent: return 0

    # add in a way to check if it passed or failed a given cut

    if mode=='fsa': printFSAEvent(row,channel)
    if mode=='isa': printISAEvent(row,treeTuple,branches,analysis)
    return 1
    
def enumerate_leps(final_state):
    '''Get the leptons available in the ntuple for a given final state'''
    out = []
    for i  in ['e', 'm', 't']:
        N = final_state.count(i)
        if N==1:
           out += i
        else:
           out += ['%s%i' % (i, n) for n in xrange(1, N+1)]
    return out

def printDetailed(row,tree,branches,analysis):
    '''Print detailed information in list form'''
    evt = branches['event']
    z = branches['z1']
    w = branches['w1']
    fs = branches['finalstate']
    strToPrint = '{run}:{lumi}:{evt}'
    strToPrint += ':{z1pt:.4f}:{z1eta:.4f}:{z1phi:.4f}:{z1iso:.4f}'
    strToPrint += ':{z2pt:.4f}:{z2eta:.4f}:{z2phi:.4f}:{z2iso:.4f}'
    strToPrint += ':{w1pt:.4f}:{w1eta:.4f}:{w1phi:.4f}:{w1iso:.4f}'
    strToPrint += ':{zdr:.4f}:{z1w1dr:.4f}:{z2w1dr:.4f}'
    strToPrint += ':{zmass:.4f}:{met:.4f}:{metphi:.4f}:{m3l:.4f}'
    print strToPrint.format(run=evt.run, lumi=evt.lumi, evt=evt.evt,
                            z1pt=z.Pt1, z1eta=z.Eta1, z1phi=z.Phi1, z1iso=z.Iso1,
                            z2pt=z.Pt2, z2eta=z.Eta2, z2phi=z.Phi2, z2iso=z.Iso2,
                            w1pt=w.Pt1, w1eta=w.Eta1, w1phi=w.Phi1, w1iso=w.Iso1,
                            zdr=z.dR, z1w1dr=w.dR1_z1_1, z2w1dr=w.dR1_z1_2,
                            zmass=z.mass, met=fs.met, metphi=fs.metPhi, m3l=fs.mass)

def printISAEvent(row,tree,branches,analysis):
    leps = ['l1','l2','l3']
    bosons = ['z1','w1']
    
    event = branches['event']
    channel = branches['channel']
    finalstate = branches['finalstate']
    print '-'*80
    print '| Event listing ISA {0:58} |'.format('')
    print '| Run: {0:7} Lumi: {1:5} Event: {2:11} {3:32} |'.format(event.run, event.lumi, event.evt, '')
    print '|{0}|'.format('-'*78)
    #print '| Channel: {0: <6} Num Vertices: {1:3} {2:45} |'.format(str(channel.channel), event.nvtx, '')
    print '| Channel: {0: <6} Num Vertices: {1:3} {2:42} |'.format('', event.nvtx, '')
    print '| Mass: {0:9.4f} MET: {1:9.4f} MET phi: {2:9.4f} {3:26} |'.format(finalstate.mass, finalstate.met, finalstate.metPhi, '')

    # print lepton info
    lep_string = '| {0: <2}: pT: {1:11.4f} eta: {2:9.4f} phi: {3:8.4f} iso: {4:7.4f} charge: {5:2} {6:3} |'
    lep_string_2 = '| Pass Loose: {0:1} Pass Tight: {1:1} {2:48} |'
    for lep in leps:
        l = branches[lep]
        #lf = branches['%sFlv'%lep]
        print '|{0}|'.format('-'*78)
        print lep_string.format(lep, l.Pt, l.Eta, l.Phi, l.Iso, l.Chg, '')
        print lep_string_2.format(l.PassLoose, l.PassTight, '')

    # print boson info
    boson_string = '| {0: <2}: Mass: {1:9.4f} Pt: {2:9.4f} {3:42} |'
    boson_1 = '| Ptl1: {0:9.4f} dRZl1: {1:6.4f} dRZl2: {2:6.4f} {3:32} |'
    boson_2 = '| Ptl1: {0:9.4f} Ptl2: {1:9.4f} dR: {2:6.4f} {3:33} |'
    for boson in bosons:
        b = branches[boson]
        #bf = branches['%sFlv'%boson]
        print '|{0}|'.format('-'*78)
        print boson_string.format(boson, b.mass, b.Pt, '')
        if boson in ['z1']:
            print boson_2.format(b.Pt1, b.Pt2, b.dR, '')
        if boson in ['w1']:
            print boson_1.format(b.Pt1, b.dR1_z1_1, b.dR1_z1_2, '')
    print '-'*80
    print ''


def printFSAEvent(row,channel):
    leps = enumerate_leps(channel)

    print '-'*80
    print '| Event listing FSA {0:58} |'.format('')
    print '| Run: {0:7} Lumi: {1:5} Event: {2:11} {3:32} |'.format(row.run, row.lumi, row.evt, '')
    print '|{0}|'.format('-'*78)
    print '| Channel: {0:5} Num Vertices: {1:3} {2:42} |'.format(channel, row.nvtx, '')

    # print lepton info
    lep_string = '| {0:2}: pT: {1:11.4f} eta: {2:9.4f} phi: {3:8.4f} iso: {4:7.4f} charge: {5:4} {6:1} |'
    elec_string = '|     CBID-Med: {0:1f} d0: {1:10.4f} dz: {2:9.4f} {3:22} |'
    muon_string = '|     ID: {0:2} d0: {1:10.4f} dz: {2:9.4f} {3:30} |'
    for lep in leps:
        iso = getattr(row,'%sRelPFIsoDBDefault' %lep) if lep[0]=='m' else getattr(row,'%sRelPFIsoRho' %lep)
        print '|{0}|'.format('-'*78)
        print lep_string.format(lep, getattr(row,'%sPt' %lep), getattr(row,'%sEta' %lep), getattr(row,'%sPhi' %lep), iso, getattr(row,'%sCharge' %lep), '')
        if lep[0]=='e':
            print elec_string.format(getattr(row,'%sCBIDMedium'%lep), getattr(row,'%sPVDXY'%lep), getattr(row,'%sPVDZ'%lep),'')
        elif lep[0]=='m':
            print muon_string.format(getattr(row,'%sPFIDTight'%lep), getattr(row,'%sPVDXY'%lep), getattr(row,'%sPVDZ'%lep),'')
    print '-'*80
    print ''

def parseEvents(eventString):
    '''Parse and event string into a list of event tuples'''
    eventString = ''.join(eventString.split())
    eventStringList = eventString.split(',')
    eventTuples = [map(int,eString.split(':')) for eString in eventStringList]
    return eventTuples

def parseEventsFile(filename):
    '''Parse and event string into a list of event tuples'''
    with open(filename,'r') as file:
        eventTuples = []
        for e in file.readlines():
            e = e.strip()
            if e[-1] == ':': e = e[:-1]
            eventTuples += [map(int,e.split(':'))]
    return eventTuples

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Print events from ntuple (FSA or DBLH)")

    parser.add_argument('files', nargs='+', help='File names w/ UNIX wildcards')
    parser.add_argument('-e','--events',nargs='?',type=str,default='',help='Comma separated list of events to print. (run:lumi:event, ...)')
    parser.add_argument('-ef','--eventsFile',nargs='?',type=str,default='',help='file with event numbers, one per line')
    parser.add_argument('-c','--cut',nargs='?',type=str,default='',help='Cut to be applied to tree')
    parser.add_argument('-ch','--channel',nargs='?',type=str,default='',help='Channel for FSA: eee, eeee, eem, ...')
    parser.add_argument('-a','--analysis',nargs='?',type=str,default='Hpp3l', help='Analysis for ISA')
    parser.add_argument('-m','--mode',nargs='?',type=str,default='isa', help='Ntuple type: isa, fsa')
    parser.add_argument('-n',nargs='?',default=1,help='Number of events to print (default 1, -1 for all)')
    parser.add_argument('-l','--list', action="store_true", help="List events passing selection")
    parser.add_argument('-d','--detailed', action="store_true", help="List events passing selection with values")
    parser.add_argument('-mc', action="store_true", help="This is MC (not data)")
    parser.add_argument('--log',nargs='?',type=str,const='INFO',default='INFO',choices=['INFO','DEBUG','WARNING','ERROR','CRITICAL'],help='Log level for logger')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')

    events = []
    if args.events:
        events = parseEvents(args.events)
        args.n = len(events)
        logging.debug('List of events with {0} entries detected'.format(args.n))
    if args.eventsFile:
        events = parseEventsFile(args.eventsFile)
        args.n = len(events)
        logging.debug('Events file with {0} entries detected'.format(args.n))
    
    files = [filename for string in args.files for filename in glob.glob(string)]
    logging.debug('Running over {0} files.'.format(len(files)))

    if args.mode=='isa':
        logging.debug('Ntuple is declared to be ISA ntuple for {0} analysis.'.format(args.analysis))
        if args.analysis == 'WZ':
            channel = 'WZ'
            final_states = ['eee','eem','emm','mmm']
            initial_states = ['z1','w1'] # in order of leptons returned in choose_objects
            object_definitions = {
                'w1': ['em','n'],
                'z1': ['em','em'],
            }
            states = [initial_states]
            alternateIds = []
            doVBF = False
        if args.analysis == 'Hpp3l':
            channel = 'Hpp3l'
            final_states = ['eee','eem','emm','mmm'] # no tau
            initial_states = ['h1','h2']
            other_states = [['z1', 'w1']]
            states = other_states + [initial_states]
            object_definitions = {
                'h1': ['em','em'],
                'h2': ['em','n'],
                'z1': ['em','em'],
                'w1': ['em','n'],
            }
            alternateIds = []
            doVBF = False

        logging.debug('Loading dummy file.')
        dummyfile = ROOT.TFile('dummy.root','recreate')
        logging.debug('Retrieving ntuple.')
        ntuple, branches = buildNtuple(object_definitions,states,channel,final_states,altIds=alternateIds,doVBF=doVBF)

    numPrinted = 0
    for file in files:
        if numPrinted >= args.n and args.n != -1:
            logging.debug('Reached max events (file)')
            break
        logging.debug('Processing file {0}.'.format(file))
        tfile = ROOT.TFile(file)
        if args.mode=='fsa':
            if not args.channel:
                print 'Indicate a channel to print'
                return 0
            fulltree = tfile.Get('%s/final/Ntuple' % args.channel)
            logging.debug('Copying tree with cut {0}.'.format(args.cut))
            tree = fulltree.CopyTree(args.cut) if args.cut else fulltree
        elif args.mode=='isa':
            fulltree = tfile.Get(args.analysis)
            dummyfile.cd()
            logging.debug('Copying tree with cut {0}.'.format(args.cut))
            tree = fulltree.CopyTree(args.cut) if args.cut else fulltree
            tree.SetBranchAddress("select",ROOT.AddressOf(branches['select'],"passTight"))
            tree.SetBranchAddress("event",ROOT.AddressOf(branches['event'],"evt"))
            tree.SetBranchAddress("channel",ROOT.AddressOf(branches['channel'],"channel"))
            tree.SetBranchAddress("genChannel",ROOT.AddressOf(branches['genChannel'],"channel"))
            tree.SetBranchAddress("finalstate",ROOT.AddressOf(branches['finalstate'],"mass"))
            tree.SetBranchAddress("l1",ROOT.AddressOf(branches['l1'],"Pt"))
            tree.SetBranchAddress("l1Flv",ROOT.AddressOf(branches['l1Flv'],"Flv"))
            tree.SetBranchAddress("l2",ROOT.AddressOf(branches['l2'],"Pt"))
            tree.SetBranchAddress("l2Flv",ROOT.AddressOf(branches['l2Flv'],"Flv"))
            tree.SetBranchAddress("l3",ROOT.AddressOf(branches['l3'],"Pt"))
            tree.SetBranchAddress("l3Flv",ROOT.AddressOf(branches['l3Flv'],"Flv"))
            tree.SetBranchAddress("z1",ROOT.AddressOf(branches['z1'],"mass"))
            tree.SetBranchAddress("z1Flv",ROOT.AddressOf(branches['z1Flv'],"Flv"))
            tree.SetBranchAddress("w1",ROOT.AddressOf(branches['w1'],"mass"))
            tree.SetBranchAddress("w1Flv",ROOT.AddressOf(branches['w1Flv'],"Flv"))
        else:
            logging.error('Unrecognized ntuple type. Valid values are isa or fsa.')
            return 0
        logging.debug('Iterate through tree')
        for row in tree:
            if args.list:
                if args.mode=='fsa':
                    print '%i:%i:%i' % (row.run, row.lumi, row.evt)
                else:
                    if args.detailed:
                        printDetailed(row,tree,branches,args.analysis)
                    else:
                        event = branches['event']
                        print '%i:%i:%i' % (event.run, event.lumi, event.evt)
            else:
                if args.mode=='isa':
                    logging.debug('Print event')
                    numPrinted += printEvent(row,tree,analysis=args.analysis,eventList=events,branches=branches,isMC=args.mc,mode=args.mode)
                    logging.debug('Event printed')
                if args.mode=='fsa':
                    numPrinted += printEvent(row,tree,channel=args.channel,eventList=events,isMC=args.mc,mode=args.mode)
            if numPrinted >= args.n and args.n != -1:
                logging.debug('Reached max events (tree)')
                break
        logging.debug('Closing file {0}.'.format(file))
        tfile.Close("R")

    if args.mode == 'isa':
        logging.debug('Closing dummy file.')
        dummyfile.Close("R")
        logging.debug('Deleting dummy file.')
        os.remove('dummy.root')

    return 0


if __name__ == "__main__":
    main()


