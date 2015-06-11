#!/usr/bin/env python

import os
import sys
import glob
import argparse
import ROOT
from analyzers.ntuples import buildNtuple

def printEvent(row,treeTuple,**kwargs):
    analysis = kwargs.pop('analysis','Hpp3l')
    eventList = kwargs.pop('eventList',[])
    mode = kwargs.pop('mode','isa')
    channel = kwargs.pop('channel','')
    isMC = kwargs.pop('isMC',False)

    printThisEvent = True
    if mode=='fsa':
        thisEvent = [row.run, row.evt] if isMC else [row.run, row.lumi, row.evt]
    elif mode=='isa':
        eventStruct = treeTuple[2]
        thisEvent = [eventStruct.run, eventStruct.evt] if isMC else [eventStruct.run, eventStruct.lumi, eventStruct.evt]
    else:
        thisEvent = []
    if eventList: # check to see if event is in eventlist
        printThisEvent = True if thisEvent in eventList else False

    if not printThisEvent: return 0

    # add in a way to check if it passed or failed a given cut

    if mode=='fsa': printFSAEvent(row,channel)
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

def printFSAEvent(row,channel):
    leps = enumerate_leps(channel)

    print '-'*80
    print '| Event listing FSA %s |' % (' '*58)
    print '| Run: %7i Lumi: %5i Event: %11i %s |' %(row.run, row.lumi, row.evt, ' '*32)
    print '|%s|' % ('-'*78)
    print '| Channel: %5s Num Vertices: %3i %s |'\
           % (channel, row.nvtx, ' '*43)

    # print lepton info
    for lep in leps:
        iso = getattr(row,'%sRelPFIsoDBDefault' %lep) if lep[0]=='m' else getattr(row,'%sRelPFIsoRho' %lep)
        print '|%s|' % ('-'*78)
        print '| %2s: pT: %11.4f eta: %9.4f phi: %8.4f iso: %7.4f charge: %2i %s |'\
              % (lep, getattr(row,'%sPt' %lep), getattr(row,'%sEta' %lep), getattr(row,'%sPhi' %lep), iso, getattr(row,'%sCharge' %lep), ' '*3)
        if lep[0]=='e':
            print '|     MVA: %10.4f d0: %10.4f dz: %9.4f %s |'\
                  % (getattr(row,'%sMVATrig'%lep), getattr(row,'%sPVDXY'%lep), getattr(row,'%sPVDZ'%lep),' '*27)
        elif lep[0]=='m':
            print '|     ID: %2i d0: %10.4f dz: %9.4f %s |'\
                  % (getattr(row,'%sPFIDTight'%lep), getattr(row,'%sPVDXY'%lep), getattr(row,'%sPVDZ'%lep),' '*30)
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
    parser.add_argument('-mc', action="store_true", help="This is MC (not data)")
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    events = []
    if args.events:
        events = parseEvents(args.events)
        args.n = len(events)
    if args.eventsFile:
        events = parseEventsFile(args.eventsFile)
        args.n = len(events)
    
    files = [filename for string in args.files for filename in glob.glob(string)]

    if args.mode=='isa':
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
        ntuple, branches = buildNtuple(object_definitions,states,channel,final_states,altIds=alternateIds,doVBF=doVBF)

    numPrinted = 0
    for file in files:
        tfile = ROOT.TFile(file)
        if args.mode=='fsa':
            if not args.channel:
                print 'Indicate a channel to print'
                return 0
            fulltree = tfile.Get('%s/final/Ntuple' % args.channel)
            tree = fulltree.CopyTree(args.cut) if args.cut else fulltree
        elif args.mode=='isa':
            fulltree = tfile.Get(args.analysis)
            tree = fulltree.CopyTree(args.cut) if args.cut else fulltree
            tree.SetBranchAddress("select",ROOT.AddressOf(branches['select'],"passTight"))
            tree.SetBranchAddress("event",ROOT.AddressOf(branches['event'],"evt"))
            tree.SetBranchAddress("finalstate",ROOT.AddressOf(branches['finalstate'],"mass"))
        else:
            print 'Unrecognized ntuple type. Valid values are isa or fsa.'
            return 0
        for row in tree:
            if args.list:
                if args.mode=='fsa':
                    print '%i:%i:%i' % (row.run, row.lumi, row.evt) if args.mc else '%i:%i:%i' % (row.run, row.lumi, row.evt)
                else:
                    event = branches['event']
                    print '%i:%i:%i' % (event.run, event.lumi, event.evt) if args.mc else '%i:%i:%i' % (event.run, event.lumi, event.evt)
        tfile.Close()


    return 0


if __name__ == "__main__":
    main()


