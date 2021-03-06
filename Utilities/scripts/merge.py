#!/usr/bin/env python
'''
A script to merge outputs of InitialStateAnalysis run

Author: Devin N. Taylor, UW-Madison
'''

import os
import sys
import glob
import pwd
import argparse
import logging

import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True

rt.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")

from InitialStateAnalysis.Analyzers.ntuples import *
from InitialStateAnalysis.Utilities.utilities import *

def parse_command_line(argv):
    parser = get_parser("Merge the output ISA ntuples")

    parser.add_argument('jobName',nargs='?',type=str,const='',help='Job Name for condor submission')
    parser.add_argument('-d','--directory',type=str,default='',help='Custom subdirectory (to keep more than one ntuple at a time)')
    args = parser.parse_args(argv)

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)

    # merge individual samples
    ntupledir = 'ntuples/%s_%sTeV_%s' % (args.analysis, args.period, args.channel)
    if args.directory: ntupledir += '/{0}'.format(args.directory)
    os.system('mkdir -p %s' % ntupledir)
    if args.jobName:
        sampledirs = ['%s/%s' % (args.jobName, name) for name in os.listdir(args.jobName)]
        for sampledir in sampledirs:
            sample = os.path.basename(sampledir)
            ntuplename = '%s/%s.root' % (ntupledir, sample)
            mergeString = 'hadd -f %s %s/*.root' % (ntuplename, sampledir)
            os.system(mergeString)

    # now merge the data samples (checking for duplicate events)
    datasets = {
        7 : ['Run2011A', 'Run2011B'],
        8 : ['Run2012A', 'Run2012B', 'Run2012C', 'Run2012D'],
        #13: ['Run2015B','Run2015C'], # 50 ns
        13: ['Run2015C','Run2015D'], # 25 ns
    }
    ntuple, branches = buildNtuple({'a':''},['a'],args.channel,[''])
    event = branches['event']
    #cutTree, cutEvents, cutsBranch = buildCutTree(['topology'])
    for dataset in datasets[args.period]:
        logger.info('Merging dataset %s' % dataset)
        # get the trees
        tchain = rt.TChain()
        #cutchain = rt.TChain()
        datafiles = glob.glob('%s/data_*_%s_*.root' % (ntupledir, dataset))
        for f in datafiles:
            logger.info(f)
            tchain.Add('%s/%s' % (f, args.channel))
            #cutchain.Add('%s/cutTree' % (f))
        # setup for iterating over tree
        events = set()
        #cutevents = set()
        numToSubtract = 0
        numEntries = tchain.GetEntries()
        #cutEntries = cutchain.GetEntries()
        tchain.SetBranchAddress('event',rt.AddressOf(event,'gen_weight'))
        #cutchain.SetBranchAddress('event',rt.AddressOf(cutEvents,'evt'))
        # clone tree
        tfile = rt.TFile('%s/data_%s.root' % (ntupledir,dataset), 'recreate')
        tree = tchain.CloneTree(0)
        #newCutTree = cutchain.CloneTree(0)
        for i in range(numEntries):
            tchain.GetEntry(i)
            eventkey = (event.run, event.lumi, event.evt)
            if eventkey in events:
                numToSubtract += 1
                continue
            events.add(eventkey)
            tree.Fill()
        #for i in range(cutEntries):
        #    cutchain.GetEntry(i)
        #    eventkey = (cutEvents.run, cutEvents.lumi, cutEvents.evt)
        #    if eventkey in cutevents: continue
        #    cutevents.add(eventkey)
        #    newCutTree.Fill()
        # and the cutflows
        # TODO: This is WRONG! only the last entry will be correct... solve later I guesss, perhaps a cutflow tree
        #cutflows = {}
        #for datafile in datafiles:
        #    file = rt.TFile(datafile)
        #    cutflowHist = file.Get('cutflow')
        #    numBins = cutflowHist.GetNbinsX()-2
        #    cutflow = []
        #    for b in range(numBins):
        #       cutflow += [cutflowHist.GetBinContent(b+1)]
        #    cutflows[datafile] = cutflow
        #cutflow = [sum(x) for x in zip(*cutflows.itervalues())]
        #cutflow = [x-numToSubtract for x in cutflow]
        #numBins = len(cutflow)
        #cutflowHist = rt.TH1F('cutflow','cutflow',numBins,0,numBins)
        #for b in range(numBins):
        #   cutflowHist.SetBinContent(b+1,cutflow[b])
        tfile.cd()
        #cutflowHist.Write()
        tfile.Write()
        tfile.Close()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
