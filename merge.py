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

sys.argv.append('-b')
import ROOT as rt
sys.argv.pop()

rt.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")

from analyzers.ntuples import *

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Run the desired analyzer on "
                                                 "FSA n-tuples")

    parser.add_argument('analysis', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to run')
    parser.add_argument('channel', type=str, choices=['WZ','TT','Hpp3l','Hpp4l','FakeRate'], help='Channel to run for given analysis')
    parser.add_argument('period', type=str, choices=['7','8','13'], help='Energy (TeV)')
    parser.add_argument('jobName',nargs='?',type=str,const='',help='Job Name for condor submission')
    args = parser.parse_args(argv)

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    # merge individual samples
    ntupledir = 'ntuples%s_%stev_%s' % (args.analysis, args.period, args.channel)
    os.system('mkdir -p %s' % ntupledir)
    sampledirs = ['%s/%s' % (args.jobName, name) for name in os.listdir(args.jobName)]
    for sampledir in sampledirs:
        sample = os.path.basename(sampledir)
        ntuplename = '%s/%s.root' % (ntupledir, sample)
        mergeString = 'hadd -f %s %s/*.root' % (ntuplename, sampledir)
        os.system(mergeString)

    # now merge the data samples (checking for duplicate events)
    datasets = {
        '7' : ['Run2011A', 'Run2011B'],
        '8' : ['Run2012A', 'Run2012B', 'Run2012A', 'Run2012D'],
        '13': []
    }
    for dataset in datasets[args.period]:
        print 'Merging dataset %s' % dataset
        # get the trees
        tchain = rt.TChain()
        datafiles = glob.glob('%s/data_*_%s_*.root' % (ntupledir, dataset))
        for f in datafiles:
            tchain.Add('%s/%s' % (f, args.analysis))
        # setup for iterating over tree
        events = set()
        numToSubtract = 0
        numEntries = tchain.GetEntries()
        ntuple, branches = buildNtuple({'a':''},['a'],args.channel,[''])
        event = branches['event']
        tchain.SetBranchAddress('event',rt.AddressOf(event,'evt'))
        # clone tree
        tfile = rt.TFile('%s/data_%s.root' % (ntupledir,dataset), 'recreate')
        tree = tchain.CloneTree(0)
        for i in range(numEntries):
            tchain.GetEntry(i)
            eventkey = (event.run, event.lumi, event.evt)
            if eventkey in events:
                numToSubtract += 1
                continue
            events.add(eventkey)
            tree.Fill()
        # and the cutflows
        # TODO: This is WRONG! only the last entry will be correct... solve later I guesss, perhaps a cutflow tree
        cutflows = {}
        for datafile in datafiles:
            file = rt.TFile(datafile)
            cutflowHist = file.Get('cutflow')
            numBins = cutflowHist.GetNbinsX()-2
            cutflow = []
            for b in range(numBins):
               cutflow += [cutflowHist.GetBinContent(b+1)]
            cutflows[datafile] = cutflow
        cutflow = [sum(x) for x in zip(*cutflows.itervalues())]
        cutflow = [x-numToSubtract for x in cutflow]
        numBins = len(cutflow)
        cutflowHist = rt.TH1F('cutflow','cutflow',numBins,0,numBins)
        for b in range(numBins):
           cutflowHist.SetBinContent(b+1,cutflow[b])
        tfile.cd()
        cutflowHist.Write()
        tfile.Write()
        tfile.Close()

    return 0


if __name__ == "__main__":
    status = main()
    sys.exit(status)
