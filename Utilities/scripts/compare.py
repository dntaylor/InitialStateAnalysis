#!/usr/bin/env python
'''
A script to compare two trees of the same structure.

Author: Devin N. Taylor, UW-Madison
'''

import os
import sys
import errno
import glob
import pwd
import subprocess
import argparse
import logging

import ROOT as rt
rt.PyConfig.IgnoreCommandLineOptions = True
rt.gROOT.SetBatch(True)

import InitialStateAnalysis.Plotters.CMS_lumi as CMS_lumi
import InitialStateAnalysis.Plotters.tdrstyle as tdrstyle

rt.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
tdrstyle.setTDRStyle()
rt.gStyle.SetOptTitle(1)

def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise

def get_leaves(obj):
    leaves = []
    objs = obj.GetListOfBranches()
    for o in objs:
        name = o.GetName()
        #logging.debug('Searching object {0}.'.format(name))
        for l in get_leaves(o):
            leaf = '{0}.{1}'.format(name,l)
            leaves += [leaf]
    if not leaves:
        for l in obj.GetListOfLeaves():
            name = l.GetName()
            #logging.debug('Adding leaf {0}'.format(name))
            leaves += [name]
    return leaves

def get_hist(tree,var,cut,label,**kwargs):
    nbins = kwargs.pop('nbins',0)
    minbin = kwargs.pop('minbin',0)
    maxbin = kwargs.pop('maxbin',0)
    color = kwargs.pop('color',0)
    histname = 'h_{0}_{1}'.format(var,label)
    if any([nbins,minbin,maxbin]):
        # bin to match others
        drawstring = '{0}>>{1}({2})'.format(var,histname,','.join([str(x) for x in [nbins,minbin,maxbin]]))
    else:
        drawstring = '{0}>>{1}(100)'.format(var,histname)
    tree.Draw(drawstring,cut,'goff')
    if not rt.gDirectory.Get(histname):
        logging.debug('No histogram for variable {0} in tree {1}.'.format(var,label))
        if any(nbins,minbin,maxbin):
            return rt.TH1F(label,label,nbins,minbin,maxbin)
        else:
            return rt.TH1F()
    hist = rt.gDirectory.Get(histname)
    hist.SetTitle(var)
    hist.SetName(label)
    hist.SetLineColor(color)
    return hist

def get_ratio(num, denom, label):
    '''Return a ratio histogram'''
    ratio = num.Clone(label)
    ratio.Sumw2()
    ratio.SetMarkerSize(0.8)
    ratio.Divide(num, denom, 1., 1., "")
    return ratio

def make_plots(*trees,**kwargs):
    cuts = kwargs.pop('cuts',[])
    labels = kwargs.pop('labels',[])
    outdir = kwargs.pop('outdir','comparisonPlots')

    varsToPlot = set()
    for tree in trees:
        varsToPlot.update(get_leaves(tree))

    if len(cuts) == 0:
        cuts = ['1'] * len(trees)
    elif len(cuts) == 1:
        cuts = cuts * len(trees)

    if len(trees) == 1:
        trees = trees * len(cuts)

    if len(labels) != len(trees):
        labels = ['Reference'] + ['Comparison_{0}'.format(x) for x in range(len(trees)-1)]

    colors = [rt.kBlack, rt.kRed, rt.kBlue, rt.kGreen, rt.kViolet, rt.kOrange]

    for var in sorted(varsToPlot):
        logging.debug('Preparing {0}'.format(var))
        canvas = rt.TCanvas(var,'var')
        canvas.SetCanvasSize(796,666)
        plotpad = rt.TPad("plotpad", "top pad", 0.0, 0.21, 1.0, 1.0)
        plotpad.SetLeftMargin(0.12)
        plotpad.SetRightMargin(0.01)
        plotpad.SetTopMargin(0.0875)
        plotpad.SetBottomMargin(28./666.)
        plotpad.SetTickx(1)
        plotpad.SetTicky(1)
        plotpad.Draw()
        ratiopad = rt.TPad("ratiopad", "bottom pad", 0.0, 0.0, 1.0, 0.21)
        ratiopad.SetTopMargin(0.)
        ratiopad.SetBottomMargin(0.5)
        ratiopad.SetLeftMargin(0.12)
        ratiopad.SetRightMargin(0.01)
        ratiopad.SetFillColor(0)
        ratiopad.SetTickx(1)
        ratiopad.SetTicky(1)
        ratiopad.Draw()
        plotpad.cd()
        hists = []
        nbins = 0
        minbin = 0
        maxbin = 0
        for t,c,l,col in zip(trees,cuts,labels,colors[:len(trees)]):
            if hists: # use first as reference
                nbins = hists[0].GetNbinsX()
                minbin = hists[0].GetXaxis().GetXmin()
                maxbin = hists[0].GetXaxis().GetXmax()
                logging.debug('Binning: {0},{1},{2}.'.format(nbins,minbin,maxbin))
                hists += [get_hist(t,var,c,l,color=col,nbins=nbins,minbin=minbin,maxbin=maxbin)]
            else:
                logging.debug('First histogram')
                hists += [get_hist(t,var,c,l,color=col)]
        ref = hists[0]
        ref.GetXaxis().SetLabelOffset(999)
        ref.Draw('hist')
        diffs = [0.] * (len(hists)-1)
        for h in range(len(hists)-1):
            hist = hists[h+1]
            hist.Draw('hist same')
            for b in range(ref.GetNbinsX()):
                diffs[h] += abs(ref.GetBinContent(b+1)-hist.GetBinContent(b+1))
        ratiopad.cd()
        ratiopad.SetGridy(0)
        ratiounity = rt.TLine(ref.GetXaxis().GetXmin(),1,ref.GetXaxis().GetXmax(),1)
        ratiounity.SetLineStyle(2)
        ratios = []
        for h in range(len(hists)-1):
            hist = hists[h+1]
            ratios += [get_ratio(hist,ref,'{0}_ratio'.format(labels[h+1]))]
            ratios[h].SetTitle("")
            ratios[h].GetYaxis().SetTitle("Ratio")
            ratios[h].SetMaximum(2.)
            ratios[h].SetMinimum(0.)
            ratios[h].SetMarkerSize(0)
            ratios[h].GetXaxis().SetLabelSize(0.19)
            ratios[h].GetXaxis().SetTitleSize(0.21)
            ratios[h].GetXaxis().SetTitleOffset(1.0)
            ratios[h].GetYaxis().SetLabelSize(0.19)
            ratios[h].GetYaxis().SetTitleSize(0.21)
            ratios[h].GetYaxis().SetTitleOffset(0.27)
            ratios[h].GetYaxis().SetNdivisions(503)
            ratios[h].Draw('hist same')
        ratiounity.Draw('same')
        plotpad.cd()
        leg = rt.TLegend(0.65,0.90-len(hists)*0.045,0.95,0.90,'','NDC')
        leg.SetTextFont(42)
        #leg.SetTextSize(0.25/len(hists))
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        for hist in hists:
            leg.AddEntry(hist,hist.GetName(),'f')
        leg.Draw()
        canvas.cd()
        if any(diffs):
            outfile = '{0}/{1}.png'.format(outdir,var)
            logging.info('Saving {0} to file.'.format(outfile))
            python_mkdir(outdir)
            canvas.Print(outfile)
        canvas.Clear()

def parse_command_line(argv):
    parser = argparse.ArgumentParser(description='Compare two ntuples')

    parser.add_argument('ntuples', nargs='+',help='Files/directories to compare (root files in directories will be merged).')
    parser.add_argument('-t','--tree',type=str,required=True,help='Name of tree in root file')
    parser.add_argument('-c','--cut',default=[],action='append',help='Cut(s) to be applied')
    parser.add_argument('-n','--name',default=[],action='append',help='Labels for each ntuple/cut pair')
    parser.add_argument('-l','--log',nargs='?',type=str,const='INFO',default='INFO',choices=['INFO','DEBUG','WARNING','ERROR','CRITICAL'],help='Log level for logger')
    args = parser.parse_args(argv)

    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')

    # first do some error checking
    numFiles = len(args.ntuples)
    numCuts = len(args.cut)
    if numFiles==1 and numCuts in [0,1]:
        logging.error('Nothing to compare (1 ntuple and 1 cut).')
        return 0
    elif numFiles>1 and numCuts>1 and numFiles!=numCuts:
        logging.error('Mismatch of multiple ntuples and cuts: {0} ntuples; {1} cuts.'.format(numFiles,numCuts))
        return 0
    else:
        logging.debug('Will run over {0} file(s) with {0} cut(s).'.format(numFiles,numCuts))

    # error check names
    numPlots = max(numFiles, numCuts)
    numLabels = len(args.name)
    if numLabels != numPlots:
        if numLabels: logging.warning('Mismatch of labels and plots: {0} labels; {1} plots. Will run without custom labels.'.format(numLabels,numPlots))
        args.name = ['Reference'] + ['Comparison_{0}'.format(x) for x in range(numPlots-1)]

    # check and merge the appropriate files
    mergeDir = '.merge_temp'
    python_mkdir(mergeDir)
    files = []
    for n in args.ntuples:
        if not os.path.exists(n):
            logging.error('{0} does not exist.'.format(n))
            return 0
        if os.path.isdir(n):
            logging.info('Merging files in {0}.'.format(n))
            outname = '{0}/{1}.root'.format(mergeDir,n.replace('/','_'))
            command = 'hadd -f {0} {1}/*.root'.format(outname,n)
            out = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]
            files += [outname]
        else:
            files += [n]

    # check tree in ntuples
    tfiles = []
    for f in files:
        tfiles += [rt.TFile(f)]
    trees = []
    for tf in tfiles:
        trees += [tf.Get(args.tree)]

    # make the plots
    make_plots(*trees,cuts=args.cut,labels=args.name)

    # cleanup merge directory
    logging.debug('Cleaning up')
    out = subprocess.Popen('rm -rf {0}'.format(mergeDir), shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT).communicate()[0]

    return 0

if __name__ == "__main__":
    main()
