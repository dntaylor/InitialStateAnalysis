#!/usr/bin/env python
import logging
import sys
import argparse
import numpy as np
import itertools

from InitialStateAnalysis.Plotters.plotUtils import getSigMap, getIntLumiMap, getChannels, getMergeDict, ZMASS, getChannelBackgrounds
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Limits.WZLimits import WZLimits
from InitialStateAnalysis.Limits.limitUtils import *
from InitialStateAnalysis.Utilities.utilities import *
from multiprocessing import Pool

def limit(analysis,region,period,chan,**kwargs):
    cut = kwargs.pop('cut','1')
    name = kwargs.pop('name','card')
    scalefactor = kwargs.pop('scalefactor','event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale')
    datacardDir = kwargs.pop('datacardDir','./datacards')
    logging.info("Processing card name {0}".format(name))

    chanCut = '{0} && channel=="{1}"'.format(cut,chan)

    limits = WZLimits(analysis,region, period, chanCut, './ntuples/%s_%iTeV_%s' % (analysis, period, region),
                    '%s/%s_%itev_%s' % (datacardDir, analysis, period, region), scalefactor=scalefactor)

    # add systematic
    bgnames = ['datadriven','ZZ','WW','TTV','VVV']
    signames = ['WZ']
    mcnames = ['WZ','ZZ','WW','TTV','VVV']

    # lumi
    # current recommendation: 12%
    lumi = {}
    for b in mcnames: lumi[b] = 1.12
    limits.add_systematics("lumi", "lnN", **lumi)

    # datadriven
    # assume 40%
    fake = {'datadriven' : 1.4}
    limits.add_systematics('fake_rate_unc','lnN',**fake)

    # eff uncertainty
    # take scale factors for leptons, propagate up and down based on statistical uncertainty
    lepvals = {
        'eee' : 1.018,
        'eem' : 1.013,
        'mme' : 1.006,
        'mmm' : 1.002,
    }
    lep = {}
    for m in mcnames: lep[m] = lepvals[chan]
    limits.add_systematics('lep_eff_unc','lnN',**lep)

    # pu
    # assume 10% uncertainty on min bias cross section, scale up and down, take largest difference in wz yield
    puvals = {
        'eee' : 1.0042,
        'eem' : 1.0015,
        'mme' : 1.0039,
        'mmm' : 1.0024,
    }
    pu = {}
    for m in mcnames: pu[m] = puvals[chan]
    limits.add_systematics('PU_unc','lnN',**pu)

    # met
    # scale all components up and down independently, add in quadrature the largest
    metvals = {
        'eee' : 1.0146, # placeholder from 8 tev
        'eem' : 1.0150,
        'mme' : 1.0159,
        'mmm' : 1.0117,
    }
    met = {}
    for m in mcnames: met[m] = metvals[chan]
    limits.add_systematics('met_unc','lnN',**met)
    
    # pdf
    # propagate pdf ucnertainties through the selection, scale up and down, take largest
    pdfvals = {
        'eee' : 1.01407, # for now just taking gen, figure it out later after new fsa
        'eem' : 1.01394,
        'mme' : 1.01399,
        'mmm' : 1.01395,
    }
    pdf = {}
    for s in signames: pdf[s] = pdfvals[chan]
    limits.add_systematics('pdf_unc','lnN',**pdf)

    # scale
    # propagate scale uncertainties through the selection, scale up and down, take largest
    scalevals = {
        'eee' : 1.04296, # again, now just taking gen, fix fsa later
        'eem' : 1.04298,
        'mme' : 1.04285,
        'mmm' : 1.04298,
    }
    scale = {}
    for s in signames: scale[s] = scalevals[chan]
    limits.add_systematics('scale_unc','lnN',**scale)

    # gen card
    limits.gen_card("{0}.txt".format(name))

def limits(analysis,region,period,**kwargs):
    for chan in ['eee','eem','mme','mmm']:
        limit(analysis,region,period,chan,name=chan,**kwargs)

def parse_command_line(argv):
    parser = get_parser("Produce datacards")

    parser.add_argument('-c','--cut',type=str,default='1',help='Cut to be applied to limits.')
    parser.add_argument('-sf','--scaleFactor',type=str,default='event.gen_weight*event.pu_weight*event.lep_scale*event.trig_scale',help='Scale factor for MC.')

    args = parser.parse_args(argv)
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    loglevel = getattr(logging,args.log)
    logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')

    limits(args.analysis,args.channel,args.period,scalefactor=args.scaleFactor,cut=args.cut)

    return 0


if __name__ == "__main__":
    main()

