#!/usr/bin/env python

import argparse
import itertools
import sys
import os
import logging
import json
import pickle

from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS

loglevel = 'INFO'
logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=loglevel, datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)

# control
doYields = True
doFakes = True
doSystematics = True

#####################
### Setup plotter ###
#####################
logger.info('Setup plotter')
analysis = 'WZ'
period = 13
channel = 'WZ'

doDataDriven = True

cut = 'select.passTight'
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 5 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 10 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 15 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 20 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 25 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS
#cut = 'z1.PassTight1 && z1.PassTight2 && w1.PassTight1 && finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && fabs(z1.mass-%f) < 30 && w1.mll_z1_1>4 && w1.mll_z1_2>4 && w1.Pt1>20. && finalstate.met>30.' % ZMASS

ntuples = 'ntuples/%s_%iTeV_%s' % (analysis,period,channel)
saves = '%s_%s_%iTeV' % (analysis,channel,period)
mergeDict = getMergeDict(period)
nl = 3
finalStates, leptons = getChannels(nl)
finalStates = ['eee','eem','mme','mmm']
s = 'WZ'
sigMap = getSigMap(nl)
channelBackground =  getChannelBackgrounds(period)
plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,datadriven=False,rootName='wz_yields')
plotter.initializeBackgroundSamples([sigMap[period][x] for x in channelBackground[channel+'datadriven']])
plotter.initializeDataSamples([sigMap[period]['data']])
intLumi = getIntLumiMap()[period]
plotter.setIntLumi(intLumi)
dataDrivenMC = channelBackground[channel+'datadriven'] + ['datadriven']
allMC = channelBackground[channel] + ['ZG','Z']

#####################
### Get FakeRates ###
#####################
if doFakes:
    logger.info('Load fake rates')
    fakesFile = 'InitialStateAnalysis/Analyzers/python/scale_factors/fakes_trigIso_dijet_13TeV.json'
    fakerates = {}
    with open(fakesFile,'r') as f:
        fakerates = json.load(f)


##############################
### Get Yields with errors ###
##############################
if doYields:
    logger.info('Get yields')
    labels = {0:'F',1:'P'}
    yields = {}
    err2 = {}
    yieldsMC = {}
    err2MC = {}
    yieldsWeightedData = {}
    err2WeightedData = {}
    yieldsWeightedMC = {}
    err2WeightedMC = {}
    yieldsWeighted = {}
    err2Weighted = {}
    nameMap = {
        0: 'z1.PassTight1',
        1: 'z1.PassTight2',
        2: 'w1.PassTight1',
    }

    
    allCuts = getPassTightDefinition(analysis,channel,period)
    myCut = cut
    for l in range(3):
        myCut = myCut.replace(nameMap[l],'1')
        myCut = myCut.replace('l{0}.PassTight'.format(l),'1')
    myCut = myCut.replace('select.passTight',allCuts)
    for chan in finalStates:
        logger.info('Get yields - {0}'.format(chan))
        yields[chan] = {}
        err2[chan] = {}
        yieldsWeighted[chan] = {}
        err2Weighted[chan] = {}
        yieldsWeightedMC[chan] = {}
        err2WeightedMC[chan] = {}
        yieldsWeightedData[chan] = {}
        err2WeightedData[chan] = {}
        for f in itertools.product('PF',repeat=3):
            thisName = ''.join(f)
            thisCut = '{0} && channel=="{1}" && fakeChannel=="{2}"'.format(myCut,chan,thisName)
            logger.info('Get yields - {0} - {1}'.format(chan,thisName))
            tempYields = plotter.getDataEntries(thisCut,doError=True)
            yields[chan][thisName] = tempYields[0]
            err2[chan][thisName] = tempYields[1]**2
            #num = '1' if thisName.count('F') in [1,3] else '-1'
            #for l in range(3):
            #    if thisName[l] == 'F': num += '*{0}'.format(fakeMap[l])
            scalefactor = 'event.fakerate'
            tempYields = plotter.getDataEntries(thisCut,customScale=scalefactor,doError=True)
            yieldsWeightedData[chan][thisName] = tempYields[0]
            err2WeightedData[chan][thisName] = tempYields[1]**2
            yieldsWeighted[chan][thisName] = tempYields[0]
            err2Weighted[chan][thisName] = tempYields[1]**2
            if thisName=='PPP':
                yieldsWeightedMC[chan][thisName] = 0.
                err2WeightedMC[chan][thisName] = 0.
            else:
                yieldsWeightedMC[chan][thisName] = 0.
                err2WeightedMC[chan][thisName] = 0.
                for bg in dataDrivenMC:
                    if bg in ['datadriven']: continue
                    mcscale = plotter.getScaleFactor() + '*' + scalefactor
                    tempYields = plotter.getNumEntries(thisCut,sigMap[period][bg],customScale=mcscale,doError=True)
                    yieldsWeightedMC[chan][thisName] += tempYields[0]
                    err2WeightedMC[chan][thisName] += tempYields[1]**2
                yieldsWeighted[chan][thisName] -= yieldsWeightedMC[chan][thisName]
                err2Weighted[chan][thisName] += err2WeightedMC[chan][thisName]

    # get the mc stuff
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,datadriven=False,rootName='wz_yields')
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in allMC])
    plotter.initializeDataSamples([sigMap[period]['data']])
    intLumi = getIntLumiMap()[period]
    plotter.setIntLumi(intLumi)

    for chan in finalStates:
        logger.info('Get mc yields - {0}'.format(chan))
        yieldsMC[chan] = {}
        err2MC[chan] = {}
        for f in itertools.product('PF',repeat=3):
            thisName = ''.join(f)
            thisCut = '{0} && channel=="{1}" && fakeChannel=="{2}"'.format(myCut,chan,thisName)
            logger.info('Get mc yields - {0} - {1}'.format(chan,thisName))
            # yields in each channel
            yieldsMC[chan][thisName] = {}
            err2MC[chan][thisName] = {}
            for b in allMC:
                tempYields = plotter.getNumEntries(thisCut, sigMap[period][b], doError=True)
                yieldsMC[chan][thisName][b] = tempYields[0]
                err2MC[chan][thisName][b] = tempYields[1]**2
            tempYields = plotter.getDataEntries(thisCut, doError=True)
            yieldsMC[chan][thisName]['data'] = tempYields[0]
            err2MC[chan][thisName]['data'] = tempYields[1]**2

    allYields = {
        'yields' : yields,
        'err2' : err2,
        'yieldsMC' : yieldsMC,
        'err2MC' : err2MC,
        'yieldsWeighted' : yieldsWeighted,
        'err2Weighted' : err2Weighted,
        'yieldsWeightedMC' : yieldsWeightedMC,
        'err2WeightedMC' : err2WeightedMC,
        'yieldsWeightedData' : yieldsWeightedData,
        'err2WeightedData' : err2WeightedData,
    }

    # save pickle file
    with open('yields.pkl','wb') as f:
        pickle.dump(allYields,f)
    # save json file
    with open('yields.json','w') as f:
        json.dump(allYields,f,sort_keys=True,indent=4)


#######################
### Get Systematics ###
#######################
if doSystematics:
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,datadriven=True,rootName='wz_systematics')
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in channelBackground[channel]])
    plotter.initializeDataSamples([sigMap[period]['data']])
    intLumi = getIntLumiMap()[period]
    plotter.setIntLumi(intLumi)
    
    def extractSystematics(default,up,down):
        syst = {}
        for chan in default:
            syst[chan] = {}
            syst[chan]['up'] = float(abs(up[chan]['data'][0] - default[chan]['data'][0]))/default[chan]['data'][0]
            syst[chan]['down'] = float(abs(default[chan]['data'][0] - down[chan]['data'][0]))/default[chan]['data'][0]
        return syst
    
    logger.info('Calculating systematics')
    systematics = {}
    
    # lepton uncertainty
    systematics['lepton'] = {
        'eee' : {'up': 0.019, 'down': 0.019},
        'eem' : {'up': 0.016, 'down': 0.016},
        'mme' : {'up': 0.016, 'down': 0.016},
        'mmm' : {'up': 0.016, 'down': 0.016},
    }

    systematics['electron'] = {
        'eee' : {'up': 0.018, 'down': 0.018},
        'eem' : {'up': 0.012, 'down': 0.012},
        'mme' : {'up': 0.005, 'down': 0.005},
        'mmm' : {'up': 0.000, 'down': 0.000},
    }

    systematics['muon'] = {
        'eee' : {'up': 0.000, 'down': 0.000},
        'eem' : {'up': 0.004, 'down': 0.004},
        'mme' : {'up': 0.011, 'down': 0.011},
        'mmm' : {'up': 0.016, 'down': 0.016},
    }

    # pileup uncertainty
    systematics['pileup'] = {
        'eee' : {'up': 0.0064, 'down': 0.0064},
        'eem' : {'up': 0.0088, 'down': 0.0088},
        'mme' : {'up': 0.0124, 'down': 0.0124},
        'mmm' : {'up': 0.0020, 'down': 0.0020},
    }

    # met uncertainty
    systematics['met'] = {
        'eee' : {'up': 0.019, 'down': 0.019},
        'eem' : {'up': 0.013, 'down': 0.013},
        'mme' : {'up': 0.033, 'down': 0.033},
        'mmm' : {'up': 0.017, 'down': 0.017},
    }

    
    # datadriven uncertainty
    # call it 40%
    systematics['datadriven'] = {
        'eee' : {'up': 0.4, 'down': 0.4},
        'eem' : {'up': 0.4, 'down': 0.4},
        'mme' : {'up': 0.4, 'down': 0.4},
        'mmm' : {'up': 0.4, 'down': 0.4},
    }

    # ZZ
    systematics['zz'] = {
        'eee' : {'up': 0.16, 'down': 0.16},
        'eem' : {'up': 0.16, 'down': 0.16},
        'mme' : {'up': 0.16, 'down': 0.16},
        'mmm' : {'up': 0.16, 'down': 0.16},
    }

    # WW
    systematics['ww'] = {
        'eee' : {'up': 0.06, 'down': 0.06},
        'eem' : {'up': 0.06, 'down': 0.06},
        'mme' : {'up': 0.06, 'down': 0.06},
        'mmm' : {'up': 0.06, 'down': 0.06},
    }

    # ZG
    systematics['zg'] = {
        'eee' : {'up': 0.06, 'down': 0.06},
        'eem' : {'up': 0.06, 'down': 0.06},
        'mme' : {'up': 0.06, 'down': 0.06},
        'mmm' : {'up': 0.06, 'down': 0.06},
    }

    # TTV
    systematics['ttv'] = {
        'eee' : {'up': 0.15, 'down': 0.15},
        'eem' : {'up': 0.15, 'down': 0.15},
        'mme' : {'up': 0.15, 'down': 0.15},
        'mmm' : {'up': 0.15, 'down': 0.15},
    }

    # VVV
    systematics['vvv'] = {
        'eee' : {'up': 0.06, 'down': 0.06},
        'eem' : {'up': 0.06, 'down': 0.06},
        'mme' : {'up': 0.06, 'down': 0.06},
        'mmm' : {'up': 0.06, 'down': 0.06},
    }



    # lumi uncertainty
    systematics['luminosity'] = {
        'eee' : {'up': 0.046, 'down': 0.046},
        'eem' : {'up': 0.046, 'down': 0.046},
        'mme' : {'up': 0.046, 'down': 0.046},
        'mmm' : {'up': 0.046, 'down': 0.046},
    }

    # pdf
    # TODO: need to propagate
    systematics['pdf'] = {
        'eee' : {'up': .01407, 'down': .01407}, # for now just taking gen, figure it out later after new fsa
        'eem' : {'up': .01394, 'down': .01394},
        'mme' : {'up': .01399, 'down': .01399},
        'mmm' : {'up': .01395, 'down': .01395},
    }

    # scale
    # TODO: need to propagate
    systematics['scale'] = {
        'eee' : {'up': .04296, 'down': .04194},
        'eem' : {'up': .04298, 'down': .04195},
        'mme' : {'up': .04285, 'down': .04159},
        'mmm' : {'up': .04298, 'down': .04196},
    }
    

    print json.dumps(systematics,sort_keys=True,indent=4)
    # save pickle file
    with open('systematics.pkl','wb') as f:
        pickle.dump(systematics,f)
    # save json file
    with open('systematics.json','w') as f:
        json.dump(systematics,f,sort_keys=True,indent=4)

