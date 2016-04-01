#!/usr/bin/env python
'''
AnalyzerBase.py

AnalyzerBase is a class that takes FSA ntuples and produces ISA ntuples.
It defines composite objects from final state objects (leptons, photons, jets, met).
Combinatorics are taken care of via user-defined minimization functions.
The output is an ntuple with useful variables for plotting and applying additional
cuts for final selections.

The user must import this class and define the preselection cuts, combinatoric
variable to be minimized, and initial state objects.

Variables and functions that must be defined in inheritor class:
    self.channel            A name for this analysis
    self.final_states       List of FSA final state strings (i.e. eemm, mj, gg, etc.)
    self.initial_states     List of ISA initial state strings (in order of objects returned in choose_objects)
                            Must match keys in object definitions.
    self.object_defintions  A dictionary of initial state objects (i.e. 'z': ['em','em'], 'w': ['em','n'], 'h': ['g','g'])
                            Allowable entries are e,m,t,j,g,n for eletron, muon, tau, jet, photon, met respectively
    choose_objects(rtrow)   A function that produces object candidates and a list minimizing variables
                            ex: return ([massdiff, -ptsum], list(leptons))
    preselection(rtrow)     A function to return a CutSequence object with ordered cuts to be applied

Author: Devin N. Taylor, UW-Madison
'''
import os
import sys
from itertools import permutations, combinations
import argparse
import datetime
import math
import logging

from scale_factors import LeptonScaleFactors, TriggerScaleFactors, ChargeIdSystematics, LeptonEfficiency, LeptonFakeRate
from pu_weights import PileupWeights
import leptonId as lepId
from ntuples import *

sys.argv.append('-b')
import ROOT as rt
sys.argv.pop()


ZMASS = 91.1876

class CutSequence(object):
    '''
    A class for defining cut orders for preselection.
    '''
    def __init__(self):
        self.cut_sequence = []
        self.results = 0

    def add(self, fun, label=''):
        self.cut_sequence.append([fun,label])

    def evaluate(self, rtrow):
        for i,cut in enumerate(self.cut_sequence):
            if not cut[0](rtrow):
                self.results = i
                return False
        self.results = i+1
        return True

    def getResults(self):
        return self.results

class CutTree(object):
    '''
    Tree for storing passing selection information.
    '''
    def __init__(self):
        # create selections
        self.labels = []
        self.selections = {}
        self.results = {}

    def add(self, fun, label):
        self.labels += [label]
        self.selections[label] = fun

    def getLabels(self):
        return self.labels

    def evaluate(self,rtrow,failed=False):
        self.results['evt'] = long(rtrow.evt)
        self.results['run'] = int(rtrow.run)
        self.results['lumi'] = int(rtrow.lumi)
        if failed:
            for label in self.selections:
                self.results[label] = False
            return False
        passAll = True
        for label in self.selections:
            cut = self.selections[label]
            self.results[label] = cut(rtrow)
            if not self.results[label]: passAll = False
        return passAll

    def getResults(self):
        return self.results

def lep_order(a, b):
    '''
    A simple function to guarantee order of leptons in FSA ntuples.
    '''
    #if len(a)==2 and len(b)==2:
    #    a_index = int(a[1])
    #    b_index = int(b[1])
    #    return a_index > b_index or a[0] > b[0]
    #return a[0] > b[0]
    try:
        a_index = int(a[1])
    except IndexError:
        a_index = 1
    try:
        b_index = int(b[1])
    except IndexError:
        b_index = 1
    
    if a[0] == b[0]:
        return a_index > b_index
    else:
        return a[0] > b[0]

def ordered(a,b):
    '''
    Return a,b in lep order.
    '''
    return [a,b] if lep_order(b,a) else [b,a]

def deltaPhi(phi0,phi1):
    result = phi0-phi1
    while result>rt.TMath.Pi():
        result -= 2*rt.TMath.Pi()
    while result<=-rt.TMath.Pi():
        result += 2*rt.TMath.Pi()
    return result

def deltaR(eta0,phi0,eta1,phi1):
    deta = eta0-eta1
    dphi = deltaPhi(phi0,phi1)
    return rt.TMath.Sqrt(deta**2+dphi**2)

class AnalyzerBase(object):
    '''
    The basic analyzer class. Inheritor classes must define
        TODO
    '''
    def __init__(self, sample_name, file_list, out_file, period, **kwargs):
        loglevel = kwargs.pop('loglevel','INFO')
        self.loglevel = getattr(logging,loglevel)
        logging.basicConfig(format='%(asctime)s.%(msecs)03d %(levelname)s %(name)s: %(message)s', level=self.loglevel, datefmt='%Y-%m-%d %H:%M:%S')
        self.sample_name = sample_name
        self.isData = ('data' in self.sample_name)
        if isinstance(file_list, basestring): # the list is a file
            self.file_names = []
            if file_list[-4:]=='root': # passed a string list of files comma seperated
                for f in file_list.split(','):
                    self.file_names += [f]
            else:
                with open(file_list,'r') as files:
                    for f in files:
                        self.file_names += [f]
        else:                                 # the list is a python list
            self.file_names = file_list
        self.out_file = out_file
        self.period = period
        self.metShift = kwargs.pop('metShift','')
        if self.metShift:
            if self.metShift not in ['ees+','ees-','mes+','mes-','tes+','tes-','ues+','ues-','jes+','jes-','jres+','jres-']:
                logging.error('{0} is not an allowed met shift'.format(self.metShift))
        if self.isData: self.metShift = '' # force no shift for data

    def __enter__(self):
        self.begin()
        return self

    def __exit__(self, type, value, traceback):
        self.finish()

    def begin(self):
        self.lepscaler = LeptonScaleFactors()
        self.trigscaler = TriggerScaleFactors()
        self.pu_weights = PileupWeights()
        self.chargeid = ChargeIdSystematics()
        self.lepeff = LeptonEfficiency()
        self.lepfake = LeptonFakeRate()

        self.file = rt.TFile(self.out_file, 'recreate')
        
        if hasattr(self,'other_states'):
            states = [self.initial_states] + self.other_states
        else:
            states = [self.initial_states]
        if not hasattr(self,'alternateIds'): self.alternateIds = []
        if not hasattr(self,'doVBF'): self.doVBF = False
        if not hasattr(self,'doMetUnc'): self.doMetUnc = False
        self.ntuple, self.branches = buildNtuple(self.object_definitions,states,self.channel,self.final_states,altIds=self.alternateIds,doVBF=self.doVBF,doMetUnc=self.doMetUnc)

        if hasattr(self,'cutTreeSelections'):
            self.cutTreeLabels = self.cutTreeSelections().getLabels()
            self.cutTree, self.eventBranch, self.cutsBranch = buildCutTree(self.cutTreeLabels)

    def analyze(self,**kwargs):
        '''
        The primary analyzer loop.
        '''
        logger = logging.getLogger(__name__)

        logger.info('%s %s Analyzing' %(self.channel, self.sample_name))
        eventMap = {}
        bestCandMap = {}
        cutflowMap = {}
        cutTreeMap = {}
        cutTreeEventsToWrite = set()
        cutTreeEventsWritten = set()
        eventsToWrite = set()
        eventsWritten = set()
        passedPreselection = set()
        numEvts = 0
        absEvts = 0
        totalWritten = 0

        # iterate over files
        for i, file_name in enumerate(self.file_names):
            self.file_name = file_name
            logger.info('%s %s Processing %i/%i files', self.channel, self.sample_name, i+1, len(self.file_names))
            sys.stdout.flush()
            if file_name.startswith('//store'): file_name = 'root://cmsxrootd.hep.wisc.edu/%s' % file_name
            if file_name.startswith('/store'): file_name = 'root://cmsxrootd.hep.wisc.edu//%s' % file_name
            if file_name.startswith('store'): file_name = 'root://cmsxrootd.hep.wisc.edu///%s' % file_name

            rtFile = rt.TFile.Open(file_name, "READ")

            # iterate over final states
            for fs in self.final_states:
                if len(self.file_names)<10: logger.info('%s %s %s' % (self.channel, self.sample_name, fs))
                tree = rtFile.Get("%s/final/Ntuple" % fs)
                #if self.period==8:
                metatree = rtFile.Get("%s/metaInfo" % fs)
                tempEvts = 0
                absTempEvts = 0
                for entry in xrange(metatree.GetEntries()):
                    metatree.GetEntry(entry)
                    absTempEvts += metatree.nevents
                    tempEvts += metatree.nevents if self.isData or self.period==8 else metatree.summedWeights # gen level processed

                self.objects = self.enumerate_objects(fs)

                # initialize event counter
                numFSEvents = 0
                totalFSEvents = tree.GetEntries()

                # iterate over each row of an fsa ntuple
                rtrow = tree
                numRows = tree.GetEntries('1')
                for r in range(numRows):
                    rtrow.GetEntry(r)
                    if numFSEvents % 10000 == 0:
                        if len(self.file_names)==1: logger.info('%s %s %s %i/%i entries' % (self.channel, self.sample_name, fs, numFSEvents, totalFSEvents))
                        sys.stdout.flush()
                    numFSEvents += 1

                    # cache to prevent excessive reads of fsa ntuple
                    self.cache = {}

                    # event number for dictionary storing
                    eventkey = (long(rtrow.evt), int(rtrow.lumi), int(rtrow.run))

                    # if we have a cutTree, do it
                    if hasattr(self,'cutTree'):
                        eventCutTree = self.cutTreeSelections()

                    # can we define the object we want?
                    candidate = self.choose_objects(rtrow)
                    if not candidate or len(candidate) < 2 or len(candidate[0])<1: # in case no objects satisfy our conditions
                        # do the cut tree stuff, this fails topology
                        if hasattr(self,'cutTree'):
                            if eventkey not in bestCandMap: # it has never passed before, store it
                                cutTreeMap[eventkey] = self.storeCutTree(rtrow,eventCutTree,failed=True)
                                cutTreeEventsToWrite.add(eventkey)
                        continue


                    # name candidate
                    self.objCand = candidate[1]

                    # now see if event is viable
                    passPreselection = self.pass_preselection(rtrow)

                    # check preselection
                    if not passPreselection:
                        if hasattr(self,'cutTree'):
                            # check combinatorics
                            if eventkey in bestCandMap:
                                bestcand = bestCandMap[eventkey]
                            else:
                                numMin = len(candidate[0])
                                bestcand = [float('inf')] * numMin
                            if self.good_to_store(rtrow,candidate[0],bestcand):
                                #bestCandMap[eventkey] = candidate[0]
                                if eventkey not in passedPreselection: # dont replace something that could pass preselection
                                    cutTreeMap[eventkey] = self.storeCutTree(rtrow,eventCutTree,failed=False)
                                    cutTreeEventsToWrite.add(eventkey)
                        continue # dont store in event tree, just cut tree

                    passedPreselection.add(eventkey)

                    # check combinatorics
                    if eventkey in bestCandMap: 
                        bestcand = bestCandMap[eventkey]
                    else:
                        numMin = len(candidate[0])
                        bestcand = [float('inf')] * numMin
                    if self.good_to_store(rtrow,candidate[0],bestcand):
                        bestCandMap[eventkey] = candidate[0]
                        ntupleRow = self.store_row(rtrow, *self.objCand)
                        eventMap[eventkey] = ntupleRow
                        eventsToWrite.add(eventkey)
                        if hasattr(self,'cutTree'):
                            cutTreeMap[eventkey] = self.storeCutTree(rtrow,eventCutTree,failed=False)
                            cutTreeEventsToWrite.add(eventkey)

            rtFile.Close("R")
            numEvts += tempEvts
            absEvts += absTempEvts

            # end of file, write the ntuples
            self.file.cd()
            for key in eventsToWrite:
                if key in eventsWritten:
                    logger.warning('%s %s Attempted to write previously written event' % (self.channel, self.sample_name))
                else:
                    self.write_row(eventMap[key])
                    self.ntuple.Fill()
            eventsWritten.update(eventsToWrite)
            eventMap = {}
            eventsToWrite = set()
            if hasattr(self,'cutTree'):
                for key in cutTreeEventsToWrite:
                    if key in cutTreeEventsWritten:
                        logger.warning('%s %s Attempted to write previously written event - cut tree' % (self.channel, self.sample_name))
                    else:
                        self.writeCutTree(cutTreeMap[key])
                        self.cutTree.Fill()
                cutTreeEventsWritten.update(cutTreeEventsToWrite)
                cutTreeMap = {}
                cutTreeEventsToWrite = set()

        logger.info('%s %s Filled Tree (%i events)' % (self.channel, self.sample_name, len(eventsWritten)))
        logger.info('%s %s Filled Cut Tree (%i events)' % (self.channel, self.sample_name, len(cutTreeEventsWritten)))

        # now we store the total processed events
        logger.info('%s %s Processed %i events' % (self.channel, self.sample_name, absEvts))
        logger.info('%s %s Weighted processed %f events' % (self.channel, self.sample_name, numEvts))

        ## and the cutflow
        cutflowVals = []
        #for val in cutflowMap.itervalues():
        #    for i in range(val+1):
        #        if len(cutflowVals)<i+1: cutflowVals.append(1)
        #        else: cutflowVals[i] += 1
        #print "%s %s %s: Cutflow: " % (str(datetime.datetime.now()), self.channel, self.sample_name), cutflowVals

        cutflowHist = rt.TH1F('cutflow','cutflow',len(cutflowVals)+1,0,len(cutflowVals)+1)
        cutflowHist.SetBinContent(1,numEvts)
        #for i in range(len(cutflowVals)):
        #    cutflowHist.SetBinContent(i+2,cutflowVals[i])
        ## rename cutflow bins if self.cutflow_labels defined
        #if hasattr(self,'cutflow_labels'):
        #    pass # TODO
        cutflowHist.Write()

    def finish(self):
        self.lepscaler.close()
        self.file.Write()
        self.file.Close()

    @staticmethod
    def enumerate_objects(final_state):
        '''Get the objects available in the ntuple for a given final state'''
        out = []
        for i  in ['e', 'm', 't', 'j', 'g']:
            N = final_state.count(i)
            if N==1:
               out += [i]
            else:
               out += ['%s%i' % (i, n) for n in xrange(1, N+1)]
        return out

    def choose_alternative_objects(self, rtorw, state):
        '''
        Dummy method for alternative object selection.
        '''
        return []

    #@staticmethod
    def good_to_store(self, rtrow, cand1, cand2):
        '''
        Iterate through minimizing variables.
        '''
        for min1, min2 in zip(cand1, cand2):
            if min1 < min2: return True
            if min1 > min2: return False
        return False

    def storeCutTree(self,rtrow,eventCutTree,failed=False):
        passAll = eventCutTree.evaluate(rtrow,failed=failed)
        results = eventCutTree.getResults()
        ntupleRow = {}
        for key in results:
            if key in ['evt','run','lumi']: # store in event branch
                ntupleRow['event.{0}'.format(key)] = results[key]
            else: # store in selections branch
                ntupleRow['selections.{0}'.format(key)] = results[key]
        return ntupleRow

    def writeCutTree(self, nrow):
        '''
        Function to write the ntuple row to the tree.
        '''
        for key,val in nrow.iteritems():
            branch, var = key.split('.')
            if branch=='event':
                setattr(self.eventBranch,var,val)
            if branch=='selections':
                setattr(self.cutsBranch,var,val)

    def store_row(self,rtrow,*objects):
        '''
        Function to return a dictionary of event values to be written to the ntuple.
        '''
        ntupleRow = {}

        ntupleRow["select.passTight"] = int(self.pass_selection(rtrow))
        ntupleRow["select.passLoose"] = int(self.pass_preselection(rtrow))
        ntupleRow["select.passDoubleMuon"] = int(rtrow.doubleMuPass if self.period==13 else rtrow.doubleMuPass or rtrow.doubleMuTrkPass)
        ntupleRow["select.passDoubleEG"] = int(rtrow.doubleEPass if self.period==13 else rtrow.doubleETightPass)
        ntupleRow["select.passMuonEG"] = int(rtrow.singleMuSingleEPass if self.period==13 else rtrow.mu17ele8isoPass)
        ntupleRow["select.passEGMuon"] = int(rtrow.singleESingleMuPass if self.period==13 else rtrow.mu8ele17isoPass)
        ntupleRow["select.passSingleMuon"] = int(rtrow.singleIsoMu20Pass or rtrow.singleIsoTkMu20Pass or rtrow.singleIsoMu27Pass) if self.period==13 else -1
        ntupleRow["select.passSingleEG"] = int(rtrow.singleE23WPLoosePass) if self.period==13 else -1
        numObjs = len(self.final_states[0])
        finalStateObjects ='emtjgn'
        allowedObjects = ''
        for fsObj in finalStateObjects:
            for fs in self.final_states:
                if fsObj in fs:
                    allowedObjects += fsObj
                    break
        numObjTypes = len(allowedObjects)
        for altId in self.alternateIds:
            ntupleRow["select.pass_%s"%altId] = int(self.ID(rtrow,*objects,**self.alternateIdMap[altId]))

        scales = self.getScales(rtrow,objects)
        
        ntupleRow["event.evt"] = long(rtrow.evt)
        ntupleRow["event.lumi"] = int(rtrow.lumi)
        ntupleRow["event.run"] = int(rtrow.run)
        ntupleRow["event.nvtx"] = int(rtrow.nvtx)
        ntupleRow["event.GenNUP"] = -1 if self.isData else int(rtrow.NUP)
        ntupleRow["event.trig_prescale"] = float(scales['trigger_prescale'])
        ntupleRow["event.lep_scale"] = float(scales['lep'])
        ntupleRow["event.lep_scale_e"] = float(scales['lepe'])
        ntupleRow["event.lep_scale_m"] = float(scales['lepm'])
        ntupleRow["event.lep_scale_up"] = float(scales['lepup'])
        ntupleRow["event.lep_scale_e_up"] = float(scales['lepeup'])
        ntupleRow["event.lep_scale_m_up"] = float(scales['lepmup'])
        ntupleRow["event.lep_scale_down"] = float(scales['lepdown'])
        ntupleRow["event.lep_scale_e_down"] = float(scales['lepedown'])
        ntupleRow["event.lep_scale_m_down"] = float(scales['lepmdown'])
        ntupleRow["event.trig_scale"] = float(scales['trig'])
        ntupleRow["event.trig_scale_up"] = float(scales['trigup'])
        ntupleRow["event.trig_scale_down"] = float(scales['trigdown'])
        ntupleRow["event.pu_weight"] = float(scales['puweight'])
        ntupleRow["event.pu_weight_up"] = float(scales['puweightup'])
        ntupleRow["event.pu_weight_down"] = float(scales['puweightdown'])
        ntupleRow["event.fakerate"] = float(scales['lepfake'])
        ntupleRow["event.fakerate_up"] = float(scales['lepfakeup'])
        ntupleRow["event.fakerate_down"] = float(scales['lepfakedown'])
        ntupleRow["event.gen_weight"] = float(scales['genweight'])
        ntupleRow["event.charge_uncertainty"] = float(scales['chargeid'])


        channelString = ''
        for x in objects: channelString += x[0]
        ntupleRow["channel.channel"] = channelString
        ntupleRow["genChannel.channel"] = self.getGenChannel(rtrow)
        ntupleRow["fakeChannel.channel"] = self.getFakeChannel(rtrow)
        if self.channel=='WZ': ntupleRow["fakeChannel_tightW.channel"] = self.getFakeChannel(rtrow,tightW=True)
        if self.channel=='WZ': ntupleRow["fakeChannel_allMedium.channel"] = self.getFakeChannel(rtrow,allMedium=True)

        ntupleRow["finalstate.mass"] = self.getObject(rtrow,'mass',*objects)
        ntupleRow["finalstate.eta"] = self.getObject(rtrow,'eta',*objects)
        ntupleRow["finalstate.phi"] = self.getObject(rtrow,'phi',*objects)
        objMet = []
        for i in objects:
            objMet += [i]
        objMet += ['met']
        ntupleRow["finalstate.mT"] = self.getObject(rtrow,'mt',*objMet)
        ntupleRow["finalstate.sT"] = self.getObject(rtrow,'st',*objects)
        ntupleRow["finalstate.hT"] = float(rtrow.Ht) if hasattr(rtrow,'Ht') else float(-1)
        ntupleRow["finalstate.met"] = self.getObject(rtrow, 'pt', 'met')
        ntupleRow["finalstate.metPhi"] = self.getObject(rtrow, 'phi', 'met')
        ntupleRow["finalstate.leadJetPt"] = float(rtrow.jet1Pt) if self.period==13 else float(-1)
        ntupleRow["finalstate.leadJetEta"] = float(rtrow.jet1Eta) if self.period==13 else float(-10)
        ntupleRow["finalstate.leadJetPhi"] = float(rtrow.jet1Phi) if self.period==13 else float(-10)
        ntupleRow["finalstate.leadJetPUMVA"] = float(rtrow.jet1PUMVA) if self.period==13 else float(-1)
        ntupleRow["finalstate.jetVeto20"] = int(rtrow.jetVeto20)
        ntupleRow["finalstate.jetVeto30"] = int(rtrow.jetVeto30)
        ntupleRow["finalstate.jetVeto40"] = int(rtrow.jetVeto40)
        ntupleRow["finalstate.bjetVeto20Loose"] = int(rtrow.bjetCISVVeto20Loose) if self.period==13 else -1
        ntupleRow["finalstate.bjetVeto30Loose"] = int(rtrow.bjetCISVVeto30Loose) if self.period==13 else -1
        ntupleRow["finalstate.bjetVeto20Medium"] = int(rtrow.bjetCISVVeto20Medium) if self.period==13 else int(rtrow.bjetCSVVeto)
        ntupleRow["finalstate.bjetVeto30Medium"] = int(rtrow.bjetCISVVeto30Medium) if self.period==13 else int(rtrow.bjetCSVVeto30)
        ntupleRow["finalstate.bjetVeto20Tight"] = int(rtrow.bjetCISVVeto20Tight) if self.period==13 else -1
        ntupleRow["finalstate.bjetVeto30Tight"] = int(rtrow.bjetCISVVeto30Tight) if self.period==13 else -1

        ntupleRow["finalstate.muonVetoTight"] = self.getObjectVetos(rtrow,'m','Tight') if self.period==13 else int(rtrow.muonVetoWZTight)
        ntupleRow["finalstate.elecVetoTight"] = self.getObjectVetos(rtrow,'e','Tight') if self.period==13 else int(rtrow.elecVetoWZTight)
        ntupleRow["finalstate.muonVetoLoose"] = self.getObjectVetos(rtrow,'m','Loose') if self.period==13 else int(rtrow.muVetoPt5IsoIdVtx)
        ntupleRow["finalstate.elecVetoLoose"] = self.getObjectVetos(rtrow,'e','Loose') if self.period==13 else int(rtrow.eVetoMVAIsoVtx)
        if self.doVBF:
            ntupleRow["finalstate.vbfMass"] = float(rtrow.vbfMass)
            ntupleRow["finalstate.vbfPt"] = float(rtrow.vbfdijetpt)
            ntupleRow["finalstate.vbfPt1"] = float(rtrow.vbfj1pt)
            ntupleRow["finalstate.vbfPt2"] = float(rtrow.vbfj2pt)
            ntupleRow["finalstate.vbfEta1"] = float(rtrow.vbfj1eta)
            ntupleRow["finalstate.vbfEta2"] = float(rtrow.vbfj2eta)
            ntupleRow["finalstate.centralJetVeto20"] = float(rtrow.vbfJetVeto20)
            ntupleRow["finalstate.centralJetVeto30"] = float(rtrow.vbfJetVeto30)

        def store_state(rtrow,ntupleRow,state,theObjects,period):
            objStart = 0
            for i in state:
                numObjects = len([ x for x in self.object_definitions[i] if x != 'n']) if theObjects else 0
                finalObjects = theObjects[objStart:objStart+numObjects]
                orderedFinalObjects = sorted(finalObjects, key = lambda x: getattr(rtrow,"%sPt" % x), reverse=True)
                if len(self.object_definitions[i]) == 1:
                    ntupleRow["%s.mass" %i] = float(-999)
                    ntupleRow["%s.Pt" %i] = float(-9)
                    ntupleRow["%s.Eta" %i] = float(-9)
                    ntupleRow["%s.Phi" %i] = float(-9)
                    ntupleRow["%s.sT" %i] = self.getObject(rtrow, "pt", finalObjects[0]) if theObjects else float(-9)
                    ntupleRow["%s.dPhi" %i] = float(-9)
                    ntupleRow["%s.dR" %i] = float(-9)
                    ntupleRow["%sFlv.Flv" %i] = finalObjects[0][0] if theObjects else 'a'
                elif 'n' == self.object_definitions[i][1]:
                    ntupleRow["%s.mass" %i] = self.getObject(rtrow, "mt", finalObjects[0], 'met') if theObjects else float(-999)
                    ntupleRow["%s.Pt" %i] = self.getObject(rtrow, "pt", finalObjects[0], 'met') if theObjects else float(-9)
                    ntupleRow["%s.sT" %i] = self.getObject(rtrow, "st", finalObjects[0], 'met') if theObjects else float(-9)
                    ntupleRow["%s.dPhi" %i] = self.getObject(rtrow, "dphi", finalObjects[0], 'met') if theObjects else float(-9)
                    ntupleRow["%sFlv.Flv" %i] = finalObjects[0][0] if theObjects else 'a'
                else:
                    finalObjOrdered = ordered(finalObjects[0], finalObjects[1]) if theObjects else []
                    ntupleRow["%s.mass" %i] = self.getObject(rtrow, "mass", finalObjOrdered[0], finalObjOrdered[1]) if theObjects else float(-999)
                    ntupleRow["%s.Pt" %i] = self.getObject(rtrow, "pt", finalObjOrdered[0], finalObjOrdered[1]) if theObjects else float(-9)
                    ntupleRow["%s.sT" %i] = self.getObject(rtrow, "st", finalObjOrdered[0], finalObjOrdered[1]) if theObjects else float(-9)
                    ntupleRow["%s.dPhi" %i] = self.getObject(rtrow, "dphi", finalObjOrdered[0], finalObjOrdered[1]) if theObjects else float(-9)
                    ntupleRow["%s.dR" %i] = self.getObject(rtrow, "dr", finalObjOrdered[0], finalObjOrdered[1]) if theObjects else float(-9)
                    ntupleRow["%sFlv.Flv" %i] = finalObjects[0][0] + finalObjects[1][0] if theObjects else 'aa'
                objCount = 0
                for obj in self.object_definitions[i]:
                    if obj=='n':
                        ntupleRow["%s.met" %i] = self.getObject(rtrow,'pt','met') if theObjects else float(-9)
                        ntupleRow["%s.metPhi" %i] = self.getObject(rtrow, 'phi', 'met') if theObjects else float(-9)

                    else:
                        objCount += 1
                        l = orderedFinalObjects[objCount-1] if theObjects else 'a'
                        ntupleRow["%s.Pt%i" % (i,objCount)] = self.getObject(rtrow, "pt", l) if theObjects else float(-9)
                        ntupleRow["%s.Eta%i" % (i,objCount)] = self.getObject(rtrow, "eta", l) if theObjects else float(-9)
                        ntupleRow["%s.Phi%i" % (i,objCount)] = self.getObject(rtrow, "phi", l) if theObjects else float(-9)
                        # TODO recalculate iso with shift
                        if theObjects:
                            if l[0]=='e': isoVar = 'RelPFIsoRho'
                            if l[0]=='m': isoVar = 'RelPFIsoDBDefault'
                            isoVal = float(getattr(rtrow, "%s%s" % (l, isoVar))) if l[0] in 'em' and theObjects else float(-9.)
                        else:
                            isoVal = float(-9)
                        ntupleRow["%s.Iso%i" % (i,objCount)] = isoVal
                        ntupleRow["%s.Dxy%i" % (i,objCount)] = float(getattr(rtrow, "%sPVDXY" % l)) if theObjects else float(-9)
                        ntupleRow["%s.Dz%i" % (i,objCount)] = float(getattr(rtrow, "%sPVDZ" % l)) if theObjects else float(-9)
                        ntupleRow["%s.SigmaIEtaIEta%i" % (i,objCount)] = float(getattr(rtrow, "%sSigmaIEtaIEta" % (l))) if l[0] in 'e' else float(-1.)
                        ntupleRow["%s.DEtaIn%i" % (i,objCount)] = float(getattr(rtrow, "%sdeltaEtaSuperClusterTrackAtVtx" % (l))) if l[0] in 'e' else float(-9.)
                        ntupleRow["%s.DPhiIn%i" % (i,objCount)] = float(getattr(rtrow, "%sdeltaPhiSuperClusterTrackAtVtx" % (l))) if l[0] in 'e' else float(-9.)
                        ntupleRow["%s.HOverE%i" % (i,objCount)] = float(getattr(rtrow, "%sHadronicOverEM" % (l))) if l[0] in 'e' else float(-1.)
                        ntupleRow["%s.OoEmOoP%i" % (i,objCount)] = float(abs((1.-getattr(rtrow, "%seSuperClusterOverP" % (l)))*1./getattr(rtrow, "%secalEnergy" % (l)))) if l[0] in 'e' else float(-1.)
                        ntupleRow["%s.TriggeringMVA%i" % (i,objCount)] = float(-9. if self.period==13 else getattr(rtrow, "%sMVATrig" % (l))) if l[0] in 'e' else float(-9.)
                        ntupleRow["%s.NonTriggeringMVA%i" % (i,objCount)] = float(getattr(rtrow, "%sMVANonTrigID" % (l)) if period==13 else getattr(rtrow, "%sMVANonTrig" % (l))) if l[0] in 'e' else float(-9.)
                        ntupleRow["%s.NormalizedChi2%i" % (i,objCount)] = float(getattr(rtrow, "%sNormTrkChi2" % (l))) if l[0] in 'm' else float(-1.)
                        ntupleRow["%s.JetPt%i" % (i,objCount)] = float(getattr(rtrow, "%sJetPt" % l)) if (theObjects and l[0] in 'emt') else float(-9.)
                        ntupleRow["%s.JetBTag%i" % (i,objCount)] = float(-9.)
                        if theObjects and l[0] in 'emt':
                            ntupleRow["%s.JetBTag%i" % (i,objCount)] = float(getattr(rtrow, "%sJetCSVBtag" % l)) if period==8 else float(getattr(rtrow, "%sJetPFCISVBtag" % l))
                        if theObjects:
                            looseScales = self.lepscaler.scale_factor(rtrow, l, lepType='Loose', period=period)
                            mediumScales = self.lepscaler.scale_factor(rtrow, l, lepType='Medium', period=period)
                            tightScales = self.lepscaler.scale_factor(rtrow, l, lepType='Tight', period=period)
                            mediumLepeff = self.lepeff.scale_factor(rtrow, l, denom='Loose', numer='Medium', period=period)
                            tightLepeff = self.lepeff.scale_factor(rtrow, l, denom='Loose', numer='Tight', period=period)
                            mediumFake = self.lepfake.scale_factor(rtrow, l, denom='Loose', numer='Medium', period=period)
                            tightFake = self.lepfake.scale_factor(rtrow, l, denom='Loose', numer='Tight', period=period)
                            mediumFakeMC = self.lepfake.scale_factor(rtrow, l, denom='Loose', numer='Medium', period=period,mc=True)
                            tightFakeMC = self.lepfake.scale_factor(rtrow, l, denom='Loose', numer='Tight', period=period,mc=True)
                        else:
                            looseScales = [-1,-1,-1]
                            mediumScales = [-1,-1,-1]
                            tightScales = [-1,-1,-1]
                            looseLepeff = [-1,-1,-1]
                            mediumLepeff = [-1,-1,-1]
                            tightLepeff = [-1,-1,-1]
                            mediumFake = [-1,-1,-1]
                            tightFake = [-1,-1,-1]
                            mediumFakeMC = [-1,-1,-1]
                            tightFakeMC = [-1,-1,-1]
                        ntupleRow["%s.LepScaleLoose%i" % (i,objCount)] = float(looseScales[0])
                        ntupleRow["%s.LepScaleMedium%i" % (i,objCount)] = float(mediumScales[0])
                        ntupleRow["%s.LepScaleTight%i" % (i,objCount)] = float(tightScales[0])
                        ntupleRow["%s.LepScaleLoose%i_up" % (i,objCount)] = float(looseScales[1])
                        ntupleRow["%s.LepScaleMedium%i_up" % (i,objCount)] = float(mediumScales[1])
                        ntupleRow["%s.LepScaleTight%i_up" % (i,objCount)] = float(tightScales[1])
                        ntupleRow["%s.LepScaleLoose%i_down" % (i,objCount)] = float(looseScales[2])
                        ntupleRow["%s.LepScaleMedium%i_down" % (i,objCount)] = float(mediumScales[2])
                        ntupleRow["%s.LepScaleTight%i_down" % (i,objCount)] = float(tightScales[2])
                        ntupleRow["%s.LepEffMedium%i" % (i,objCount)] = float(mediumLepeff[0])
                        ntupleRow["%s.LepEffMedium_up%s" % (i,objCount)] = float(mediumLepeff[1])
                        ntupleRow["%s.LepEffMedium_down%s" % (i,objCount)] = float(mediumLepeff[2])
                        ntupleRow["%s.LepEffTight%i" % (i,objCount)] = float(tightLepeff[0])
                        ntupleRow["%s.LepEffTight_up%s" % (i,objCount)] = float(tightLepeff[1])
                        ntupleRow["%s.LepEffTight_down%s" % (i,objCount)] = float(tightLepeff[2])
                        ntupleRow["%s.LepFakeMedium%i" % (i,objCount)] = float(mediumFake[0])
                        ntupleRow["%s.LepFakeMedium_up%i" % (i,objCount)] = float(mediumFake[1])
                        ntupleRow["%s.LepFakeMedium_down%i" % (i,objCount)] = float(mediumFake[2])
                        ntupleRow["%s.LepFakeTight%i" % (i,objCount)] = float(tightFake[0])
                        ntupleRow["%s.LepFakeTight_up%i" % (i,objCount)] = float(tightFake[1])
                        ntupleRow["%s.LepFakeTight_down%i" % (i,objCount)] = float(tightFake[2])
                        ntupleRow["%s.LepFakeMCMedium%i" % (i,objCount)] = float(mediumFakeMC[0])
                        ntupleRow["%s.LepFakeMCMedium_up%i" % (i,objCount)] = float(mediumFakeMC[1])
                        ntupleRow["%s.LepFakeMCMedium_down%i" % (i,objCount)] = float(mediumFakeMC[2])
                        ntupleRow["%s.LepFakeMCTight%i" % (i,objCount)] = float(tightFakeMC[0])
                        ntupleRow["%s.LepFakeMCTight_up%i" % (i,objCount)] = float(tightFakeMC[1])
                        ntupleRow["%s.LepFakeMCTight_down%i" % (i,objCount)] = float(tightFakeMC[2])
                        ntupleRow["%s.Chg%i" % (i,objCount)] = float(getattr(rtrow, "%sCharge" % l)) if theObjects else float(-9)
                        ntupleRow["%s.PassLoose%i" % (i,objCount)] = float(self.ID(rtrow,l,**self.getIdArgs('Loose'))) if theObjects else float(-9)
                        ntupleRow["%s.PassMedium%i" % (i,objCount)] = float(self.ID(rtrow,l,**self.getIdArgs('Medium'))) if theObjects else float(-9)
                        ntupleRow["%s.PassTight%i" % (i,objCount)] = float(self.ID(rtrow,l,**self.getIdArgs('Tight'))) if theObjects else float(-9)
                        ntupleRow["%s.GenIsPrompt%i" % (i,objCount)] = -2000
                        ntupleRow["%s.GenPdgId%i" % (i,objCount)] = -2000
                        ntupleRow["%s.GenPatPdgId%i" % (i,objCount)] = -2000
                        ntupleRow["%s.MotherGenPdgId%i" % (i,objCount)] = -2000
                        if not self.isData and theObjects:
                            ntupleRow["%s.GenIsPrompt%i" % (i,objCount)] = float(getattr(rtrow, "%sGenPrompt" % l))
                            ntupleRow["%s.GenPdgId%i" % (i,objCount)] = float(getattr(rtrow, "%sGenPdgId" % l))
                            ntupleRow["%s.GenPatPdgId%i" % (i,objCount)] = float(getattr(rtrow, "%sGenParticle" % l))
                            ntupleRow["%s.MotherGenPdgId%i" % (i,objCount)] = float(getattr(rtrow, "%sGenMotherPdgId" % l))
                        ntupleRow["%s.ChargeConsistent%i" % (i,objCount)] = -1
                        if theObjects and period==8:
                            if l[0]=='e':
                                ntupleRow["%s.ChargeConsistent%i" % (i,objCount)] = int(getattr(rtrow,'%sChargeIdTight' %l))
                        ntupleRow["%s.ExpectedMissingInnerHits%i" % (i,objCount)] = int(getattr(rtrow,'%sMissingHits' %l)) if l[0]=='e' else int(-1)
                        ntupleRow["%s.PassConversionVeto%i" % (i,objCount)] = int(getattr(rtrow,'%sPassesConversionVeto' %l) if period==13 else not getattr(rtrow, "%sHasConversion" %l)) if l[0]=='e' else int(-1)
                        ntupleRow["%s.IsGlobalMuon%i" % (i,objCount)] = int(getattr(rtrow,'%sIsGlobal' %l)) if l[0]=='m' else int(-1)
                        ntupleRow["%s.IsPFMuon%i" % (i,objCount)] = int(getattr(rtrow,'%sIsPFMuon' %l)) if l[0]=='m' and period==13 else int(-1)
                        ntupleRow["%s.IsTrackerMuon%i" % (i,objCount)] = int(getattr(rtrow,'%sIsTracker' %l)) if l[0]=='m' else int(-1)
                        ntupleRow["%s.ValidMuonHits%i" % (i,objCount)] = int(getattr(rtrow,'%sMuonHits' %l)) if l[0]=='m' else int(-1)
                        ntupleRow["%s.MatchedStations%i" % (i,objCount)] = int(getattr(rtrow,'%sMatchedStations' %l)) if l[0]=='m' else int(-1)
                        ntupleRow["%s.ValidPixelHits%i" % (i,objCount)] = int(getattr(rtrow,'%sPixHits' %l)) if l[0]=='m' else int(-1)
                        ntupleRow["%s.TrackerLayers%i" % (i,objCount)] = int(getattr(rtrow,'%sTkLayersWithMeasurement' %l)) if l[0]=='m' else int(-1)
                        # manually add w z deltaRs
                        if i=='w1' and len(theObjects)==3:
                            oZ1 = ordered(theObjects[0],theObjects[2]) if theObjects else []
                            oZ2 = ordered(theObjects[1],theObjects[2]) if theObjects else []
                            ntupleRow["w1.dR1_z1_1"] = self.getObject(rtrow,"dr",oZ1[0],oZ1[1]) if theObjects else float(-9)
                            ntupleRow["w1.dR1_z1_2"] = self.getObject(rtrow,"dr",oZ2[0],oZ2[1]) if theObjects else float(-9)
                            ntupleRow["w1.mll_z1_1"] = self.getObject(rtrow,"mass",oZ1[0],oZ1[1]) if theObjects else float(-9)
                            ntupleRow["w1.mll_z1_2"] = self.getObject(rtrow,"mass",oZ2[0],oZ2[1]) if theObjects else float(-9)
                        if i=='w1' and theObjects and self.period==13:
                            lEta = self.getObject(rtrow,'eta',theObjects[-1])
                            lPhi = self.getObject(rtrow,'phi',theObjects[-1])
                            # TODO implement jets
                            jEta = rtrow.jet1Eta
                            jPhi = rtrow.jet1Phi
                            dr = deltaR(lEta,lPhi,jEta,jPhi)
                            ntupleRow["w1.dR1_leadJet"] = float(dr)
                        # do alternate IDs
                        for altId in self.alternateIds:
                            ntupleRow["%s.pass_%s_%i"%(i,altId,objCount)] = int(self.ID(rtrow,l,**self.alternateIdMap[altId]) if theObjects else float(-9))
                objStart += numObjects


        # initial state objects
        store_state(rtrow,ntupleRow,self.initial_states,objects,self.period)

        # alternative state objects
        if hasattr(self,'other_states'):
            for state in self.other_states:
                store_state(rtrow,ntupleRow,state,self.choose_alternative_objects(rtrow,state),self.period)
                

        # final state objects
        lepCount = 0
        jetCount = 0
        phoCount = 0
        orderedAllObjects = sorted(objects, key = lambda x: getattr(rtrow,"%sPt" % x), reverse=True)
        for obj in orderedAllObjects:
            if obj[0] in 'emt':
                charName = 'l'
                lepCount += 1
                objCount = lepCount
            if obj[0] == 'j':
                charName = 'j'
                jetCount += 1
                objCount = jetCount
            if obj[0] == 'g':
                charName = 'g'
                phoCount += 1
                objCount = phoCount
            ntupleRow["%s%i.Pt" % (charName,objCount)] = self.getObject(rtrow, "pt", obj)
            ntupleRow["%s%i.Eta" % (charName,objCount)] = self.getObject(rtrow, "eta", obj)
            ntupleRow["%s%i.Phi" % (charName,objCount)] = self.getObject(rtrow, "phi", obj)
            # TODO isolation implementation with shifts
            if obj[0]=='e': isoVar = 'RelPFIsoRho'
            if obj[0]=='m': isoVar = 'RelPFIsoDBDefault'
            ntupleRow["%s%i.Iso" % (charName,objCount)] = float(getattr(rtrow, "%s%s" % (obj, isoVar))) if obj[0] in 'em' else float(-1.)
            ntupleRow["%s%i.SigmaIEtaIEta" % (charName,objCount)] = float(getattr(rtrow, "%sSigmaIEtaIEta" % (obj))) if obj[0] in 'e' else float(-1.)
            ntupleRow["%s%i.DEtaIn" % (charName,objCount)] = float(getattr(rtrow, "%sdeltaEtaSuperClusterTrackAtVtx" % (obj))) if obj[0] in 'e' else float(-9.)
            ntupleRow["%s%i.DPhiIn" % (charName,objCount)] = float(getattr(rtrow, "%sdeltaPhiSuperClusterTrackAtVtx" % (obj))) if obj[0] in 'e' else float(-9.)
            ntupleRow["%s%i.HOverE" % (charName,objCount)] = float(getattr(rtrow, "%sHadronicOverEM" % (obj))) if obj[0] in 'e' else float(-1.)
            ntupleRow["%s%i.OoEmOoP" % (charName,objCount)] = float(abs((1.-getattr(rtrow, "%seSuperClusterOverP" % (obj)))*1./getattr(rtrow, "%secalEnergy" % (obj)))) if obj[0] in 'e' else float(-1.)
            ntupleRow["%s%i.TriggeringMVA" % (charName,objCount)] = float(-9. if self.period==13 else getattr(rtrow, "%sMVATrig" % (obj))) if obj[0] in 'e' else float(-9.)
            ntupleRow["%s%i.NonTriggeringMVA" % (charName,objCount)] = float(getattr(rtrow, "%sMVANonTrigID" % (obj)) if self.period==13 else getattr(rtrow, "%sMVANonTrig" % (obj))) if obj[0] in 'e' else float(-9.)
            ntupleRow["%s%i.NormalizedChi2" % (charName,objCount)] = float(getattr(rtrow, "%sNormTrkChi2" % (obj))) if obj[0] in 'm' else float(-1.)
            ntupleRow["%s%i.JetPt" % (charName,objCount)] = float(getattr(rtrow, "%sJetPt" % obj)) if obj[0] in 'emt' else float(-1.)
            if obj[0] in 'emt':
                ntupleRow["%s%i.JetBTag" % (charName,objCount)] = float(getattr(rtrow, "%sJetCSVBtag" % obj)) if self.period==8 else float(getattr(rtrow, "%sJetPFCISVBtag" % obj))
                ntupleRow["%s%i.Dxy" % (charName,objCount)] = float(getattr(rtrow, "%sPVDXY" % obj))
                ntupleRow["%s%i.Dz" % (charName,objCount)] = float(getattr(rtrow, "%sPVDZ" % obj))
            looseScales = self.lepscaler.scale_factor(rtrow, obj, lepType='Loose', period=self.period)
            mediumScales = self.lepscaler.scale_factor(rtrow, obj, lepType='Medium', period=self.period)
            tightScales = self.lepscaler.scale_factor(rtrow, obj, lepType='Tight', period=self.period)
            mediumLepeff = self.lepeff.scale_factor(rtrow, obj, denom='Loose', numer='Medium', period=self.period)
            tightLepeff = self.lepeff.scale_factor(rtrow, obj, denom='Loose', numer='Tight', period=self.period)
            mediumFake = self.lepfake.scale_factor(rtrow, obj, denom='Loose', numer='Medium', period=self.period)
            tightFake = self.lepfake.scale_factor(rtrow, obj, denom='Loose', numer='Tight', period=self.period)
            mediumFakeMC = self.lepfake.scale_factor(rtrow, obj, denom='Loose', numer='Medium', period=self.period, mc=True)
            tightFakeMC = self.lepfake.scale_factor(rtrow, obj, denom='Loose', numer='Tight', period=self.period, mc=True)
            ntupleRow["%s%i.LepScaleLoose" % (charName,objCount)] = float(looseScales[0])
            ntupleRow["%s%i.LepScaleMedium" % (charName,objCount)] = float(mediumScales[0])
            ntupleRow["%s%i.LepScaleTight" % (charName,objCount)] = float(tightScales[0])
            ntupleRow["%s%i.LepScaleLoose_up" % (charName,objCount)] = float(looseScales[1])
            ntupleRow["%s%i.LepScaleMedium_up" % (charName,objCount)] = float(mediumScales[1])
            ntupleRow["%s%i.LepScaleTight_up" % (charName,objCount)] = float(tightScales[1])
            ntupleRow["%s%i.LepScaleLoose_down" % (charName,objCount)] = float(looseScales[2])
            ntupleRow["%s%i.LepScaleMedium_down" % (charName,objCount)] = float(mediumScales[2])
            ntupleRow["%s%i.LepScaleTight_down" % (charName,objCount)] = float(tightScales[2])
            ntupleRow["%s%i.LepEffMedium" % (charName,objCount)] = float(mediumLepeff[0])
            ntupleRow["%s%i.LepEffMedium_up" % (charName,objCount)] = float(mediumLepeff[1])
            ntupleRow["%s%i.LepEffMedium_down" % (charName,objCount)] = float(mediumLepeff[2])
            ntupleRow["%s%i.LepEffTight" % (charName,objCount)] = float(tightLepeff[0])
            ntupleRow["%s%i.LepEffTight_up" % (charName,objCount)] = float(tightLepeff[1])
            ntupleRow["%s%i.LepEffTight_down" % (charName,objCount)] = float(tightLepeff[2])
            ntupleRow["%s%i.LepFakeMedium" % (charName,objCount)] = float(mediumFake[0])
            ntupleRow["%s%i.LepFakeMedium_up" % (charName,objCount)] = float(mediumFake[1])
            ntupleRow["%s%i.LepFakeMedium_down" % (charName,objCount)] = float(mediumFake[2])
            ntupleRow["%s%i.LepFakeTight" % (charName,objCount)] = float(tightFake[0])
            ntupleRow["%s%i.LepFakeTight_up" % (charName,objCount)] = float(tightFake[1])
            ntupleRow["%s%i.LepFakeTight_down" % (charName,objCount)] = float(tightFake[2])
            ntupleRow["%s%i.LepFakeMCMedium" % (charName,objCount)] = float(mediumFakeMC[0])
            ntupleRow["%s%i.LepFakeMCMedium_up" % (charName,objCount)] = float(mediumFakeMC[1])
            ntupleRow["%s%i.LepFakeMCMedium_down" % (charName,objCount)] = float(mediumFakeMC[2])
            ntupleRow["%s%i.LepFakeMCTight" % (charName,objCount)] = float(tightFakeMC[0])
            ntupleRow["%s%i.LepFakeMCTight_up" % (charName,objCount)] = float(tightFakeMC[1])
            ntupleRow["%s%i.LepFakeMCTight_down" % (charName,objCount)] = float(tightFakeMC[2])
            ntupleRow["%s%i.Chg" % (charName,objCount)] = float(getattr(rtrow, "%sCharge" % obj))
            ntupleRow["%s%i.PassLoose" % (charName,objCount)] = float(self.ID(rtrow,obj,**self.getIdArgs('Loose')))
            ntupleRow["%s%i.PassMedium" % (charName,objCount)] = float(self.ID(rtrow,obj,**self.getIdArgs('Medium')))
            ntupleRow["%s%i.PassTight" % (charName,objCount)] = float(self.ID(rtrow,obj,**self.getIdArgs('Tight')))
            ntupleRow["%s%iFlv.Flv" % (charName,objCount)] = obj[0]
            ntupleRow["%s%i.GenIsPrompt" % (charName,objCount)] = -2000
            ntupleRow["%s%i.GenPdgId" % (charName,objCount)] = -2000
            ntupleRow["%s%i.GenPatPdgId" % (charName,objCount)] = -2000
            ntupleRow["%s%i.MotherGenPdgId" % (charName,objCount)] = -2000
            if not self.isData and obj[0] in 'emt':
                ntupleRow["%s%i.GenIsPrompt" % (charName,objCount)] = float(getattr(rtrow, "%sGenPrompt" % obj))
                ntupleRow["%s%i.GenPdgId" % (charName,objCount)] = float(getattr(rtrow, "%sGenPdgId" % obj))
                ntupleRow["%s%i.GenPatPdgId" % (charName,objCount)] = float(getattr(rtrow, "%sGenParticle" % obj))
                ntupleRow["%s%i.MotherGenPdgId" % (charName,objCount)] = float(getattr(rtrow, "%sGenMotherPdgId" % obj))
            ntupleRow["%s%i.ChargeConsistent" % (charName,objCount)] = int(getattr(rtrow,'%sChargeIdTight' %obj)) if obj[0]=='e' and self.period==8 else int(-1)
            ntupleRow["%s%i.ExpectedMissingInnerHits" % (charName,objCount)] = int(getattr(rtrow,'%sMissingHits' %obj)) if obj[0]=='e' else int(-1)
            ntupleRow["%s%i.PassConversionVeto" % (charName,objCount)] = int(getattr(rtrow,'%sPassesConversionVeto' %obj) if self.period==13 else not getattr(rtrow,'%sHasConversion' %obj)) if obj[0]=='e' else int(-1)
            ntupleRow["%s%i.IsGlobalMuon" % (charName,objCount)] = int(getattr(rtrow,'%sIsGlobal' %obj)) if obj[0]=='m' else int(-1)
            ntupleRow["%s%i.IsPFMuon" % (charName,objCount)] = int(getattr(rtrow,'%sIsPFMuon' %obj)) if obj[0]=='m'  and self.period==13 else int(-1)
            ntupleRow["%s%i.IsTrackerMuon" % (charName,objCount)] = int(getattr(rtrow,'%sIsTracker' %obj)) if obj[0]=='m' else int(-1)
            ntupleRow["%s%i.ValidMuonHits" % (charName,objCount)] = int(getattr(rtrow,'%sMuonHits' %obj)) if obj[0]=='m' else int(-1)
            ntupleRow["%s%i.MatchedStations" % (charName,objCount)] = int(getattr(rtrow,'%sMatchedStations' %obj)) if obj[0]=='m' else int(-1)
            ntupleRow["%s%i.ValidPixelHits" % (charName,objCount)] = int(getattr(rtrow,'%sPixHits' %obj)) if obj[0]=='m' else int(-1)
            ntupleRow["%s%i.TrackerLayers" % (charName,objCount)] = int(getattr(rtrow,'%sTkLayersWithMeasurement' %obj)) if obj[0]=='m' else int(-1)

        return ntupleRow

    def write_row(self, nrow):
        '''
        Function to write the ntuple row to the tree.
        '''
        for key,val in nrow.iteritems():
            branch, var = key.split('.')
            setattr(self.branches[branch],var,val)

    def pass_preselection(self, rtrow):
        '''
        Wrapper for preselection defined by user.
        '''
        if 'preselection' in self.cache: return self.cache['preselection']
        cuts = self.preselection(rtrow)
        cutResults = cuts.evaluate(rtrow)
        self.cache['cutflow'] = cuts
        self.cache['preselection'] = cutResults
        return cutResults

    def pass_selection(self,rtrow):
        '''
        Wrapper for the selection defined by the user (tight selection whereas preselection
        is the loose selection for fake rate method).
        '''
        if 'selection' in self.cache: return self.cache['selection']
        cuts = self.selection(rtrow)
        cutResults = cuts.evaluate(rtrow)
        self.cache['selection'] = cutResults
        return cutResults

    def npass(self,rtrow,numObjects,**kwargs):
        '''
        Get the number that pass the full selection for fake rate method.
        numObjects is dictionary with structure:
          numObjects = {
            'e': ne,
            'm': nm,
            't': nt,
            ...
          }
        All keys are optional, if a key is not present the number is assumed to be zero
        and no check will be performed. keys = ['e', 'm', 't', 'j', 'g']
        kwargs are the arguments for the ID definition.
        '''
        passLoose = self.pass_preselection(rtrow)
        if not passLoose: return False
        for obj, num in numObjects.iteritems():
            if self.numPassID(rtrow,obj,**kwargs)!=num: return False
        return True

    def numPassID(self,rtrow,flav,**kwargs):
        '''
        Uses self.ID() to count number of a given object that pass the ID.
        '''
        num = 0
        for obj in self.objects:
            if obj[0]==flav:
                num += self.ID(rtrow,obj,**kwargs)
        return num

    def ID(self,rtrow,*objects,**kwargs):
        '''
        An ID accessor method.
        '''
        idDef = kwargs.pop('idDef',{})
        isoCut = kwargs.pop('isoCut',{})
        for obj in objects:
            if obj[0] not in idDef: continue
            type = idDef[obj[0]]
            if 'ID_%s_%s' %(type,obj) in self.cache:
                if not self.cache['ID_%s_%s'%(type,obj)]: return False
            else:
                result = lepId.lep_id(rtrow,self.period,obj,idType=idDef[obj[0]],metShift=self.metShift)
                self.cache['ID_%s_%s'%(type,obj)] = result
                if not result: return False
        # TODO support iso cut with shift
        if isoCut:
            for obj in objects:
                if obj[0] not in isoCut: continue
                if isoCut[obj[0]]<=0.: continue
                if obj[0] == 'e':
                    isotype = "RelPFIsoRho"
                if obj[0] == 'm':
                    isotype = "RelPFIsoDBDefaultR04"
                if obj[0] in 'tjgn': continue # no iso cut on tau
                if getattr(rtrow, '%s%s' %(obj,isotype)) > isoCut[obj[0]]: return False
        return True

    def getScales(self,rtrow,objects,**lepargs):
        '''Return the scale factors in a dictionary'''
        scales = {
            'lep'             : 1,
            'lepe'            : 1,
            'lepm'            : 1,
            'lepup'           : 1,
            'lepeup'          : 1,
            'lepmup'          : 1,
            'lepdown'         : 1,
            'lepedown'        : 1,
            'lepmdown'        : 1,
            'lepfake'         : 0,
            'lepfakeup'       : 0,
            'lepfakedown'     : 0,
            'trig'            : 1,
            'trigup'          : 1,
            'trigdown'        : 1,
            'puweight'        : 1,
            'puweightup'      : 1,
            'puweightdown'    : 1,
            'genweight'       : 1,
            'chargeid'        : 1,
            'trigger_prescale': 1,
        }
        if self.period==8: # TODO: move when we have numbers for 13 tev
            chargeid  = self.chargeid.systematic(rtrow, *objects, period=self.period)
            scales['chargeid'] = chargeid
        trigeff_data = self.trigscaler.scale_factor(rtrow, *objects, period=self.period, useData=True, metShift=self.metShift)
        trigeffup_data = self.trigscaler.scale_factor(rtrow, *objects, period=self.period, shiftUp=True, useData=True, metShift=self.metShift)
        trigeffdown_data = self.trigscaler.scale_factor(rtrow, *objects, period=self.period, shiftDown=True, useData=True, metShift=self.metShift)
        trigeff_mc = self.trigscaler.scale_factor(rtrow, *objects, period=self.period, useData=False, metShift=self.metShift)
        trigeffup_mc = self.trigscaler.scale_factor(rtrow, *objects, period=self.period, shiftUp=True, useData=False, metShift=self.metShift)
        trigeffdown_mc = self.trigscaler.scale_factor(rtrow, *objects, period=self.period, shiftDown=True, useData=False, metShift=self.metShift)
        trigscale = trigeff_data/trigeff_mc if trigeff_mc else 1.
        trigscaleup = trigeffup_data/trigeffup_mc if trigeffup_mc else 1.
        trigscaledown = trigeffdown_data/trigeffdown_mc if trigeffdown_mc else 1.
        scales['trig']     = trigscale
        scales['trigup']   = trigscaleup
        scales['trigdown'] = trigscaledown
        puweight  = self.pu_weights.weight(rtrow, period=self.period)
        scales['puweight'] = puweight[0]
        scales['puweightup'] = puweight[1]
        scales['puweightdown'] = puweight[2]
        # do different based on category
        lepscales = [1.,1.,1.]
        lepescales = [1.,1.,1.]
        lepmscales = [1.,1.,1.]
        wl = self.objCand[-1]
        for obj in objects:
            passLoose = self.ID(rtrow,obj,**self.getIdArgs('Loose'))
            passMedium = self.ID(rtrow,obj,**self.getIdArgs('Medium'))
            passTight = self.ID(rtrow,obj,**self.getIdArgs('Tight'))
            if self.tightW: 
                lepType = 'Medium' if passMedium else 'Loose'
                if obj==wl and passTight: lepType = 'Tight'
                ls = self.lepscaler.scale_factor(rtrow, obj, period=self.period, metShift=self.metShift, lepType=lepType **lepargs)
            if not self.tightW: ls = self.lepscaler.scale_factor(rtrow, obj, period=self.period, metShift=self.metShift, lepType='Tight' if passTight else 'Loose', **lepargs)
            lepscales[0] *= ls[0]
            lepscales[1] *= ls[1]
            lepscales[2] *= ls[2]
            if obj[0]=='e':
                lepescales[0] *= ls[0]
                lepescales[1] *= ls[1]
                lepescales[2] *= ls[2]
            else:
                lepescales[0] *= ls[0]
                lepescales[1] *= ls[0] # dont shift
                lepescales[2] *= ls[0]
            if obj[0]=='m':
                lepmscales[0] *= ls[0]
                lepmscales[1] *= ls[1]
                lepmscales[2] *= ls[2]
            else:
                lepmscales[0] *= ls[0] # dont shift
                lepmscales[1] *= ls[0]
                lepmscales[2] *= ls[0]

        #lepscales = self.lepscaler.scale_factor(rtrow, *objects, period=self.period, **lepargs)
        scales['lep']      = lepscales[0]
        scales['lepe']     = lepescales[0]
        scales['lepm']     = lepmscales[0]
        scales['lepup']    = lepscales[1]
        scales['lepeup']   = lepescales[1]
        scales['lepmup']   = lepmscales[1]
        scales['lepdown']  = lepscales[2]
        scales['lepedown'] = lepescales[2]
        scales['lepmdown'] = lepmscales[2]
        passtight = []
        for obj in objects:
            if self.tightW: passtight += [1 if self.ID(rtrow,obj,**self.getIdArgs('Tight' if obj==wl else 'Medium')) else 0]
            if not self.tightW: passtight += [1 if self.ID(rtrow,obj,**self.getIdArgs('Tight')) else 0]
        totalfake = [1.,1.,1.]
        sign = 1. if sum(passtight) in [3,2,0] else -1.
        for obj in objects:
            if self.tightW: lepfake = self.lepfake.scale_factor(rtrow, obj, period=self.period, metShift=self.metShift, denom='Loose', numer='Tight' if obj==wl else 'Medium', **lepargs)
            if not self.tightW: lepfake = self.lepfake.scale_factor(rtrow, obj, period=self.period, metShift=self.metShift, denom='Loose', numer='Tight', **lepargs)
            if not passtight[objects.index(obj)]:
                totalfake[0] *= lepfake[0]
                totalfake[1] *= lepfake[1]
                totalfake[2] *= lepfake[2]
        totalfake = [x*sign for x in totalfake]
        scales['lepfake'] = totalfake[0]
        scales['lepfakeup'] = totalfake[1]
        scales['lepfakedown'] = totalfake[2]
        scales['trigger_prescale'] = self.getTriggerPrescale(rtrow)
        genweight = rtrow.GenWeight if hasattr(rtrow,'GenWeight') else 1.
        scales['genweight'] = genweight
        return scales

    def getTriggerPrescale(self,rtrow):
        '''Default trigger prescale'''
        return 1

    def getGenChannel(self,rtrow):
        '''Dummy return gen channel string'''
        return 'a'

    def returnTrue(self,rtrow):
        return True

    def getFakeChannel(self,rtrow,**kwargs):
        chan = ''
        for obj in self.objCand:
            chan += 'P' if self.ID(rtrow,obj,**self.getIdArgs('Tight')) else 'F'
        return chan

    def getObject(self,rtrow,var,*objs,**kwargs):
        '''Get modified object'''
        key = 'getObject' + var
        for obj in objs: key += obj
        if key in self.cache: return self.cache[key]
        if len(objs)==1 and objs[0]=='met': # get met with shift
            if var=='eta': return 0.
            metVar = 'type1_pfMet'
            #shift = kwargs.pop('shift','')
            shiftStrings = {
                'ees+': 'ElectronEnUp',
                'ees-': 'ElectronEnDown',
                'mes+': 'MuonEnUp',
                'mes-': 'MuonEnDown',
                'tes+': 'TauEnUp',
                'tes-': 'TauEnDown',
                'ues+': 'UnclusteredEnUp',
                'ues-': 'UnclusteredEnDown',
                'jes+': 'JetEnUp',
                'jes-': 'JetEnDown',
                'jres+': 'JetResUp',
                'jres-': 'JetResDown',
            }
            shift = shiftStrings[self.metShift] if self.metShift else '' # shift MC only
            varName = {'pt': 'Pt', 'phi': 'Phi'}
            varName2 = {'pt': 'Et', 'phi': 'Phi'}
            shiftString = '{0}_shifted{1}_{2}'.format(metVar,varName[var],shift) if shift else '{0}{1}'.format(metVar,varName2[var])
            val = getattr(rtrow,shiftString)
        elif len(objs)==1 and objs[0] in self.object_definitions: # composite object, must pass 'objects' kwarg
            obj = objs[0]
            masses = {'e':0.511e-3, 'm':0.1056, 't':1.776, 'j':0}
            if self.object_definitions[obj] == 1:
                val = -9.
            elif self.object_definitions[obj][1] == 'n': # obj + met
                o = kwargs.get('objects',[])
                if not o: val = -9
                shift = self.metShift
                val = self.getObject(rtrow,var,o[0],'met',**kwargs)
            else: # assumes two lepton constituents
                o = kwargs.get('objects',[])
                if len(o) != 2: val = -9
                o = ordered(o[0],o[1])
                val = self.getObject(rtrow,var,*o,**kwargs)

        elif len(objs)==1 and objs[0][0] in ['e','m','t']: # lepton
            obj = objs[0]
            shift = self.metShift
            shiftString = ''
            masses = {'e':0.511e-3, 'm':0.1056, 't':1.776, 'j':0}
            if shift: # shift MC only
                shiftMap = {
                    'ees+' : '_ElectronEnUp',
                    'ees-' : '_ElectronEnDown',
                    'mes+' : '_MuonEnUp',
                    'mes-' : '_MuonEnDown',
                }
                if obj[0] == 'e' and shift in shiftMap and shift[0] == 'e':
                    shiftString = shiftMap[shift]
                if obj[0] == 'm' and shift in shiftMap and shift[0] == 'm':
                    shiftString = shiftMap[shift]
            # pt
            if var=='pt' or var=='st': val = getattr(rtrow,'{0}Pt{1}'.format(obj,shiftString))
            # eta
            if var=='eta': val = getattr(rtrow,'{0}Eta{1}'.format(obj,shiftString))
            # phi
            if var=='phi': val = getattr(rtrow,'{0}Phi{1}'.format(obj,shiftString))
            # mass
            if var=='mass': val = masses[objs[0]]

        elif len(objs)>1: # multiple objects
            vecs = []
            pts = []
            etas = []
            phis = []
            masses = {'e':0.511e-3, 'm':0.1056, 't':1.776, 'j':0}
            for obj in objs:
                pt = self.getObject(rtrow,'pt',obj,**kwargs)
                eta = self.getObject(rtrow,'eta',obj,**kwargs)
                phi = self.getObject(rtrow,'phi',obj,**kwargs)
                if obj[0] in masses and obj!='met':
                    mass = masses[obj[0]]
                else:
                    mass = 0.
                vec = rt.TLorentzVector()
                vec.SetPtEtaPhiM(pt,eta,phi,mass)
                vecs += [vec]
                pts += [pt]
                etas += [eta]
                phis += [phi]
            vec = rt.TLorentzVector()
            for v in vecs: vec += v
            if var=='pt':
                val = vec.Pt()
            if var=='mass' or var=='m':
                val = vec.M()
            if var=='mt':
                val = 0.
                #if val < 20:
                #    print '{0}: mt = {1}'.format(', '.join(objs),val)
                #    print '    pts: {0}'.format(', '.join(['{0}: {1}'.format(o,p) for o,p in zip(objs,pts)]))
                #    print '    etas: {0}'.format(', '.join(['{0}: {1}'.format(o,e) for o,e in zip(objs,etas)]))
                #    print '    phis: {0}'.format(', '.join(['{0}: {1}'.format(o,p) for o,p in zip(objs,phis)]))
                #    mt = math.sqrt(2*pts[0]*pts[1]*(1-math.cos(deltaPhi(phis[0],phis[1]))))
                #    print '    mt by hand: {0}'.format(mt)
                #    print '    mt from fsa: {0}'.format(getattr(rtrow,'{0}MtToPfMet_type1'.format(objs[0])))
                if len(objs)>1 and objs[-1]=='met':
                    #val = math.sqrt(2*pts[0]*pts[1]*(1-math.cos(deltaPhi(phis[0],phis[1]))))
                    #val = getattr(rtrow,'{0}MtToPfMet_type1'.format(objs[0]))
                    candVec = rt.TLorentzVector()
                    for v in vecs[:-1]: candVec += v
                    metVec = vecs[-1]
                    val = math.sqrt(abs((candVec.Et()+metVec.Et())**2 - (vec.Pt())**2))
            if var=='eta':
                val = vec.Eta()
            if var=='phi':
                val = vec.Phi()
            if var=='dr' and len(objs)==2:
                val = deltaR(etas[0],phis[0],etas[1],phis[1])
            if var=='dr' and len(objs)!=2:
                print 'WHAT ARE YOU DOING!'
            if var=='dphi' and len(objs)==2:
                val = deltaPhi(phis[0],phis[1])
            if var=='dphi' and len(objs)!=2:
                print 'WHAT ARE YOU DOING!'
            if var=='st':
                val = sum(pts)

        else:
            val = -9
        self.cache[key] = val
        return val

    def getObjectVetos(self,rtrow,flv,vetoType):
        key = 'getObjectVetos'+flv+vetoType
        if key in self.cache: return self.cache[key]
        vetoMap = {
            'e' : {
                'Tight': 'eVetoTight',
                'Loose': 'eVetoLoose',
            },
            'm' : {
                'Tight': 'muVetoMedium',
                'Loose': 'muVetoLoose',
            },
        }
        vetoString = vetoMap[flv][vetoType]
        if self.metShift: 
            if self.metShift=='mes+' and flv=='m': vetoString += '_mesUp'
            if self.metShift=='mes-' and flv=='m': vetoString += '_mesDown'
            if self.metShift=='ees+' and flv=='e': vetoString += '_eesUp'
            if self.metShift=='ees-' and flv=='e': vetoString += '_eesDown'
        val = getattr(rtrow,vetoString)
        self.cache[key] = val
        return val
