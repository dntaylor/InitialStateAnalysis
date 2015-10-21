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
        self.selections = {}
        self.results = {}

    def add(self, fun, label):
        self.selections[label] = fun

    def evaluate(self,rtrow):
        self.results['evt'] = rtrow.evt
        self.results['run'] = rtrow.run
        self.results['lumi'] = rtrow.lumi
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
    while result>rt.TMath.Pi(): result -= 2*rt.TMath.Pi()
    while result<=rt.TMath.Pi(): result += 2*rt.TMath.Pi()
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
            with open(file_list,'r') as files:
                for f in files:
                    self.file_names += [f]
        else:                                 # the list is a python list
            self.file_names = file_list
        self.out_file = out_file
        self.period = period

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
        self.ntuple, self.branches = buildNtuple(self.object_definitions,states,self.channel,self.final_states,altIds=self.alternateIds,doVBF=self.doVBF)

        if hasattr(self,'cutTreeLabels'):
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
        eventsToWrite = set()
        eventsWritten = set()
        numEvts = 0
        totalWritten = 0

        # iterate over files
        for i, file_name in enumerate(self.file_names):
            self.file_name = file_name
            logger.info('%s %s Processing %i/%i files', self.channel, self.sample_name, i+1, len(self.file_names))
            sys.stdout.flush()
            if file_name.startswith('/store'): file_name = 'root://cmsxrootd.hep.wisc.edu//%s' % file_name

            rtFile = rt.TFile.Open(file_name, "READ")

            # iterate over final states
            for fs in self.final_states:
                if len(self.file_names)<10: logger.info('%s %s %s' % (self.channel, self.sample_name, fs))
                tree = rtFile.Get("%s/final/Ntuple" % fs)
                #if self.period==8:
                metatree = rtFile.Get("%s/metaInfo" % fs)
                tempEvts = 0
                for entry in xrange(metatree.GetEntries()):
                    metatree.GetEntry(entry)
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
                    eventkey = (rtrow.evt, rtrow.lumi, rtrow.run)

                    # if we have a cutTree, do it
                    if hasattr(self,'cutTree'):
                        eventCutTree = 0

                    # can we define the object we want?
                    candidate = self.choose_objects(rtrow)
                    if not candidate: # in case no objects satisfy our conditions
                        continue

                    # name candidate
                    self.objCand = candidate[1]

                    # now see if event is viable
                    passPreselection = self.pass_preselection(rtrow)

                    # check preselection
                    if not passPreselection:
                        continue

                    # check combinatorics
                    if eventkey in bestCandMap: 
                        bestcand = bestCandMap[eventkey]
                    else:
                        numMin = len(self.objCand)
                        bestcand = [float('inf')] * numMin
                    if self.good_to_store(rtrow,candidate[0],bestcand):
                        bestCandMap[eventkey] = candidate[0]
                        ntupleRow = self.store_row(rtrow, *self.objCand)
                        eventMap[eventkey] = ntupleRow
                        eventsToWrite.add(eventkey)

            rtFile.Close("R")
            numEvts += tempEvts

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

        logger.info('%s %s Filled Tree (%i events)' % (self.channel, self.sample_name, len(eventsWritten)))

        # now we store the total processed events
        logger.info('%s %s Processed %i events' % (self.channel, self.sample_name, numEvts))

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
        numObjs = len(self.final_states[0])
        finalStateObjects ='emtjgn'
        allowedObjects = ''
        for fsObj in finalStateObjects:
            for fs in self.final_states:
                if fsObj in fs:
                    allowedObjects += fsObj
                    break
        numObjTypes = len(allowedObjects)
        #for prompts in list(product(range(numObjs+1),repeat=numObjTypes)):
        #    if sum(prompts) > numObjs: continue
        #    promptString = ''.join([str(x) for x in prompts])
        #    promptDict = {}
        #    for o in range(len(allowedObjects)):
        #        promptDict[allowedObjects[o]] = prompts[o]
        #    ntupleRow["select.pass_%s"%promptString] = int(self.npass(rtrow,promptDict,**self.getIdArgs('Tight')))
        for altId in self.alternateIds:
            ntupleRow["select.pass_%s"%altId] = int(self.ID(rtrow,*self.objects,**self.alternateIdMap[altId]))

        scales = self.getScales(rtrow,objects,**self.lepargs)
        
        ntupleRow["event.evt"] = int(rtrow.evt)
        ntupleRow["event.lumi"] = int(rtrow.lumi)
        ntupleRow["event.run"] = int(rtrow.run)
        ntupleRow["event.nvtx"] = int(rtrow.nvtx)
        ntupleRow["event.GenNUP"] = -1 if self.isData else int(rtrow.NUP)
        ntupleRow["event.trig_prescale"] = int(scales['trigger_prescale'])
        ntupleRow["event.lep_scale"] = float(scales['lep'])
        ntupleRow["event.lep_scale_up"] = float(scales['lepup'])
        ntupleRow["event.lep_scale_down"] = float(scales['lepdown'])
        ntupleRow["event.trig_scale"] = float(scales['trig'])
        ntupleRow["event.pu_weight"] = float(scales['puweight'][0])
        ntupleRow["event.pu_weight_up"] = float(scales['puweight'][1])
        ntupleRow["event.pu_weight_down"] = float(scales['puweight'][2])
        ntupleRow["event.gen_weight"] = float(scales['genweight'])
        ntupleRow["event.charge_uncertainty"] = float(scales['chargeid'])

        channelString = ''
        for x in objects: channelString += x[0]
        ntupleRow["channel.channel"] = channelString
        ntupleRow["genChannel.channel"] = self.getGenChannel(rtrow)

        def getMT(rtrow,period,*objects):
            masses = {'e':0.511e-3, 'm':0.1056, 't':1.776, 'j':0}
            metVar = 'type1_pfMet'
            mtVar = 'PFMET' if period==13 else 'PfMet_Ty1'
            visP4 = rt.TLorentzVector()
            for l in objects:
                p4 = rt.TLorentzVector()
                pt = getattr(rtrow, "%sPt" % l)
                eta = getattr(rtrow, "%sEta" % l)
                phi = getattr(rtrow, "%sPhi" % l)
                mass = masses[l[0]]
                p4.SetPtEtaPhiM(pt,eta,phi,mass)
                visP4 += p4
            ptMet = getattr(rtrow, "%sEt" % metVar)
            phiMet = getattr(rtrow, "%sPhi" % metVar)
            pxMet = ptMet*rt.TMath.Cos(phiMet)
            pyMet = ptMet*rt.TMath.Sin(phiMet)
            mt = rt.TMath.Sqrt(visP4.M2() + 2*(ptMet*visP4.Et() - (pxMet*visP4.Px() + pyMet*visP4.Py())))
            return mt

        ntupleRow["finalstate.mass"] = float(rtrow.Mass)
        ntupleRow["finalstate.eta"] = float(rtrow.Eta) if self.period==13 else float(-9)
        ntupleRow["finalstate.phi"] = float(rtrow.Phi) if self.period==13 else float(-9)
        ntupleRow["finalstate.mT"] = float(getMT(rtrow,self.period,*objects))
        ntupleRow["finalstate.sT"] = float(sum([getattr(rtrow, "%sPt" % x) for x in objects]))
        ntupleRow["finalstate.hT"] = float(rtrow.Ht) if hasattr(rtrow,'Ht') else float(-1)
        metVar = 'type1_pfMet'
        ntupleRow["finalstate.met"] = float(getattr(rtrow, '%sEt' %metVar))
        ntupleRow["finalstate.metPhi"] = float(getattr(rtrow,'%sPhi' %metVar))
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
        ntupleRow["finalstate.muonVetoTight"] = int(rtrow.muVetoTight) if self.period==13 else int(rtrow.muonVetoWZTight)
        ntupleRow["finalstate.elecVetoTight"] = int(rtrow.eVetoTight) if self.period==13 else int(rtrow.elecVetoWZTight)
        ntupleRow["finalstate.muonVetoLoose"] = int(rtrow.muVeto) if self.period==13 else int(rtrow.muVetoPt5IsoIdVtx)
        ntupleRow["finalstate.elecVetoLoose"] = int(rtrow.eVeto) if self.period==13 else int(rtrow.eVetoMVAIsoVtx)
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
            masses = {'e':0.511e-3, 'm':0.1056, 't':1.776, 'j':0}
            objStart = 0
            metVar = 'type1_pfMet'
            mtVar = 'PfMet_Ty1' if period==13 else 'PfMet_Ty1'
            for i in state:
                numObjects = len([ x for x in self.object_definitions[i] if x != 'n']) if theObjects else 0
                finalObjects = theObjects[objStart:objStart+numObjects]
                orderedFinalObjects = sorted(finalObjects, key = lambda x: getattr(rtrow,"%sPt" % x), reverse=True)
                if len(self.object_definitions[i]) == 1:
                    ntupleRow["%s.mass" %i] = float(-999)
                    ntupleRow["%s.Pt" %i] = float(-9)
                    ntupleRow["%s.Eta" %i] = float(-9)
                    ntupleRow["%s.Phi" %i] = float(-9)
                    ntupleRow["%s.sT" %i] = float(getattr(rtrow, "%sPt" % finalObjects[0])) if theObjects else float(-9)
                    ntupleRow["%s.dPhi" %i] = float(-9)
                    ntupleRow["%s.dR" %i] = float(-9)
                    ntupleRow["%sFlv.Flv" %i] = finalObjects[0][0] if theObjects else 'a'
                elif 'n' == self.object_definitions[i][1]:
                    if theObjects: # get lorentz vectors
                        pt1 = getattr(rtrow, "%sPt" % finalObjects[0])
                        eta1 = getattr(rtrow, "%sEta" % finalObjects[0])
                        phi1 = getattr(rtrow, "%sPhi" % finalObjects[0])
                        mass1 = masses[finalObjects[0][0]]
                        px1 = pt1*rt.TMath.Cos(phi1)
                        py1 = pt1*rt.TMath.Sin(phi1)
                        ptMet = getattr(rtrow, "%sEt" % metVar)
                        phiMet = getattr(rtrow, "%sPhi" % metVar)
                        pxMet = ptMet*rt.TMath.Cos(phiMet)
                        pyMet = ptMet*rt.TMath.Sin(phiMet)
                        wpt = rt.TMath.Sqrt((px1+pxMet)**2 + (py1+pyMet)**2)
                    else:
                        wpt = -9
                    ntupleRow["%s.mass" %i] = float(getattr(rtrow, "%sMtTo%s" % (finalObjects[0], mtVar))) if theObjects else float(-999)
                    ntupleRow["%s.Pt" %i] = float(wpt)
                    ntupleRow["%s.sT" %i] = float(getattr(rtrow, "%sPt" % finalObjects[0]) + getattr(rtrow, '%sEt' %metVar)) if theObjects else float(-9)
                    ntupleRow["%s.dPhi" %i] = float(getattr(rtrow, "%sToMETDPhi" % finalObjects[0])) if theObjects else float(-9)
                    ntupleRow["%sFlv.Flv" %i] = finalObjects[0][0] if theObjects else 'a'
                else:
                    finalObjOrdered = ordered(finalObjects[0], finalObjects[1]) if theObjects else []
                    #if theObjects: # get lorentz vectors
                    #    vec1 = rt.TLorentzVector()
                    #    pt1 = getattr(rtrow, "%sPt" % finalObjOrdered[0])
                    #    eta1 = getattr(rtrow, "%sEta" % finalObjOrdered[0])
                    #    phi1 = getattr(rtrow, "%sPhi" % finalObjOrdered[0])
                    #    mass1 = masses[finalObjOrdered[0][0]]
                    #    vec1.SetPtEtaPhiM(pt1,eta1,phi1,mass1)
                    #    vec2 = rt.TLorentzVector()
                    #    pt2 = getattr(rtrow, "%sPt" % finalObjOrdered[1])
                    #    eta2 = getattr(rtrow, "%sEta" % finalObjOrdered[1])
                    #    phi2 = getattr(rtrow, "%sPhi" % finalObjOrdered[1])
                    #    mass2 = masses[finalObjOrdered[1][0]]
                    #    vec2.SetPtEtaPhiM(pt2,eta2,phi2,mass2)
                    #    vec3 = vec1+vec2
                    #else:
                    #    vec3 = 0
                    ntupleRow["%s.mass" %i] = float(getattr(rtrow, "%s_%s_Mass" % (finalObjOrdered[0], finalObjOrdered[1]))) if theObjects else float(-999)
                    ntupleRow["%s.Pt" %i] = float(getattr(rtrow, "%s_%s_Pt" % (finalObjOrdered[0], finalObjOrdered[1]))) if theObjects else float(-9)
                    ntupleRow["%s.sT" %i]   = float(sum([getattr(rtrow, "%sPt" % x) for x in finalObjects])) if theObjects else float(-9)
                    ntupleRow["%s.dPhi" %i] = float(getattr(rtrow, "%s_%s_DPhi" % (finalObjOrdered[0], finalObjOrdered[1]))) if theObjects else float(-9)
                    ntupleRow["%s.dR" %i] = float(getattr(rtrow, "%s_%s_DR" % (finalObjOrdered[0], finalObjOrdered[1]))) if theObjects else float(-9)
                    ntupleRow["%sFlv.Flv" %i] = finalObjects[0][0] + finalObjects[1][0] if theObjects else 'aa'
                objCount = 0
                for obj in self.object_definitions[i]:
                    if obj=='n':
                        ntupleRow["%s.met" %i] = float(getattr(rtrow,'%sEt' %metVar)) if theObjects else float(-9)
                        ntupleRow["%s.metPhi" %i] = float(getattr(rtrow, '%sPhi' %metVar)) if theObjects else float(-9)
                    else:
                        objCount += 1
                        l = orderedFinalObjects[objCount-1] if theObjects else 'a'
                        ntupleRow["%s.Pt%i" % (i,objCount)] = float(getattr(rtrow, "%sPt" % l)) if theObjects else float(-9)
                        ntupleRow["%s.Eta%i" % (i,objCount)] = float(getattr(rtrow, "%sEta" % l)) if theObjects else float(-9)
                        ntupleRow["%s.Phi%i" % (i,objCount)] = float(getattr(rtrow, "%sPhi" % l)) if theObjects else float(-9)
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
                            looseScales = self.lepscaler.scale_factor(rtrow, l, loose=True, period=period)
                            tightScales = self.lepscaler.scale_factor(rtrow, l, loose=False, period=period)
                            looseEff = self.lepeff.scale_factor(rtrow, l, loose=True, period=period)
                            tightEff = self.lepeff.scale_factor(rtrow, l, loose=False, period=period)
                            tightFake = self.lepfake.scale_factor(rtrow, l, loose=False, period=period)
                        else:
                            looseScales = [-1,-1,-1]
                            tightScales = [-1,-1,-1]
                            looseEff = [-1,-1,-1]
                            tightEff = [-1,-1,-1]
                            tightFake = [-1,-1,-1]
                        ntupleRow["%s.LepScaleLoose%i" % (i,objCount)] = float(looseScales[0])
                        ntupleRow["%s.LepScaleTight%i" % (i,objCount)] = float(tightScales[0])
                        ntupleRow["%s.LepScaleLoose%i_up" % (i,objCount)] = float(looseScales[1])
                        ntupleRow["%s.LepScaleTight%i_up" % (i,objCount)] = float(tightScales[1])
                        ntupleRow["%s.LepScaleLoose%i_down" % (i,objCount)] = float(looseScales[2])
                        ntupleRow["%s.LepScaleTight%i_down" % (i,objCount)] = float(tightScales[2])
                        ntupleRow["%s.LepEffLoose%i" % (i,objCount)] = float(looseEff[0])
                        ntupleRow["%s.LepEffTight%i" % (i,objCount)] = float(tightEff[0])
                        ntupleRow["%s.LepFake%i" % (i,objCount)] = float(tightFake[0])
                        ntupleRow["%s.Chg%i" % (i,objCount)] = float(getattr(rtrow, "%sCharge" % l)) if theObjects else float(-9)
                        ntupleRow["%s.PassLoose%i" % (i,objCount)] = float(self.ID(rtrow,l,**self.getIdArgs('Loose'))) if theObjects else float(-9)
                        ntupleRow["%s.PassTight%i" % (i,objCount)] = float(self.ID(rtrow,l,**self.getIdArgs('Tight'))) if theObjects else float(-9)
                        ntupleRow["%s.GenPdgId%i" % (i,objCount)] = -2000
                        ntupleRow["%s.MotherGenPdgId%i" % (i,objCount)] = -2000
                        if not self.isData and theObjects:
                            ntupleRow["%s.GenPdgId%i" % (i,objCount)] = float(getattr(rtrow, "%sGenPdgId" % l))
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
                            ntupleRow["w1.dR1_z1_1"] = float(getattr(rtrow,"%s_%s_DR" % (oZ1[0],oZ1[1]))) if theObjects else float(-9)
                            ntupleRow["w1.dR1_z1_2"] = float(getattr(rtrow,"%s_%s_DR" % (oZ2[0],oZ2[1]))) if theObjects else float(-9)
                            ntupleRow["w1.mll_z1_1"] = float(getattr(rtrow,"%s_%s_Mass" % (oZ1[0],oZ1[1]))) if theObjects else float(-9)
                            ntupleRow["w1.mll_z1_2"] = float(getattr(rtrow,"%s_%s_Mass" % (oZ2[0],oZ2[1]))) if theObjects else float(-9)
                        if i=='w1' and theObjects and self.period==13:
                            lEta = getattr(rtrow,'%sEta' %theObjects[-1])
                            lPhi = getattr(rtrow,'%sPhi' %theObjects[-1])
                            jEta = rtrow.jet1Eta
                            jPhi = rtrow.jet1Phi
                            dr = math.sqrt((lEta-jEta)**2 + (lPhi-jPhi)**2)
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
            ntupleRow["%s%i.Pt" % (charName,objCount)] = float(getattr(rtrow, "%sPt" % obj))
            ntupleRow["%s%i.Eta" % (charName,objCount)] = float(getattr(rtrow, "%sEta" % obj))
            ntupleRow["%s%i.Phi" % (charName,objCount)] = float(getattr(rtrow, "%sPhi" % obj))
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
            looseScales = self.lepscaler.scale_factor(rtrow, obj, loose=True, period=self.period)
            tightScales = self.lepscaler.scale_factor(rtrow, obj, loose=False, period=self.period)
            looseEff = self.lepeff.scale_factor(rtrow, obj, loose=True, period=self.period)
            tightEff = self.lepeff.scale_factor(rtrow, obj, loose=False, period=self.period)
            tightFake = self.lepfake.scale_factor(rtrow, obj, loose=False, period=self.period)
            ntupleRow["%s%i.LepScaleLoose" % (charName,objCount)] = float(looseScales[0])
            ntupleRow["%s%i.LepScaleTight" % (charName,objCount)] = float(tightScales[0])
            ntupleRow["%s%i.LepScaleLoose_up" % (charName,objCount)] = float(looseScales[1])
            ntupleRow["%s%i.LepScaleTight_up" % (charName,objCount)] = float(tightScales[1])
            ntupleRow["%s%i.LepScaleLoose_down" % (charName,objCount)] = float(looseScales[2])
            ntupleRow["%s%i.LepScaleTight_down" % (charName,objCount)] = float(tightScales[2])
            ntupleRow["%s%i.LepEffLoose" % (charName,objCount)] = float(looseEff[0])
            ntupleRow["%s%i.LepEffTight" % (charName,objCount)] = float(tightEff[0])
            ntupleRow["%s%i.LepFake" % (charName,objCount)] = float(tightFake[0])
            ntupleRow["%s%i.Chg" % (charName,objCount)] = float(getattr(rtrow, "%sCharge" % obj))
            ntupleRow["%s%i.PassLoose" % (charName,objCount)] = float(self.ID(rtrow,obj,**self.getIdArgs('Loose')))
            ntupleRow["%s%i.PassTight" % (charName,objCount)] = float(self.ID(rtrow,obj,**self.getIdArgs('Tight')))
            ntupleRow["%s%iFlv.Flv" % (charName,objCount)] = obj[0]
            ntupleRow["%s%i.GenPdgId" % (charName,objCount)] = -2000
            ntupleRow["%s%i.MotherGenPdgId" % (charName,objCount)] = -2000
            if not self.isData and obj[0] in 'emt':
                ntupleRow["%s%i.GenPdgId" % (charName,objCount)] = float(getattr(rtrow, "%sGenPdgId" % obj))
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
                result = lepId.lep_id(rtrow,self.period,obj,idType=idDef[obj[0]])
                self.cache['ID_%s_%s'%(type,obj)] = result
                if not result: return False
        if isoCut:
            for obj in objects:
                if obj[0] not in isoCut: continue
                if obj[0] == 'e':
                    isotype = "RelPFIsoRho"
                if obj[0] == 'm':
                    isotype = "RelPFIsoDBDefault"
                if obj[0] in 'tjgn': continue # no iso cut on tau
                if getattr(rtrow, '%s%s' %(obj,isotype)) > isoCut[obj[0]]: return False
        return True

    def getScales(self,rtrow,objects,**lepargs):
        '''Return the scale factors in a dictionary'''
        scales = {
            'lep'             : 1,
            'lepup'           : 1,
            'lepdown'         : 1,
            'lepeff'          : 1,
            'lepeffup'        : 1,
            'lepeffdown'      : 1,
            'lepfake'         : 0,
            'lepfakeup'       : 0,
            'lepfakedown'     : 0,
            'trig'            : 1,
            'puweight'        : [1,1,1],
            'genweight'       : 1,
            'chargeid'        : 1,
            'trigger_prescale': 1,
        }
        if self.period==8: # TODO: move when we have numbers for 13 tev
            trigscale = self.trigscaler.scale_factor(rtrow, *objects, period=self.period)
            chargeid  = self.chargeid.systematic(rtrow, *objects, period=self.period)
            scales['trig']     = trigscale
            scales['chargeid'] = chargeid
        puweight  = self.pu_weights.weight(rtrow, period=self.period)
        scales['puweight'] = puweight
        lepscales = self.lepscaler.scale_factor(rtrow, *objects, period=self.period, **lepargs)
        scales['lep']      = lepscales[0]
        scales['lepup']    = lepscales[1]
        scales['lepdown']  = lepscales[2]
        if self.period==13:
            lepeff = self.lepeff.scale_factor(rtrow, *objects, period=self.period, **lepargs)
            scales['lepeff'] = lepeff[0]
            scales['lepeffup'] = lepeff[1]
            scales['lepeffdown'] = lepeff[2]
            lepfake = self.lepfake.scale_factor(rtrow, *objects, period=self.period, **lepargs)
            scales['lepfake'] = lepfake[0]
            scales['lepfakeup'] = lepfake[1]
            scales['lepfakedown'] = lepfake[2]
        scales['trigger_prescale'] = self.getTriggerPrescale(rtrow)
        genweight = rtrow.GenWeight if hasattr(rtrow,'GenWeight') else 1.
        scales['genweight']= genweight
        return scales

    def getTriggerPrescale(self,rtrow):
        '''Default trigger prescale'''
        return 1

    def getGenChannel(self,rtrow):
        '''Dummy return gen channel string'''
        return 'a'

    def getObject(self,rtrow,obj,var,**kwargs):
        '''Get modified object'''
        if obj=='met':
            metVar = 'type1_pfMet'
            shift = kwargs.pop('metShift','')
            shiftMap = {
              'jres': 'JetRes',
              'jes' : 'JetEn',
              'mes' : 'MuonEn',
              'ees' : 'ElectronEn',
              'tes' : 'TauEn',
              'ues' : 'UnclusteredEn',
              'pes' : 'PhotonEn',
            }
            change = {'+': 'Up', '-': 'Down'}
            varName = {'pt': 'Pt', 'phi': 'Phi'}
            varName2 = {'pt': 'Et', 'phi': 'Phi'}
            shiftString = '{0}_shifted{1}_{2}{3}'.format(metVar,varNam[var],shiftMap[shift[:-1]],change[-1]) if shift else '{0}{1}'.format(metVar,varName2[var])
            return getattr(rtrow,shiftString)
        elif obj in self.object_definitions: # composite
            masses = {'e':0.511e-3, 'm':0.1056, 't':1.776, 'j':0}
            if self.object_definitions[obj] == 1:
                return -9.
            elif self.object_definitions[obj][1] == 'n':
                o = kwargs.get('objects',[])
                if not o: return -9
                pt1 = self.getObject(rtrow,o[0],'pt',**kwargs)
                eta1 = self.getObject(rtrow,o[0],'eta',**kwargs)
                phi1 = self.getObject(rtrow,o[0],'phi',**kwargs)
                mass1 = masses[o[0][0]]
                vec1 = rt.TLorentzVector()
                vec1.SetPtEtaPhiM(pt1,eta1,phi1,mass1)
                px1 = vec1.Px()
                py1 = vec1.Py()
                ptMet = self.getObject(rtrow,'met','pt',**kwargs)
                phiMet = self.getObject(rtrow,'met','phi',**kwargs)
                vecMet = rt.TLorentzVector()
                vecMet.SetPtEtaPhiM(ptMet,0.,phiMet,0.)
                pxMet = vecMet.Px()
                pyMet = vecMet.Py
                vecW = vec1 + vecMet
                if var=='pt':
                    return vecW.Pt()
                if var=='mass' or var=='mt' or var=='m':
                    return vecW.Mt()
                if var=='dphi':
                    return deltaPhi(phi1,phiMet)
                if var=='st':
                    return pt1 + ptMet
            else: # assumes two lepton constituents
                o = kwargs.get('objects',[])
                if len(o) != 2: return 0.
                o = ordered(o[0],o[1])
                vec1 = rt.TLorentzVector()
                pt1 = self.getObject(rtrow,o[0],'pt',**kwargs)
                eta1 = self.getObject(rtrow,o[0],'eta',**kwargs)
                phi1 = self.getObject(rtrow,o[0],'phi',**kwargs)
                mass1 = masses[o[0][0]]
                vec1.SetPtEtaPhiM(pt1,eta1,phi1,mass1)
                vec2 = rt.TLorentzVector()
                pt2 = self.getObject(rtrow,o[1],'pt',**kwargs)
                eta2 = self.getObject(rtrow,o[1],'eta',**kwargs)
                phi2 = self.getObject(rtrow,o[1],'phi',**kwargs)
                mass2 = masses[o[1][0]]
                vec2.SetPtEtaPhiM(pt2,eta2,phi2,mass2)
                vec3 = vec1+vec2
                if var=='pt':
                    return vec3.Pt()
                if var=='mass' or var=='m':
                    return vec3.M()
                if var=='mt':
                    return vec3.Mt()
                if var=='dphi':
                    return deltaPhi(phi1,phi2)
                if var=='dr':
                    return deltaR(eta1,phi1,eta2,phi2)
                if var=='st':
                    return pt1 + pt2


        elif obj[0] in ['e','m','t']: # lepton
            # pt
            if var=='pt': return getattr(rtrow,'{0}Pt'.format(obj))
            # eta
            if var=='eta': return getattr(rtrow,'{0}Eta'.format(obj))
            # phi
            if var=='phi': return getattr(rtrow,'{0}Phi'.format(obj))
