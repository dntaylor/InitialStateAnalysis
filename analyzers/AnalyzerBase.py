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

#from scale_factors import LeptonScaleFactors
#from pu_weights import PileupWeights
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

    def add(self, fun):
        self.cut_sequence.append(fun)

    def evaluate(self, rtrow):
        for i,cut in enumerate(self.cut_sequence):
            if not cut(rtrow): return False, i
        return True, i+1

def lep_order(a, b):
    '''
    A simple function to guarantee order of leptons in FSA ntuples.
    '''
    if len(a)==2 and len(b)==2:
        a_index = int(a[1])
        b_index = int(b[1])
        return a_index > b_index or a[0] > b[0]
    return a[0] > b[0]

def ordered(a,b):
    '''
    Return a,b in lep order.
    '''
    return [a,b] if lep_order(b,a) else [b,a]

class AnalyzerBase(object):
    '''
    The basic analyzer class. Inheritor classes must define
        TODO
    '''
    def __init__(self, sample_location, out_file, period):
        self.sample_name = os.path.basename(sample_location)
        self.file_names = os.listdir(sample_location)
        self.out_file = out_file
        self.sample_location = sample_location
        self.period = period

    def __enter__(self):
        self.begin()
        return self

    def __exit__(self, type, value, traceback):
        self.finish()

    def begin(self):
        #self.lepscaler = LeptonScaleFactors()
        #self.pu_weights = PileupWeights()

        self.file = rt.TFile(self.out_file, 'recreate')
        
        if hasattr(self,'other_states'):
            states = [self.initial_states] + self.other_states
        else:
            states = [self.initial_states]
        if not hasattr(self,'alternateIds'): self.alternateIds = []
        if not hasattr(self,'doVBF'): self.doVBF = False
        self.ntuple, self.branches = buildNtuple(self.object_definitions,states,self.channel,self.final_states,altIds=self.alternateIds,doVBF=self.doVBF)

    def analyze(self,**kwargs):
        '''
        The primary analyzer loop.
        '''
        self.eventMap = {}
        self.bestCandMap = {}
        self.cutflowMap = {}
        eventsToWrite = set()
        eventsWritten = set()
        numEvts = 0
        totalWritten = 0

        # iterate over files
        for i, file_name in enumerate(self.file_names):
            print "%s %s %s: Processing %i/%i files" % (str(datetime.datetime.now()), self.channel, self.sample_name, i+1, len(self.file_names))

            file_path = os.path.join(self.sample_location, file_name)
            rtFile = rt.TFile(file_path, "READ")

            # iterate over final states
            for fs in self.final_states:
                if len(self.file_names)<10: print "%s %s: %s" % (self.channel, self.sample_name, fs)
                tree = rtFile.Get("%s/final/Ntuple" % fs)
                metatree = rtFile.Get("%s/eventCount" % fs)
                tempEvts = metatree.GetEntries()

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
                        if len(self.file_names)==1: print "%s %s %s: %s %i/%i entries" % (str(datetime.datetime.now()), self.channel, self.sample_name, fs, numFSEvents, totalFSEvents)
                    numFSEvents += 1

                    # cache to prevent excessive reads of fsa ntuple
                    self.cache = {}

                    # now see if event is viable
                    passPreselection = self.pass_preselection(rtrow)
                    eventkey = (rtrow.evt, rtrow.lumi, rtrow.run)
                    old = self.cutflowMap[eventkey] if eventkey in self.cutflowMap else -1
                    if self.num>old:
                        self.cutflowMap[eventkey] = self.num
                    if not passPreselection:
                        continue

                    # can we define the object we want?
                    candidate = self.choose_objects(rtrow)
                    if not candidate: continue # in case no objects satisfy our requirements
                    if self.num+1>old:
                        self.cutflowMap[eventkey] = self.num+1

                    # check combinatorics
                    if eventkey in self.bestCandMap: 
                        bestcand = self.bestCandMap[eventkey]
                    else:
                        numMin = len(candidate[0])
                        bestcand = [float('inf')] * numMin
                    if self.good_to_store(rtrow,candidate[0],bestcand):
                        self.bestCandMap[eventkey] = candidate[0]
                        if self.num+2>old:
                            self.cutflowMap[eventkey] = self.num+2
                        ntupleRow = self.store_row(rtrow, *candidate[1])
                        self.eventMap[eventkey] = ntupleRow
                        eventsToWrite.add(eventkey)

            rtFile.Close("R")
            numEvts += tempEvts

            # end of file, write the ntuples
            self.file.cd()
            for key in eventsToWrite:
                if key in eventsWritten:
                    print "%s %s %s: Error: attempted to write previously written event" % (str(datetime.datetime.now()), self.channel, self.sample_name)
                else:
                    self.write_row(self.eventMap[key])
                    self.ntuple.Fill()
            eventsWritten.update(eventsToWrite)
            self.eventMap = {}
            eventsToWrite = set()

        # now we store all events that are kept
        print "%s %s %s: Filling Tree" % (str(datetime.datetime.now()), self.channel, self.sample_name)
        #self.file.cd()
        #for key in eventsToWrite:
        #    self.write_row(self.eventMap[key])
        #    self.ntuple.Fill()
        print "%s %s %s: Filled Tree (%i events)" % (str(datetime.datetime.now()), self.channel, self.sample_name, len(eventsWritten))

        # now we store the total processed events
        print "%s %s %s: Processed %i events" % (str(datetime.datetime.now()), self.channel, self.sample_name, numEvts)

        # and the cutflow
        cutflowVals = []
        for val in self.cutflowMap.itervalues():
            for i in range(val+1):
                if len(cutflowVals)<i+1: cutflowVals.append(1)
                else: cutflowVals[i] += 1
        print "%s %s %s: Cutflow: " % (str(datetime.datetime.now()), self.channel, self.sample_name), cutflowVals

        cutflowHist = rt.TH1F('cutflow','cutflow',len(cutflowVals)+1,0,len(cutflowVals)+1)
        cutflowHist.SetBinContent(1,numEvts)
        for i in range(len(cutflowVals)):
            cutflowHist.SetBinContent(i+2,cutflowVals[i])
        # rename cutflow bins if self.cutflow_labels defined
        if hasattr(self,'cutflow_labels'):
            pass # TODO
        cutflowHist.Write()

    def finish(self):
        self.file.Write()
        self.file.Close()

    @staticmethod
    def enumerate_objects(final_state):
        '''Get the objects available in the ntuple for a given final state'''
        out = []
        for i  in ['e', 'm', 't', 'j', 'g']:
            N = final_state.count(i)
            if N==1:
               out += i
            else:
               out += ['%s%i' % (i, n) for n in xrange(1, N+1)]
        return out

    def choose_alternative_objects(self, rtorw, state):
        '''
        Dummy method for alternative object selection.
        '''
        return []

    @staticmethod
    def good_to_store(rtrow, cand1, cand2):
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
        numObjs = len(self.final_states[0])
        finalStateObjects ='emtjgn'
        allowedObjects = ''
        for fsObj in finalStateObjects:
            for fs in self.final_states:
                if fsObj in fs:
                    allowedObjects += fsObj
                    break
        numObjTypes = len(allowedObjects)
        for prompts in list(product(range(numObjs+1),repeat=numObjTypes)):
            if sum(prompts) > numObjs: continue
            promptString = ''.join([str(x) for x in prompts])
            promptDict = {}
            for o in range(len(allowedObjects)):
                promptDict[allowedObjects[o]] = prompts[o]
            ntupleRow["select.pass_%s"%promptString] = int(self.npass(rtrow,promptDict,**self.getIdArgs('Tight')))
        for altId in self.alternateIds:
            ntupleRow["select.pass_%s"%altId] = int(self.ID(rtrow,*self.objects,**self.alternateIdMap[altId]))

        
        ntupleRow["event.evt"] = int(rtrow.evt)
        ntupleRow["event.lumi"] = int(rtrow.lumi)
        ntupleRow["event.run"] = int(rtrow.run)
        ntupleRow["event.nvtx"] = int(rtrow.nvtx)
        #ntupleRow["event.lep_scale"] = float(self.lepscaler.scale_factor(rtrow, *objects, tight=True))
        #ntupleRow["event.pu_weight"] = float(self.pu_weights.weight(rtrow))
        ntupleRow["event.lep_scale"] = float(1.)
        ntupleRow["event.pu_weight"] = float(1.)

        channelString = ''
        for x in objects: channelString += x[0]
        ntupleRow["channel.channel"] = channelString

        ntupleRow["finalstate.mass"] = float(rtrow.Mass)
        ntupleRow["finalstate.sT"] = float(sum([getattr(rtrow, "%sPt" % x) for x in objects]))
        ntupleRow["finalstate.met"] = float(rtrow.pfMetEt)
        ntupleRow["finalstate.metPhi"] = float(rtrow.pfMetPhi)
        ntupleRow["finalstate.jetVeto20"] = int(rtrow.jetVeto20)
        ntupleRow["finalstate.jetVeto30"] = int(rtrow.jetVeto30)
        ntupleRow["finalstate.jetVeto40"] = int(rtrow.jetVeto40)
        ntupleRow["finalstate.bjetVeto20Loose"] = int(rtrow.bjetCISVVeto20Loose)
        ntupleRow["finalstate.bjetVeto30Loose"] = int(rtrow.bjetCISVVeto30Loose)
        ntupleRow["finalstate.bjetVeto20Medium"] = int(rtrow.bjetCISVVeto20Medium)
        ntupleRow["finalstate.bjetVeto30Medium"] = int(rtrow.bjetCISVVeto30Medium)
        ntupleRow["finalstate.bjetVeto20Tight"] = int(rtrow.bjetCISVVeto20Tight)
        ntupleRow["finalstate.bjetVeto30Tight"] = int(rtrow.bjetCISVVeto30Tight)
        ntupleRow["finalstate.muonVeto5"] = int(rtrow.muVetoPt5IsoIdVtx)
        ntupleRow["finalstate.muonVeto10Loose"] = int(rtrow.muGlbIsoVetoPt10)
        ntupleRow["finalstate.muonVeto15"] = int(rtrow.muVetoPt15IsoIdVtx)
        ntupleRow["finalstate.elecVeto10"] = int(rtrow.eVetoMVAIsoVtx)
        if self.doVBF:
            ntupleRow["finalstate.vbfMass"] = float(rtrow.vbfMass)
            ntupleRow["finalstate.vbfPt"] = float(rtrow.vbfdijetpt)
            ntupleRow["finalstate.vbfPt1"] = float(rtrow.vbfj1pt)
            ntupleRow["finalstate.vbfPt2"] = float(rtrow.vbfj2pt)
            ntupleRow["finalstate.vbfEta1"] = float(rtrow.vbfj1eta)
            ntupleRow["finalstate.vbfEta2"] = float(rtrow.vbfj2eta)
            ntupleRow["finalstate.centralJetVeto20"] = float(rtrow.vbfJetVeto20)
            ntupleRow["finalstate.centralJetVeto30"] = float(rtrow.vbfJetVeto30)

        def store_state(rtrow,ntupleRow,state,theObjects):
            objStart = 0
            for i in state:
                numObjects = len([ x for x in self.object_definitions[i] if x != 'n']) if theObjects else 0
                finalObjects = theObjects[objStart:objStart+numObjects]
                orderedFinalObjects = sorted(finalObjects, key = lambda x: getattr(rtrow,"%sPt" % x))
                if len(self.object_definitions[i]) == 1:
                    ntupleRow["%s.mass" %i] = float(-9)
                    ntupleRow["%s.sT" %i] = float(getattr(rtrow, "%sPt" % finalObjects[0])) if theObjects else float(-9)
                    ntupleRow["%s.dPhi" %i] = float(-9)
                    ntupleRow["%s.dR" %i] = float(-9)
                    ntupleRow["%sFlv.Flv" %i] = finalObjects[0][0] if theObjects else 'a'
                elif 'n' == self.object_definitions[i][1]:
                    ntupleRow["%s.mass" %i] = float(getattr(rtrow, "%sMtToPFMET" % finalObjects[0])) if theObjects else float(-9)
                    ntupleRow["%s.sT" %i] = float(getattr(rtrow, "%sPt" % finalObjects[0]) + rtrow.pfMetEt) if theObjects else float(-9)
                    ntupleRow["%s.dPhi" %i] = float(getattr(rtrow, "%sToMETDPhi" % finalObjects[0])) if theObjects else float(-9)
                    ntupleRow["%sFlv.Flv" %i] = finalObjects[0][0] if theObjects else 'a'
                else:
                    finalObjOrdered = ordered(finalObjects[0], finalObjects[1])
                    ntupleRow["%s.mass" %i] = float(getattr(rtrow, "%s_%s_Mass" % (finalObjOrdered[0], finalObjOrdered[1]))) if theObjects else float(-9)
                    ntupleRow["%s.sT" %i]   = float(sum([getattr(rtrow, "%sPt" % x) for x in finalObjects])) if theObjects else float(-9)
                    ntupleRow["%s.dPhi" %i] = float(getattr(rtrow, "%s_%s_DPhi" % (finalObjOrdered[0], finalObjOrdered[1]))) if theObjects else float(-9)
                    ntupleRow["%s.dR" %i] = float(getattr(rtrow, "%s_%s_DR" % (finalObjOrdered[0], finalObjOrdered[1]))) if theObjects else float(-9)
                    ntupleRow["%sFlv.Flv" %i] = finalObjects[0][0] + finalObjects[1][0] if theObjects else 'aa'
                objCount = 0
                for obj in self.object_definitions[i]:
                    if obj=='n':
                        ntupleRow["%s.met" %i] = float(rtrow.pfMetEt) if theObjects else float(-9)
                        ntupleRow["%s.metPhi" %i] = float(rtrow.pfMetPhi) if theObjects else float(-9)
                    else:
                        objCount += 1
                        ntupleRow["%s.Pt%i" % (i,objCount)] = float(getattr(rtrow, "%sPt" % orderedFinalObjects[objCount-1])) if theObjects else float(-9)
                        ntupleRow["%s.Eta%i" % (i,objCount)] = float(getattr(rtrow, "%sEta" % orderedFinalObjects[objCount-1])) if theObjects else float(-9)
                        ntupleRow["%s.Phi%i" % (i,objCount)] = float(getattr(rtrow, "%sPhi" % orderedFinalObjects[objCount-1])) if theObjects else float(-9)
                        if orderedFinalObjects[objCount-1][0]=='e': isoVar = 'RelPFIsoRho'
                        if orderedFinalObjects[objCount-1][0]=='m': isoVar = 'RelPFIsoDBDefault'
                        isoVal = float(getattr(rtrow, "%s%s" % (orderedFinalObjects[objCount-1], isoVar))) if orderedFinalObjects[objCount-1][0] in 'em' and theObjects else float(-9.)
                        ntupleRow["%s.Iso%i" % (i,objCount)] = isoVal
                        ntupleRow["%s.Chg%i" % (i,objCount)] = float(getattr(rtrow, "%sCharge" % orderedFinalObjects[objCount-1])) if theObjects else float(-9)
                        ntupleRow["%s.PassTight%i" % (i,objCount)] = float(self.ID(rtrow,orderedFinalObjects[objCount-1],**self.getIdArgs('Tight'))) if theObjects else float(-9)
                        # manually add w z deltaRs
                        if i=='w1':
                            oZ1 = ordered(theObjects[0],theObjects[2])
                            oZ2 = ordered(theObjects[1],theObjects[2])
                            ntupleRow["w1.dR1_z1_1"] = float(getattr(rtrow,"%s_%s_DR" % (oZ1[0],oZ1[1])))
                            ntupleRow["w1.dR1_z1_2"] = float(getattr(rtrow,"%s_%s_DR" % (oZ2[0],oZ2[1])))
                            ntupleRow["w1.mll_z1_1"] = float(getattr(rtrow,"%s_%s_Mass" % (oZ1[0],oZ1[1])))
                            ntupleRow["w1.mll_z1_2"] = float(getattr(rtrow,"%s_%s_Mass" % (oZ2[0],oZ2[1])))
                        # do alternate IDs
                        for altId in self.alternateIds:
                            ntupleRow["%s.pass_%s_%i"%(i,altId,objCount)] = int(self.ID(rtrow,orderedFinalObjects[objCount-1],**self.alternateIdMap[altId]) if theObjects else float(-9))
                objStart += numObjects


        # initial state objects
        store_state(rtrow,ntupleRow,self.initial_states,objects)

        # alternative state objects
        if hasattr(self,'other_states'):
            for state in self.other_states:
                store_state(rtrow,ntupleRow,state,self.choose_alternative_objects(rtrow,state))
                

        # final state objects
        lepCount = 0
        jetCount = 0
        phoCount = 0
        orderedAllObjects = sorted(objects, key = lambda x: getattr(rtrow,"%sPt" % x))
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
            ntupleRow["%s%i.Chg" % (charName,objCount)] = float(getattr(rtrow, "%sCharge" % obj))
            ntupleRow["%s%i.PassTight" % (charName,objCount)] = float(self.ID(rtrow,obj,**self.getIdArgs('Tight')))
            ntupleRow["%s%iFlv.Flv" % (charName,objCount)] = obj[0]

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
        cutResults,self.num = cuts.evaluate(rtrow)
        self.cache['preselection'] = cutResults
        return cutResults

    def pass_selection(self,rtrow):
        '''
        Wrapper for the selection defined by the user (tight selection whereas preselection
        is the loose selection for fake rate method).
        '''
        if 'selection' in self.cache: return self.cache['selection']
        cuts = self.selection(rtrow)
        cutResults,self.num = cuts.evaluate(rtrow)
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
