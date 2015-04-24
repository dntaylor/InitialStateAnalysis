'''
Utilities for building ntuples used in ISA.

Author: Devin N. Taylor, UW-Madison
'''
from itertools import product

import ROOT as rt
from array import array

def buildNtuple(object_definitions,states,channelName,final_states,**kwargs):
    '''
    A function to build an initial state ntuple for AnalyzerBase.py
    '''
    alternateIds = kwargs.pop('altIds',[])
    doVBF = kwargs.pop('doVBF',False)

    finalStateObjects = 'emtjgn'
    structureDict = {}
    structOrder = []

    # define selection bools and fake rate things
    numObjs = len(final_states[0])
    allowedObjects = ''
    for fsObj in finalStateObjects:
        for fs in final_states:
            if fsObj in fs:
                allowedObjects += fsObj
                break
    numObjTypes = len(allowedObjects)
    strToProcess = "struct structSelect_t {"
    strForBranch = ""
    strToProcess += "Int_t passTight;"
    strForBranch += "passTight/I:"
    strToProcess += "Int_t passLoose;"
    strForBranch += "passLoose:"
    for prompts in list(product(range(numObjs+1),repeat=numObjTypes)):
        if sum(prompts) > numObjs: continue
        promptString = ''.join([str(x) for x in prompts])
        strToProcess += "Int_t pass_%s;" % promptString
        strForBranch += "pass_%s:" % promptString
    for altId in alternateIds:
        strToProcess += "Int_t pass_%s;" % altId
        strForBranch += "pass_%s:" % altId
    strToProcess += "};"
    strForBranch = strForBranch[:-1] # remove trailing :
    rt.gROOT.ProcessLine(strToProcess)
    selectStruct = rt.structSelect_t()
    structureDict['select'] = [selectStruct, selectStruct, strForBranch]
    structOrder += ['select']    

    # define common root classes
    rt.gROOT.ProcessLine(
    "struct structEvent_t {\
       Int_t   evt;\
       Int_t   run;\
       Int_t   lumi;\
       Int_t   nvtx;\
       Float_t lep_scale;\
       Float_t pu_weight;\
    };");
    eventStruct = rt.structEvent_t()
    structureDict['event'] = [eventStruct, eventStruct,'evt/I:run:lumi:nvtx:lep_scale/F:pu_weight']
    structOrder += ['event']

    rt.gROOT.ProcessLine(
    "struct structChannel_t {\
       Char_t  channel[9];\
    };");
    channelStruct = rt.structChannel_t()
    structureDict['channel'] = [channelStruct, rt.AddressOf(channelStruct,'channel'),'channel/C']
    structOrder += ['channel']

    fsStrToProcess = "struct structFinalState_t {\
       Float_t mass;\
       Float_t sT;\
       Float_t met;\
       Float_t metPhi;"
    fsStrForBranch = "mass/F:sT:met:metPhi:"

    if doVBF:
        fsStrToProcess += "Float_t vbfMass;\
                           Float_t vbfPt;\
                           Float_t vbfPt1;\
                           Float_t vbfPt2;\
                           Float_t vbfEta1;\
                           Float_t vbfEta2;"
        fsStrForBranch += "vbfMass:vbfPt:vbfPt1:vbfPt2:vbfEta1:vbfEta2:"

    fsStrToProcess += "Int_t   jetVeto20;\
       Int_t   jetVeto30;\
       Int_t   jetVeto40;\
       Int_t   bjetVeto20;\
       Int_t   bjetVeto30;\
       Int_t   muonVeto5;\
       Int_t   muonVeto10Loose;\
       Int_t   muonVeto15;\
       Int_t   elecVeto10;"
    fsStrForBranch += "jetVeto20/I:jetVeto30:jetVeto40:bjetVeto20:bjetVeto30:muonVeto5:muonVeto10Loose:muonVeto15:elecVeto10:"

    if doVBF:
        fsStrToProcess += "Int_t   centralJetVeto20;\
                           Int_t   centralJetVeto30;"
        fsStrForBranch += "centralJetVeto20:centralJetVeto30:"

    fsStrToProcess += "};"
    fsStrForBranch = fsStrForBranch[:-1]
    rt.gROOT.ProcessLine(fsStrToProcess);
    finalStateStruct = rt.structFinalState_t()
    structureDict['finalstate'] = [finalStateStruct, finalStateStruct, fsStrForBranch]
    structOrder += ['finalstate']

    rt.gROOT.ProcessLine(
    "struct structObject_t {\
       Float_t Pt;\
       Float_t Eta;\
       Float_t Phi;\
       Float_t Iso;\
       Int_t   Chg;\
       Int_t   PassTight;\
    };");
    rt.gROOT.ProcessLine(
    "struct structObjChar_t {\
       Char_t  Flv[2];\
    };");
    lepCount = 0
    jetCount = 0
    phoCount = 0
    for key in states[0]:
        val = object_definitions[key]
        for obj in val:
            if obj=='n': continue
            else:
                objStruct = rt.structObject_t()
                flvStruct = rt.structObjChar_t()
                if obj in 'emt': 
                    charName = 'l'
                    lepCount += 1
                    objCount = lepCount
                if obj == 'j':
                    charName = 'j'
                    jetCount += 1
                    objCount = jetCount
                if obj == 'g':
                    charName = 'g'
                    phoCount += 1
                    objCount = phoCount
                structureDict['%s%i' % (charName, objCount)] = [objStruct, objStruct, 'Pt/F:Eta:Phi:Iso:Chg/I:PassTight']
                structureDict['%s%iFlv' % (charName, objCount)] = [flvStruct, rt.AddressOf(flvStruct,'Flv'),'Flv/C']
                structOrder += ['%s%i' % (charName, objCount)]
                structOrder += ['%s%iFlv' % (charName, objCount)]

    # define objects for each initial state
    for state in states:
        for key in state:
            val = object_definitions[key]
            strForBranch = ""
            strToProcess = "struct struct%s_t {" % key.upper()
            strForBranch += "mass/F:sT:dPhi:"
            strToProcess += "\
                Float_t mass;\
                Float_t sT;\
                Float_t dPhi;"
            if 'n' not in val:
                strToProcess += "Float_t dR;"
                strForBranch += "dR:"
            objCount = 0
            for obj in val:
                if obj == 'n':
                    strForBranch += "met:metPhi:"
                    strToProcess += "\
                        Float_t met;\
                        Float_t metPhi;"
                else:
                    objCount += 1
                    strForBranch += "Pt%i:Eta%i:Phi%i:Iso%i:" % (objCount, objCount, objCount, objCount)
                    strToProcess += "\
                        Float_t Pt%i;\
                        Float_t Eta%i;\
                        Float_t Phi%i;\
                        Float_t Iso%i;" % (objCount, objCount, objCount, objCount)
                    # do the deltaRs
                    #for s in states:
                    #    for k in state:
                    #        if k==key and s==state: continue # dont dR the same object
                    #        v = object_definitions[k]
                    #        oCount = 0
                    #        for o in v:
                    #            if o == 'n': continue
                    #            else:
                    #                oCount += 1
                    #                strForBranch += "dR%i_%s_%i:" % (objCount,k,oCount)
                    #                strToProcess += "Float_t dR%i_%s_%i;" % (objCount,k,oCount)
                    # manually add the W Z deltaRs for now
                    if key == 'w1':
                        strForBranch += "dR1_z1_1:dR1_z1_2:mll_z1_1:mll_z1_2:"
                        strToProcess += "Float_t dR1_z1_1; Float_t dR1_z1_2; Float_t mll_z1_1; Float_t mll_z1_2;"
            # do the alt IDs
            objCount = 0
            for obj in val:
                if obj == 'n': continue
                else:
                    objCount += 1
                    strForBranch += "Chg%i/I:PassTight%i:" % (objCount, objCount) if objCount == 1 else\
                                    "Chg%i:PassTight%i:" % (objCount, objCount)
                    strToProcess += "\
                        Int_t   Chg%i;\
                        Int_t   PassTight%i;" % (objCount, objCount)
                    for altId in alternateIds:
                        strToProcess += "Int_t pass_%s_%i;" % (altId, objCount)
                        strForBranch += "pass_%s_%i:" % (altId, objCount)
            strForBranch = strForBranch[:-1] # remove trailing :
            strToProcess += "\
                };"
            rt.gROOT.ProcessLine(strToProcess)
            initialStruct = getattr(rt,"struct%s_t" % key.upper())()
            structureDict[key] = [initialStruct, initialStruct, strForBranch]
            structOrder += [key]

    rt.gROOT.ProcessLine(
    "struct structInitialChar_t {\
       Char_t  Flv[3];\
    };");
    for key in object_definitions:
        initialFlvStruct = rt.structInitialChar_t()
        structureDict['%sFlv' % key] = [initialFlvStruct,rt.AddressOf(initialFlvStruct,'Flv'),'Flv/C']
        structOrder += ['%sFlv' % key]

    # now create the tree
    tree = rt.TTree(channelName,channelName)
    allBranches = {}
    for key in structOrder:
        val = structureDict[key]
        tree.Branch(key,val[1],val[2])
        allBranches[key] = val[0]

    return (tree, allBranches)

    

