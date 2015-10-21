'''
Utilities for building ntuples used in ISA.

Author: Devin N. Taylor, UW-Madison
'''
from itertools import product

import ROOT as rt
from array import array

def buildCutTree(cutlabels,**kwargs):
    eventLeafs = ['evt','run','lumi']
    eventBranchLineToProcess = "struct eventBranch_t {" + " ".join(["Int_t {0};".format(x) for x in eventLeafs]) + "}"
    cutBranchLineToProcess = "struct cutBranch_t {" + " ".join(["Int_t {0};".format(x) for x in cutlabels]) + "}"
    eventBranchStrForBranch = '{0}/I:'.format(eventLeafs[0]) + ':'.join(eventLeafs[1:])
    cutBranchStrForBranch = '{0}/I:'.format(cutlabels[0]) + ':'.join(cutlabels[1:])
    if not hasattr(rt,"eventBranch_t"): rt.gROOT.ProcessLine(eventBranchLineToProcess)
    if not hasattr(rt,"cutBranch_t"): rt.gROOT.ProcessLine(cutBranchLineToProcess)
    eventBranchStruct = rt.eventBranch_t()
    cutBranchStruct = rt.cutBranch_t()
    tree = rt.TTree('cutTree','cutTree')
    tree.Branch('event',eventBranchStruct,eventBranchStrForBranch)
    tree.Branch('selections',cutBranchStruct,cutBranchStrForBranch)
    return (tree, eventBranchStruct, cutBranchStruct)


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
    for trig in ['passDoubleMuon','passDoubleEG','passMuonEG','passEGMuon']:
        strToProcess += "Int_t {0};".format(trig)
        strForBranch += "{0}:".format(trig)
    #for prompts in list(product(range(numObjs+1),repeat=numObjTypes)):
    #    if sum(prompts) > numObjs: continue
    #    promptString = ''.join([str(x) for x in prompts])
    #    strToProcess += "Int_t pass_%s;" % promptString
    #    strForBranch += "pass_%s:" % promptString
    for altId in alternateIds:
        strToProcess += "Int_t pass_%s;" % altId
        strForBranch += "pass_%s:" % altId
    strToProcess += "};"
    strForBranch = strForBranch[:-1] # remove trailing :
    if not hasattr(rt,'structSelect_t'): rt.gROOT.ProcessLine(strToProcess)
    selectStruct = rt.structSelect_t()
    structureDict['select'] = [selectStruct, selectStruct, strForBranch]
    structOrder += ['select']    

    # define common root classes
    if not hasattr(rt,'structEvent_t'):
        rt.gROOT.ProcessLine(
        "struct structEvent_t {\
           Int_t   evt;\
           Int_t   run;\
           Int_t   lumi;\
           Int_t   nvtx;\
           Int_t   GenNUP;\
           Int_t   trig_prescale;\
           Float_t lep_scale;\
           Float_t lep_scale_up;\
           Float_t lep_scale_down;\
           Float_t trig_scale;\
           Float_t pu_weight;\
           Float_t pu_weight_up;\
           Float_t pu_weight_down;\
           Float_t gen_weight;\
           Float_t charge_uncertainty;\
        };");
    eventStruct = rt.structEvent_t()
    structureDict['event'] = [eventStruct, eventStruct,'evt/I:run:lumi:nvtx:GenNUP:trig_prescale:lep_scale/F:lep_scale_up:lep_scale_down:trig_scale:pu_weight:pu_weight_up:pu_weight_down:gen_weight:charge_uncertainty']
    structOrder += ['event']

    if not hasattr(rt,'structChannel_t'):
        rt.gROOT.ProcessLine(
        "struct structChannel_t {\
           Char_t  channel[9];\
        };");
    channelStruct = rt.structChannel_t()
    structureDict['channel'] = [channelStruct, rt.AddressOf(channelStruct,'channel'),'channel/C']
    structOrder += ['channel']
    genChannelStruct = rt.structChannel_t()
    structureDict['genChannel'] = [genChannelStruct, rt.AddressOf(genChannelStruct,'channel'),'channel/C']
    structOrder += ['genChannel']

    fsStrToProcess = "struct structFinalState_t {\
       Float_t mass;\
       Float_t mT;\
       Float_t sT;\
       Float_t hT;\
       Float_t met;\
       Float_t metPhi;\
       Float_t leadJetPt;\
       Float_t leadJetEta;\
       Float_t leadJetPhi;\
       Float_t leadJetPUMVA;"
    fsStrForBranch = "mass/F:mT:sT:hT:met:metPhi:leadJetPt:leadJetEta:leadJetPhi:leadJetPUMVA:"

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
       Int_t   bjetVeto20Loose;\
       Int_t   bjetVeto30Loose;\
       Int_t   bjetVeto20Medium;\
       Int_t   bjetVeto30Medium;\
       Int_t   bjetVeto20Tight;\
       Int_t   bjetVeto30Tight;\
       Int_t   muonVetoTight;\
       Int_t   elecVetoTight;\
       Int_t   muonVetoLoose;\
       Int_t   elecVetoLoose;"
    fsStrForBranch += "jetVeto20/I:jetVeto30:jetVeto40:bjetVeto20Loose:bjetVeto30Loose:bjetVeto20Medium:bjetVeto30Medium:bjetVeto20Tight:bjetVeto30Tight:muonVetoTight:elecVetoTight:muonVetoLoose:elecVetoLoose:"

    if doVBF:
        fsStrToProcess += "Int_t   centralJetVeto20;\
                           Int_t   centralJetVeto30;"
        fsStrForBranch += "centralJetVeto20:centralJetVeto30:"

    fsStrToProcess += "};"
    fsStrForBranch = fsStrForBranch[:-1]
    if not hasattr(rt,'structFinalState_t'):
        rt.gROOT.ProcessLine(fsStrToProcess);
    finalStateStruct = rt.structFinalState_t()
    structureDict['finalstate'] = [finalStateStruct, finalStateStruct, fsStrForBranch]
    structOrder += ['finalstate']

    if not hasattr(rt,'structObject_t'):
        rt.gROOT.ProcessLine(
        "struct structObject_t {\
           Float_t Pt;\
           Float_t Eta;\
           Float_t Phi;\
           Float_t Iso;\
           Float_t Dxy;\
           Float_t Dz;\
           Float_t SigmaIEtaIEta;\
           Float_t DEtaIn;\
           Float_t DPhiIn;\
           Float_t HOverE;\
           Float_t OoEmOoP;\
           Float_t TriggeringMVA;\
           Float_t NonTriggeringMVA;\
           Float_t NormalizedChi2;\
           Float_t JetPt;\
           Float_t JetBTag;\
           Float_t LepScaleLoose;\
           Float_t LepScaleTight;\
           Float_t LepScaleLoose_up;\
           Float_t LepScaleTight_up;\
           Float_t LepScaleLoose_down;\
           Float_t LepScaleTight_down;\
           Float_t LepEffLoose;\
           Float_t LepEffTight;\
           Float_t LepFake;\
           Int_t   Chg;\
           Int_t   PassLoose;\
           Int_t   PassTight;\
           Int_t   GenPdgId;\
           Int_t   MotherGenPdgId;\
           Int_t   ChargeConsistent;\
           Int_t   ExpectedMissingInnerHits;\
           Int_t   PassConversionVeto;\
           Int_t   IsGlobalMuon;\
           Int_t   IsPFMuon;\
           Int_t   IsTrackerMuon;\
           Int_t   ValidMuonHits;\
           Int_t   MatchedStations;\
           Int_t   ValidPixelHits;\
           Int_t   TrackerLayers;\
        };");
    if not hasattr(rt,'structObjChar_t'):
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
                structureDict['%s%i' % (charName, objCount)] = [objStruct, objStruct, 'Pt/F:Eta:Phi:Iso:Dxy:Dz:SigmaIEtaIEta:DEtaIn:DPhiIn:HOverE:OoEmOoP:TriggeringMVA:NonTriggeringMVA:NormalizedChi2:JetPt:JetBTag:LepScaleLoose:LepScaleTight:LepScaleLoose_up:LepScaleTight_up:LepScaleLoose_down:LepScaleTight_down:LepEffLoose:LepEffTight:LepFake:Chg/I:PassLoose:PassTight:GenPdgId:MotherGenPdgId:ChargeConsistent:ExpectedMissingInnerHits:PassConversionVeto:IsGlobalMuon:IsPFMuon:IsTrackerMuon:ValidMuonHits:MatchedStations:ValidPixelHits:TrackerLayers']
                structureDict['%s%iFlv' % (charName, objCount)] = [flvStruct, rt.AddressOf(flvStruct,'Flv'),'Flv/C']
                structOrder += ['%s%i' % (charName, objCount)]
                structOrder += ['%s%iFlv' % (charName, objCount)]

    # define objects for each initial state
    for state in states:
        for key in state:
            val = object_definitions[key]
            strForBranch = ""
            strToProcess = "struct struct%s_t {" % (key.upper())
            strForBranch += "mass/F:Pt:sT:dPhi:"
            strToProcess += "\
                Float_t mass;\
                Float_t Pt;\
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
                    strForBranch += "Pt{0}:Eta{0}:Phi{0}:Iso{0}:Dxy{0}:Dz{0}:SigmaIEtaIEta{0}:DEtaIn{0}:DPhiIn{0}:HOverE{0}:OoEmOoP{0}:TriggeringMVA{0}:NonTriggeringMVA{0}:NormalizedChi2{0}:JetPt{0}:JetBTag{0}:LepScaleLoose{0}:LepScaleTight{0}:LepScaleLoose{0}_up:LepScaleTight{0}_up:LepScaleLoose{0}_down:LepScaleTight{0}_down:LepEffLoose{0}:LepEffTight{0}:LepFake{0}:".format(objCount)
                    strToProcess += "\
                        Float_t Pt{0};\
                        Float_t Eta{0};\
                        Float_t Phi{0};\
                        Float_t Iso{0};\
                        Float_t Dxy{0};\
                        Float_t Dz{0};\
                        Float_t SigmaIEtaIEta{0};\
                        Float_t DEtaIn{0};\
                        Float_t DPhiIn{0};\
                        Float_t HOverE{0};\
                        Float_t OoEmOoP{0};\
                        Float_t TriggeringMVA{0};\
                        Float_t NonTriggeringMVA{0};\
                        Float_t NormalizedChi2{0};\
                        Float_t JetPt{0};\
                        Float_t JetBTag{0};\
                        Float_t LepScaleLoose{0};\
                        Float_t LepScaleTight{0};\
                        Float_t LepScaleLoose{0}_up;\
                        Float_t LepScaleTight{0}_up;\
                        Float_t LepScaleLoose{0}_down;\
                        Float_t LepScaleTight{0}_down;\
                        Float_t LepEffLoose{0};\
                        Float_t LepEffTight{0};\
                        Float_t LepFake{0};".format(objCount)
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
                        strForBranch += "dR1_z1_1:dR1_z1_2:mll_z1_1:mll_z1_2:dR1_leadJet:"
                        strToProcess += "Float_t dR1_z1_1; Float_t dR1_z1_2; Float_t mll_z1_1; Float_t mll_z1_2; Float_t dR1_leadJet;"
            # do the alt IDs
            objCount = 0
            for obj in val:
                if obj == 'n': continue
                else:
                    objCount += 1
                    strForBranch += "Chg{0}/I:PassLoose{0}:PassTight{0}:GenPdgId{0}:MotherGenPdgId{0}:ChargeConsistent{0}:ExpectedMissingInnerHits{0}:PassConversionVeto{0}:IsGlobalMuon{0}:IsPFMuon{0}:IsTrackerMuon{0}:ValidMuonHits{0}:MatchedStations{0}:ValidPixelHits{0}:TrackerLayers{0}:".format(objCount)
                    strToProcess += "\
                        Int_t   Chg{0};\
                        Int_t   PassLoose{0};\
                        Int_t   PassTight{0};\
                        Int_t   GenPdgId{0};\
                        Int_t   MotherGenPdgId{0};\
                        Int_t   ChargeConsistent{0};\
                        Int_t   ExpectedMissingInnerHits{0};\
                        Int_t   PassConversionVeto{0};\
                        Int_t   IsGlobalMuon{0};\
                        Int_t   IsPFMuon{0};\
                        Int_t   IsTrackerMuon{0};\
                        Int_t   ValidMuonHits{0};\
                        Int_t   MatchedStations{0};\
                        Int_t   ValidPixelHits{0};\
                        Int_t   TrackerLayers{0};".format(objCount)
                    for altId in alternateIds:
                        strToProcess += "Int_t pass_{0}_{1};".format(altId, objCount)
                        strForBranch += "pass_{0}_{1}:".format(altId, objCount)
            strForBranch = strForBranch[:-1] # remove trailing :
            strToProcess += "\
                };"
            if not hasattr(rt,'struct%s_t' % (key.upper())):
                rt.gROOT.ProcessLine(strToProcess)
            initialStruct = getattr(rt,"struct{0}_t".format(key.upper()))()
            structureDict[key] = [initialStruct, initialStruct, strForBranch]
            structOrder += [key]

    if not hasattr(rt,'structInitialChar_t'):
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

    

