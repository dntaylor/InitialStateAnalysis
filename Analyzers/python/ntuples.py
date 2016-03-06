'''
Utilities for building ntuples used in ISA.

Author: Devin N. Taylor, UW-Madison
'''
from itertools import product

import ROOT as rt
from array import array

def buildCutTree(cutlabels,**kwargs):
    eventLeafs = ['evt','run','lumi']
    eventBranchLineToProcess = "struct eventBranch_t {ULong64_t evt; Int_t run; Int_t lumi;}"
    cutBranchLineToProcess = "struct cutBranch_t {" + " ".join(["Int_t {0};".format(x) for x in cutlabels]) + "}"
    eventBranchStrForBranch = 'evt/l:run/I:lumi'
    cutBranchStrForBranch = '{0}/I:'.format(cutlabels[0]) + ':'.join(cutlabels[1:])
    if not hasattr(rt,"eventBranch_t"): rt.gROOT.ProcessLine(eventBranchLineToProcess)
    if not hasattr(rt,"cutBranch_t"): rt.gROOT.ProcessLine(cutBranchLineToProcess)
    eventBranchStruct = rt.eventBranch_t()
    cutBranchStruct = rt.cutBranch_t()
    tree = rt.TTree('cutTree','cutTree')
    tree.Branch('event',eventBranchStruct,eventBranchStrForBranch)
    tree.Branch('selections',cutBranchStruct,cutBranchStrForBranch)
    return (tree, eventBranchStruct, cutBranchStruct)

def getStruct(name,varDict):
    strToProcess = 'struct {0}'.format(name) + ' {'
    strForBranch = ''
    typeMap = {'I': 'Int_t', 'F': 'Float_t', 'l': 'ULong64_t', 'C': 'Char_t'}
    for t in ['F','I','l']:
        if t not in varDict: continue
        for v,var in enumerate(varDict[t]):
            strToProcess += '{0} {1};'.format(typeMap[t],var)
            if v==0:
                strForBranch += '{0}/{1}:'.format(var,t)
            else:
                strForBranch += '{0}:'.format(var)
    strToProcess += '};'
    strForBranch = strForBranch[:-1]
    if not hasattr(rt,name): rt.gROOT.ProcessLine(strToProcess)
    struct = getattr(rt,name)()
    return struct, strForBranch

def buildNtuple(object_definitions,states,channelName,final_states,**kwargs):
    '''
    A function to build an initial state ntuple for AnalyzerBase.py
    '''
    alternateIds = kwargs.pop('altIds',[])
    doVBF = kwargs.pop('doVBF',False)
    doMetUnc = kwargs.pop('doMetUnc',False)

    finalStateObjects = 'emtjgn'
    structureDict = {}
    structOrder = []

    # count objects
    nl = 0
    nj = 0
    np = 0
    for key in states[0]:
        val = object_definitions[key]
        for obj in val:
            if obj=='n': continue
            else:
                if obj in 'emt':
                    nl += 1
                if obj == 'j':
                    nj += 1
                if obj == 'g':
                    np += 1


    # some common things
    baseMetVars = ['met','metPhi','mT']
    metVars = []
    metVars += baseMetVars
    baseLepMetVars = ['mass','Pt','sT','dPhi']
    lepMetVars = []
    if doMetUnc:
        for s in ['JetRes','JetEn','MuonEn','ElectronEn','TauEn','UnclusteredEn','PhotonEn']:
            for d in ['Up','Down']:
                for m in baseMetVars:
                    metVars += ['{0}{1}{2}'.format(m,s,d)]
                for m in baseLepMetVars:
                    lepMetVars += ['{0}{1}{2}'.format(m,s,d)]

    lepFloats = ['Pt', 'Eta', 'Phi', 'Iso', 'Dxy', 'Dz', 'SigmaIEtaIEta', 'DEtaIn', 'DPhiIn', 'HOverE', 'OoEmOoP', 'TriggeringMVA', 'NonTriggeringMVA', 'NormalizedChi2', 'JetPt', 'JetBTag']
    for t in ['','_up','_down']:
        for l in ['Loose','Medium','Tight']:
            lepFloats += ['LepScale{0}{1}'.format(l,t)]
            if l=='Loose': continue
            lepFloats += ['LepEff{0}{1}'.format(l,t)]
            lepFloats += ['LepFake{0}{1}'.format(l,t)]
    lepInts = ['Chg', 'PassLoose', 'PassMedium', 'PassTight', 'GenIsPrompt', 'GenPdgId', 'GenPatPdgId', 'MotherGenPdgId', 'ChargeConsistent', 'ExpectedMissingInnerHits', 'PassConversionVeto', 'IsGlobalMuon', 'IsPFMuon', 'IsTrackerMuon', 'ValidMuonHits', 'MatchedStations', 'ValidPixelHits', 'TrackerLayers']



    # define selection bools
    numObjs = len(final_states[0])
    allowedObjects = ''
    for fsObj in finalStateObjects:
        for fs in final_states:
            if fsObj in fs:
                allowedObjects += fsObj
                break
    numObjTypes = len(allowedObjects)
    
    selectVars = {}
    selectVars['I'] = ['passTight','passLoose']
    for trig in ['DoubleMuon','DoubleEG','MuonEG','EGMuon','SingleMuon','SingleEG']:
        selectVars['I'] += ['pass{0}'.format(trig)]
    for altId in alternateIds:
        selectVars['I'] += ['pass_{0}'.format(altId)]
    selectName = 'structSelect_t'
    selectStruct, selectStrForBranch = getStruct(selectName,selectVars)

    structureDict['select'] = [selectStruct, selectStruct, selectStrForBranch]
    structOrder += ['select']    

    # define common root classes
    eventVars = {}
    eventVars['I'] = ['run','lumi','nvtx','GenNUP']
    eventVars['l'] = ['evt'] # ULong64
    eventVars['F'] = ['gen_weight','charge_uncertainty','trig_prescale']
    scaleVars = ['lep_scale','lep_scale_e','lep_scale_m','pu_weight','fakerate','trig_scale']
    for v in scaleVars:
        for t in ['','_up','_down']:
            eventVars['F'] += ['{0}{1}'.format(v,t)]
    eventName = 'structEvent_t'
    eventStruct, eventStrForBranch = getStruct(eventName,eventVars)
    structOrder += ['event']
    structureDict['event'] = [eventStruct,eventStruct,eventStrForBranch]

    # add channels
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
    fakeChannelStruct = rt.structChannel_t()
    structureDict['fakeChannel'] = [fakeChannelStruct, rt.AddressOf(fakeChannelStruct,'channel'),'channel/C']
    structOrder += ['fakeChannel']
    if channelName=='WZ':
        fakeChannelTightWStruct = rt.structChannel_t()
        structureDict['fakeChannel_tightW'] = [fakeChannelTightWStruct, rt.AddressOf(fakeChannelTightWStruct,'channel'),'channel/C']
        structOrder += ['fakeChannel_tightW']
        fakeChannelAllMediumStruct = rt.structChannel_t()
        structureDict['fakeChannel_allMedium'] = [fakeChannelAllMediumStruct, rt.AddressOf(fakeChannelAllMediumStruct,'channel'),'channel/C']
        structOrder += ['fakeChannel_allMedium']

    # add final state structs
    fsVars = {}
    fsVars['F'] = ['mass','eta','phi','sT']
    fsVars['F'] += metVars
    fsVars['F'] += ['leadJetPt','leadJetEta','leadJetPhi','leadJetPUMVA']
    fsVars['I'] = ['jetVeto20','jetVeto30','jetVeto40']
    fsVars['I'] += ['elecVetoLoose','elecVetoTight','muonVetoLoose','muonVetoTight']
    for pt in ['20','30']:
        for l in ['Loose','Medium','Tight']:
            fsVars['I'] += ['bjetVeto{0}{1}'.format(pt,l)]
    if doVBF:
        fsVars['F'] += ['vbfMass','vbfPt','vbfPt1','vbfPt2','vbfEta1','vbfEta2']
        fsVars['I'] += ['centralJetVeto20','centralJetVeto30']
    fsName = 'structFinalState_t'
    fsStruct, fsStrForBranch = getStruct(fsName,fsVars)
    structureDict['finalstate'] = [fsStruct, fsStruct, fsStrForBranch]
    structOrder += ['finalstate']

    # add object struct
    objVars = {}
    objVars['F'] = lepFloats
    objVars['I'] = lepInts
    objName = 'structObject_t'

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
                objStruct, objStrForBranch = getStruct(objName,objVars)
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
                structureDict['%s%i' % (charName, objCount)] = [objStruct, objStruct, objStrForBranch]
                structureDict['%s%iFlv' % (charName, objCount)] = [flvStruct, rt.AddressOf(flvStruct,'Flv'),'Flv/C']
                structOrder += ['%s%i' % (charName, objCount)]
                structOrder += ['%s%iFlv' % (charName, objCount)]

    # define objects for each initial state
    for state in states:
        for key in state:
            val = object_definitions[key]
            stateVars = {'I':[],'F':[]}
            stateName = 'struct{0}_t'.format(key.upper())
            stateVars['F'] += ['mass','Pt','sT','dPhi']
            if 'n' not in val:
                stateVars['F'] += ['dR']
            objCount = 0
            for obj in val:
                if obj == 'n':
                    #stateVars['F'] += metVars
                    stateVars['F'] += lepMetVars
                else:
                    objCount += 1
                    for v in lepFloats:
                        stateVars['F'] += ['{0}{1}'.format(v,objCount)]
                    for v in lepInts:
                        stateVars['I'] += ['{0}{1}'.format(v,objCount)]
                    # manually add the W Z deltaRs for now
                    if key == 'w1':
                        stateVars['F'] += ['dR1_z1_1','dR1_z1_2','mll_z1_1','mll_z1_2','dR1_leadJet']
                    for altId in alternateIds:
                        stateVars['I'] += ['pass_{0}_{1}'.format(altId, objCount)]
            initialStruct, initialStrForBranch = getStruct(stateName,stateVars)
            structureDict[key] = [initialStruct, initialStruct, initialStrForBranch]
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

    

