import itertools
import os
import sys
import errno
import argparse

_3L_MASSES = [170, 200, 250, 300, 350, 400, 450, 500, 600, 700]
_4L_MASSES = [110, 130, 150, 170, 200, 250, 300, 350, 400, 450, 500, 600, 700]

ZMASS = 91.1876

def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise

def defineCutFlowMap(region,channels,mass):
    # define regions (based on number of taus in higgs candidate)
    regionMap = { 'Hpp3l' : {}, 'Hpp4l' : {}, 'WZ' : {} }
    regionMap['Hpp3l'][0] = {
        'st' : 'finalstate.sT>1.1*%f+60.' %mass,
        'zveto' : 'fabs(z1.mass-%f)>80.' %ZMASS,
        'met' : None,
        'dphi' : 'hN.dPhi<%f/600.+1.95' %mass,
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp3l'][1] = {
        'st' : 'finalstate.sT>0.85*%f+125.' %mass,
        'zveto' : 'fabs(z1.mass-%f)>80.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'hN.dPhi<%f/200.+1.15' %mass,
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp3l'][2] = {
        'st' : '(finalstate.sT>%f-10||finalstate.sT>200.)' %mass,
        'zveto' : 'fabs(z1.mass-%f)>50.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'hN.dPhi<2.1',
        'mass' : 'hN.mass>0.5*%f-20.&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp4l'][0] = {
        'st' : 'finalstate.sT>0.6*%f+130.' %mass,
        'zveto' : None,
        'dphi' : None,
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp4l'][1] = {
        'st' : '(finalstate.sT>%f+100.||finalstate.sT>400.)' %mass,
        'zveto' : 'fabs(z1.mass-%f)>10.&&fabs(z2.mass-%f)>10.' %(ZMASS,ZMASS),
        'dphi' : None,
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp4l'][2] = {
        'st' : 'finalstate.sT>120.',
        'zveto' : 'fabs(z1.mass-%f)>50.&&fabs(z2.mass-%f)>50.' %(ZMASS,ZMASS),
        'dphi' : 'hN.dPhi<2.5',
        'mass' : None
    }
    regionMap['WZ'][0] = {
        'mass' : 'finalstate.mass>100',
        'zpt' : '(z1.Pt1>20.&z1.Pt2>10.)',
        'zmass' : 'fabs(z1.mass-%f)<20.' % ZMASS,
        'bveto' : 'finalstate.bjetVeto30Medium==0',
        'wdr' : 'w1.dR1_z1_1>0.1 & w1.dR1_z1_2>0.1',
        'wpt' : 'w1.Pt1>20.',
        'met' : 'w1.met>30.',
    }
    # define cutmap to be returned
    cutMap = { 'cuts' : [], 'labels': [], 'labels_simple': [], 'preselection': [] }
    if region == 'Hpp3l':
        cutMap['labels'] = ['Preselection','s_{T}','Z Veto','E_{T}^{miss}','#Delta#phi','Mass window']
        cutMap['labels_simple'] = ['Preselection','sT','Z Veto','MET','dPhi','Mass window']
        cutMap['preselection'] = ['All events', 'Three Lepton', 'Trigger', 'Fiducial',\
                                  'Trigger threshold', 'Lepton ID', 'Isolation', 'QCD Rejection',\
                                  'Lepton Charge', '4th Lepton Veto']
        cuts = { 'pre' : '',
                 'st' : '',
                 'zveto' : '',
                 'met' : '',
                 'dphi' : '',
                 'mass' : '' }
    elif region == 'Hpp4l':
        cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta#phi','Mass window']
        cutMap['labels_simple'] = ['Preselection','sT','Z Veto','dPhi','Mass window']
        cutMap['preselection'] = ['All events', 'Four Lepton', 'Trigger', 'Fiducial',\
                                  'Trigger threshold', 'Lepton ID', 'Isolation', 'QCD Rejection',\
                                  'Lepton Charge']
        cuts = { 'pre' : '',
                 'st' : '',
                 'zveto' : '',
                 'dphi' : '',
                 'mass' : '' }
    elif region == 'WZ' or region == 'TT':
        cutMap['labels'] = ['Preselection (ID)', 'Z lepton p_{T}', 'Z window', 'b-jet Veto', 'W lepton p_{T}', 'E_{T}^{miss}']
        cutMap['labels_simple'] = ['Presel (ID)', 'Z lep pt', 'Z window', 'bjet Veto', 'W lep pt', 'MET']
        cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        cutMap['cuts'] = ['1', regionMap['WZ'][0]['zpt'], regionMap['WZ'][0]['zmass'], regionMap['WZ'][0]['bveto'],\
                          regionMap['WZ'][0]['wpt'], regionMap['WZ'][0]['met']]
        # this is 8 tev, add a period check later
        cutMap['labels'] = ['Preselection (ID)', 'Z lepton p_{T}', 'Z window', 'Mass 3l', 'W #DeltaR to Z', 'W lepton p_{T}', 'E_{T}^{miss}']
        cutMap['labels_simple'] = ['Presel (ID)', 'Z lep pt', 'Z window', 'mass3l', 'W dR', 'W lep pt', 'MET']
        cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        cutMap['cuts'] = ['1', regionMap['WZ'][0]['zpt'], regionMap['WZ'][0]['zmass'], regionMap['WZ'][0]['mass'],\
                          regionMap['WZ'][0]['wdr'], regionMap['WZ'][0]['wpt'], regionMap['WZ'][0]['met']]
    else:
        cutMap['cuts'] = '1'
        cutMap['labels'] = ['%s Full Selection' %region]
        cutMap['labels_simple'] = [region]
    if region not in ['Hpp3l','Hpp4l']: return cutMap
    usedLepPairs = []
    for channel in channels:
        lepPairs = [channel[:2]]
        #if region=='Hpp4l': lepPairs += [channel[2:]]
        if lepPairs in usedLepPairs: continue
        for cut in cuts:
            if cuts[cut] and cuts[cut][-2:]!='||': cuts[cut] += '||'
        hNum = 0
        tempCut = {}
        for lepPair in lepPairs:
            hNum += 1
            numTau = lepPair.count('t')
            for cut in regionMap[region][numTau]:
                thisCut = regionMap[region][numTau][cut]
                if thisCut is not None:
                    if cut in tempCut:
                        tempCut[cut] += '&&'
                    else:
                        tempCut[cut] = ''
                    tempCut[cut] += 'h%iFlv=="%s"&&%s' % (hNum, lepPair, thisCut.replace('N',str(hNum)))
                else:
                    if cut in tempCut:
                        tempCut[cut] += '&&'
                    else:
                        tempCut[cut] = ''
                    tempCut[cut] += 'h%iFlv=="%s"' % (hNum, lepPair)
        for cut in cuts:
            if cut in tempCut:
                cuts[cut] += '(%s)' % tempCut[cut]
        usedLepPairs += [lepPairs]
    for cut in cuts:
        if not cuts[cut]:
            cuts[cut] = '('+'||'.join(['h1Flv=="%s"' %x[0] for x in usedLepPairs])+')'
        else:
            if cuts[cut][-2:]=='||': cuts[cut] = cuts[cut][:-2]
            cuts[cut] = '(%s)' % cuts[cut]

    if region == 'Hpp3l':
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['met'], cuts['dphi'], cuts['mass']]
    if region == 'Hpp4l':
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dphi'], cuts['mass']]
    return cutMap

def getChannels(numLeptons,**kwargs):
    '''Get channels for a given region.'''
    runTau = kwargs.pop('runTau',False)
    leptons = ['l%i' %(x+1) for x in range(numLeptons)]
    lepTypes = 'emt' if runTau else 'em'
    lepPairs = [x[0]+x[1] for x in itertools.combinations_with_replacement(lepTypes,2)]
    if numLeptons == 3: channels = [x[0]+x[1] for x in itertools.product(lepPairs,lepTypes)]
    else: channels = [x[0]+x[1] for x in itertools.product(lepPairs,lepPairs)]
    return channels,leptons

def getMergeDict(period):
    '''Return a dictionary of samples to merge in plotting.'''
    sampleMergeDict = {}
    # 13 tev
    if period==13:
        sampleMergeDict['SingleTop'] = ['TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola',\
                                        'TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola',\
                                        'TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola',\
                                        'TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola',\
                                        'T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola',\
                                        'Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola']
        sampleMergeDict['Diboson']   = ['WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola',\
                                        'ZZTo4L_Tune4C_13TeV-powheg-pythia8']
        sampleMergeDict['WZJets']    = ['WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola']
        sampleMergeDict['ZZJets']    = ['ZZTo4L_Tune4C_13TeV-powheg-pythia8']
        sampleMergeDict['TTJets']    = ['TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola']
        sampleMergeDict['ZJets']     = ['DYJetsToLL_M-50_13TeV-madgraph-pythia8']
        sampleMergeDict['WJets']     = ['WJetsToLNu_13TeV-madgraph-pythia8-tauola']
        sampleMergeDict['TTVJets']   = ['TTWJets_Tune4C_13TeV-madgraph-tauola',\
                                        'TTZJets_Tune4C_13TeV-madgraph-tauola']
    # 8 TeV sample aliases
    if period==8:
        sampleMergeDict['WWJets']    = ['WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola']
        sampleMergeDict['WZJets']    = ['WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola',\
                                        'WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola']
        sampleMergeDict['Diboson']   = ['WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola',\
                                        'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola',\
                                        'WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola',\
                                        'ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola']
        #sampleMergeDict['ZZJets']    = ['ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola']
        sampleMergeDict['ZZJets']    = ['ZZTo2e2mu_8TeV-powheg-pythia6',\
                                        'ZZTo2e2tau_8TeV-powheg-pythia6',\
                                        'ZZTo2mu2tau_8TeV-powheg-pythia6',\
                                        'ZZTo4e_8TeV-powheg-pythia6',\
                                        'ZZTo4mu_8TeV-powheg-pythia6',\
                                        'ZZTo4tau_8TeV-powheg-pythia6']
        sampleMergeDict['SingleTop'] = ['T_s-channel_TuneZ2star_8TeV-powheg-tauola',\
                                        'T_t-channel_TuneZ2star_8TeV-powheg-tauola',\
                                        'T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola',\
                                        'Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola',\
                                        'Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola',\
                                        'Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola']
        sampleMergeDict['TTJets']    = ['TTJetsFullLepMGDecays', 'TTJetsSemiLepMGDecays']
        sampleMergeDict['ZJets']     = ['Z1jets_M50',\
                                        'Z2jets_M50_S10',\
                                        'Z3jets_M50',\
                                        'Z4jets_M50']
        #sampleMergeDict['ZJets']     = ['DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball']
        sampleMergeDict['data']      = ['data_Run2012A',\
                                        'data_Run2012B',\
                                        'data_Run2012C',\
                                        'data_Run2012D']
        sampleMergeDict['TTVJets']   = ['TTZJets',\
                                        'TTWJets',\
                                        'TTWWJets',\
                                        'TTGJets']
        sampleMergeDict['VVVJets']   = ['ZZZNoGstarJets',\
                                        'WWZNoGstarJets',\
                                        'WWWJets']
    return sampleMergeDict


def getSigMap(numLeptons,mass):
    '''Return a signal map for a given running period'''
    sigMap = {
        8 : {
             'ZZ'  : 'ZZJets',
             'WZ'  : 'WZJets',
             'WW'  : 'WWJets',
             'Z'   : 'ZJets',
             'DB'  : 'Diboson',
             'TT'  : 'TTJets',
             'T'   : 'SingleTop',
             'TTV' : 'TTVJets',
             'TTZ' : 'TTZJets',
             'TTW' : 'TTWJets',
             'VVV' : 'VVVJets',
             'ZZZ' : 'ZZZNoGstarJets',
             'WWZ' : 'WWZNoGstarJets',
             'WWW' : 'WWWJets',
             'Sig' : 'HPlusPlusHMinusHTo3L_M-%i_8TeV-calchep-pythia6' % mass\
                      if numLeptons==3 else 'HPlusPlusHMinusMinusHTo4L_M-%i_8TeV-pythia6' % mass,
             'data': 'data'
        },
        13 : {
             'ZZ'  : 'ZZJets',
             'WZ'  : 'WZJets',
             'Z'   : 'ZJets',
             'W'   : 'WJets',
             'DB'  : 'Diboson',
             'TT'  : 'TTJets',
             'T'   : 'SingleTop',
             'TTV' : 'TTVJets',
             'TTZ' : 'TTZJets_Tune4C_13TeV-madgraph-tauola',
             'TTW' : 'TTWJets_Tune4C_13TeV-madgraph-tauola',
             'Sig' : 'DBLH_m500',
        }
    }
    return sigMap

def getIntLumiMap():
    '''Get map of integrated luminosity to scale MC'''
    intLumiMap = {
        7 : 4900,
        8 : 19700,
        #13: 15000
        13: 1000
    }
    return intLumiMap

def getChannelStringsCuts(region,channels):
    channelCharMap = {'e':'e', 'm':'#mu', 't':'#tau'}
    channelStrings = []
    for channel in channels:
        channelString = ''
        for c in channel:
            channelString += channelCharMap[c]
        channelStrings += [channelString]
    channelCuts = ['channel=="%s"' % x for x in channels]
    channelsWZ = [['ee','e'],['ee','m'],['mm','e'],['mm','m']]
    channelStringsWZ = ['(ee)e','(ee)#mu','(#mu#mu)e','(#mu#mu)#mu']
    channelCutsWZ = ['z1Flv=="%s"&&w1Flv=="%s"' %(x[0],x[1]) for x in channelsWZ]
    plotChannelCuts = channelCuts
    if region in ['WZ']: plotChannelCuts = channelCutsWZ
    plotChannelStrings = channelStrings
    if region in ['WZ']: plotChannelStrings = channelStringsWZ
    return plotChannelStrings, plotChannelCuts
