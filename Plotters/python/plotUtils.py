import itertools
import os
import sys
import errno
import argparse

_3L_MASSES = [170, 200, 250, 300, 350, 400, 450, 500, 600, 700]
_4L_MASSES = [130, 150, 170, 200, 250, 300, 350, 400, 450, 500, 600, 700]

ZMASS = 91.1876

def python_mkdir(dir):
    '''A function to make a unix directory as well as subdirectories'''
    try:
        os.makedirs(dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dir):
            pass
        else: raise

def getChannelSidebandSignalRegion(region,channel,**kwargs):
    mass = kwargs.pop('mass',500) # for higgs
    regionMap = {'Hpp3l':{}, 'Hpp4l':{}}
    regionMap['Hpp3l'][0] = {
        'sr' : 'hN.mass>0.9*%f && hN.mass<1.1*%f' %(mass,mass),
        'sb' : '((hN.mass<0.9*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
    }
    regionMap['Hpp3l'][1] = {
        'sr' : 'hN.mass>0.5*%f && hN.mass<1.1*%f' %(mass,mass),
        'sb' : '((hN.mass<0.5*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
    }
    regionMap['Hpp3l'][2] = {
        'sr' : 'hN.mass>0.5*%f && hN.mass<1.1*%f' %(mass,mass),
        'sb' : '((hN.mass<0.5*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
    }
    regionMap['Hpp4l'][0] = {
        'sr' : 'hN.mass>0.9*%f && hN.mass<1.1*%f' %(mass,mass),
        #'sb' : '((hN.mass<0.9*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
        'sb' : '(hN.mass>12. && hN.mass<800. && !(hN.mass>0.9*%f && hN.mass<1.1*%f))' %(mass,mass)
    }
    regionMap['Hpp4l'][1] = {
        'sr' : 'hN.mass>0.5*%f && hN.mass<1.1*%f' %(mass,mass),
        #'sb' : '((hN.mass<0.5*%f && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
        'sb' : '(hN.mass>12. && hN.mass<800. && !(hN.mass>0.5*%f && hN.mass<1.1*%f))' %(mass,mass)
    }
    regionMap['Hpp3l'][2] = {
        'sr' : 'hN.mass>0.5*%f && hN.mass<1.1*%f' %(mass,mass),
        #'sb' : '((hN.mass<0.5*%f-20. && hN.mass>12.) ||  (hN.mass>1.1*%f && hN.mass<800.))' %(mass,mass)
        'sb' : '(hN.mass>12. && hN.mass<800. && !(hN.mass>0.5*%f && hN.mass<1.1*%f))' %(mass,mass)
    }
    if region == 'Hpp3l':
        numTaus = channel[:2].count('t')
        theMap = regionMap['Hpp3l'][numTaus]
        cutMap = {'srcut': theMap['sr'].replace('hN','h1'), 'sbcut': theMap['sb'].replace('hN','h1')}
    if region == 'Hpp4l':
        h1Taus = channel[:2].count('t')
        h2Taus = channel[2:].count('t')
        h1Map = regionMap['Hpp4l'][h1Taus]
        h2Map = regionMap['Hpp4l'][h2Taus]
        cutMap = {
            'srcut': h1Map['sr'].replace('hN','h1') + ' && ' + h2Map['sr'].replace('hN','h2'),
            'sbcut': h1Map['sb'].replace('hN','h1') + ' && ' + h2Map['sb'].replace('hN','h2'),
        }
    return cutMap


def getChannelCutFlowMap(region,channel,**kwargs):
    mass = kwargs.pop('mass',500) # for higgs
    regionMap = { 'Hpp3l' : {}, 'Hpp4l' : {}, 'WZ' : {}, 'Z' : {}, 'TT' : {}, }
    regionMap['Hpp3l'][0] = {
        #'st' : 'finalstate.sT>1.1*%f+60.' %mass,
        'st' : 'finalstate.sT>1.07*%f+45.' %mass,
        'zveto' : 'abs(z1.mass-%f)>80.' %ZMASS,
        'met' : None,
        'dphi' : 'abs(hN.dPhi)<%f/600.+1.95' %mass,
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass,
        'dr' : 'hN.dR<%f/380.+2.06'%mass if mass<400 else 'hN.dR<%f/1200.+2.77'%mass,
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp3l'][1] = {
        #'st' : 'finalstate.sT>0.85*%f+125.' %mass,
        'st' : 'finalstate.sT>0.72*%f+50.' %mass,
        'zveto' : 'abs(z1.mass-%f)>80.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'abs(hN.dPhi)<%f/200.+1.15' %mass,
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'dr' : 'hN.dR<%f/380.+1.96'%mass if mass<400 else 'hN.dR<%f/1000.+2.6'%mass,
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp3l'][2] = {
        #'st' : '(finalstate.sT>%f-10||finalstate.sT>200.)' %mass,
        'st' : 'finalstate.sT>0.44*%f+65' %mass,
        'zveto' : 'abs(z1.mass-%f)>50.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'abs(hN.dPhi)<2.1',
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'dr' : 'hN.dR<%f/380.+1.86'%mass if mass<400 else 'hN.dR<%f/750.+2.37'%mass,
        #'mass' : 'hN.mass>0.5*%f-20.&&hN.mass<1.1*%f' %(mass,mass),
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp4l'][0] = {
        'st' : 'finalstate.sT>0.6*%f+130.' %mass,
        'zveto' : None,
        'dphi' : None,
        'dr' : None,
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp4l'][1] = {
        'st' : '(finalstate.sT>%f+100.||finalstate.sT>400.)' %mass,
        'zveto' : 'abs(z1.mass-%f)>10.&&abs(z2.mass-%f)>10.' %(ZMASS,ZMASS),
        'dphi' : None,
        'dr' : None,
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass),
    }
    regionMap['Hpp4l'][2] = {
        'st' : 'finalstate.sT>120.',
        'zveto' : 'abs(z1.mass-%f)>50.&&abs(z2.mass-%f)>50.' %(ZMASS,ZMASS),
        'dphi' : 'abs(hN.dPhi)<2.5',
        'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'mass' : 'hN.mass>0.5*%f-20.&&hN.mass<1.1*%f' %(mass,mass),
    }

    cutMap = {}
    if region == 'Hpp3l':
        numTaus = channel[:2].count('t')
        theMap = regionMap['Hpp3l'][numTaus]
        cuts = {}
        for cut in theMap:
            if theMap[cut]:
                cuts[cut] = theMap[cut].replace('hN','h1')
            else:
                cuts[cut] = '1'
        cuts['pre'] = '1'
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['met'], cuts['dphi'], cuts['mass']]
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dphi'], cuts['mass']]
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dr'], cuts['mass']]
        cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta R','Mass window']
    if region == 'Hpp4l':
        h1Taus = channel[:2].count('t')
        h2Taus = channel[2:].count('t')
        h1Map = regionMap['Hpp4l'][h1Taus]
        h2Map = regionMap['Hpp4l'][h2Taus]
        cuts = {}
        for cut in h1Map:
            thisCut = []
            if h1Map[cut]:
                thisCut += [h1Map[cut].replace('hN','h1')]
            if h2Map[cut]:
                thisCut += [h2Map[cut].replace('hN','h2')]
            if thisCut:
                cuts[cut] = ' && '.join(thisCut)
            else:
                cuts[cut] = '1'
        cuts['pre'] = '1'
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dr'], cuts['mass']]
        cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta R','Mass window']
    return cutMap


def defineCutFlowMap(region,channels,mass):
    # define regions (based on number of taus in higgs candidate)
    regionMap = { 'Hpp3l' : {}, 'Hpp4l' : {}, 'WZ' : {}, 'Z' : {}, 'TT' : {}, }
    regionMap['Hpp3l'][0] = {
        #'st' : 'finalstate.sT>1.1*%f+60.' %mass,
        'st' : 'finalstate.sT>1.07*%f+45.' %mass,
        'zveto' : 'fabs(z1.mass-%f)>80.' %ZMASS,
        'met' : None,
        'dphi' : 'fabs(hN.dPhi)<%f/600.+1.95' %mass,
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass,
        'dr' : 'hN.dR<%f/380.+2.06'%mass if mass<400 else 'hN.dR<%f/1200.+2.77'%mass,
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp3l'][1] = {
        #'st' : 'finalstate.sT>0.85*%f+125.' %mass,
        'st' : 'finalstate.sT>0.72*%f+50.' %mass,
        'zveto' : 'fabs(z1.mass-%f)>80.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'fabs(hN.dPhi)<%f/200.+1.15' %mass,
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'dr' : 'hN.dR<%f/380.+1.96'%mass if mass<400 else 'hN.dR<%f/1000.+2.6'%mass,
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp3l'][2] = {
        #'st' : '(finalstate.sT>%f-10||finalstate.sT>200.)' %mass,
        'st' : 'finalstate.sT>0.44*%f+65' %mass,
        'zveto' : 'fabs(z1.mass-%f)>50.' %ZMASS,
        'met' : 'finalstate.met>20.',
        'dphi' : 'fabs(hN.dPhi)<2.1',
        #'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'dr' : 'hN.dR<%f/380.+1.86'%mass if mass<400 else 'hN.dR<%f/750.+2.37'%mass,
        #'mass' : 'hN.mass>0.5*%f-20.&&hN.mass<1.1*%f' %(mass,mass)
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp4l'][0] = {
        'st' : 'finalstate.sT>0.6*%f+130.' %mass,
        'zveto' : None,
        'dphi' : None,
        'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'mass' : 'hN.mass>0.9*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp4l'][1] = {
        'st' : '(finalstate.sT>%f+100.||finalstate.sT>400.)' %mass,
        'zveto' : 'fabs(z1.mass-%f)>10.&&fabs(z2.mass-%f)>10.' %(ZMASS,ZMASS),
        'dphi' : None,
        'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
        'mass' : 'hN.mass>0.5*%f&&hN.mass<1.1*%f' %(mass,mass)
    }
    regionMap['Hpp4l'][2] = {
        'st' : 'finalstate.sT>120.',
        'zveto' : 'fabs(z1.mass-%f)>50.&&fabs(z2.mass-%f)>50.' %(ZMASS,ZMASS),
        'dphi' : 'fabs(hN.dPhi)<2.5',
        'dr' : 'hN.dR<%f/1400.+2.43' %mass, # TODO optimize
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
    regionMap['Z'][0] = {
        'zpt' : '(z1.Pt1>20.&z1.Pt2>10.)',
        'zmass' : 'fabs(z1.mass-%f)<20.' % ZMASS,
        'metveto' : 'met<30.',
        'bveto' : 'finalstate.bjetVeto30Medium==0',
    }
    regionMap['TT'][0] = {
        'zveto' : 'fabs(z1.mass-%f)>20.' % ZMASS,
        'met': 'met>30.',
        'jets' : 'finalstate.jetVeto30>1',
        'bjets' : 'finalstate.bjetVeto30Medium>0',
    }
    # define cutmap to be returned
    cutMap = { 'cuts' : [], 'labels': [], 'labels_simple': [], 'preselection': [] }
    if region == 'Hpp3l':
        #cutMap['labels'] = ['Preselection','s_{T}','Z Veto','E_{T}^{miss}','#Delta#phi','Mass window']
        #cutMap['labels_simple'] = ['Preselection','sT','Z Veto','MET','dPhi','Mass window']
        #cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta#phi','Mass window']
        cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta R','Mass window']
        #cutMap['labels_simple'] = ['Preselection','sT','Z Veto','dPhi','Mass window']
        cutMap['labels_simple'] = ['Preselection','sT','Z Veto','dR','Mass window']
        cutMap['preselection'] = ['All events', 'Three Lepton', 'Trigger', 'Fiducial',\
                                  'Trigger threshold', 'Lepton ID', 'Isolation', 'QCD Rejection',\
                                  'Lepton Charge', '4th Lepton Veto']
        cuts = { 'pre' : '',
                 'st' : '',
                 'zveto' : '',
                 'met' : '',
                 'dphi' : '',
                 'dr' : '',
                 'mass' : '' }
    elif region == 'Hpp4l':
        #cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta#phi','Mass window']
        #cutMap['labels'] = ['Preselection','s_{T}','Z Veto','#Delta R','Mass window']
        cutMap['labels'] = ['Preselection','s_{T}','Mass window']
        #cutMap['labels_simple'] = ['Preselection','sT','Z Veto','dPhi','Mass window']
        #cutMap['labels_simple'] = ['Preselection','sT','Z Veto','dR','Mass window']
        cutMap['labels_simple'] = ['Preselection','sT','Mass window']
        cutMap['preselection'] = ['All events', 'Four Lepton', 'Trigger', 'Fiducial',\
                                  'Trigger threshold', 'Lepton ID', 'Isolation', 'QCD Rejection',\
                                  'Lepton Charge']
        cuts = { 'pre' : '',
                 'st' : '',
                 'zveto' : '',
                 'dphi' : '',
                 'dr' : '',
                 'mass' : '' }
    elif region == 'WZ':
        cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        cutMap['labels'] = ['Preselection (ID)',\
                            'Z lepton p_{T}',\
                            'Z window',\
                            'Mass 3l',\
                            #'b-jet Veto',\
                            'W #DeltaR to Z',\
                            'W lepton p_{T}',\
                            'E_{T}^{miss}'\
                           ]
        cutMap['labels_simple'] = ['Presel (ID)',\
                                   'Z lep pt',\
                                   'Z window',\
                                   'mass3l',\
                                   #'bjet Veto',\
                                   'W dR',\
                                   'W lep pt',\
                                   'MET'\
                                  ]
        cutMap['cuts'] = ['1',\
                          regionMap['WZ'][0]['zpt'],\
                          regionMap['WZ'][0]['zmass'],\
                          regionMap['WZ'][0]['mass'],\
                          #regionMap['WZ'][0]['bveto'],\
                          regionMap['WZ'][0]['wdr'],\
                          regionMap['WZ'][0]['wpt'],\
                          regionMap['WZ'][0]['met']\
                         ]
        # this is 8 tev, add a period check later
        #cutMap['labels'] = ['Preselection (ID)', 'Z lepton p_{T}', 'Z window', 'Mass 3l', 'W #DeltaR to Z', 'W lepton p_{T}', 'E_{T}^{miss}']
        #cutMap['labels_simple'] = ['Presel (ID)', 'Z lep pt', 'Z window', 'mass3l', 'W dR', 'W lep pt', 'MET']
        #cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        #cutMap['cuts'] = ['1', regionMap['WZ'][0]['zpt'], regionMap['WZ'][0]['zmass'], regionMap['WZ'][0]['mass'],\
        #                  regionMap['WZ'][0]['wdr'], regionMap['WZ'][0]['wpt'], regionMap['WZ'][0]['met']]
    elif region=='Z':
        cutMap['labels'] = ['Preselection (ID)', 'Z lepton p_{T}', 'Z window', 'met veto', 'b veto']
        cutMap['labels_simple'] = ['Presel (ID)', 'Z lep pt', 'Z window', 'met veto', 'b veto']
        cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        cutMap['cuts'] = ['1', regionMap['Z'][0]['zpt'], regionMap['Z'][0]['zmass'], regionMap['Z'][0]['metveto'], regionMap['Z'][0]['bveto']]
    elif region=='TT':
        cutMap['labels'] = ['Preselection (ID)', 'Z veto', 'met', 'jets', 'b jets']
        cutMap['labels_simple'] = ['Presel (ID)', 'Z veto', 'met', 'jets', 'b jets']
        cutMap['preselection'] = ['All events','Three lepton','Trigger','Fiducial','4th lepton veto']
        cutMap['cuts'] = ['1', regionMap['TT'][0]['zveto'], regionMap['TT'][0]['met'], regionMap['TT'][0]['jets'], regionMap['TT'][0]['bjets']]
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
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['met'], cuts['dphi'], cuts['mass']]
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dphi'], cuts['mass']]
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dr'], cuts['mass']]
    if region == 'Hpp4l':
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dphi'], cuts['mass']]
        #cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['zveto'], cuts['dr'], cuts['mass']]
        cutMap['cuts'] = [cuts['pre'], cuts['st'], cuts['mass']]
    return cutMap

def getChannels(numLeptons,**kwargs):
    '''Get channels for a given region.'''
    runTau = kwargs.pop('runTau',False)
    leptons = ['l%i' %(x+1) for x in range(numLeptons)]
    lepTypes = 'emt' if runTau else 'em'
    lepPairs = [x[0]+x[1] for x in itertools.product(lepTypes,repeat=2)]
    if numLeptons == 2:
        channels = [x[0]+x[1] for x in itertools.product(lepTypes,lepTypes)]
    elif numLeptons == 3:
        channels = [x[0]+x[1] for x in itertools.product(lepPairs,lepTypes)]
    else:
        channels = [x[0]+x[1] for x in itertools.product(lepPairs,lepPairs)]
    return channels,leptons

def getMergeDict(period):
    '''Return a dictionary of samples to merge in plotting.'''
    sampleMergeDict = {}
    # 13 tev
    if period==13:
        # phys14 samples
        #sampleMergeDict['SingleTop'] = {
        #     'TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola': '1',
        #     'TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola': '1',
        #     'TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola'   : '1',
        #     'TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola'   : '1',
        #     'T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola'          : '1',
        #     'Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola'       : '1',
        #}
        #sampleMergeDict['Diboson']   = {
        #    'WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola': '1',
        #    'ZZTo4L_Tune4C_13TeV-powheg-pythia8'       : '1',
        #}
        #sampleMergeDict['WZJets']    = {
        #    'WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola': '1',
        #}
        #sampleMergeDict['ZZJets']    = {
        #    'ZZTo4L_Tune4C_13TeV-powheg-pythia8': '1',
        #}
        #sampleMergeDict['TTJets']    = {
        #    'TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola': '1',
        #}
        #sampleMergeDict['ZJets']     = {
        #    'DYJetsToLL_M-50_13TeV-madgraph-pythia8': '1',
        #}
        #sampleMergeDict['WJets']     = {
        #    'WJetsToLNu_13TeV-madgraph-pythia8-tauola': '1',
        #}
        #sampleMergeDict['TTVJets']   = {
        #    'TTWJets_Tune4C_13TeV-madgraph-tauola': '1',
        #    'TTZJets_Tune4C_13TeV-madgraph-tauola': '1',
        #}
        # RunIISpring15Dr
        sampleMergeDict['SingleTop'] = {
            'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1'       : '1',
            #'ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1'       : '1',
            #'ST_t-channel_5f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1'       : '1',
            'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1' : '1',
            'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1'     : '1',
            #'ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'  : '1',
            'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'     : '1',
            #'ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'      : '1',
            'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'         : '1',
        }
        sampleMergeDict['TTJets'] = {
            'TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : '1',
            #'TTTo2L2Nu_13TeV-powheg'                         : '1',
            #'TT_TuneCUETP8M1_13TeV-powheg-pythia8'           : '1',
        }
        sampleMergeDict['ZJets'] = {
            'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : '1',
            'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'     : '1',
        }
        sampleMergeDict['WJets'] = {
            'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : '1',
        }
        sampleMergeDict['WWJets'] = {
            'WWTo2L2Nu_13TeV-powheg'        : '1',
            'WWTo4Q_13TeV-powheg'           : '1',
            'WWToLNuQQ_13TeV-powheg'        : '1',
            #'WW_TuneCUETP8M1_13TeV-pythia8' : '1',
        }
        sampleMergeDict['WZJets'] = {
            #'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8' : '1',
            'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8'     : '1',
            #'WZ_TuneCUETP8M1_13TeV-pythia8'                  : '1',
        }
        sampleMergeDict['ZZJets'] = {
            'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8' : '1',
            'ZZTo4L_13TeV_powheg_pythia8'                 : '1',
            'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8'   : '1',
            #'ZZ_TuneCUETP8M1_13TeV-pythia8'               : '1',
        }
        sampleMergeDict['TTVJets'] = {
        }
        sampleMergeDict['QCD'] = {
            #'QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    : '1',
            #'QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    : '1',
            #'QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    : '1',
            #'QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    : '1',
            #'QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8' : '1',
            #'QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8' : '1',

            #'QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8' : '1',

            #'QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   : '1',
            #'QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8'  : '1',
            #'QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8' : '1',
            #'QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8' : '1',
            #'QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8' : '1',

            'QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8'      : '1',
            'QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8'     : '1',
            'QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8'     : '1',
            'QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8'     : '1',
            'QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8'     : '1',
            'QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8'    : '1',
            'QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8'   : '1',
            'QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8'   : '1',
            'QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8'   : '1',
            'QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8'   : '1',
            'QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8'  : '1',
            'QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8' : '1',
            'QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8' : '1',
            'QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8' : '1',
            'QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8' : '1',
            'QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8'  : '1',
        }
        sampleMergeDict['data']      = {
            'data_Run2015B': '1',
        }
    # 8 TeV sample aliases
    if period==8:
        sampleMergeDict['WWJets']    = {
            'WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola': '1',
            #'WW_TuneZ2star_8TeV_pythia6_tauola': '1',
        }
        sampleMergeDict['WZJets']    = {
            'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola': '1',
            'WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola'    : '1',
            #'WZ_TuneZ2star_8TeV_pythia6_tauola': '1',
        }
        sampleMergeDict['ZZJets']    = {
            'ZZTo2e2mu_8TeV-powheg-pythia6'                : '1',
            'ZZTo2e2tau_8TeV-powheg-pythia6'               : '1',
            'ZZTo2mu2tau_8TeV-powheg-pythia6'               : '1',
            'ZZTo4e_8TeV-powheg-pythia6'                    : '1',
            'ZZTo4mu_8TeV-powheg-pythia6'                   : '1',
            'ZZTo4tau_8TeV-powheg-pythia6'                  : '1',
            #'ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola'    : '1',
            'ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola'  : '1',
            'ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola' : '1',
            #'ZZ_TuneZ2star_8TeV_pythia6_tauola'             : '1',
            'GluGluToZZTo4L_8TeV-gg2zz-pythia6'             : '1',
            'GluGluToZZTo2L2L_TuneZ2star_8TeV-gg2zz-pythia6': '1',
        }
        sampleMergeDict['HZZ']       = {
            'GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6': '1',
        }
        sampleMergeDict['WGJets']    = {
            'WGstarToLNu2E_TuneZ2star_8TeV-madgraph-tauola' : '1',
            'WGstarToLNu2Mu_TuneZ2star_8TeV-madgraph-tauola': '1',
        }
        sampleMergeDict['SingleTop'] = {
            'T_s-channel_TuneZ2star_8TeV-powheg-tauola'       : '1',
            'T_t-channel_TuneZ2star_8TeV-powheg-tauola'       : '1',
            'T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'   : '1',
            'Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola'    : '1',
            'Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola'    : '1',
            'Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola': '1',
        }
        sampleMergeDict['TTJets']    = {
            'TTJetsFullLepMGDecays': '1',
            'TTJetsSemiLepMGDecays': '1',
        }
        sampleMergeDict['ZJets']     = {
            'DYJetsToLL_M-10To50filter_8TeV-madgraph'         : '1',
            'Z1jets_M50'                                      : '1',
            'Z2jets_M50_S10'                                  : '1',
            'Z3jets_M50'                                      : '1',
            'Z4jets_M50'                                      : '1',
            'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball': 'event.GenNUP==5',
            #'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball': '1',
            #'DYJetsToLL_M-10To50filter_8TeV-madgraph_filtered'         : '1',
            #'Z1jets_M50_filtered'                                      : '1',
            #'Z2jets_M50_S10_filtered'                                  : '1',
            #'Z3jets_M50_filtered'                                      : '1',
            #'Z4jets_M50_filtered'                                      : '1',
            #'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_filtered': 'event.GenNUP==5',
            #'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_filtered': '1',
        }
        sampleMergeDict['WJets']     = {
            'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-v1': '1',
            'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-v2': '1',
        }
        sampleMergeDict['ZG']        = {
            'ZGToLLG_8TeV-madgraph': '1',
        }
        sampleMergeDict['data']      = {
            'data_Run2012A': '1',
            'data_Run2012B': '1',
            'data_Run2012C': '1',
            'data_Run2012D': '1',
            #'data_DoubleElectron_Run2012A_22Jan2013_v1': '1',
            #'data_DoubleElectron_Run2012B_22Jan2013_v1': '1',
            #'data_DoubleElectron_Run2012C_22Jan2013_v1': '1',
            #'data_DoubleElectron_Run2012D_22Jan2013_v1': '1',
            #'data_DoubleMu_Run2012A_22Jan2013_v1': '1',
            #'data_DoubleMuParked_Run2012B_22Jan2013_v1': '1',
            #'data_DoubleMuParked_Run2012C_22Jan2013_v1': '1',
            #'data_DoubleMuParked_Run2012D_22Jan2013_v1': '1',
            #'data_MuEG_Run2012A_22Jan2013_v1': '1',
            #'data_MuEG_Run2012B_22Jan2013_v1': '1',
            #'data_MuEG_Run2012C_22Jan2013_v1': '1',
            #'data_MuEG_Run2012D_22Jan2013_v1': '1',
        }
        sampleMergeDict['TTVJets']   = {
            'TTZJets' : '1',
            'TTWJets' : '1',
            'TTWWJets': '1',
            'TTGJets' : '1',
        }
        sampleMergeDict['VVVJets']   = {
            'ZZZNoGstarJets': '1',
            'WZZNoGstarJets': '1',
            'WWZNoGstarJets': '1',
            'WWWJets'       : '1',
        }
    return sampleMergeDict


def getSigMap(numLeptons,mass):
    '''Return a signal map for a given running period'''
    sigMap = {
        8 : {
             'ZZ'  : 'ZZJets',
             'WZ'  : 'WZJets',
             'WW'  : 'WWJets',
             'WG'  : 'WGJets',
             'Z'   : 'ZJets',
             'Zlow': 'DYJetsToLL_M-10To50filter_8TeV-madgraph',
             'W'   : 'WJets',
             'ZG'  : 'ZG',
             'HZZ' : 'HZZ',
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
             'WW'  : 'WWJets',
             'Z'   : 'ZJets',
             'W'   : 'WJets',
             'TT'  : 'TTJets',
             'T'   : 'SingleTop',
             'TTV' : 'TTVJets',
             'QCD' : 'QCD',
             'Sig' : 'DBLH_m500',
             'data': 'data'
        }
    }
    if numLeptons==3:
        for m in _3L_MASSES:
            sigMap[8][m] = 'HPlusPlusHMinusHTo3L_M-%i_8TeV-calchep-pythia6' % m
    if numLeptons==4:
        for m in _4L_MASSES:
            sigMap[8][m] = 'HPlusPlusHMinusMinusHTo4L_M-%i_8TeV-pythia6' % m
    return sigMap

def getIntLumiMap():
    '''Get map of integrated luminosity to scale MC'''
    intLumiMap = {
        7 : 4900,
        8 : 19700,
        #13: 1000,
        13: 40.03,
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
    channelsZ = [['ee'],['mm']]
    channelStringsZ = ['ee','#mu#mu']
    channelCutsZ = ['z1Flv=="%s"' %(x[0]) for x in channelsZ]
    plotChannelCuts = channelCuts
    if region in ['WZ']: plotChannelCuts = channelCutsWZ
    if region in ['Z']: plotChannelCuts = channelCutsZ
    plotChannelStrings = channelStrings
    if region in ['WZ']: plotChannelStrings = channelStringsWZ
    if region in ['Z']: plotChannelStrings = channelStringsZ
    return plotChannelStrings, plotChannelCuts

def getGenChannelStringsCuts(region,genChannels):
    channelCharMap = {'e':'e', 'm':'#mu', 't':'#tau'}
    channelStrings = []
    channelCuts = []
    for channel in genChannels:
        channelString = ''
        for c in channel:
            channelString += channelCharMap[c]
        channelStrings += [channelString]
        theGenChannels = [channel+c for c in ['e','m','t']]
        genCut = '(' + '||'.join(['genChannel=="%s"' %x for x in theGenChannels]) + ')'
        channelCuts += [genCut]
    return channelStrings, channelCuts
            
    

def getChannelBackgrounds(runPeriod):
    channelBackground = {
        'Hpp2l'   : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'TT'      : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Z'       : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'W'       : ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'FakeRate': ['T', 'TT', 'TTV', 'W', 'Z', 'VVV', 'ZZ', 'WW', 'WZ'],
        'Hpp3l'   : ['TT', 'TTV', 'Z', 'VVV', 'ZZ', 'WZ'],
        'WZ'      : ['TT', 'TTV', 'Z', 'VVV', 'ZZ', 'WZ'],
        'LowMass' : ['TT', 'TTV', 'Z', 'VVV', 'ZZ', 'WZ'],
        'Hpp4l'   : ['TT', 'TTV', 'Z', 'VVV', 'ZZ', 'WZ'],
    }   
    if runPeriod==13:
        channelBackground = {
            'WZ'      : ['T', 'TT', 'TTV', 'Z', 'WW', 'ZZ', 'WZ'],
            'W'       : ['T', 'TT', 'TTV', 'W', 'Z', 'WW', 'ZZ', 'WZ'],
            'FakeRate': ['T', 'TT', 'TTV', 'W', 'Z', 'WW', 'ZZ', 'WZ'],
            'Z'       : ['T', 'TT', 'TTV', 'Z', 'WW', 'ZZ', 'WZ'],
            'TT'      : ['T', 'TT', 'TTV', 'Z', 'WW', 'ZZ', 'WZ'],
            'Hpp2l'   : ['T', 'TT', 'TTV', 'Z', 'WW', 'ZZ', 'WZ'],
            'Hpp3l'   : ['TT', 'TTV', 'Z', 'ZZ', 'WZ'],
            'LowMass' : ['TT', 'TTV', 'Z', 'ZZ', 'WZ'],
            'Hpp4l'   : ['TT', 'Z', 'TTV', 'ZZ', 'WZ']
        }
    return channelBackground
