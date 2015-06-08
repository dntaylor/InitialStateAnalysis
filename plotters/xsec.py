'''
Table of cross sections to be used within the FSA framework with the PlotterBase class.
'''

NB = 1.e3
PB = 1.
FB = 1.e-3

BR_WJJ = 0.676
BR_WLNU = 0.324

Z_XSEC = 3503.7 * PB
#TT_XSEC = 252.89 * PB
TT_XSEC = 241.5* PB

xsecs = { 7 : {},
          8 : {},
          13: {} }
xsecs[13] = {
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    'WJetsToLNu_13TeV-madgraph-pythia8-tauola' : 20508.9 * PB,
    'DYJetsToLL_M-50_13TeV-madgraph-pythia8' : 6025.2 * PB,
    
    # https://twiki.cern.ch/wiki/bin/view/LHCPhysics/TtbarNNLO
    'TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola' : 831.76 * PB,

    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopSigma
    'TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola' : 7.20 * PB, 
    'TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola' : 136.02 * PB, 
    'T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola' : 35.6 * PB,
    'TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola' : 4.16 * PB,
    'TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola' : 80.95 * PB,
    'Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola' : 35.6 * PB,

    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
    # Z -> ll = 0.101 * inclusive
    # w -> lv = 0.324 * inclusive

    # from McM
    'WZJetsTo3LNu_Tune4C_13TeV-madgraph-tauola' : 1.634 * PB,
    'ZZTo4L_Tune4C_13TeV-powheg-pythia8' : 1.218 * PB,

    'TTWJets_Tune4C_13TeV-madgraph-tauola' : 1.152 * PB,
    'TTZJets_Tune4C_13TeV-madgraph-tauola' : 2.232 * PB,

    'DBLH_m500' : 0.001652 * PB, # from pythia
}

xsecs[8] = {
    'HPlusPlusHMinusMinusHTo4L_M-110_8TeV-pythia6': 352.49 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-130_8TeV-pythia6': 186.21 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-150_8TeV-pythia6': 106.55 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-170_8TeV-pythia6': 64.641 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-200_8TeV-pythia6': 33.209 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-250_8TeV-pythia6': 12.724 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-300_8TeV-pythia6': 5.5458 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-350_8TeV-pythia6': 2.6413 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-400_8TeV-pythia6': 1.3414 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-450_8TeV-pythia6': 0.71531 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-500_8TeV-pythia6': 0.39604 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-600_8TeV-pythia6': 0.13271 * FB,
    'HPlusPlusHMinusMinusHTo4L_M-700_8TeV-pythia6': 0.48382e-1 * FB,

    # from mario, but the 4l dont match the reference... http://cms.hep.kbfi.ee/~mario/dblh/xsinfo.py
    'HPlusPlusHMinusHTo3L_M-170_8TeV-calchep-pythia6': 114.838 * FB,
    'HPlusPlusHMinusHTo3L_M-200_8TeV-calchep-pythia6': 51.59 * FB,
    'HPlusPlusHMinusHTo3L_M-250_8TeV-calchep-pythia6': 21.172 * FB,
    'HPlusPlusHMinusHTo3L_M-300_8TeV-calchep-pythia6': 8.308 * FB,
    'HPlusPlusHMinusHTo3L_M-350_8TeV-calchep-pythia6': 4.154 * FB,
    'HPlusPlusHMinusHTo3L_M-400_8TeV-calchep-pythia6': 1.9832 * FB,
    'HPlusPlusHMinusHTo3L_M-450_8TeV-calchep-pythia6': 1.02912 * FB,
    'HPlusPlusHMinusHTo3L_M-500_8TeV-calchep-pythia6': 0.55476 * FB,
    'HPlusPlusHMinusHTo3L_M-600_8TeV-calchep-pythia6': 0.17286 * FB,
    'HPlusPlusHMinusHTo3L_M-700_8TeV-calchep-pythia6': 0.05762 * FB,

    'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-v1' : 36257.2 * PB,
    'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-v2' : 36257.2 * PB,

    'DYJetsToLL_M-10To50filter_8TeV-madgraph':          915 * PB,
    'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball': 3503.7 * PB,

    'Z1jets_M50':     Z_XSEC * 0.190169492,
    'Z2jets_M50_S10': Z_XSEC * 0.061355932,
    'Z3jets_M50':     Z_XSEC * 0.017322034,
    'Z4jets_M50':     Z_XSEC * 0.007810169,

    'TTJetsFullLepMGDecays': TT_XSEC * BR_WLNU * BR_WLNU,
    'TTJetsSemiLepMGDecays': TT_XSEC * BR_WJJ * BR_WLNU * 2,

    'TTGJets' : 2.166 * PB,
    'TTWJets' : 0.2057 * PB,
    'TTWWJets' : 0.002 * PB,
    'TTZJets' : 0.232 * PB,

    'T_s-channel_TuneZ2star_8TeV-powheg-tauola':        3.79 * PB,
    'T_t-channel_TuneZ2star_8TeV-powheg-tauola':        56.4 * PB,
    'T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola':    11.1 * PB,
    'Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola':     1.56 * PB,
    'Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola':     30.7 * PB,
    'Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola': 11.1 * PB,

    'ZZTo4e_8TeV-powheg-pythia6':      0.07691 * PB,
    'ZZTo4mu_8TeV-powheg-pythia6':     0.07691 * PB,
    'ZZTo4tau_8TeV-powheg-pythia6':    0.07691 * PB,
    'ZZTo2e2mu_8TeV-powheg-pythia6':   0.1767 * PB,
    'ZZTo2e2tau_8TeV-powheg-pythia6':  0.1767 * PB,
    'ZZTo2mu2tau_8TeV-powheg-pythia6': 0.1767 * PB,

    'ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola':   0.1769 * PB,
    'ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola': 2.4487 * PB,

    'ZGToLLG_8TeV-madgraph': 132.6 * PB,

    'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola': 1.09 * 2.207 * PB, # From CMS WZ analysis data/theo = 1.09
    'WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola':     1.09 * 1.058 * PB, # From CMS WZ analysis data/theo = 1.09

    'WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola' : 54.838*(0.1075+0.1057+0.1125)*(0.1075+0.1057+0.1125) * PB,

    'ZZZNoGstarJets' : 0.0192 * PB,
    'WWZNoGstarJets' : 0.0633 * PB,
    'WZZNoGstarJets' : 0.01968 * PB,
    'WWWJets' : 0.08217 * PB,
    'WWGJets' : 1.44 * PB,
}
