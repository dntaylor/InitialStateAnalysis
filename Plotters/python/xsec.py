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
    # these are the phys14 samples
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

    # here we have all the samples for RunIISpring15DR (all from mcm unless otherwise noted)
    'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'          : 18610   * PB,
    'DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_filtered' : 18610   * PB,
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'              :  6025.2 * PB,
    'DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_filtered'     :  6025.2 * PB,

    'QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    : 1273190000      * PB,
    'QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    :  558528000      * PB,
    'QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    :  139803000      * PB,
    'QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'    :   19222500      * PB,
    'QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'   :    2758420      * PB,
    'QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  :     469797      * PB,
    'QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  :     117989      * PB,
    'QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  :       7820.25   * PB,
    'QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  :        645.528  * PB,
    'QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8'  :        187.109  * PB,
    'QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8' :         32.3486 * PB,
    'QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8' :         10.4305 * PB,

    'QCD_Pt-20toInf_MuEnrichedPt15_TuneCUETP8M1_13TeV_pythia8' : 720648000 * PB,

    'QCD_Pt-15to20_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   : 1279000000 * PB,
    'QCD_Pt-20to30_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   :  557600000 * PB,
    'QCD_Pt-30to50_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   :  136000000 * PB,
    'QCD_Pt-50to80_EMEnriched_TuneCUETP8M1_13TeV_pythia8'   :   19800000 * PB,
    'QCD_Pt-80to120_EMEnriched_TuneCUETP8M1_13TeV_pythia8'  :    2800000 * PB,
    'QCD_Pt-120to170_EMEnriched_TuneCUETP8M1_13TeV_pythia8' :     477000 * PB,
    'QCD_Pt-170to300_EMEnriched_TuneCUETP8M1_13TeV_pythia8' :     114000 * PB,
    'QCD_Pt-300toInf_EMEnriched_TuneCUETP8M1_13TeV_pythia8' :       9000 * PB,
    
    'QCD_Pt_5to10_TuneCUETP8M1_13TeV_pythia8'      : 61018300000           * PB,
    'QCD_Pt_10to15_TuneCUETP8M1_13TeV_pythia8'     :  5887580000           * PB,
    'QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8'     :  1837410000           * PB,
    'QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8'     :   140932000           * PB,
    'QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8'     :    19204300           * PB,
    'QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8'    :     2762530           * PB,
    'QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8'   :      471100           * PB,
    'QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8'   :      117276           * PB,
    'QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8'   :        7823           * PB,
    'QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8'   :         648.2         * PB,
    'QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8'  :          32.293       * PB,
    'QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8' :           9.4183      * PB,
    'QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8' :           0.84265     * PB,
    'QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8' :           0.114943    * PB,
    'QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8' :           0.00682981  * PB,
    'QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8'  :           0.000165445 * PB,
    
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/SingleTopSigma
    'ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1'       :  11.36 * PB,
    'ST_t-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1'       :  70.69 * PB,
    'ST_t-channel_5f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1'       :  70.69 * PB,
    'ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1' :  80.95 * PB,
    'ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1'     : 136.02 * PB,
    'ST_tW_antitop_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'  :  35.6  * PB,
    'ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'     :  35.6  * PB,
    'ST_tW_top_5f_DS_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'      :  35.6  * PB,
    'ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1'         :  35.6  * PB,

    # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    'TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : 831.76 * PB,
    'TTTo2L2Nu_13TeV-powheg'                         :  87.31 * PB,
    'TT_TuneCUETP8M1_13TeV-powheg-pythia8'           : 831.76 * PB,

    'TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8'         : 0.2529 * PB,
    'TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8'                  : 0.5297 * PB,
    'TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8' : 0.2043 * PB,
    'TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8'  : 0.4062 * PB,
    'tZq_ll_4f_13TeV-amcatnlo-pythia8_TuneCUETP8M1'                : 0.0758 * PB,

    'WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' : 61526.7 * PB,

    'ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'          : 117.864 * PB,
    'ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_filtered' : 117.864 * PB,
    'WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8'         : 489     * PB,

    'WWTo2L2Nu_13TeV-powheg'        : 10.481 * PB,
    'WWTo4Q_13TeV-powheg'           : 45.20  * PB,
    'WWToLNuQQ_13TeV-powheg'        : 43.53  * PB,
    'WW_TuneCUETP8M1_13TeV-pythia8' : 63.21  * PB,

    'WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8' : 10.71    * PB,
    'WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8'   :  3.05    * PB,
    'WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8'    :  5.60    * PB,
    'WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8'     :  4.42965 * PB,
    'WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8' :  4.71    * PB,
    'WZ_TuneCUETP8M1_13TeV-pythia8'                  : 47.13    * PB,

    'ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8'   :  6.842 * PB,
    'ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8':  4.04  * PB,
    'ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8' :  3.28  * PB,
    'ZZTo2L2Nu_13TeV_powheg_pythia8'              :  0.564 * PB,
    'ZZTo4L_13TeV_powheg_pythia8'                 :  1.256 * PB * 1.1,
    'ZZTo4L_13TeV-amcatnloFXFX-pythia8'           :  1.212 * PB * 1.1,
    'ZZ_TuneCUETP8M1_13TeV-pythia8'               : 16.523 * PB,

    'GluGluToZZTo2e2mu_BackgroundOnly_13TeV_MCFM'   : 0.003194 * PB * 1.7,
    'GluGluToZZTo2e2tau_BackgroundOnly_13TeV_MCFM'  : 0.003194 * PB * 1.7,
    'GluGluToZZTo2mu2tau_BackgroundOnly_13TeV_MCFM' : 0.003194 * PB * 1.7,
    'GluGluToZZTo4e_BackgroundOnly_13TeV_MCFM'      : 0.001586 * PB * 1.7,
    'GluGluToZZTo4mu_BackgroundOnly_13TeV_MCFM'     : 0.001586 * PB * 1.7,
    'GluGluToZZTo4tau_BackgroundOnly_13TeV_MCFM'    : 0.001586 * PB * 1.7,
    
    'WWW_TuneCUETP8M1_13TeV-amcatnlo-pythia8' : 0.1651  * PB,
    'WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8' : 0.1651  * PB,
    'WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8' : 0.05565 * PB,
    'ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8' : 0.01398 * PB,

}

xsecs[8] = {
    'HPlusPlusHMinusMinusHTo4L_M-110_8TeV-pythia6': 352.49     * FB,
    'HPlusPlusHMinusMinusHTo4L_M-130_8TeV-pythia6': 186.21     * FB,
    'HPlusPlusHMinusMinusHTo4L_M-150_8TeV-pythia6': 106.55     * FB,
    'HPlusPlusHMinusMinusHTo4L_M-170_8TeV-pythia6':  64.641    * FB,
    'HPlusPlusHMinusMinusHTo4L_M-200_8TeV-pythia6':  33.209    * FB,
    'HPlusPlusHMinusMinusHTo4L_M-250_8TeV-pythia6':  12.724    * FB,
    'HPlusPlusHMinusMinusHTo4L_M-300_8TeV-pythia6':   5.5458   * FB,
    'HPlusPlusHMinusMinusHTo4L_M-350_8TeV-pythia6':   2.6413   * FB,
    'HPlusPlusHMinusMinusHTo4L_M-400_8TeV-pythia6':   1.3414   * FB,
    'HPlusPlusHMinusMinusHTo4L_M-450_8TeV-pythia6':   0.71531  * FB,
    'HPlusPlusHMinusMinusHTo4L_M-500_8TeV-pythia6':   0.39604  * FB,
    'HPlusPlusHMinusMinusHTo4L_M-600_8TeV-pythia6':   0.13271  * FB,
    'HPlusPlusHMinusMinusHTo4L_M-700_8TeV-pythia6':   0.048382 * FB,

    # from mario, but the 4l dont match the reference... http://cms.hep.kbfi.ee/~mario/dblh/xsinfo.py
    # these were 7 tev NNLO
    #'HPlusPlusHMinusHTo3L_M-170_8TeV-calchep-pythia6': 114.838 * FB,
    #'HPlusPlusHMinusHTo3L_M-200_8TeV-calchep-pythia6': 51.59 * FB,
    #'HPlusPlusHMinusHTo3L_M-250_8TeV-calchep-pythia6': 21.172 * FB,
    #'HPlusPlusHMinusHTo3L_M-300_8TeV-calchep-pythia6': 8.308 * FB,
    #'HPlusPlusHMinusHTo3L_M-350_8TeV-calchep-pythia6': 4.154 * FB,
    #'HPlusPlusHMinusHTo3L_M-400_8TeV-calchep-pythia6': 1.9832 * FB,
    #'HPlusPlusHMinusHTo3L_M-450_8TeV-calchep-pythia6': 1.02912 * FB,
    #'HPlusPlusHMinusHTo3L_M-500_8TeV-calchep-pythia6': 0.55476 * FB,
    #'HPlusPlusHMinusHTo3L_M-600_8TeV-calchep-pythia6': 0.17286 * FB,
    #'HPlusPlusHMinusHTo3L_M-700_8TeV-calchep-pythia6': 0.05762 * FB,

    # calculated from 8 tev lo using pp k factor
    'HPlusPlusHMinusHTo3L_M-170_8TeV-calchep-pythia6': 97.73    * 1.2651 * FB,
    'HPlusPlusHMinusHTo3L_M-200_8TeV-calchep-pythia6': 50.63    * 1.2632 * FB,
    'HPlusPlusHMinusHTo3L_M-250_8TeV-calchep-pythia6': 19.79    * 1.256 * FB, # interpolated
    'HPlusPlusHMinusHTo3L_M-300_8TeV-calchep-pythia6':  8.836   * 1.2504 * FB,
    'HPlusPlusHMinusHTo3L_M-350_8TeV-calchep-pythia6':  4.292   * 1.247 * FB, # interpolated
    'HPlusPlusHMinusHTo3L_M-400_8TeV-calchep-pythia6':  2.224   * 1.2401 * FB,
    'HPlusPlusHMinusHTo3L_M-450_8TeV-calchep-pythia6':  1.206   * 1.237 * FB, # interpolated
    'HPlusPlusHMinusHTo3L_M-500_8TeV-calchep-pythia6':  0.6773  * 1.2337 * FB,
    'HPlusPlusHMinusHTo3L_M-600_8TeV-calchep-pythia6':  0.2302  * 1.223 * FB, # interpolated
    'HPlusPlusHMinusHTo3L_M-700_8TeV-calchep-pythia6':  0.08403 * 1.212 * FB, # interpolated

    'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-v1' : 36257.2 * PB,
    'WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-v2' : 36257.2 * PB,

    'W1JetsToLNu_TuneZ2Star_8TeV-madgraph' : 6662.8 * PB,
    'W2JetsToLNu_TuneZ2Star_8TeV-madgraph' : 2159.2 * PB,
    'W3JetsToLNu_TuneZ2Star_8TeV-madgraph' : 640.4 * PB,
    'W4JetsToLNu_TuneZ2Star_8TeV-madgraph' : 264. * PB,
    'WbbJetsToLNu_Massive_TuneZ2star_8TeV-madgraph-pythia6_tauola' : 138.9 * PB,

    'WGstarToLNu2E_TuneZ2star_8TeV-madgraph-tauola' : 5.87 * PB,
    'WGstarToLNu2Mu_TuneZ2star_8TeV-madgraph-tauola' : 1.91 * PB,

    'DYJetsToLL_M-10To50filter_8TeV-madgraph':          915 * PB,
    'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball': Z_XSEC,
    'DYJetsToLL_M-10To50filter_8TeV-madgraph_filtered':          915 * PB,
    'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_filtered': Z_XSEC,

    'Z1jets_M50':     Z_XSEC * 0.190169492,
    'Z2jets_M50_S10': Z_XSEC * 0.061355932,
    'Z3jets_M50':     Z_XSEC * 0.017322034,
    'Z4jets_M50':     Z_XSEC * 0.007810169,
    'Z1jets_M50_filtered':     Z_XSEC * 0.190169492,
    'Z2jets_M50_S10_filtered': Z_XSEC * 0.061355932,
    'Z3jets_M50_filtered':     Z_XSEC * 0.017322034,
    'Z4jets_M50_filtered':     Z_XSEC * 0.007810169,

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

    'GluGluToZZTo4L_8TeV-gg2zz-pythia6': 4.8 * FB,
    'GluGluToZZTo2L2L_TuneZ2star_8TeV-gg2zz-pythia6': 12.03 * FB,

    'GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6': 5.38752 * FB,

    'ZZTo4e_8TeV-powheg-pythia6':      0.07691 * PB,
    'ZZTo4mu_8TeV-powheg-pythia6':     0.07691 * PB,
    'ZZTo4tau_8TeV-powheg-pythia6':    0.07691 * PB,
    'ZZTo2e2mu_8TeV-powheg-pythia6':   0.1767 * PB,
    'ZZTo2e2tau_8TeV-powheg-pythia6':  0.1767 * PB,
    'ZZTo2mu2tau_8TeV-powheg-pythia6': 0.1767 * PB,

    'ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola':   0.1769 * PB,
    'ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola': 2.4487 * PB,
    'ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola': 0.36 * PB,

    'ZZ_TuneZ2star_8TeV_pythia6_tauola' : 8.2 * PB,
    'WZ_TuneZ2star_8TeV_pythia6_tauola' : 33.6 * PB,
    'WW_TuneZ2star_8TeV_pythia6_tauola' : 56.0 * PB,

    'ZGToLLG_8TeV-madgraph': 132.6 * PB,
    'ZGToLLG_8TeV-madgraph_filtered': 132.6 * PB,

    'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola': 1.09 * 2.207 * PB, # From CMS WZ analysis data/theo = 1.09
    'WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola':     1.09 * 1.058 * PB, # From CMS WZ analysis data/theo = 1.09

    'WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola' : 54.838*(0.1075+0.1057+0.1125)*(0.1075+0.1057+0.1125) * PB,

    'ZZZNoGstarJets' : 0.0192 * PB,
    'WWZNoGstarJets' : 0.0633 * PB,
    'WZZNoGstarJets' : 0.01968 * PB,
    'WWWJets' : 0.08217 * PB,
    'WWGJets' : 1.44 * PB,
}
