

def getSystUncertaintyMap(analysis,region,period,mode):
    if analysis=='WZ' and region=='WZ' and period==13:
        return getWZUncertainty(mode,'all')
    #if analysis=='Hpp3l' and period==8:
    #    return getHpp3lUncertainty(mode)
    #if analysis=='Hpp4l' and period==8:
    #    return getHpp4lUncertainty(mode)
    return {}

def getHpp3lUncertainty(mode):
    unc = {}
    if 'HPlusPlusH' in mode:
        unc['sig_mc_err'] = 0.15
    if mode in ['WZJets','ZZJets','ZG','TTVJets','VVVJets','TTJets','ZJets','ZG'] or 'HPlusPlusH' in mode:
        unc['lumi'] = 0.026
        unc['lep'] = 0.03
    return unc

def getHpp4lUncertainty(mode):
    unc = {}
    if 'HPlusPlusH' in mode:
        unc['sig_mc_err'] = 0.15
    if mode in ['WZJets','ZZJets','ZG','TTVJets','VVVJets','TTJets','ZJets','ZG'] or 'HPlusPlusH' in mode:
        unc['lumi'] = 0.026
        unc['lep'] = 0.04
    return unc

def getWZUncertainty(mode,chan):
    # add on flat uncertainty
    unc = {
        'dd_e'  : 0.054,
        'dd_m'  : 0.039,
        'btag'  : 0.021,
        'met'   : 0.02,
        'eff_e' : 0.019,
        'eff_m' : 0.015,
        'pileup': 0.008,
        'zz'    : 0.004,
        'pdf'   : 0.01,
        'lumi'  : 0.027,
    }
    # old way
    #unc = {}
    #if mode=='datadriven':
    #    unc['fake_rate_unc'] = 0.3
    #if mode=='ZZJets':
    #    unc['zz_xsec'] = 0.16
    #if mode=='ZG':
    #    unc['zg_xsec_theory'] = 0.06
    #if mode=='TTVJets':
    #    unc['ttv_xsec_theory'] = 0.15
    #    unc['btag'] = 0.07
    #if mode=='VVVJets':
    #    unc['vvv_xsec_theory'] = 0.06
    #if mode in ['WZJets','ZZJets','ZG','TTVJets','VVVJets']:
    #    unc['lumi'] = 0.027
    #    unc['lep_eff'] = 0.03 # total uncertainty, maybe break into electron muon and channels later
    #    unc['pu_unc'] = 0.01
    #    unc['met_unc'] = 0.02
    #if mode=='WZJets':
    #    unc['pdf_unc'] = 0.01
    #    unc['btag'] = 0.01
    return unc
