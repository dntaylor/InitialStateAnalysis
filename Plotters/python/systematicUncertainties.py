

def getSystUncertaintyMap(analysis,region,period,mode):
    if analysis=='WZ' and region=='WZ' and period==13:
        return getWZUncertainty(mode)
    return {}

def getWZUncertainty(mode):
    unc = {}
    if mode=='datadriven':
        unc['fake_rate_unc'] = 0.4
    if mode=='ZZJets':
        unc['zz_xsec'] = 0.16
    if mode=='ZG':
        unc['zg_xsec_theory'] = 0.06
    if mode=='TTVJets':
        unc['ttv_xsec_theory'] = 0.15
    if mode=='VVVJets':
        unc['vvv_xsec_theory'] = 0.06
    if mode in ['WZJets','ZZJets','ZG','TTVJets','VVVJets']:
        unc['lumi'] = 0.046
        unc['lep_eff'] = 0.016 # total uncertainty, maybe break into electron muon and channels later
        unc['pu_unc'] = 0.01
        unc['met_unc'] = 0.02
    if mode=='WZJets':
        unc['scale_unc'] = 0.043
        unc['pdf_unc'] = 0.014
    return unc
