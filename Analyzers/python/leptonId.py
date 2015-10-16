'''
Lepton ID's available in ISA

Author: Devin N. Taylor, UW-Madison
'''

import sys

sys.argv.append('b')
import ROOT as rt
sys.argv.pop()


def lep_id(rtrow, period, *lep, **kwargs):
    idType = kwargs.get('idType','')

    if idType:
        for l in lep:
            if l[0]=='e': lep_method = 'elec_id'
            if l[0]=='m': lep_method = 'muon_id'
            if l[0]=='t': lep_method = 'tau_id'
            if not eval('%s(rtrow,l,period,idType)'%lep_method): return False
        
    return True

def elec_id(rtrow, l, period, idType):
    if idType=='NonTrig':
        if not _elec_mva_nontriggering(rtrow, l, period): return False
    if idType=='Trig':
        if not _elec_mva_triggering(rtrow, l, period): return False
    if idType=='Veto':
        if not getattr(rtrow, '%sCBIDVeto' % l): return False
    if idType=='Loose':
        if not getattr(rtrow, '%sCBIDLoose' % l): return False
    if idType=='Medium':
        if not getattr(rtrow, '%sCBIDMedium' % l): return False
    if idType=='Tight':
        if not getattr(rtrow, '%sCBIDTight' % l): return False
    if idType=='VetoNoIso':
        if not getattr(rtrow, '%sCBIDVetoNoIso' % l): return False
    if idType=='LooseNoIso':
        if not getattr(rtrow, '%sCBIDLooseNoIso' % l): return False
    if idType=='MediumNoIso':
        if not getattr(rtrow, '%sCBIDMediumNoIso' % l): return False
    if idType=='TightNoIso':
        if not getattr(rtrow, '%sCBIDTightNoIso' % l): return False
    if idType=='ZZLoose':
        if not _elec_zz_loose(rtrow,l,period): return False
    if idType=='ZZTight':
        if not _elec_zz_tight(rtrow,l,period): return False
    if idType=='WZLoose':
        if not elec_WZ_loose(rtrow,l,period): return False
    if idType=='WZTight':
        if not elec_WZ_tight(rtrow,l,period): return False
    if idType=='4l':
        if not elec_4l_id(rtrow,l,period): return False
    return True

def muon_id(rtrow, l, period, idType):
    if idType=='Tight':
        if not getattr(rtrow,'%sPFIDTight'%l): return False
    if idType=='Loose':
        if not getattr(rtrow,'%sPFIDLoose'%l): return False
    if idType=='ZZLoose':
        if not _muon_zz_loose(rtrow,l,period): return False
    if idType=='ZZTight':
        if not _muon_zz_tight(rtrow,l,period): return False
    if idType=='WZLoose':
        if not muon_WZ_loose(rtrow,l,period): return False
    if idType=='WZTight':
        if not muon_WZ_tight(rtrow,l,period): return False
    if idType=='4l':
        if not muon_4l_id(rtrow,l,period): return False
    return True

def tau_id(rtrow, l, period, idType):
    if not getattr(rtrow, "%sDecayModeFinding" %l): return False  # really should be old DM, but not available in PHYS14 right now, all miniAOD pass
    if not getattr(rtrow, "%sAgainstElectronMediumMVA5" % l): return False
    if not getattr(rtrow, "%sAgainstMuonTight3" % l): return False
    if idType=='Loose':
        if not getattr(rtrow, "%sByLooseCombinedIsolationDeltaBetaCorr3Hits" %l): return False
    if idType=='Medium':
        if not getattr(rtrow, "%sByMediumCombinedIsolationDeltaBetaCorr3Hits" %l): return False
    if idType=='Tight':
        if not getattr(rtrow, "%sByTightCombinedIsolationDeltaBetaCorr3Hits" %l):return False
    return True

def _muon_zz_loose(rtrow, l, period):
    if getattr(rtrow, "%sPt" % l) < 5: return False
    if abs(getattr(rtrow, "%sEta" % l)) > 2.4: return False
    if abs(getattr(rtrow, '%sPVDZ' % l)) > 1.: return False
    if abs(getattr(rtrow, '%sPVDXY' % l)) > 0.5:return False
    if abs(getattr(rtrow, '%sSIP3D' % l)) > 4.:return False
    if getattr(rtrow, '%sBestTrackType' % l)==2: return False
    isGlobal = getattr(rtrow, '%sIsGlobal' % l)
    isTracker = getattr(rtrow, '%sIsTracker' % l)
    matchedStations = getattr(rtrow, '%sMatchedStations' % l)
    return isGlobal or (isTracker and matchedStations>0)

def _muon_zz_tight(rtrow, l, period):
    if not _muon_zz_loose(rtrow,l,period): return False
    return getattr(rtrow,'%sIsPFMuon' %l)

def _elec_zz_loose(rtrow, l, period):
    if getattr(rtrow, "%sPt" % l) < 7: return False
    if abs(getattr(rtrow, "%sEta" % l)) > 2.5: return False
    if abs(getattr(rtrow, '%sPVDZ' % l)) > 1.: return False
    if abs(getattr(rtrow, '%sPVDXY' % l)) > 0.5:return False
    if getattr(rtrow,'%sMissingHits'%l) > 1: return False
    if abs(getattr(rtrow, '%sSIP3D' % l)) > 4.:return False
    return True

def _elec_zz_tight(rtrow, l, period):
    if not _elec_zz_loose(rtrow,l,period): return False
    return _elec_mva_nontriggering_zz(rtrow,l,period)

def _elec_mva_nontriggering_zz(rtrow, l, period):
    pt = getattr(rtrow, "%sPt" % l)
    eta = abs(getattr(rtrow, "%sSCEta" % l))
    mva = getattr(rtrow, "%sMVANonTrigID" % l) if period == 13 else getattr(rtrow, "%sMVANonTrig" % l)

    if 5.0 < pt < 10.0:
        return (eta < 0.8 and mva > -0.265) or (0.8 < eta < 1.479 and mva > -0.556) or (1.479 < eta and mva > -0.551)
    elif 10.0 < pt:
        return (eta < 0.8 and mva > -0.072) or (0.8 < eta < 1.479 and mva > -0.286) or (1.479 < eta and mva > -0.267)
    else:
        return False

def _elec_mva_nontriggering(rtrow, l, period):
    pt = getattr(rtrow, "%sPt" % l)
    eta = abs(getattr(rtrow, "%sSCEta" % l))
    mva = getattr(rtrow, "%sMVANonTrigID" % l) if period == 13 else getattr(rtrow, "%sMVANonTrig" % l)

    if 5.0 < pt < 10.0:
        return (eta < 0.8 and mva > 0.47) or (0.8 < eta < 1.479 and mva > 0.004) or (1.479 < eta and mva > 0.295)

    elif 10.0 < pt:
        return (eta < 0.8 and mva > -0.34) or (0.8 < eta < 1.479 and mva > -0.65) or (1.479 < eta and mva > 0.6)

    else:
        return False

def _elec_mva_triggering(rtrow, l, period):
    pt = getattr(rtrow, "%sPt" % l)
    eta = abs(getattr(rtrow, "%sSCEta" % l))
    mva = getattr(rtrow, "%sMVATrigID" % l) if period == 13 else getattr(rtrow, "%sMVATrig" % l)

    if 10.0 < pt < 20.0:
        return (eta < 0.8 and mva > 0.00) or (0.8 < eta < 1.479 and mva > 0.10) or (1.479 < eta and mva > 0.62)

    elif 20.0 < pt:
        return (eta < 0.8 and mva > 0.94) or (0.8 < eta < 1.479 and mva > 0.85) or (1.479 < eta and mva > 0.92)

    else:
        return False

########################
### Old WZ 8 TeV IDs ###
########################

def elec_WZ_loose(rtrow, l, period):
    pt = getattr(rtrow, "%sPt" % l)
    eta = abs(getattr(rtrow, "%sEta" % l))
    sceta = abs(getattr(rtrow, "%sSCEta" % l))
    sieie = getattr(rtrow, "%sSigmaIEtaIEta" % l)
    dphi = getattr(rtrow, "%sdeltaPhiSuperClusterTrackAtVtx" %l)
    deta = getattr(rtrow, "%sdeltaEtaSuperClusterTrackAtVtx" %l)
    hoe = getattr(rtrow, "%sHadronicOverEM" %l)
    eiso = getattr(rtrow, "%sEcalIsoDR03" %l)
    hiso = getattr(rtrow, "%sHcalIsoDR03" %l)
    tiso = getattr(rtrow, "%sTrkIsoDR03" %l)
    conv = getattr(rtrow, "%sHasConversion" %l)
    misshits = getattr(rtrow, "%sMissingHits" %l)

    passid = True
    if pt < 10: passid = False
    if eta > 2.5: passid = False
    if sceta < 1.479:
        if sieie > 0.01: passid = False
        if dphi > 0.15: passid = False
        if deta > 0.007: passid = False
        if hoe > 0.12: passid = False
        if max(eiso-1,0)/pt > 0.2: passid = False
    if sceta >= 1.479:
        if sieie > 0.03: passid = False
        if dphi > 0.1: passid = False
        if deta > 0.009: passid = False
        if hoe > 0.1: passid = False
        if eiso/pt > 0.2: passid = False
    if tiso/pt > 0.2: passid = False
    if hiso/pt > 0.2: passid = False
    if conv: passid = False
    if misshits: passid = False

    return passid

def elec_WZ_tight(rtrow, l, period):
    d0 = getattr(rtrow, "%sPVDXY" %l)
    dz = getattr(rtrow, "%sPVDZ" %l)
    mva = _elec_mva_triggering(rtrow, l, period)
    reliso = getattr(rtrow, "%sRelPFIsoRho" %l)

    wzloose = elec_WZ_loose(rtrow, l, period)

    #chgId = getattr(rtrow,'%sChargeIdTight' %l)
    chgId = 1.

    return reliso < 0.15 and d0 < 0.02 and dz < 0.1 and wzloose and mva and chgId

def muon_WZ_loose(rtrow, l, period):
    pt = getattr(rtrow, "%sPt" % l)
    eta = abs(getattr(rtrow, "%sEta" % l))
    #d0 = getattr(rtrow, "%sPVDXY" %l)
    #dz = getattr(rtrow, "%sPVDZ" %l)
    tightid = getattr(rtrow, '%sPFIDTight' %l)
    reliso = getattr(rtrow, "%sRelPFIsoDBDefault" %l)

    passid = True
    if pt < 10: passid = False
    if eta > 2.4: passid = False
    if not tightid: passid = False
    if reliso > 0.2: passid = False

    return passid

def muon_WZ_tight(rtrow, l, period):
    loose = muon_WZ_loose(rtrow, l, period)
    reliso = getattr(rtrow, "%sRelPFIsoDBDefault" %l)

    passid = loose
    if reliso > 0.12: passid = False

    return passid

####################
### 8 TeV 4l IDs ###
####################

def muon_4l_id(rtrow,l,period):
    dz = getattr(rtrow, "%sPVDZ" % l) < 1.0
    dxy = getattr(rtrow, "%sPVDXY" % l) < 0.5
    sip = getattr(rtrow, "%sIP3DS" % l) < 4.0
    mu_type = getattr(rtrow, "%sIsTracker" % l) or getattr(rtrow, "%sIsGlobal" % l)
    return all([dz, dxy, sip, mu_type])

def elec_4l_id(rtrow, l, period):
    dz = getattr(rtrow, "%sPVDZ" % l) < 1.0
    dxy = getattr(rtrow, "%sPVDXY" % l) < 0.5
    sip = getattr(rtrow, "%sIP3DS" % l) < 4.0
    nhit = getattr(rtrow, "%sMissingHits" % l) <= 1
    mva = _elec_mva_nontriggering(rtrow, l, period)
    return all([dz, dxy, sip, nhit, mva])
