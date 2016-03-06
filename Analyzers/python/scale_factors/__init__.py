import sys
import os
import glob
import pickle
import json
import csv
from operator import itemgetter, attrgetter

sys.argv.append('-b')
import ROOT as rt
sys.argv.pop()

class ChargeIdSystematics(object):

    def __init__(self):
        self.missid = {
            'EB': {
                10 : 0.00350,
                20 : 0.00331,
                30 : 0.00319,
                40 : 0.00289,
                60 : 0.00581,
                100: 0.00758,
            },
            'EE': {
                10 : 0.02700,
                20 : 0.02717,
                30 : 0.01937,
                40 : 0.01664,
                60 : 0.04597,
                100: 0.06232,
            }
        }



    def systematic(self, rtrow, *leps, **kwargs):
        val = 1.
        for l in leps:
            if l[0]=='e':
                pt = getattr(rtrow,'%sPt' %l)
                eta = getattr(rtrow,'%sEta' %l)
                pt_reg = [x for x in [10,20,30,40,60,100,200] if pt>=x][0]
                eta_reg = 'EB' if abs(eta)<1.479 else 'EE'
                val *= (1.-self.missid[eta_reg][pt_reg])
        return 1+(1-val)
                

class TriggerScaleFactors(object):

    def __init__(self):
        self.ww_scales = self.init_ww_scales()
        # WZ 13 TeV
        with open(os.path.join(os.path.dirname(__file__),'muons_13TeV.json'),'r') as ef:
            self.muons_13TeV = json.load(ef)
        with open(os.path.join(os.path.dirname(__file__),'electrons_13TeV.json'),'r') as ef:
            self.electrons_13TeV = json.load(ef)
        with open(os.path.join(os.path.dirname(__file__),'electronTrigger_13TeV.json'),'r') as ef:
            self.electronTrigger_13TeV = json.load(ef)
        with open(os.path.join(os.path.dirname(__file__),'trackerMuonDZ_13TeV.json'),'r') as ef:
            self.trackerMuonDZ_13TeV = json.load(ef)
        with open(os.path.join(os.path.dirname(__file__),'globalMuonDZ_13TeV.json'),'r') as ef:
            self.globalMuonDZ_13TeV = json.load(ef)
        # HWW 13 TeV
        with open(os.path.join(os.path.dirname(__file__),'HLT_Ele23_WPLoose.txt'),'r') as ef:
            self.hww_singleEle23 = self.init_hww_scales_e(ef)
        with open(os.path.join(os.path.dirname(__file__),'SingleMu_IsoTkMu20_Run2015D_25ns_PTvsETA_HWW.txt'),'r') as ef:
            self.hww_singleMu20 = self.init_hww_scales_m(ef)

    def init_hww_scales_e(self,ef):
        result = []
        for line in ef:
            etalow, etahigh, ptlow, pthigh, eff, err = line.split()
            element = {}
            element["data"]     = float(eff)
            element["data_err"] = float(err)
            element["pt_lo"]    = float(ptlow)
            element["pt_hi"]    = float(pthigh)
            element["eta_lo"]   = float(etalow)
            element["eta_hi"]   = float(etahigh)
            result += [element]
        return result

    def init_hww_scales_m(self,ef):
        result = []
        for line in ef:
            etalow, etahigh, ptlow, pthigh, eff, errup, errdown = line.split()
            element = {}
            element["data"]     = float(eff)
            element["data_err"] = (float(errup)+float(errdown))/2.
            element["pt_lo"]    = float(ptlow)
            element["pt_hi"]    = float(pthigh)
            element["eta_lo"]   = float(etalow)
            element["eta_hi"]   = float(etahigh)
            result += [element]
        return result

    def init_ww_scales(self):
        scales = {}
        with open(os.path.join(os.path.dirname(__file__),'WW_140416_TriggerEfficiencies.txt'),'r') as file:
            for line in file.readlines()[1:]:
                line.rstrip()
                leg, etalow, etahigh, ptlow, pthigh, eff, errdown, errup = line.split()
                if leg not in scales:
                    scales[leg] = []
                scales[leg].append([float(etalow),float(etahigh),float(ptlow),float(pthigh),float(eff),float(errdown),float(errup)])
        return scales

    def scale_factor(self, rtrow, *lep_list, **kwargs):
        shift = kwargs.get('metShift','')
        def getObjPt(l):
            ptString = '{0}Pt'.format(l)
            if l[0]=='m' and shift=='mes+': ptString += '_MuonEnUp'
            if l[0]=='m' and shift=='mes-': ptString += '_MuonEnDown'
            if l[0]=='e' and shift=='ees+': ptString += '_ElectronEnUp'
            if l[0]=='e' and shift=='ees-': ptString += '_ElectronEnDown'
            return getattr(rtrow,ptString)
        def getObjEta(l):
            if l[0]=='m': return getattr(rtrow,'{0}Eta'.format(l))
            if l[0]=='e': return getattr(rtrow,'{0}SCEta'.format(l))
        lep_objs = [(x, getObjPt(x), getObjEta(x)) for x in lep_list]
        lep_ord = sorted(lep_objs, key=itemgetter(1), reverse=True)
        period = kwargs.pop('period',8)
        if period==8:
            if len(lep_ord)==3:
                eff = 1-(\
                        (1-self.double_lead_eff(*lep_ord[0])) * (1-self.double_lead_eff(*lep_ord[1])) * (1-self.double_lead_eff(*lep_ord[2]))\
                        + self.double_lead_eff(*lep_ord[0]) * (1-self.double_trail_eff(*lep_ord[1])) * (1-self.double_trail_eff(*lep_ord[2]))\
                        + self.double_lead_eff(*lep_ord[1]) * (1-self.double_trail_eff(*lep_ord[2])) * (1-self.double_trail_eff(*lep_ord[0]))\
                        + self.double_lead_eff(*lep_ord[2]) * (1-self.double_trail_eff(*lep_ord[0])) * (1-self.double_trail_eff(*lep_ord[1]))\
                        )
            elif len(lep_ord) == 4:
                eff = 1 - (\
                        (1 - self.double_lead_eff(*lep_ord[0])) * \
                        (1 - self.double_lead_eff(*lep_ord[1])) * \
                        (1 - self.double_lead_eff(*lep_ord[2])) * \
                        (1 - self.double_lead_eff(*lep_ord[3])) \

                      + self.double_lead_eff(*lep_ord[0]) * \
                        (1 - self.double_trail_eff(*lep_ord[1])) * \
                        (1 - self.double_trail_eff(*lep_ord[2])) * \
                        (1 - self.double_trail_eff(*lep_ord[3])) \

                      + self.double_lead_eff(*lep_ord[1]) * \
                        (1 - self.double_trail_eff(*lep_ord[2])) * \
                        (1 - self.double_trail_eff(*lep_ord[3])) * \
                        (1 - self.double_trail_eff(*lep_ord[0])) \

                      + self.double_lead_eff(*lep_ord[2]) * \
                        (1 - self.double_trail_eff(*lep_ord[3])) * \
                        (1 - self.double_trail_eff(*lep_ord[0])) * \
                        (1 - self.double_trail_eff(*lep_ord[1])) \

                      + self.double_lead_eff(*lep_ord[3]) * \
                        (1 - self.double_trail_eff(*lep_ord[0])) * \
                        (1 - self.double_trail_eff(*lep_ord[1])) * \
                        (1 - self.double_trail_eff(*lep_ord[2])) \
                        )
            elif len(lep_ord)==2:
                eff = 1-(\
                        (1-self.double_lead_eff(*lep_ord[0])) * (1-self.double_lead_eff(*lep_ord[1]))\
                        + self.double_lead_eff(*lep_ord[0]) * (1-self.double_trail_eff(*lep_ord[1]))\
                        + self.double_lead_eff(*lep_ord[1]) * (1-self.double_trail_eff(*lep_ord[0]))\
                        )
                #eff = self.double_lead_eff(*lep_ord[0]) * self.double_trail_eff(*lep_ord[1])
            else:
                eff = 1
        else:
            if len(lep_ord)==3:
                # single * double efficiency
                eff = 1-(\
                        # none pass single
                        ((1-self.single_eff_13(*lep_ord[0],**kwargs)) * (1-self.single_eff_13(*lep_ord[1],**kwargs)) * (1-self.single_eff_13(*lep_ord[2],**kwargs)))\
                        # none pass lead
                        *((1-self.double_lead_eff_13(*lep_ord[0],**kwargs)) * (1-self.double_lead_eff_13(*lep_ord[1],**kwargs)) * (1-self.double_lead_eff_13(*lep_ord[2],**kwargs))\
                        # one pass lead but none pass trail
                        + self.double_lead_eff_13(*lep_ord[0],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[1],**kwargs)) * (1-self.double_trail_eff_13(*lep_ord[2],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[1],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[2],**kwargs)) * (1-self.double_trail_eff_13(*lep_ord[0],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[2],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[0],**kwargs)) * (1-self.double_trail_eff_13(*lep_ord[1],**kwargs))\
                        # one pass lead, one pass trail, one fail trail, fail dz
                        + self.double_lead_eff_13(*lep_ord[0],**kwargs) * self.double_trail_eff_13(*lep_ord[1],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[2])) * (1-self.double_dz_eff_13(lep_ord[0][0],*lep_ord[1],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[0],**kwargs) * self.double_trail_eff_13(*lep_ord[2],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[1])) * (1-self.double_dz_eff_13(lep_ord[0][0],*lep_ord[2],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[1],**kwargs) * self.double_trail_eff_13(*lep_ord[2],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[0])) * (1-self.double_dz_eff_13(lep_ord[1][0],*lep_ord[2],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[1],**kwargs) * self.double_trail_eff_13(*lep_ord[0],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[2])) * (1-self.double_dz_eff_13(lep_ord[1][0],*lep_ord[0],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[2],**kwargs) * self.double_trail_eff_13(*lep_ord[0],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[1])) * (1-self.double_dz_eff_13(lep_ord[2][0],*lep_ord[0],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[2],**kwargs) * self.double_trail_eff_13(*lep_ord[1],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[0])) * (1-self.double_dz_eff_13(lep_ord[2][0],*lep_ord[1],**kwargs))\
                        # one pass lead, two pass trail, both fail dz
                        + self.double_lead_eff_13(*lep_ord[0],**kwargs) * self.double_trail_eff_13(*lep_ord[1],**kwargs) * self.double_trail_eff_13(*lep_ord[2])\
                          * (1-self.double_dz_eff_13(lep_ord[0][0],*lep_ord[1],**kwargs)) * (1-self.double_dz_eff_13(lep_ord[0][0],*lep_ord[2],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[1],**kwargs) * self.double_trail_eff_13(*lep_ord[2],**kwargs) * self.double_trail_eff_13(*lep_ord[0])\
                          * (1-self.double_dz_eff_13(lep_ord[1][0],*lep_ord[2],**kwargs)) * (1-self.double_dz_eff_13(lep_ord[1][0],*lep_ord[0],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[2],**kwargs) * self.double_trail_eff_13(*lep_ord[0],**kwargs) * self.double_trail_eff_13(*lep_ord[1])\
                          * (1-self.double_dz_eff_13(lep_ord[2][0],*lep_ord[0],**kwargs)) * (1-self.double_dz_eff_13(lep_ord[2][0],*lep_ord[1],**kwargs)))\
                        )
                #eff = 1. # override, i dont have single lepton efficiencies
            elif len(lep_ord)==2: # for or of double lepton triggers
                eff = 1-(\
                        # none pass lead
                        (1-self.double_lead_eff_13(*lep_ord[0],**kwargs)) * (1-self.double_lead_eff_13(*lep_ord[1],**kwargs))\
                        # one pass lead but none pass trail
                        + self.double_lead_eff_13(*lep_ord[0],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[1],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[1],**kwargs) * (1-self.double_trail_eff_13(*lep_ord[0],**kwargs))\
                        # one pass lead, one pass trail, fail dz
                        + self.double_lead_eff_13(*lep_ord[0],**kwargs) * self.double_trail_eff_13(*lep_ord[1],**kwargs) * (1-self.double_dz_eff_13(lep_ord[0][0],*lep_ord[1],**kwargs))\
                        + self.double_lead_eff_13(*lep_ord[1],**kwargs) * self.double_trail_eff_13(*lep_ord[0],**kwargs) * (1-self.double_dz_eff_13(lep_ord[1][0],*lep_ord[0],**kwargs))\
                        )
            elif len(lep_ord)==1: # just get the trail leg one
                eff = self.double_trail_eff_13(*lep_ord[0],**kwargs) if lep_ord[0][0]=='e' or lep_ord[0][1]<20 else self.double_lead_eff_13(*lep_ord[0],**kwargs)
            else:
                eff = 1.
        return eff

    def single_eff_13(self,l,pt,eta,**kwargs):
        if l[0]=='e': return self.single_e_13(pt,eta,**kwargs)
        if l[0]=='m': return self.single_m_13(pt,eta,**kwargs)
        return 1.

    def double_lead_eff_13(self,l,pt,eta,**kwargs):
        if l[0]=='e': return self.double_lead_e_13(pt,eta,**kwargs)
        if l[0]=='m': return self.double_lead_m_13(pt,eta,**kwargs)
        return 1.

    def double_trail_eff_13(self,l,pt,eta,**kwargs):
        if l[0]=='e': return self.double_trail_e_13(pt,eta,**kwargs)
        if l[0]=='m': return self.double_trail_m_13(pt,eta,**kwargs)
        return 1.

    def double_dz_eff_13(self,l0,l1,pt,eta,**kwargs):
        if l0[0]=='e' and l1[0]=='e': return self.double_dz_e_13(pt,eta,**kwargs)
        if l0[0]=='m' and l1[0]=='m': return self.double_dz_m_13(pt,eta,**kwargs)
        return 1.

    # electron
    def single_e_13(self,pt,eta,**kwargs):
        shiftUp = kwargs.pop('shiftUp',False)
        shiftDown = kwargs.pop('shiftDown',False)
        val,err = self.get_eff_err_mod(pt,eta,self.hww_singleEle23,**kwargs)
        if shiftUp: return val+err
        if shiftDown: return val-err
        return val

    def double_lead_e_13(self,pt,eta,**kwargs):
        shiftUp = kwargs.pop('shiftUp',False)
        shiftDown = kwargs.pop('shiftDown',False)
        val,err = self.get_eff_err(pt,eta,self.electronTrigger_13TeV,'passingHLTEle17Ele12Leg1',**kwargs)
        if shiftUp: return val+err
        if shiftDown: return val-err
        return val

    def double_trail_e_13(self,pt,eta,**kwargs):
        shiftUp = kwargs.pop('shiftUp',False)
        shiftDown = kwargs.pop('shiftDown',False)
        val,err = self.get_eff_err(pt,eta,self.electronTrigger_13TeV,'passingHLTEle17Ele12Leg2',**kwargs)
        if shiftUp: return val+err
        if shiftDown: return val-err
        return val

    def double_dz_e_13(self,pt,eta,**kwargs):
        shiftUp = kwargs.pop('shiftUp',False)
        shiftDown = kwargs.pop('shiftDown',False)
        val,err = self.get_eff_err(pt,eta,self.electronTrigger_13TeV,'passingHLTDZFilter',**kwargs)
        if shiftUp: return val+err
        if shiftDown: return val-err
        return val

    # muon
    def single_m_13(self,pt,eta,**kwargs):
        shiftUp = kwargs.pop('shiftUp',False)
        shiftDown = kwargs.pop('shiftDown',False)
        val,err = self.get_eff_err_mod(pt,eta,self.hww_singleMu20,**kwargs)
        if shiftUp: return val+err
        if shiftDown: return val-err
        return val

    def double_lead_m_13(self,pt,eta,**kwargs):
        shiftUp = kwargs.pop('shiftUp',False)
        shiftDown = kwargs.pop('shiftDown',False)
        val,err = self.get_eff_err(pt,eta,self.muons_13TeV,'passingMu17',**kwargs)
        if shiftUp: return val+err
        if shiftDown: return val-err
        return val

    def double_trail_m_13(self,pt,eta,**kwargs):
        shiftUp = kwargs.pop('shiftUp',False)
        shiftDown = kwargs.pop('shiftDown',False)
        valMu,errMu = self.get_eff_err(pt,eta,self.muons_13TeV,'passingMu8',**kwargs)
        valTkMu,errTkMu = self.get_eff_err(pt,eta,self.muons_13TeV,'passingTkMu8',**kwargs)
        val = (valMu+valTkMu)/2.
        valUp = (valMu+errMu+valTkMu+errTkMu)/2.
        valDown = (valMu-errMu+valTkMu-errTkMu)/2.
        if shiftUp: return valUp
        if shiftDown: return valDown
        return val

    def double_dz_m_13(self,pt,eta,**kwargs):
        shiftUp = kwargs.pop('shiftUp',False)
        shiftDown = kwargs.pop('shiftDown',False)
        valMu,errMu = self.get_eff_err(pt,eta,self.globalMuonDZ_13TeV,'passingDZ',**kwargs)
        valTkMu,errTkMu = self.get_eff_err(pt,eta,self.trackerMuonDZ_13TeV,'passingDZ',**kwargs)
        val = (valMu+valTkMu)/2.
        valUp = (valMu+errMu+valTkMu+errTkMu)/2.
        valDown = (valMu-errMu+valTkMu-errTkMu)/2.
        if shiftUp: return valUp
        if shiftDown: return valDown
        return val

    def get_eff_err(self,pt,eta,effDict,effKey,**kwargs):
        useData = kwargs.pop('useData',True)
        efflist = effDict[effKey]
        for eff in efflist:
            ptlow = eff['pt_lo']
            pthi = eff['pt_hi']
            etalow = eff['abseta_lo']
            etahi = eff['abseta_hi']
            if pt>=ptlow and pt<pthi and abs(eta)>=etalow and abs(eta)<etahi:
                val = eff['data'] if useData else eff['mc']
                err = eff['data_err'] if useData else eff['mc_err']
                return val, err
        return 1.0, 0.0

    def get_eff_err_mod(self,pt,eta,efflist,**kwargs):
        useData = kwargs.pop('useData',True)
        for eff in efflist:
            ptlow = eff['pt_lo']
            pthi = eff['pt_hi']
            etalow = eff['eta_lo']
            etahi = eff['eta_hi']
            if pt>=ptlow and pt<pthi and eta>=etalow and eta<etahi:
                val = eff['data']
                err = eff['data_err']
                return val, err
        return 1.0, 0.0


    # 8TeV stuff
    def single_eff(self,l,pt,eta):
        if l[0]=='e': return self.single_e(pt,eta)
        if l[0]=='m': return self.single_m(pt,eta)
        return 1.

    def double_lead_eff(self,l,pt,eta):
        if l[0]=='e': return self.double_lead_e(pt,eta)
        if l[0]=='m': return self.double_lead_m(pt,eta)
        return 1.

    def double_trail_eff(self,l,pt,eta):
        if l[0]=='e': return self.double_trail_e(pt,eta)
        if l[0]=='m': return self.double_trail_m(pt,eta)
        return 1.

    def get_eff(self,leg,pt,eta):
        for etalow, etahigh, ptlow, pthigh, eff, errdown, errup in self.ww_scales[leg]:
            if eta>=etalow and eta<etahigh and pt>=ptlow and pt<pthigh:
                return eff
        return 1.

    def single_e(self,pt,eta):
        return self.get_eff('SingleEl',pt,eta)
        pts = [10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,35.,40.,50.,7000.]
        if eta < 1.5:
           effs = [0.,0.,0.,0.,0.,0.,0.05,0.62,0.9,0.91,0.91,0.91]
        elif eta < 2.5:
           effs = [0.,0.,0.,0.,0.,0.02,0.11,0.43,0.70,0.74,0.74,0.74]
        else: return 1.0
        for p in range(len(pts)-1):
            if pt<pts[p+1]: return effs[p]
        return 1.0

    def single_m(self,pt,eta):
        return self.get_eff('SingleMu',pt,eta)
        pts = [10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,35.,40.,50.,7000.]
        if eta < 0.8:
            effs = [0.0, 0.0, 0.0, 0.0, 0.01, 0.5, 0.9, 0.91, 0.91, 0.91, 0.94, 0.94]
        elif eta < 1.2:
            effs = [0.0, 0.0, 0.0, 0.0, 0.01, 0.46, 0.82, 0.83, 0.83, 0.83, 0.85, 0.85]
        elif eta < 2.1:
            effs = [0.0, 0.0, 0.0, 0.0, 0.01, 0.48, 0.79, 0.8, 0.8, 0.8, 0.82, 0.82]
        elif eta < 2.5:
            effs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        else: return 1.0
        for p in range(len(pts)-1):
            if pt<pts[p+1]: return effs[p]
        return 1.0

    def double_lead_e(self,pt,eta):
        return self.get_eff('DoubleElLead',pt,eta)
        pts = [10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,35.,40.,50.,7000.]
        if eta < 1.5:
            effs = [0.0, 0.0, 0.0, 0.75, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97, 0.97]
        elif eta < 2.5:
            effs = [0.0, 0.01, 0.16, 0.7, 0.95, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98, 0.98]
        else: return 1.0
        for p in range(len(pts)-1):
            if pt<pts[p+1]: return effs[p]
        return 1.0

    def double_trail_e(self,pt,eta):
        return self.get_eff('DoubleElTrail',pt,eta)
        pts = [10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,35.,40.,50.,7000.]
        if eta < 1.5:
            effs = [0.77, 0.84, 0.88, 0.91, 0.93, 0.94, 0.95, 0.95, 0.96, 0.96, 0.96, 0.96]
        elif eta < 2.5:
            effs = [0.56, 0.67, 0.76, 0.83, 0.88, 0.91, 0.93, 0.94, 0.95, 0.96, 0.96, 0.96]
        else: return 1.0
        for p in range(len(pts)-1):
            if pt<pts[p+1]: return effs[p]
        return 1.0

    def double_lead_m(self,pt,eta):
        return self.get_eff('DoubleMuLead',pt,eta)
        pts = [10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,35.,40.,50.,7000.]
        if eta < 1.2:
            effs = [0.0, 0.0, 0.1, 0.93, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94, 0.94]
        elif eta < 2.1:
            effs = [0.0, 0.0, 0.1, 0.89, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90]
        elif eta < 2.5:
            effs = [0.0, 0.0, 0.17, 0.84, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86, 0.86]
        else: return 1.0
        for p in range(len(pts)-1):
            if pt<pts[p+1]: return effs[p]
        return 1.0

    def double_trail_m(self,pt,eta):
        return self.get_eff('DoubleMuTrail',pt,eta)
        pts = [10.,12.5,15.,17.5,20.,22.5,25.,27.5,30.,35.,40.,50.,7000.]
        if eta < 1.2:
            effs = [0.94, 0.95, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96, 0.96]
        elif eta < 2.1:
            effs = [0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92, 0.92]
        elif eta < 2.5:
            effs = [0.88, 0.89, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90, 0.90]
        else: return 1.0
        for p in range(len(pts)-1):
            if pt<pts[p+1]: return effs[p]
        return 1.0

class LeptonEfficiency(object):

    def __init__(self):
        # WZ 13TeV
        with open(os.path.join(os.path.dirname(__file__),'muons_13TeV.json'),'r') as mf:
            self.m_id_dict_13tev = json.load(mf)
        with open(os.path.join(os.path.dirname(__file__),'electrons_13TeV.json'),'r') as ef:
            self.e_id_dict_13tev = json.load(ef)
        # MIT efficiencies
        efilename = os.path.join(os.path.dirname(__file__),'SingleElectron_efficiencies_electronTnP.root')
        self.eEffTfile = rt.TFile.Open(efilename,'READ')
        self.eEffHist_loose = self.eEffTfile.Get('eff_Loose_ele')
        self.eEffHist_medium = self.eEffTfile.Get('eff_Medium_ele')
        self.eEffHist_tight = self.eEffTfile.Get('eff_Tight_ele')
        mfilename = os.path.join(os.path.dirname(__file__),'SingleMuon_efficiencies_muonTnP.root')
        self.mEffTfile = rt.TFile.Open(mfilename,'READ')
        self.mEffHist_loose = self.mEffTfile.Get('eff_Loose_mu')
        self.mEffHist_medium = self.mEffTfile.Get('eff_Medium_mu')
        self.mEffHist_tight = self.mEffTfile.Get('eff_Medium_mu')


    def close(self):
        self.eEffTfile.Close()
        self.mEffTfile.Close()

    def scale_factor(self, row, *lep_list, **kwargs):
        period = kwargs.pop('period',13)
        shift = kwargs.pop('metShift','')
        numer = kwargs.pop('numer','Tight')
        denom = kwargs.pop('denom','Loose')
        out = []
        for l in lep_list:
            lep_type = l[0]

            if lep_type == 'm':
                out += [self.getMuonEfficiency(l,row,shift,numer,denom)]
            elif lep_type == 'e':
                out += [self.getElectronEfficiency(l,row,shift,numer,denom)]
            elif lep_type == 't':
                out += [[1,1,1]] # TODO
            else:
                raise TypeError("Lepton type %s not recognized" % lep_type)

        final = [1,1,1]
        for o in out:
            final[0] *= o[0]
            final[1] *= o[1]
            final[2] *= o[2]

        return final

    def getMuonEfficiency(self,l,row,shift,numer,denom):
        if shift in ['mes+','mes-']:
            shiftName = 'MuonEnUp' if shift=='mes+' else 'MuonEnDown'
            pt = getattr(row,'{0}Pt_{1}'.format(l,shiftName))
        else:
            pt = getattr(row,'{0}Pt'.format(l))
        eta = abs(getattr(row,'{0}Eta'.format(l)))
        hists = {
            'Loose' : self.mEffHist_loose,
            'Medium' : self.mEffHist_medium,
            'Tight' : self.mEffHist_tight,
        }
        nid = self.get_eff_err_new(pt,eta,hists[numer])
        did = self.get_eff_err_new(pt,eta,hists[denom])
        default = nid[0]/did[0]
        up = (nid[0]+nid[1])/(did[0]+did[1])
        down = (nid[0]-nid[1])/(did[0]-did[1])
        # loose to tight efficiency = tight/loose
        #tid  = self.get_eff_err(pt,eta,self.m_id_dict_13tev,'passingIDWZTight')
        #tiso = self.get_eff_err(pt,eta,self.m_id_dict_13tev,'passingIsoWZTight_passingIDWZTight')
        #lid  = self.get_eff_err(pt,eta,self.m_id_dict_13tev,'passingIDWZLoose')
        #liso = self.get_eff_err(pt,eta,self.m_id_dict_13tev,'passingIsoWZLoose_passingIDWZLoose')
        #default = (tid[0]*tiso[0])/(lid[0]*liso[0])
        #up = ((tid[0]+tid[1])*(tiso[0]+tiso[1]))/((lid[0]+lid[1])*(liso[0]+liso[1]))
        #down = ((tid[0]-tid[1])*(tiso[0]-tiso[1]))/((lid[0]-lid[1])*(liso[0]-liso[1]))
        return [default, up, down]

    def getElectronEfficiency(self,l,row,shift,numer,denom):
        if shift in ['ees+','ees-']:
            shiftName = 'ElectronEnUp' if shift=='ees+' else 'ElectronEnDown'
            pt = getattr(row,'{0}Pt_{1}'.format(l,shiftName))
        else:
            pt = getattr(row,'{0}Pt'.format(l))
        eta = abs(getattr(row,'{0}SCEta'.format(l)))
        hists = {
            'Loose' : self.eEffHist_loose,
            'Medium' : self.eEffHist_medium,
            'Tight' : self.eEffHist_tight,
        }
        nid = self.get_eff_err_new(pt,eta,hists[numer])
        did = self.get_eff_err_new(pt,eta,hists[denom])
        default = nid[0]/did[0]
        up = (nid[0]+nid[1])/(did[0]+did[1])
        down = (nid[0]-nid[1])/(did[0]-did[1])
        # loose to tight efficiency = tight/loose
        #tid = self.get_eff_err(pt,eta,self.e_id_dict_13tev,'passingMedium')
        #lid = self.get_eff_err(pt,eta,self.e_id_dict_13tev,'passingLoose')
        #default = (tid[0])/(lid[0])
        #up = (tid[0]+tid[1])/(lid[0]+lid[1])
        #down = (tid[0]-tid[1])/(lid[0]-lid[1])
        return [default, up, down]

    def get_eff_err_new(self,pt,eta,hist,**kwargs):
        if eta>2.4: eta = 2.39 # TODO: stupid hack for electrons!!!
        if pt<10: pt = 10.1
        val = hist.GetBinContent(hist.FindBin(eta,pt))
        err = hist.GetBinError(hist.FindBin(eta,pt))
        return val,err

    def get_eff_err(self,pt,eta,effDict,effKey,**kwargs):
        useData = kwargs.pop('useData',True)
        efflist = effDict[effKey]
        for eff in efflist:
            ptlow = eff['pt_lo']
            pthi = eff['pt_hi']
            etalow = eff['abseta_lo']
            etahi = eff['abseta_hi']
            if pt>=ptlow and pt<pthi and eta>=etalow and eta<etahi:
                val = eff['data'] if useData else eff['mc']
                err = eff['data_err'] if useData else eff['mc_err']
                return val, err
        return 1.0, 0.0


class LeptonFakeRate(object):

    def __init__(self):
        # WZ 13TeV
        #with open(os.path.join(os.path.dirname(__file__),'fakes.json'),'r') as f:
        #with open(os.path.join(os.path.dirname(__file__),'fakes_trigIso_dijet_13TeV.json'),'r') as f:
        #with open(os.path.join(os.path.dirname(__file__),'fakes_trigIso_13TeV.json'),'r') as f:
        with open(os.path.join(os.path.dirname(__file__),'fakes_veryTight_dijet_13TeV.json'),'r') as f:
            self.fake_dict_13tev = json.load(f)
        # WZ 8TeV
        with open(os.path.join(os.path.dirname(__file__),'fakes_8TeV.json'),'r') as f:
            self.fake_dict_8tev = json.load(f)
        # spain fake rates
        efilename = os.path.join(os.path.dirname(__file__),'EGFR_RunII_25ns_jet35_08Jan.root')
        self.efrTfile = rt.TFile.Open(efilename,'READ')
        self.eFakeHist = self.efrTfile.Get('FR_pT_eta_EWKcorr')
        mfilename = os.path.join(os.path.dirname(__file__),'MuFR_RunII_25ns_jet20_08Jan.root')
        self.mfrTfile = rt.TFile.Open(mfilename,'READ')
        self.mFakeHist = self.mfrTfile.Get('FR_pT_eta_EWKcorr')
        # my fakerates
        filename = os.path.join(os.path.dirname(__file__),'fakes_dijet_13TeV.root')
        self.fakeTfile = rt.TFile.Open(filename,'READ')
        self.eTightFakeHist = self.fakeTfile.Get('FakeRateProbeElecTight')
        self.eMediumFakeHist = self.fakeTfile.Get('FakeRateProbeElecMedium')
        self.mMediumFakeHist = self.fakeTfile.Get('FakeRateProbeMuonMedium')

    def close(self):
        self.efrTfile.Close()
        self.mfrTfile.Close()
        self.fakeTfile.Close()

    def scale_factor(self, row, *lep_list, **kwargs):
        period = kwargs.pop('period',13)
        shift = kwargs.pop('metShift','')
        numer = kwargs.pop('numer','Tight')
        denom = kwargs.pop('denom','Loose')
        out = []
        for l in lep_list:
            lep_type = l[0]

            if lep_type == 'm':
                out += [self.m_wz_fake_veto(row,l,period,shift,numer,denom)]
            elif lep_type == 'e':
                out += [self.e_wz_fake_veto(row,l,period,shift,numer,denom)]
            elif lep_type == 't':
                out += [[0,0,0]] # TODO
            else:
                raise TypeError("Lepton type %s not recognized" % lep_type)

        final = [1,1,1]
        for o in out:
            final[0] *= o[0]
            final[1] *= o[1]
            final[2] *= o[2]
        return final

    def get_fake_err(self,row,l,pt,eta,numer,denom):
        # TODO: different numerator and denominator
        #if numer=='Medium' and denom=='Loose':
        #    if l[0]=='e': return 0.152, 0.010
        #hist = self.eFakeHist if l[0]=='e' else self.mFakeHist
        if l[0]=='m':
            hist = self.mMediumFakeHist
        if l[0]=='e':
            if numer=='Medium':
                hist = self.eMediumFakeHist
            if numer=='Tight':
                hist = self.eTightFakeHist
        fake = hist.GetBinContent(hist.FindBin(pt,eta))
        err = hist.GetBinError(hist.FindBin(pt,eta))
        return fake,err
        #fakelist = self.fake_dict_13tev[fakename] if period == 13 else self.fake_dict_8tev[fakename]
        #for fakedict in fakelist:
        #    ptlow   = fakedict['pt_low']
        #    pthigh  = fakedict['pt_high']
        #    etalow  = fakedict['eta_low']
        #    etahigh = fakedict['eta_high']
        #    if pt>=ptlow and pt<pthigh and eta>=etalow and eta<etahigh:
        #        fake = fakedict['fakerate']
        #        err = fakedict['error']
        #        return fake, err
        #return 0.0, 0.0

    def e_wz_fake_veto(self, row, l, period,shift,numer,denom):
        #fake, err = self.get_fake_err(row,l,'ZTightProbeElecTight',period,shift)
        #fake, err = self.get_fake_err(row,l,"FakeRateProbeElecTight",period,shift)
        pt = getattr(row, "%sPt" % l)
        eta = abs(getattr(row, "%sSCEta" % l))
        if shift=='ees+': pt = getattr(row, "%sPt_ElectronEnUp" % l)
        if shift=='ees-': pt = getattr(row, "%sPt_ElectronEnDown" % l)
        fake, err = self.get_fake_err(row,l,pt,eta,numer,denom)
        default = fake
        up = min([fake+err,1.])
        down = max([fake-err,0.])
        return [default, up, down]

    def m_wz_fake_veto(self, row, l, period,shift,numer,denom):
        #fake, err = self.get_fake_err(row,l,'ZTightProbeMuonTight',period,shift)
        #fake, err = self.get_fake_err(row,l,"FakeRateProbeMuonTight",period,shift)
        pt = getattr(row, "%sPt" % l)
        eta = abs(getattr(row, "%sEta" % l))
        if shift=='mes+': pt = getattr(row, "%sPt_MuonEnUp" % l)
        if shift=='mes-': pt = getattr(row, "%sPt_MuonEnDown" % l)
        fake, err = self.get_fake_err(row,l,pt,eta,numer,denom)
        default = fake
        up = min([fake+err,1.])
        down = max([fake-err,0.])
        return [default, up, down]


class LeptonScaleFactors(object):

    def __init__(self):
        # WZ 8TeV
        with open(os.path.join(os.path.dirname(__file__),'MuonEfficiencies_Run2012ReReco_53X.pkl'),'r') as mf:
            self.m_id_dict = pickle.load(mf)
        
        with open(os.path.join(os.path.dirname(__file__),'MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl'),'r') as mf:
            self.m_iso_dict = pickle.load(mf)

        path = os.path.join(os.path.dirname(__file__), 'MuonScaleFactors_2011_2012.root')
        self.m_rtfile = rt.TFile(path, 'READ')
        self.m_hist = self.m_rtfile.Get("TH2D_ALL_2012")

        # WZ 13TeV
        with open(os.path.join(os.path.dirname(__file__),'muons_13TeV.json'),'r') as mf:
            self.m_id_dict_13tev = json.load(mf)
        with open(os.path.join(os.path.dirname(__file__),'electrons_13TeV.json'),'r') as ef:
            self.e_id_dict_13tev = json.load(ef)

        # 4l 8TeV
        path = os.path.join(os.path.dirname(__file__),'CombinedMethod_ScaleFactors_RecoIdIsoSip.root')
        self.e_rtfile = rt.TFile(path, 'READ')
        self.e_hist = self.e_rtfile.Get("h_electronScaleFactor_RecoIdIsoSip")

        # MIT scale factors
        efilename = os.path.join(os.path.dirname(__file__),'scalefactors_ele.root')
        self.eScaleTfile = rt.TFile.Open(efilename,'READ')
        self.eScaleHist_loose = self.eScaleTfile.Get('unfactorized_scalefactors_Loose_ele')
        self.eScaleHist_medium = self.eScaleTfile.Get('unfactorized_scalefactors_Medium_ele')
        self.eScaleHist_tight = self.eScaleTfile.Get('unfactorized_scalefactors_Tight_ele')
        mfilename = os.path.join(os.path.dirname(__file__),'scalefactors_mu.root')
        self.mScaleTfile = rt.TFile.Open(mfilename,'READ')
        self.mScaleHist_loose = self.mScaleTfile.Get('unfactorized_scalefactors_Loose_mu')
        self.mScaleHist_medium = self.mScaleTfile.Get('unfactorized_scalefactors_Medium_mu')
        self.mScaleHist_tight = self.mScaleTfile.Get('unfactorized_scalefactors_Medium_mu')


    def close(self):
        self.e_rtfile.Close()
        self.m_rtfile.Close()
        self.mScaleTfile.Close()
        self.eScaleTfile.Close()

    def scale_factor(self, row, *lep_list, **kwargs):
        period = kwargs.pop('period',8)
        shift = kwargs.pop('metShift','')
        lepType = kwargs.pop('lepType','Tight')
        out = 1.0
        out = []
        for l in lep_list:
            lep_type = l[0]

            if lep_type == 'm':
                if period==8:
                    out += [self.m_4l_scale(row,l)] if lepType=='Loose' else [self.m_tight_scale(row, l)]
                if period==13:
                    thisScale = self.m_wz_scale(row,l,shift,lepType)
                    out += [thisScale]
            elif lep_type == 'e':
                if period==8:
                    out += [self.e_4l_scale(row,l)] if lepType=='Loose' else [self.e_ww_scale(row, l)]
                if period==13:
                    thisScale = self.e_wz_scale(row,l,shift,lepType)
                    out += [thisScale]
            elif lep_type == 't':
                out += [[1,1,1]] # TODO
            else:
                raise TypeError("Lepton type %s not recognized" % lep_type)

        final = [1,1,1]
        for o in out:
            final[0] *= o[0]
            final[1] *= o[1]
            final[2] *= o[2]

        return final

    def get_scale_err(self,row,l,idname,shift):
        pt = getattr(row, "%sPt" % l)
        if l[0]=='e' and shift=='ees+': pt = getattr(row, "%sPt_ElectronEnUp" % l)
        if l[0]=='e' and shift=='ees-': pt = getattr(row, "%sPt_ElectronEnDown" % l)
        if l[0]=='m' and shift=='mes+': pt = getattr(row, "%sPt_MuonEnUp" % l)
        if l[0]=='m' and shift=='mes-': pt = getattr(row, "%sPt_MuonEnDown" % l)
        eta = abs(getattr(row, "%sSCEta" % l)) if l[0]=='e' else abs(getattr(row, "%sEta" % l))
        ldict = getattr(self,'{0}_id_dict_13tev'.format(l[0]))
        idlist = ldict[idname]
        for iddict in idlist:
            ptlow = iddict['pt_lo']
            pthi = iddict['pt_hi']
            etalow = iddict['abseta_lo']
            etahi = iddict['abseta_hi']
            if pt>=ptlow and pt<pthi and eta>=etalow and eta<etahi:
                scale = iddict['ratio']
                err = iddict['ratio_err']
                return scale, err
        return 1.0, 0.0

    def get_scale_err_new(self,row,l,hist,shift):
        pt = getattr(row, "%sPt" % l)
        if l[0]=='e' and shift=='ees+': pt = getattr(row, "%sPt_ElectronEnUp" % l)
        if l[0]=='e' and shift=='ees-': pt = getattr(row, "%sPt_ElectronEnDown" % l)
        if l[0]=='m' and shift=='mes+': pt = getattr(row, "%sPt_MuonEnUp" % l)
        if l[0]=='m' and shift=='mes-': pt = getattr(row, "%sPt_MuonEnDown" % l)
        eta = abs(getattr(row, "%sSCEta" % l)) if l[0]=='e' else abs(getattr(row, "%sEta" % l))
        if eta>2.4: eta = 2.39 # TODO: stupid hack for MIT electrons
        val = hist.GetBinContent(hist.FindBin(eta,pt))
        err = hist.GetBinError(hist.FindBin(eta,pt))
        return val,err

    def m_wz_scale(self, row, l, shift, lepType):
        hists = {
            'Loose' : self.mScaleHist_loose,
            'Medium' : self.mScaleHist_medium,
            'Tight' : self.mScaleHist_tight,
        }
        scale,err = self.get_scale_err_new(row,l,hists[lepType],shift)
        return [scale, scale+err, scale-err]

        #idString = ''
        #isoString = ''
        #if lepType=='Loose':
        #    idString = 'passingIDWZLoose'
        #    isoString = 'passingIsoWZLoose_passingIDWZLoose'
        #if lepType=='Medium':
        #    idString = 'passingIDWZTight'
        #    isoString = 'passingIsoWZTight_passingIDWZTight'
        #if lepType=='Tight':
        #    idString = 'passingIDWZTight'
        #    isoString = 'passingIsoWZTight_passingIDWZTight'
        ## id
        #idscale, iderr = self.get_scale_err(row,l,idString,shift) if idString else [1.,0.]
        ## iso
        #isoscale, isoerr = self.get_scale_err(row,l,isoString,shift) if isoString else [1.,0.]
        #return [idscale*isoscale, (idscale+iderr)*(isoscale+isoerr), (idscale-iderr)*(isoscale-isoerr)]

    def e_wz_scale(self, row, l,shift, lepType):
        hists = {
            'Loose' : self.eScaleHist_loose,
            'Medium' : self.eScaleHist_medium,
            'Tight' : self.eScaleHist_tight,
        }
        scale,err = self.get_scale_err_new(row,l,hists[lepType],shift)
        return [scale, scale+err, scale-err]

        #idString = ''
        #isoString = ''
        #if lepType=='Loose':
        #    idString = 'passingLoose'
        #if lepType=='Medium':
        #    idString = 'passingMedium'
        #if lepType=='Tight':
        #    idString = 'passingTight'
        ## id
        #idscale, iderr = self.get_scale_err(row,l,idString,shift) if idString else [1.,0.]
        ## iso
        #isoscale, isoerr = self.get_scale_err(row,l,isoString,shift) if isoString else [1.,0.]
        #return [idscale*isoscale, (idscale+iderr)*(isoscale+isoerr), (idscale-iderr)*(isoscale-isoerr)]

    def e_ww_scale(self, row, l):
        pt = getattr(row, "%sPt" % l)
        eta = abs(getattr(row, "%sSCEta" % l))
        pts = [10, 15, 20, 30, 40, 50, 200]
        if eta < 0.8:
            effs = [0.662, 0.901, 0.943, 0.961, 0.976, 0.974]
            errs = [0.019, 0.009, 0.003, 0.001, 0.001, 0.001]
        elif eta < 1.4442:
            effs = [0.730, 0.942, 0.950, 0.944, 0.967, 0.970]
            errs = [0.027, 0.014, 0.003, 0.151, 0.092, 0.001]
        elif eta < 1.556:
            effs = [0.808, 0.857, 0.917, 0.964, 0.954, 0.986]
            errs = [0.094, 0.066, 0.010, 0.005, 0.004, 0.009]
        elif eta < 2.0:
            effs = [0.606, 0.834, 0.922, 0.925, 0.961, 0.963]
            errs = [0.037, 0.019, 0.005, 0.002, 0.001, 0.168]
        elif eta < 2.5:
            effs = [0.644, 0.759, 0.972, 0.981, 0.982, 0.970]
            errs = [0.032, 0.017, 0.006, 0.005, 0.002, 0.003]
        else:
            return [1., 1, 1]
        for p in range(len(pts)-1):
            if pt<pts[p+1]: return [effs[p], effs[p] + errs[p], effs[p] - errs[p]]
        return [1.0, 1, 1]

    def getTightScale(self,pt,eta):
        tight = self.m_id_dict['Tight']
        iso = self.m_iso_dict['combRelIsoPF04dBeta<012_Tight']
        if eta < 0.9:
           etaName = 'ptabseta<0.9'
        elif eta < 1.2:
           etaName = 'ptabseta0.9-1.2'
        elif eta < 2.1:
           etaName = 'ptabseta1.2-2.1'
        else:
           etaName = 'ptabseta2.1-2.4'

        if pt < 20: ptName = '10_20'
        elif pt < 25: ptName = '20_25'
        elif pt < 30: ptName = '25_30'
        elif pt < 35: ptName = '30_35'
        elif pt < 40: ptName = '35_40'
        elif pt < 50: ptName = '40_50'
        elif pt < 60: ptName = '50_60'
        elif pt < 90: ptName = '60_90'
        elif pt < 140: ptName = '90_140'
        else: ptName = '140_300'

        id_sf = tight[etaName][ptName]['data/mc']['efficiency_ratio']
        id_sf_up = id_sf + tight[etaName][ptName]['data/mc']['err_hi']
        id_sf_down = id_sf - tight[etaName][ptName]['data/mc']['err_low']

        iso_sf = iso[etaName][ptName]['data/mc']['efficiency_ratio']
        iso_sf_up = iso_sf + iso[etaName][ptName]['data/mc']['err_hi']
        iso_sf_down = iso_sf - iso[etaName][ptName]['data/mc']['err_low']

        return [id_sf * iso_sf, id_sf_up * iso_sf_up, id_sf_down * iso_sf_down]

    def m_tight_scale(self, row, l):
        pt = getattr(row, "%sPt" % l)
        eta = abs(getattr(row, "%sEta" % l))
        corr = self.getTightScale(pt,eta)
        return corr

    def m_ww_scale(self, row, l):
        pt = getattr(row, "%sPt" % l)
        eta = abs(getattr(row, "%sEta" % l))
        pts = [10., 15., 20, 25, 30, 50, 150]
        if eta < 0.9:
            effs = [0.992, 0.961, 0.982, 1.000, 0.993, 0.994]
            errs = [0.012, 0.005, 0.002, 0.001, 0.001, 0.001]
        elif eta < 1.2:
            effs = [0.971, 0.951, 0.982, 0.993, 0.991, 0.991]
            errs = [0.358, 0.005, 0.001, 0.002, 0.001, 0.001]
        elif eta < 2.5:
            effs = [1.002, 0.995, 1.020, 1.019, 1.002, 1.005]
            errs = [0.005, 0.003, 0.002, 0.001, 0.001, 0.001]
        else:
            return [1,1,1]
        for p in range(len(pts)-1):
            if pt<pts[p+1]: return [effs[p], effs[p]+errs[p], effs[p]-errs[p]]
        return [1.0,1,1]

    def e_4l_scale(self, row, l):
        pt = getattr(row, "%sPt" % l)
        eta = getattr(row, "%sEta" % l)
        global_bin = self.e_hist.FindBin(pt, eta)
        scl = self.e_hist.GetBinContent(global_bin)
        err = self.e_hist.GetBinError(global_bin)

        if scl < 0.1:
            scl = 1.0
            err = 0.02

        return [scl, scl + err, scl - err]

    def m_4l_scale(self, row, l):
        pt = getattr(row, "%sPt" % l)
        eta = getattr(row, "%sEta" % l)
        global_bin = self.m_hist.FindBin(pt, eta)
        scl = self.m_hist.GetBinContent(global_bin)
        err = self.m_hist.GetBinError(global_bin)

        if scl < 0.1:
            scl = 1.0
            err = 0.005

        return [scl, scl + err, scl - err]

class BjetScaleFactors(object):

    def __init__(self):
        # WZ 13TeV
        with open(os.path.join(os.path.dirname(__file__),'CSVv2_13TeV_74X.csv'),'rb') as csv:
            self.btag_scales = csv.reader(csv,delimiter=',')

    def close(self):
        return

    def scale_factor(self, row, numjets, numbjets, **kwargs):
        workingPoint = kwargs.pop('workingPoint','tight')
        period = kwargs.pop('period',13)
        shift = kwargs.pop('metShift','')

        opPoints = {
            'loose': 0,
            'medium': 1,
            'tight': 2,
        }
        shiftMap = {
            'jes+': 'JetEnUp',
            'jes-': 'JetEnDown',
        }

        scale = [1.0, 1.0, 1.0]
        for i in range(numjets):
            jetPtName = 'jet{0}Pt'.format(i+1)
            jetEtaName = 'jet{0}Eta'.format(i+1)
            if shift in shiftMap:
                jetPtName += '_{0}'.format(shiftMap[shift])
                jetEtaName += '_{0}'.format(shiftMap[shift])
            if not hasattr(row,jetPtName) or not hasattr(row,jetEtaName):
                continue
            pt = getattr(row,jetPtName)
            eta = getattr(row,jetEtaName)
            if pt < 30: continue
            thisScale = self.bjet_scale(opPoints[workingPoint],pt,eta)
            scale[0] *= thisScale[0]
            scale[1] *= thisScale[1]
            scale[2] *= thisScale[2]

    def bjet_scale(self,opPoint,pt,eta):
        thisScale = [1.,1.,1.]
        for s in self.btag_scales[1:]:
            OperatingPoint, measurementType, sysType, jetFlavor, etaMin, etaMax, ptMin, ptMax, discrMin, discrMax, formula = s
            if int(OperatingPoint)==opPoint and measurementType=='mujets' and int(jetFlavor)==0:
                if eta>float(etaMin) and eta<float(etaMax) and pt>float(ptMin) and pt<float(ptMax):
                    thisFormula = formula.replace(x,str(pt))
                    val = eval(thisFormula)
                    if sysType=='central': thisScale[0] = val if False else 1-val
                    if sysType=='up': thisScale[1] = val if False else 1-val
                    if sysType=='down': thisScale[2] = val if False else 1-val
        return thisScale
