import sys
import os
import glob
import pickle
import json
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
        lep_objs = [(x, getattr(rtrow,"%sPt"%x), abs(getattr(rtrow,"%sEta"%x))) for x in lep_list]
        lep_ord = sorted(lep_objs, key=itemgetter(1), reverse=True)
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
        return eff

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
        with open(os.path.join(os.path.dirname(__file__),'muon.json'),'r') as mf:
            self.m_id_dict_13tev = json.load(mf)
        with open(os.path.join(os.path.dirname(__file__),'electron.json'),'r') as ef:
            self.e_id_dict_13tev = json.load(ef)

        # 4l 8TeV
        path = os.path.join(os.path.dirname(__file__),'CombinedMethod_ScaleFactors_RecoIdIsoSip.root')
        self.e_rtfile = rt.TFile(path, 'READ')
        self.e_hist = self.e_rtfile.Get("h_electronScaleFactor_RecoIdIsoSip")

    def close(self):
        self.e_rtfile.Close()
        self.m_rtfile.Close()

    def scale_factor(self, row, *lep_list, **kwargs):
        tight = kwargs.pop('tight',False)
        loose = kwargs.pop('loose',False)
        period = kwargs.pop('period',8)
        out = 1.0
        out = []
        for l in lep_list:
            lep_type = l[0]

            if lep_type == 'm':
                if period==8:
                    out += [self.m_4l_scale(row,l)] if loose else [self.m_tight_scale(row, l)]
                if period==13:
                    out += [self.m_wz_scale_loose(row,l)] if loose else [self.m_wz_scale_tight(row,l)]
            elif lep_type == 'e':
                if period==8:
                    out += [self.e_4l_scale(row,l)] if loose else [self.e_ww_scale(row, l)]
                if period==13:
                    out += [self.e_wz_scale_loose(row,l)] if loose else [self.e_wz_scale_tight(row,l)]
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

    def get_scale_err(self,row,l,idname):
        pt = getattr(row, "%sPt" % l)
        eta = abs(getattr(row, "%sSCEta" % l)) if l[0]=='e' else abs(getattr(row, "%sEta" % l))
        ldict = getattr(self,'{0}_id_dict_13tev'.format(l[0]))
        idlist = ldict[idname]
        pre = 'supercluster_' if l[0]=='e' else ''
        for iddict in idlist:
            ptlow = iddict['{0}et_lo'.format(pre)]
            pthi = iddict['{0}et_hi'.format(pre)]
            etalow = iddict['{0}eta_lo'.format(pre)]
            etahi = iddict['{0}eta_hi'.format(pre)]
            if pt>ptlow and pt<pthi and eta>etalow and eta<etahi:
                scale = iddict['ratio']
                err = iddict['ratio_err']
                return scale, err
        return 1.0, 0.0

    def m_wz_scale_loose(self, row, l):
        # id
        idscale, iderr = self.get_scale_err(row,l,'MuonWZIDLoose')
        # iso
        isoscale, isoerr = self.get_scale_err(row,l,'MuonWZIsolationLooseFromIDLoose')
        return [idscale*isoscale, (idscale+iderr)*(isoscale+isoerr), (idcale-iderr)*(isoscale-isoerr)]

    def m_wz_scale_tight(self, row, l):
        # id
        idscale, iderr = self.get_scale_err(row,l,'MuonWZIDTight')
        # iso
        isoscale, isoerr = self.get_scale_err(row,l,'MuonWZIsolationTightFromIDTight')
        return [idscale*isoscale, (idscale+iderr)*(isoscale+isoerr), (idcale-iderr)*(isoscale-isoerr)]

    def e_wz_scale_loose(self, row, l):
        # id
        idscale, iderr = self.get_scale_err(row,l,'cutBasedElectronID-Spring15-25ns-V1-standalone-loose')
        return [idscale, idscale+iderr, idcale-iderr]

    def e_wz_scale_tight(self, row, l):
        # id
        idscale, iderr = self.get_scale_err(row,l,'cutBasedElectronID-Spring15-25ns-V1-standalone-medium')
        return [idscale, idscale+iderr, idcale-iderr]


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
