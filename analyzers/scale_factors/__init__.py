import sys
import os
import glob
import pickle
from operator import itemgetter, attrgetter

sys.argv.append('-b')
import ROOT as rt
sys.argv.pop()

class TriggerScaleFactors(object):

    def __init__(self):
        pass

    def scale_factor(self, rtrow, *lep_list, **kwargs):
        lep_objs = [(x, getattr(rtrow,"%sPt"%x), abs(getattr(rtrow,"%sEta"%x))) for x in lep_list]
        lep_ord = sorted(lep_objs, key=itemgetter(1), reverse=True)
        eff = 1-(\
                (1-self.double_lead_eff(*lep_ord[0])) * (1-self.double_lead_eff(*lep_ord[1])) * (1-self.double_lead_eff(*lep_ord[2]))\
                + self.double_lead_eff(*lep_ord[0]) * (1-self.double_trail_eff(*lep_ord[1])) * (1-self.double_trail_eff(*lep_ord[2]))\
                + self.double_lead_eff(*lep_ord[1]) * (1-self.double_trail_eff(*lep_ord[2])) * (1-self.double_trail_eff(*lep_ord[0]))\
                + self.double_lead_eff(*lep_ord[2]) * (1-self.double_trail_eff(*lep_ord[0])) * (1-self.double_trail_eff(*lep_ord[1]))\
                )
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

    def single_e(self,pt,eta):
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
        m_id_file = open(os.path.join(os.path.dirname(__file__),'MuonEfficiencies_Run2012ReReco_53X.pkl'),'r')
        self.m_id_dict = pickle.load(m_id_file)
        m_iso_file = open(os.path.join(os.path.dirname(__file__),'MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl'),'r')
        self.m_iso_dict = pickle.load(m_iso_file)


    def scale_factor(self, row, *lep_list, **kwargs):
        tight = kwargs.pop('tight',False)
        out = 1.0
        out = []
        for l in lep_list:
            lep_type = l[0]

            if lep_type == 'm':
                out += [self.m_tight_scale(row, l)]
            elif lep_type == 'e':
                out += [self.e_ww_scale(row, l)]
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

