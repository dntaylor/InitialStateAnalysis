'''
Data styles to be used with the PlotterBase class.
'''

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import copy

dataStyles = {}


# Some colors
Gray         = ROOT.TColor.GetColor('#B8B8B8')
GrayA        = ROOT.TColor.GetColor('#C8C8C8')
Purple       = ROOT.TColor.GetColor('#AD33FF')
PurpleA      = ROOT.TColor.GetColor('#7924B2')
Yellow       = ROOT.TColor.GetColor('#FFFF00')
YellowA      = ROOT.TColor.GetColor('#FFCC26')
Gold         = ROOT.TColor.GetColor('#CCCC00')
GoldA        = ROOT.TColor.GetColor('#D1D119')
DarkYellow   = ROOT.TColor.GetColor('#FFCC00')
DarkYellowA  = ROOT.TColor.GetColor('#E6B800')
Orange       = ROOT.TColor.GetColor('#DC7612')
OrangeA      = ROOT.TColor.GetColor('#BD3200')
Blue         = ROOT.TColor.GetColor('#107FC9')
BlueA        = ROOT.TColor.GetColor('#0E4EAD')
Navy         = ROOT.TColor.GetColor('#003399')
NavyA        = ROOT.TColor.GetColor('#00297A')
Steel        = ROOT.TColor.GetColor('#9999FF')
SteelA       = ROOT.TColor.GetColor('#B8B8FF')
DarkRed      = ROOT.TColor.GetColor('#A30000')
DarkRedA     = ROOT.TColor.GetColor('#8F0000')
Red          = ROOT.TColor.GetColor('#F01800')
RedA         = ROOT.TColor.GetColor('#780000')
Green        = ROOT.TColor.GetColor('#36802D')
GreenA       = ROOT.TColor.GetColor('#234D20')
BlueGreen    = ROOT.TColor.GetColor('#00CC99')
BlueGreenA   = ROOT.TColor.GetColor('#00A37A')
LightGreen   = ROOT.TColor.GetColor('#66FF99')
LightGreenA  = ROOT.TColor.GetColor('#52CC7A')
Lime         = ROOT.TColor.GetColor('#9ED54C')
LimeA        = ROOT.TColor.GetColor('#59A80F')
Aqua         = ROOT.TColor.GetColor('#66FFFF')
AquaA        = ROOT.TColor.GetColor('#52CCCC')

dataStyles['data'] = {
    'legendstyle' : 'pe',
    'drawstyle' : 'pe',
    'name' : "Observed",
}

dataStyles['datadriven'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : GrayA,
    'fillcolor' : Gray,
    #'name' : "Data driven",
    'name' : "Non-Prompt",
    'fillstyle': 1001,
}

dataStyles['ZZJets'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : BlueA,
    'fillcolor' : Blue,
    'name' : "ZZ",
    'fillstyle': 1001,
}

dataStyles['ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola'] = copy.deepcopy(dataStyles['ZZJets'])
dataStyles['Diboson'] = copy.deepcopy(dataStyles['ZZJets'])
dataStyles['Diboson']['name'] = 'Diboson'

dataStyles['HZZ'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : SteelA,
    'fillcolor' : Steel,
    'name' : "HZZ",
    'fillstyle': 1001,
}


dataStyles['ZG'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : DarkRedA,
    'fillcolor' : DarkRed,
    'name' : "Z#gamma",
    'fillstyle': 1001,
}

dataStyles['ZGFiltered'] = copy.deepcopy(dataStyles['ZG'])
dataStyles['ZGToLLG_8TeV-madgraph'] = copy.deepcopy(dataStyles['ZG'])

dataStyles['WG'] = copy.deepcopy(dataStyles['ZG'])
dataStyles['WG']['name'] = 'W#gamma'
dataStyles['WG']['fillcolor'] = DarkYellow
dataStyles['WG']['linecolor'] = DarkYellowA


dataStyles['WZJets'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : PurpleA,
    'fillcolor' : Purple,
    'name' : "WZ",
    'fillstyle': 1001,
}

dataStyles['WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola'] = copy.deepcopy(dataStyles['WZJets'])
dataStyles['WWJets'] = copy.deepcopy(dataStyles['WZJets'])
dataStyles['WWJets']['name'] = 'WW'
dataStyles['WWJets']['linecolor'] = AquaA
dataStyles['WWJets']['fillcolor'] = Aqua

dataStyles['ZZZNoGstarJets'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : NavyA,
    'fillcolor' : Navy,
    'name' : "ZZZ+Jets",
    'fillstyle': 1001,
}

dataStyles['VVVJets'] = copy.deepcopy(dataStyles['ZZZNoGstarJets'])
dataStyles['VVVJets']['name'] = "VVV"

dataStyles['WWZNoGstarJets'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : PurpleA,
    'fillcolor' : Purple,
    'name' : "WWZ+Jets",
    'fillstyle': 1001,
}

dataStyles['WWWJets'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : AquaA,
    'fillcolor' : Aqua,
    'name' : "WWW+Jets",
    'fillstyle': 1001,
}

dataStyles['WWGJets'] = copy.deepcopy(dataStyles['WWWJets'])
dataStyles['WWGJets']['name'] = "WW#gamma+Jets"

dataStyles['ZJets'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'fillcolor' : Yellow,
    'linecolor' : YellowA,
    'name' : "Z+Jets",
    'fillstyle': 1001,
}

dataStyles['ZJetsFiltered'] = copy.deepcopy(dataStyles['ZJets'])
dataStyles['DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'] = copy.deepcopy(dataStyles['ZJets'])
dataStyles['DYJetsToLL_M-10To50filter_8TeV-madgraph'] = copy.deepcopy(dataStyles['ZJets'])
dataStyles['DYJetsToLL_M-10To50filter_8TeV-madgraph']['name'] = "Z Low Mass"
dataStyles['DYJetsToLL_M-10To50filter_8TeV-madgraph']['fillcolor'] = Gold
dataStyles['DYJetsToLL_M-10To50filter_8TeV-madgraph']['linecolor'] = GoldA
dataStyles['Z1jets_M50'] = copy.deepcopy(dataStyles['ZJets'])
dataStyles['Z2jets_M50_S10'] = copy.deepcopy(dataStyles['ZJets'])
dataStyles['Z3jets_M50'] = copy.deepcopy(dataStyles['ZJets'])
dataStyles['Z4jets_M50'] = copy.deepcopy(dataStyles['ZJets'])
dataStyles['WJets'] = copy.deepcopy(dataStyles['ZJets'])
dataStyles['WJets']['name'] = "W+Jets"
dataStyles['WJets']['fillcolor'] = DarkYellow
dataStyles['WJets']['linecolor'] = DarkYellowA

dataStyles['TTJets'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : GreenA,
    'fillcolor' : Green,
    'name' : "t#bar{t}+Jets",
    'fillstyle': 1001,
}

dataStyles['TTZJets'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'linecolor' : BlueGreenA,
    'fillcolor' : BlueGreen,
    'name' : "t#bar{t}Z+Jets",
    'fillstyle': 1001,
}

dataStyles['TTZJets_Tune4C_13TeV-madgraph-tauola'] = copy.deepcopy(dataStyles['TTZJets'])
dataStyles['TTWJets'] = copy.deepcopy(dataStyles['TTZJets'])
dataStyles['TTWJets']['name'] = "t#bar{t}W+Jets"
dataStyles['TTWJets']['linecolor'] = LightGreenA
dataStyles['TTWJets']['fillcolor'] = LightGreen
dataStyles['TTWJets_Tune4C_13TeV-madgraph-tauola'] = copy.deepcopy(dataStyles['TTWJets'])
dataStyles['TTWWJets'] = copy.deepcopy(dataStyles['TTZJets'])
dataStyles['TTWWJets']['name'] = "t#bar{t}WW+Jets"
dataStyles['TTGJets'] = copy.deepcopy(dataStyles['TTZJets'])
dataStyles['TTGJets']['name'] = "t#bar{t}#gamma+Jets"
dataStyles['TTVJets'] = copy.deepcopy(dataStyles['TTZJets'])
dataStyles['TTVJets']['name'] = "t#bar{t}V"

dataStyles['TToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'fillcolor' : Lime,
    'linecolor' : LimeA,
    'name' : "t_{s}",
    'fillstyle': 1001,
}

dataStyles['TToLeptons_t-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'fillcolor' : Lime,
    'linecolor' : LimeA,
    'name' : "t_{t}",
    'fillstyle': 1001,
}

dataStyles['T_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'fillcolor' : Lime,
    'linecolor' : LimeA,
    'name' : "t_{t}W",
    'fillstyle': 1001,
}

dataStyles['TBarToLeptons_s-channel-CSA14_Tune4C_13TeV-aMCatNLO-tauola'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'fillcolor' : Lime,
    'linecolor' : LimeA,
    'name' : "#bar{t}_s",
    'fillstyle': 1001,
}

dataStyles['TBarToLeptons_t-channel_Tune4C_CSA14_13TeV-aMCatNLO-tauola'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'fillcolor' : Lime,
    'linecolor' : LimeA,
    'name' : "#bar{t}_{t}",
    'fillstyle': 1001,
}

dataStyles['Tbar_tW-channel-DR_Tune4C_13TeV-CSA14-powheg-tauola'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'fillcolor' : Lime,
    'linecolor' : LimeA,
    'name' : "#bar{t}_{t}W",
    'fillstyle': 1001,
}

dataStyles['SingleTop'] = {
    'legendstyle' : 'f',
    'drawstyle' : 'hist',
    'fillcolor' : Lime,
    'linecolor' : LimeA,
    'name' : "Single Top",
    'fillstyle': 1001,
}

for mass in [110, 130, 150, 170, 200, 250, 300, 350, 400, 450, 500, 600, 700]:
    dataStyles['DBLH_m%i' % mass] = {
        'legendstyle' : 'f',
        'drawstyle' : 'hist',
        'linecolor' : OrangeA,
        'fillcolor' : Orange,
        'name' : "#Phi^{++}#Phi^{--}#rightarrow4l (%i GeV)" % mass,
        'fillstyle': 1001,
    }
    dataStyles['HPlusPlusHMinusMinusHTo4L_M-%i_8TeV-pythia6' % mass] = copy.deepcopy(dataStyles['DBLH_m%i' % mass])
    dataStyles['HPlusPlusHMinusHTo3L_M-%i_8TeV-calchep-pythia6' % mass] = copy.deepcopy(dataStyles['DBLH_m%i' % mass])
    dataStyles['HPlusPlusHMinusHTo3L_M-%i_8TeV-calchep-pythia6' % mass]['name'] = "#Phi^{++}#Phi^{-}#rightarrow3l (%i GeV)" % mass

