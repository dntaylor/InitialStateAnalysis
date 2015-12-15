from copy import deepcopy

def singlePlot(plotMethod,args,kwargs):
    '''Make a single plot'''
    plotMethod(*args,**kwargs)

def multiplePlots(plotter,method,plotMap):
    '''
    Make multiple plots in series
    plotMap is of the form:
        {
            'plotkey' : {
                'args' : [...],
                'kwargs' : {...},
            },
            ...
        }
    '''
    plotMethod = getattr(plotter,method)
    for p in plotMap:
        singlePlot(plotMethod,plotMap[p]['args'],plotMap[p]['kwargs'])

def multiplePlotsParallel(plotters,methods,plotMaps):
    '''Parallelize the plotting, plotters must be unique,
       all lists must be same length'''
    jobs = []
    for p,m,pm in zip(plotters,methods,plotMaps):
        jobs += [(p,m,pm)]
    # TODO: parallelize here
    for j in jobs:
        multiplePlots(*j)


leptonParams = {
    'Pt'               : { 'args': ['{obj}.{pre}Pt{post}',               [20,0,200],            '{savedir}/{name}/Pt'],              'kwargs': {'xaxis': 'p_{{T}}^{{{pretty}}} (GeV)', 'yaxis': 'Events/10 GeV', 'overflow': True}},
    'Iso'              : { 'args': ['{obj}.{pre}Iso{post}',              [50,0,.5],             '{savedir}/{name}/Iso'],             'kwargs': {'xaxis': '\\rm{{Rel. Iso.}} ({pretty})', 'overflow': True}},
    'Eta'              : { 'args': ['{obj}.{pre}Eta{post}',              [30,-3.0,3.0],         '{savedir}/{name}/Eta'],             'kwargs': {'xaxis': '\\eta^{{{pretty}}}', 'legendpos': 43, 'numcol': 3}},
    'Phi'              : { 'args': ['{obj}.{pre}Phi{post}',              [30,-3.14159,3.14159], '{savedir}/{name}/Eta'],             'kwargs': {'xaxis': '\\phi^{{{pretty}}}', 'legendpos': 43, 'numcol': 3}},
    'Dxy'              : { 'args': ['{obj}.{pre}Dxy{post}',              [50,-0.1,0.1],         '{savedir}/{name}/Dxy'],             'kwargs': {'xaxis': 'd_{{0}}^{{{pretty}}}', 'legendpos': 43}},
    'Dz'               : { 'args': ['{obj}.{pre}Dz{post}',               [50,-0.1,0.1],         '{savedir}/{name}/Dz'],              'kwargs': {'xaxis': 'd_{{Z}}^{{{pretty}}}', 'legendpos': 43}},
    #'SigmaIEtaIEta'    : { 'args': ['{obj}.{pre}SigmaIEtaIEta{post}',    [40,0.,0.08],          '{savedir}/{name}/SigmaIEtaIEta'],   'kwargs': {'xaxis': '\\sigma_{{i\\eta i\\eta}}^{{{pretty}}}', 'legendpos': 43}},
    #'DEtaIn'           : { 'args': ['{obj}.{pre}DEtaIn{post}',           [50,-0.1,0.1],         '{savedir}/{name}/DEtaIn'],          'kwargs': {'xaxis': '\\Delta\\eta_{{In}}^{{{pretty}}}', 'legendpos': 43}},
    #'DPhiIn'           : { 'args': ['{obj}.{pre}DPhiIn{post}',           [50,-0.1,0.1],         '{savedir}/{name}/DPhiIn'],          'kwargs': {'xaxis': '\\Delta\\phi_{{In}}^{{{pretty}}}', 'legendpos': 43}},
    #'HOverE'           : { 'args': ['{obj}.{pre}HOverE{post}',           [50,0.,0.1],           '{savedir}/{name}/HOverE'],          'kwargs': {'xaxis': '\\frac{{H}}{{E}}^{{{pretty}}}', 'legendpos': 43}},
    #'OoEmOoP'          : { 'args': ['{obj}.{pre}OoEmOoP{post}',          [50,0.,0.1],           '{savedir}/{name}/OoEmOoP'],         'kwargs': {'xaxis': '|\\frac{{1}}{{E}}-\\frac{{1}}{{P}}| ({pretty})', 'legendpos': 43}},
    'TriggeringMVA'    : { 'args': ['{obj}.{pre}TriggeringMVA{post}',    [50,-1.,1.],           '{savedir}/{name}/TriggeringMVA'],    'kwargs': {'xaxis': 'MVA Score ({pretty})', 'legendpos': 43}},
    'NonTriggeringMVA' : { 'args': ['{obj}.{pre}NonTriggeringMVA{post}', [50,-1.,1.],           '{savedir}/{name}/NonTriggeringMVA'], 'kwargs': {'xaxis': 'MVA Score ({pretty})', 'legendpos': 43}},
}


def addLeptonPlotParams(savedir,obj,pre,post,name,pretty,doDetailed):
    '''Add the lepton object plots'''
    plotParams = {}
    basicVars = ['Pt','Iso','Eta','Phi']
    for v in leptonParams:
        if not doDetailed and v not in basicVars: continue
        vname = '{0}{1}'.format(name,v)
        plotParams[vname] = {
            'args'  : [leptonParams[v]['args'][0].format(obj=obj,pre=pre,post=post), leptonParams[v]['args'][1], leptonParams[v]['args'][2].format(savedir=savedir,name=name)],
            'kwargs': deepcopy(leptonParams[v]['kwargs']),
        }
        plotParams[vname]['kwargs']['xaxis'] = leptonParams[v]['kwargs']['xaxis'].format(pretty=pretty)
    return plotParams

dileptonParams = {
    'Mass': { 'args': ['{obj}.mass', [24,0,600],'{savedir}/{name}/Mass'], 'kwargs': {'yaxis': 'Events/25 GeV', 'xaxis': 'm_{{{pretty}}} (GeV)', 'overflow': True} },
    'DPhi': { 'args': ['{obj}.dPhi', [32,0,3.2],'{savedir}/{name}/dPhi'], 'kwargs': {'yaxis': 'Events/0.1 rad', 'xaxis': '\\Delta\\phi_{{{pretty}}} (rad)'} },
    'Pt'  : { 'args': ['{obj}.Pt',   [40,0,400],'{savedir}/{name}/Pt'],   'kwargs': {'yaxis': 'Events/10 GeV', 'xaxis': 'p_{{T}}^{{{pretty}}} (GeV)', 'overflow': True} },
    'DR'  : { 'args': ['{obj}.dR',   [60,0,6],  '{savedir}/{name}/dR'],   'kwargs': {'xaxis': '\\Delta R({pretty})'} },
}

def addDileptonPlotParams(savedir,obj,name,pretty,doDetailed):
    '''Add dilepton object plots'''
    plotParams = {}
    for v in dileptonParams:
        vname = '{0}{1}'.format(name,v)
        plotParams[vname] = {
            'args'  : [dileptonParams[v]['args'][0].format(obj=obj), dileptonParams[v]['args'][1], dileptonParams[v]['args'][2].format(savedir=savedir,name=name)],
            'kwargs': deepcopy(dileptonParams[v]['kwargs'],
        }
        plotParams[vname]['kwargs']['xaxis'] = dileptonParams[v]['kwargs']['xaxis'].format(pretty=pretty)
    return plotParams

leptonMetParams = {
    'Pt'  : { 'args': ['{obj}.Pt',  [40,0,400],'{savedir}/{name}/Pt'],   'kwargs' : { 'yaxis': 'Events/10 GeV', 'xaxis': 'p_{{T}}^{{{pretty}}} (GeV)', 'overflow': True} },
    'Mass': { 'args': ['{obj}.mass',[40,0,400],'{savedir}/{name}/Mass'], 'kwargs' : { 'yaxis': 'Events/10 GeV', 'xaxis': 'm_{{T}}^{{{pretty}}} (GeV)', 'overflow': True} },
    'DPhi': { 'args': ['{obj}.dPhi',[32,0,3.2],'{savedir}/{name}/dPhi'], 'kwargs' : { 'yaxis': 'Events/0.1 rad', 'xaxis': '\\Delta\\phi({pretty} lepton, E_{{T}}^{{miss}}) (rad)'} },
}

def addLeptonMetPlotParams(savedir,obj,name,pretty,doDetailed):
    '''Add lepton+met object plots'''
    plotParams = {}
    for v in leptonMetParams:
        vname = '{0}{1}'.format(name,v)
        plotParams[vname] = {
            'args'  : [leptonMetParams[v]['args'][0].format(obj=obj), leptonMetParams[v]['args'][1], leptonMetParams[v]['args'][2].format(savedir=savedir,name=name)],
            'kwargs': deepcopy(leptonMetParams[v]['kwargs'],
        }
        plotParams[vname]['kwargs']['xaxis'] = leptonMetParams[v]['kwargs']['xaxis'].format(pretty=pretty)
    return plotParams


def addEventParams(savedir,doDetailed):
    eventParams = {
        'st'                   : { 'args': ['finalstate.sT',                                    [40,0,1000],          '{0}sT'.format(savedir)],                    'kwargs': {'xaxis': 'S_{T} (GeV)', 'yaxis': 'Events/25 GeV'} },
        'st_log'               : { 'args': ['finalstate.sT',                                    [40,0,1000],          '{0}sT_log'.format(savedir)],                'kwargs': {'xaxis': 'S_{T} (GeV)', 'yaxis': 'Events/25 GeV', 'logy': 1} },
        'numElectronsLoose'    : { 'args': ['finalstate.elecVetoLoose',                         [8,0,8],              '{0}numElectronsLoose'.format(savedir)],     'kwargs': {'xaxis': 'Number of Electrons (p_{T}>10 GeV)'} },
        'numMuonsLoose'        : { 'args': ['finalstate.muonVetoLoose',                         [8,0,8],              '{0}numMuonsLoose'.format(savedir)],         'kwargs': {'xaxis': 'Number of Muons (p_{T}>10 GeV)'} },
        'numLeptonsLoose'      : { 'args': ['finalstate.elecVetoLoose+finalstate.muonVetoLoose',[8,0,8],              '{0}numLeptonsLoose'.format(savedir)],       'kwargs': {'xaxis': 'Number of Leptons (p_{T}>10 GeV)'} },
        'numJets30'            : { 'args': ['finalstate.jetVeto30',                             [8,0,8],              '{0}numJets30'.format(savedir)],             'kwargs': {'xaxis': 'Number of Jets (p_{T}>30 GeV)'} },
        'numBJets30Medium'     : { 'args': ['finalstate.bjetVeto30Medium',                      [8,0,8],              '{0}numBJets30Medium'.format(savedir)],      'kwargs': {'xaxis': 'Number of b Jets (p_{T}>30 GeV)'} },
        'numBJets30Tight'      : { 'args': ['finalstate.bjetVeto30Tight',                       [8,0,8],              '{0}numBJets30Tight'.format(savedir)],       'kwargs': {'xaxis': 'Number of b Jets (p_{T}>30 GeV)'} },
        'met'                  : { 'args': ['finalstate.met',                                   [20,0,200],           '{0}met'.format(savedir)],                   'kwargs': {'yaxis': 'Events/10 GeV', 'xaxis': 'E_{T}^{miss} (GeV)', 'overflow': True} },
        'mass'                 : { 'args': ['finalstate.mass',                                  [25,0,500],           '{0}mass'.format(savedir)],                  'kwargs': {'yaxis': 'Events/20 GeV', 'xaxis': 'm_{3\\ell} (GeV)', 'overflow': True} },
        'mass_zoom'            : { 'args': ['finalstate.mass',                                  [100,0,500],          '{0}mass_zoom'.format(savedir)],             'kwargs': {'yaxis': 'Events/5 GeV', 'xaxis': 'm_{3\\ell} (GeV)', 'overflow': True} },
        'mT'                   : { 'args': ['finalstate.mT',                                    [50,0,1000],          '{0}mT'.format(savedir)],                    'kwargs': {'yaxis': 'Events/20 GeV', 'xaxis': 'm_T^{3\\ell+MET} (GeV)', 'overflow': True} },
        'puVertices'           : { 'args': ['event.nvtx',                                       [50,0,50],            '{0}puVertices'.format(savedir)],            'kwargs': {'xaxis': 'Number PU Vertices'} },
        'puVertices_noreweight': { 'args': ['event.nvtx',                                       [50,0,50],            '{0}puVertices_noreweight'.format(savedir)], 'kwargs': {'xaxis': 'Number PU Vertices', 'scalefactor': 'event.gen_weight*event.lep_scale*event.trig_scale'} },
        'JetPt'                : { 'args': ['finalstate.leadJetPt',                             [30,0,300],           '{0}JetPt'.format(savedir)],                 'kwargs': {'yaxis': 'Events/10 GeV', 'xaxis': 'p_{T}^{jet} (GeV)', 'overflow': True} },
        'JetEta'               : { 'args': ['finalstate.leadJetEta',                            [50,-5.0,5.0],        '{0}JetEta'.format(savedir)],                'kwargs': {'xaxis': '\\eta^{jet}'} },
        'JetPhi'               : { 'args': ['finalstate.leadJetPhi',                            [30,-3.14159,3.14159],'{0}JetPhi'.format(savedir)],                'kwargs': {'xaxis': '\\phi^{jet}'} },
    }
    return eventParams

def buildPlotParams(analysis,region,period,savedir,nl,doDetailed):
    if savedir: savedir = savedir + '/'
    plotParams = {}

    # each lepton
    leptonArgs = []
    if doDetailed:
        for i in range(nl):
            leptonArgs += [(savedir,'l{0}'.format(i+1),'','','l{0}'.format(i+),'\\ell{0}'.format(l+1),doDetailed)]
    if analysis in ['Hpp4l']:
        leptonArgs += [(savedir,'h1','','1','hpp/Leading','\\Phi^{++} \\text{Leading Lepton}',doDetailed)]
        leptonArgs += [(savedir,'h1','','2','hpp/SubLeading','\\Phi^{++} \\text{Subleading Lepton}',doDetailed)]
        leptonArgs += [(savedir,'h2','','1','hmm/Leading','\\Phi^{--} \\text{Leading Lepton}',doDetailed)]
        leptonArgs += [(savedir,'h2','','2','hmm/SubLeading','\\Phi^{--} \\text{Subleading Lepton}',doDetailed)]
    if analysis in ['Hpp3l']:
        leptonArgs += [(savedir,'h1','','1','hpp/Leading','\\Phi^{\\pm\\pm} \\text{Leading Lepton}',doDetailed)]
        leptonArgs += [(savedir,'h1','','2','hpp/SubLeading','\\Phi^{\\pm\\pm} \\text{Subleading Lepton}',doDetailed)]
        leptonArgs += [(savedir,'h2','','1','hm/Lepton','\\Phi^{\\mp} \\text{Lepton}',doDetailed)]
    if analysis in ['WZ','Hpp3l']:
        leptonArgs += [(savedir,'z1','','1','z1/Leading','Z Leading Lepton',doDetailed)]
        leptonArgs += [(savedir,'z1','','2','z1/SubLeading','Z Subleading Lepton',doDetailed)]
    if analysis in ['ZZ','Hpp4l']:
        leptonArgs += [(savedir,'z1','','1','z1/Leading','Z1 Leading Lepton',doDetailed)]
        leptonArgs += [(savedir,'z1','','2','z1/SubLeading','Z1 Subleading Lepton',doDetailed)]
        leptonArgs += [(savedir,'z2','','1','z2/Leading','Z2 Leading Lepton',doDetailed)]
        leptonArgs += [(savedir,'z2','','2','z2/SubLeading','Z2 Subleading Lepton',doDetailed)]
    if analysis in ['WZ','Hpp3l']:
        leptonArgs += [(savedir,'w1','','1','w1/Lepton','W Lepton',doDetailed)]

    for la in leptonArgs:
        plotParams.update(addLeptonPlotParams(*la))

    # dilepton objects
    dileptonArgs = []
    if analysis in ['Hpp4l']:
        dileptonArgs += [(savedir,'h1','hpp','\\ell^{+}\\ell^{+}',doDetailed]
        dileptonArgs += [(savedir,'h2','hmm','\\ell^{-}\\ell^{-}',doDetailed]
    if analysis in ['Hpp3l']:
        dileptonArgs += [(savedir,'h1','hpp','\\ell^{\\pm}\\ell^{\\pm}',doDetailed]
    if analysis in ['Hpp4l','ZZ','WZ','Hpp3l']:
        dileptonArgs += [(savedir,'z1','z1','\\ell^{+}\\ell^{-}',doDetailed]
    if analysis in ['Hpp4l','ZZ']:
        dileptonArgs += [(savedir,'z2','z2','\\ell^{+}\\ell^{-}',doDetailed]

    for da in dileptonArgs:
        plotParams.update(addDileptonPlotParams(*da))

    # add leptonMet objects
    leptonMetArgs = []
    if analysis in ['Hpp3l']:
        leptonMetArgs += [(savedir,'h2','hm','\\ell^{\\mp}',doDetailed]
    if analysis in ['Hpp3l','WZ']:
        leptonMetArgs += [(savedir,'w1','w1','W',doDetailed]

    for lma in leptonMetArgs:
        plotParams.update(addLeptonMetPlotParams(*lma))

    # add event params
    plotParams.update(addEventParams)

    # customization

    # customize stuff
    if 'z1Mass' in plotParams:
        plotParams['z1Mass']['args'][1] = [13,58.5,123.5]
        plotParams['z1Mass']['kwargs']['yaxis'] = 'Events/5 GeV'
    if 'z2Mass' in plotParams:
        plotParams['z2Mass']['args'][1] = [13,58.5,123.5]
        plotParams['z2Mass']['kwargs']['yaxis'] = 'Events/5 GeV'


    return plotParams

def getDefaultPlotterArgs():
    return {
        'cut'         : '',
        'xaxis'       : '',
        'yaxis'       : 'Events',
        'xrange'      : [],
        'xmin'        : None,
        'xmax'        : None,
        'ymin'        : None,
        'ymax'        : None,
        'overflow'    : False,
        'underflow'   : False,
        'nostack'     : False,
        'normalize'   : False,
        'blinder'     : [],
        'boxes'       : [],
        'logy'        : 0,
        'logx'        : 0,
        'nobg'        : 0,
        'lumitext'    : 11,
        'legendpos'   : 33,
        'numcol'      : 1,
        'signalscale' : 1,
        'isprelim'    : 1,
        'scalefactor' : '',
        'yscale'      : 1.25,
    }


