#!/usr/bin/env python

from InitialStateAnalysis.Plotters.plotUtils import _3L_MASSES, _4L_MASSES
import os
import glob

# a script to generate the needed tables for the hpp analysis note

tableHeader = '\
\\begin{table}[h]\n\
\\begin{footnotesize}\n\
    \\centering\n\
'

tabularHeader = '\
    \\begin{tabular}{STRUCTURE}\n\
    \\hline\n\
'

tabularFooter = '\
    \\end{tabular}\n\
'

caption = '\
    \\caption{CAPTION}\n\
    \\label{tab:LABEL}\n\
'

tableFooter = '\
\\end{footnotesize}\n\
\\end{table}\n\
'

tableStrings = {'ee100' : '100\\% decay to electrons',
                'em100' : '100\\% decay to $e\\mu$ pairs',
                'mm100' : '100\\% decay to muons',
                'et100' : '100\\% decay to $e\\tau$ pairs',
                'mt100' : '100\\% decay to $\\mu\\tau$ pairs',
                'tt100' : '100\\% decay to $\\tau\\tau$ pairs',
		'BP1' : 'benchmark point 1',
                'BP2' : 'benchmark point 2',
                'BP3' : 'benchmark point 3',
                'BP4' : 'benchmark point 4'}

fsLabels = {'Hpp3l' : '$\\ell\\ell\\ell$    ',
            'Hpp4l' : '$\\ell\\ell\\ell\\ell$'}

tables = {}

# | Signal                      | Background Estimation        | Data     |
# | Mass Chan | AP      PP      | N_SB alpha BG Est. | MC Est. | Observed |
# | x    x    | x +/- y x +/- y | N    a     x +/- y | x +/- y | x +/- y  |
# ...
labels0 = '        \\multicolumn{4}{|l|}{Signal}                                                & \\multicolumn{4}{|l|}{Background Estimation}                               & Data                    \\\\ \\hline\n'
labels1 = '        Mass & Final state      & AP                      & PP                      & $N_{SB}$  & $\\alpha$  & Data-driven             & MC Estimate             & Observed                \\\\ \\hline\n'
labels0 = '        \\multicolumn{3}{|l|}{Signal}                                 & \\multicolumn{2}{|l|}{Background Estimation}            & Data                      \\\\ \\hline\n'
labels1 = '        Mass & AP                        & PP                        & Data-driven                & MC Estimate               & Observed                  \\\\ \\hline\n'

rows = {
    'pp'  : '        {mass:<4} & {recochan:16} & ---                     & {pp:9.3f} \\pm {ppe:9.3f} & {nsb:9d} & {a:9.3g} & {bg:9.3g} \\pm {bge:9.3g} & {mc:9.3g} \\pm {mce:9.3f} & {d:9d} \\pm {de:9.3}',
    'ap'  : '        {mass:<4} & {recochan:16} & {ap:9.3f} \\pm {ape:9.3f} & ---                       & {nsb:9d} & {a:9.3g} & {bg:9.3g} \\pm {bge:9.3g} & {mc:9.3g} \\pm {mce:9.3f} & {d:9d} \\pm {de:9.3}',
    'both': '        {mass:<4} & {recochan:16} & {ap:9.3f} \\pm {ape:9.3f} & {pp:9.3f} \\pm {ppe:9.3f} & {nsb:9d} & {a:9.3g} & {bg:9.3g} \\pm {bge:9.3g} & {mc:9.3g} \\pm {mce:9.3f} & {d:9d} \\pm {de:9.3}',
}

rows = {
    'pp'  : '        {mass:<4} & ---                       & ${pp:9.3f} \\pm {ppe:9.3f}$ &  ${bg:9.3g} \\pm {bge:9.3g}$ & ${mc:9.3g} \\pm {mce:9.3f}$ & ${d:9d} \\pm {de:9.3}$',
    'ap'  : '        {mass:<4} & ${ap:9.3f} \\pm {ape:9.3f}$ & ---                         &  ${bg:9.3g} \\pm {bge:9.3g}$ & ${mc:9.3g} \\pm {mce:9.3f}$ & ${d:9d} \\pm {de:9.3}$',
    'both': '        {mass:<4} & ${ap:9.3f} \\pm {ape:9.3f}$ & ${pp:9.3f} \\pm {ppe:9.3f}$ &  ${bg:9.3g} \\pm {bge:9.3g}$ & ${mc:9.3g} \\pm {mce:9.3f}$ & ${d:9d} \\pm {de:9.3}$',
}

valuenames = {
    'pp'  : 'datacards/{analysis}_8tev_{channel}/{bp}/{mass}/values_{bp}_{recochan}.txt',
    'ap'  : 'datacards/{analysis}_8tev_{channel}/{bp}/{mass}/values_{bp}_{recochan}.txt',
    'both': 'datacards/{analysis}_8tev_{channel}/{bp}/{mass}/values_{bp}_{recochan}_APandPP.txt',
}
alphavaluenames = {
    'pp'  : 'datacards/{analysis}_8tev_{channel}/{bp}/{mass}/alphavalues_{bp}_{recochan}.txt',
    'ap'  : 'datacards/{analysis}_8tev_{channel}/{bp}/{mass}/alphavalues_{bp}_{recochan}.txt',
    'both': 'datacards/{analysis}_8tev_{channel}/{bp}/{mass}/alphavalues_{bp}_{recochan}_APandPP.txt',
}

channels = {
    'Hpp3l': ['eee','eem','eme','emm','mme','mmm'],
    'Hpp4l': ['eeee','eeem','eemm','emem','emmm','mmmm'],
}

names = {
    'e' : 'e',
    'm' : '\\mu',
    't' : '\\tau',
}

for analysis,channel in [('Hpp3l','Hpp3l'),('Hpp4l','Hpp4l')]:
    for bp in ['ee100', 'em100', 'mm100', 'et100', 'mt100', 'BP1', 'BP2', 'BP3', 'BP4']:
        replacements = {'STRUCTURE' : '|l|ll|ll|l|',
                        'CAPTION' : 'Background estimation for %s-lepton final states for %s.' % ('three' if analysis=='Hpp3l' else 'four',tableStrings[bp]),
                        'LABEL' : 'bg%s_%s' % (bp,analysis)}
        
        elementstring = ''
        for mass in _4L_MASSES:
            massblockstring = ''
            valdict = {}
            valfile = valuenames['both'] if analysis=='Hpp3l' else valuenames['pp']
            alphafile = alphavaluenames['both'] if analysis=='Hpp3l' else alphavaluenames['pp']
            recochans = channels[analysis]
            for recochan in recochans:
                vfname = valfile.format(analysis=analysis,channel=channel,bp=bp,mass=mass,recochan=recochan)
                afname = alphafile.format(analysis=analysis,channel=channel,bp=bp,mass=mass,recochan=recochan)
                if os.path.isfile(vfname) and os.path.isfile(afname):
                    #massblockstring += ' \n'
                    valdict[recochan] = {}
                    valdict[recochan]['mass'] = mass
                    valdict[recochan]['recochan'] = '$' + ''.join([names[x] for x in recochan]) + '$'
                    with open(vfname,'r') as f:
                        text = f.read()
                        vals = [float(x) for x in text.split(':')]
                        valdict[recochan]['mc'] = vals[0]
                        valdict[recochan]['mce'] = vals[1]
                        valdict[recochan]['d'] = int(vals[4])
                        valdict[recochan]['de'] = vals[5]
                        valdict[recochan]['pp'] = vals[7]
                        valdict[recochan]['ppe'] = vals[8]
                        if analysis in ['Hpp3l'] and mass in _3L_MASSES:
                            valdict[recochan]['ap'] = vals[9]
                            valdict[recochan]['ape'] = vals[10]
                    with open(afname,'r') as f:
                        atext = f.read()
                        avals = [float(x) for x in atext.split(':')]
                        nSB,eSB,nSR,eSR,alpha,nSBData,eSBData,nBGSR,eBGSR,nSRData,eSRData = avals
                        valdict[recochan]['nsb'] = int(nSBData)
                        valdict[recochan]['a'] = alpha
                        valdict[recochan]['bg'] = nBGSR
                        valdict[recochan]['bge'] = eBGSR
                    #rowstring = rows['both'] if analysis=='Hpp3l' and mass in _3L_MASSES else rows['pp']
                    #recochanstring = rowstring.format(**valdict[recochan])
                    #massblockstring += recochanstring
                    #massblockstring += ' \\\\'
            chans = valdict.keys()
            valdict['all'] = {}
            valdict['all']['mass'] = mass
            valnames =  ['mc','mce','d','de','bg','bge','pp','ppe']
            if analysis in ['Hpp3l'] and mass in _3L_MASSES: valnames += ['ap','ape']
            for key in valnames:
                if 'e' in key:
                    valdict['all'][key] = sum([valdict[x][key]**2 for x in chans])**0.5
                else:
                    valdict['all'][key] = sum([valdict[x][key] for x in chans])
            #massblockstring += ' \\hline \n'
            rowstring = rows['both'] if analysis=='Hpp3l' and mass in _3L_MASSES else rows['pp']
            recochanstring = rowstring.format(**valdict['all'])
            massblockstring += recochanstring
            massblockstring += ' \\\\ \n' if mass != _4L_MASSES[-1] else ' \\\\ \\hline \n'
            elementstring += massblockstring

        tables['bg%s_%s' % (bp,analysis)] = tableHeader + tabularHeader + labels0 + labels1 + elementstring + tabularFooter + caption + tableFooter
        for key, val in replacements.iteritems():
            tables['bg%s_%s' % (bp,analysis)] = tables['bg%s_%s' % (bp,analysis)].replace(key,val)

        print tables['bg%s_%s' % (bp,analysis)]

