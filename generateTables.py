#!/usr/bin/env python

from plotters.plotUtils import _3L_MASSES, _4L_MASSES
import os

# a script to generate the needed tables for the hpp analysis note

tableHeader = '\
\\begin{table}\n\
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
                'em100' : '100\\% decay to $e\mu$ pairs',
                'mm100' : '100\\% decay to muons',
		'BP1' : 'benchmark point 1',
                'BP2' : 'benchmark point 2',
                'BP3' : 'benchmark point 3',
                'BP4' : 'benchmark point 4'}

fsLabels = {'Hpp3l' : '$\\ell\\ell\\ell$    ',
            'Hpp4l' : '$\\ell\\ell\\ell\\ell$'}

tables = {}

for mode in ['ee100', 'em100', 'mm100', 'BP1', 'BP2', 'BP3', 'BP4']:
    # Background estimation
    
    # 7 tables, 3 for 100% decay, 4 for BP
    
    # Mass, final state, MC estimate, sideband estimate, observation, pair production, associated production
    replacements = {'STRUCTURE' : '|l|lll|l|',
                    #'CAPTION' : 'Background estimation for %s.' % tableStrings[mode],
                    'CAPTION' : 'Background estimation for three-lepton final states for %s.' % tableStrings[mode],
                    'LABEL' : 'bg%s' % mode}
    #labels = '        Mass    & Final state        & MC estimate     & Sideband estimate & Observation              & Pair prod.      & Assoc. prod.      \\\\ \\hline\n'
    labels = '        Mass    & MC estimate       & Sideband estimate   & Signal            & Observation                \\\\ \\hline\n'
    elementStrings = ''
    for mass in _3L_MASSES:
        #for fs in ['Hpp3l', 'Hpp4l']:
        for fs in [('Hpp3l','Hpp3l')]:
            #rowString = '        %i GeV & %s & ' % (mass,fsLabels[fs])
            rowString = '        %i GeV & ' % mass
            filename = 'datacards/%s_8tev_%s/%s/%i/values.txt' % (fs[0],fs[1],mode,mass)
            if os.path.isfile(filename):
                file = open(filename)
                text = file.read()
                file.close()
                vals = [float(x) for x in text.split(':')]
            else:
                vals = [0]*11
            rowString += '$%6.3g \\pm %6.3g$ & $%6.3g \\pm %6.3g$   & $%6.3g \\pm %6.3g$ & $%6.3g \\pm %6.3g \\pm %6.3g$' % (vals[0], vals[1], vals[2], vals[3], vals[9], vals[10], vals[4], vals[5], vals[6])
            #if fs=='4l':
            #    rowString += ' \\\\ \\hline'
            #else:
            #    rowString += ' \\\\'
            rowString += ' \\\\'
            elementStrings += '%s\n' % rowString
    elementStrings += '        \\hline\n'

    tables['bg%s' % mode] = tableHeader + tabularHeader + labels + elementStrings + tabularFooter + caption + tableFooter
    for key, val in replacements.iteritems():
        tables['bg%s' % mode] = tables['bg%s' % mode].replace(key,val)

    print tables['bg%s' % mode]

