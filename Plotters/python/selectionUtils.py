
selections = {
    # (analysis, region, period) : {},
    ('WZ','WZ',13) : {
        # mode : {},
        'SMP-15-006' : {
            # cutname : {},
            'zmass' : {
                'logic' : 'AND',
                'cuts' : [
                    'z1.mass > 60.',
                    'z1.mass < 120.',
                ],
            },
            'mass3l' : {
                'cuts' : 'finalstate.mass>100.',
            },
            'zpt' : {
                'logic' : 'AND',
                'cuts' : [
                    'z1.Pt1>20.',
                    'z2.Pt1>10.',
                ],
            },
            'wmll' : {
                'logic' : 'AND',
                'cuts' : [
                    'w1.mll_z1_1>4',
                    'w1.mll_z1_2>4',
                ],
            },
            'wpt' : {
                'cuts' : 'w1.Pt1>20.',
            },
            'met' : {
                'cuts' : 'finalstate.met>30.',
            },
        },
    },
}


logicMap = {
    'AND' : '&&',
    'OR'  : '||',
}

def getSelection(analysis,region,period,mode):
    '''Get the analysis selections for a given working point.'''
    if (analysis,region,period) not in selections:
        print 'Invalid analysis'
        return '1'
    modes = selections[(analysis,region,period)]
    if mode not in modes:
        print 'Invalid mode {0}'.format(mode)
        return '1'
    allcuts = modes[mode]
    cutstrings = []
    for cut in sorted(allcuts.keys()):
        if isinstance(allcuts[cut]['cuts'], basestring):
            cutstrings += ['('+allcuts[cut]['cuts']+')']
        else:
            logic = logicMap[allcuts[cut]['logic'].upper()]
            ls = ' {0} '.format(logic)
            cutstrings += ['('+ls.join(allcuts[cut]['cuts'])+')']
    fullCut = ' && '.join(cutstrings)
    return fullCut
