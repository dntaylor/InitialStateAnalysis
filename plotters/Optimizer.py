'''
Optimizer.py

Optimize selections.

Author: Devin N. Taylor, UW-Madison
'''
import math
import ROOT

class FunctionOptimizer( ROOT.TPyMultiGenFunction ):
    '''Multi-variable optimizer'''
    def __init__(self,args):
        self.plotter = args[0]
        self.sigSelection = args[1]
        self.bgSelection = args[2]
        self.cuts = args[3]
        self.formula = args[4]
        self.maximize = args[5]
        ROOT.TPyMultiGenFunction.__init__(self,self)

    def NDim(self): return len(self.cuts)

    def DoEval(self,args):
        sigSel = self.sigSelection
        bgSel = self.bgSelection
        for cut,arg in zip(self.cuts,args):
            sigSel += ' & %s %f' % (cut,arg)
            bgSel += ' & %s %f' % (cut,arg)
        sig = self.plotter.getSignalEntries(sigSel)
        bg = self.plotter.getBackgroundEntries(bgSel)
        #result = self.formula.EvalPar(sig,bg) if bg else sig
        result = sig/math.sqrt(bg) if bg else 0.
        print 'Current args:', [args[x] for x in range(len(self.cuts))], 'result:', result
        return -result if self.maximize else result

class Optimizer(object):
    '''Class to optimize a series of cuts'''
    def __init__(self,plotter):
        self.plotter = plotter
        self.cuts = []
        self.formula = ROOT.TFormula('f','[0]/sqrt([1])')

    def setSelection(self,selection,**kwargs):
        '''Set preselection for both background and signal.'''
        signalSelection = kwargs.pop('signalSelection','1')
        backgroundSelection = kwargs.pop('backgroundSelection','1')
        self.selection = selection
        self.signalSelection = signalSelection
        self.backgroundSelection = backgroundSelection

    def addCut(self,name,func,min,max,step):
        '''Add a cut to optimize'''
        self.cuts += [{'name': name, 'func': func, 'min': min, 'max': max, 'step': step}]

    def setSignificance(self,formula):
        '''Set the significance formula, variables are SIG and BG.'''
        formula.replace('SIG','[0]')
        formula.replace('BG','[1]')
        self.formula = ROOT.TFormula('f',formula)
        if self.formula.Compile()!=0: print 'Formula failed to compile: %s' % formula

    def optimize(self):
        '''Optimize selections'''
        sigSel = self.selection + ' & ' + self.signalSelection
        bgSel = self.selection + ' & ' + self.backgroundSelection
        # do the simple N-D optimization with minuit simplex method
        if False: # broken
            mini = ROOT.Math.Factory.CreateMinimizer("Minuit2", "Simplex")
            mini.SetMaxFunctionCalls(1000)
            mini.SetMaxIterations(1000)
            mini.SetTolerance(0.0001)
            mini.SetPrintLevel(1)
            args = [self.plotter,sigSel,bgSel,[cut['func'] for cut in self.cuts],self.formula,True]
            miniFunct = FunctionOptimizer(args)
            mini.SetFunction(miniFunct)
            for v,cut in enumerate(self.cuts):
                mini.SetLimitedVariable(v,cut['name'],float(cut['max']-cut['min'])/2.,10.*cut['step'],cut['min'],cut['max'])
            mini.Minimize()
        # get efficiencies
        for cut in self.cuts:
            print cut
            cutRange = [cut['min']+x*cut['step'] for x in range(int(cut['max']-cut['min']/cut['step']))]
            effs = [self.plotter.getSignalEntries('%s & %s %f' %(sigSel, cut['func'], cutVal))/self.plotter.getSignalEntries(sigSel) for cutVal in cutRange]
            print effs
            effCutoffs = [0.7,0.8,0.9,0.95,0.99,0.999]
            if effs[0]>effs[-1]: # reversed list
                for c in effCutoffs:
                    highEffs = [i for i,e in enumerate(effs) if e<c]
                    if highEffs: 
                        val = cutRange[highEffs[0]]
                    else:
                        val = -1
                    print '%s: Eff=%f: %f' % (cut['name'],c,val)
            else:
                for c in effCutoffs:
                    highEffs = [i for i,e in enumerate(effs) if e>c]
                    if highEffs: 
                        val = cutRange[highEffs[0]]
                    else:
                        val = -1
                    print '%s: Eff=%f: %f' % (cut['name'],c,val)


