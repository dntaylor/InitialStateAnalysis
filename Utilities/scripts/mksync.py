#!/usr/bin/env python

from InitialStateAnalysis.Plotters.Plotter import Plotter
from InitialStateAnalysis.Plotters.plotUtils import *
from InitialStateAnalysis.Plotters.plotUtils import ZMASS
import argparse
import itertools
import sys

def sync(analysis,channel,period,**kwargs):
    '''Print sync information to file.'''
    runTau = kwargs.pop('runTau',False)
    blind = kwargs.pop('blind',True)
    doEfficiency = kwargs.pop('doEfficiency',False)
    doCutflow = kwargs.pop('doCutflow',False)
    doYields = kwargs.pop('doYields',False)
    doCounts = kwargs.pop('doCounts',False)
    doCorrelation = kwargs.pop('doCorrelation',False)
    doCategories = kwargs.pop('doCategories',False)
    doChargeId = kwargs.pop('doChargeId',False)
    do4l = kwargs.pop('do4l',False)
    doFinalStates = kwargs.pop('doFinalStates',False)
    doDataDriven = kwargs.pop('doDataDriven',True)
    cut = kwargs.pop('cut','1')
    bp = kwargs.pop('bp','')
    mass = kwargs.pop('mass',500)

    print ''
    print '%s:%s:%iTeV' % (analysis, channel, period)
    print 'Selection to be used:'
    print cut
    if analysis in ['Hpp3l','Hpp4l']:
        print 'Mass %i' % mass
    print ''
    sys.stdout.flush()

    # sync on WZ sample
    # sync on channels: eee, eem, emm, mmm
    ntuples = 'ntuples/%s_%sTeV_%s' % (analysis,period,channel)
    saves = '%s_%s_%sTeV' % (analysis,channel,period)
    mergeDict = getMergeDict(period)
    nl = 3 if analysis == 'WZ' or analysis == 'Hpp3l' else 4
    finalStates, leptons = getChannels(nl)
    if analysis in ['WZ']: finalStates = ['eee','eem','mme','mmm']
    s = 'Sig'
    if analysis in ['WZ']: s = 'WZ'
    if do4l: s = 'SigPP'
    sigMap = getSigMap(nl,mass)
    channelBackground =  getChannelBackgrounds(period)
    plotter = Plotter(channel,ntupleDir=ntuples,saveDir=saves,period=period,mergeDict=mergeDict,dataDriven=doDataDriven)
    plotter.initializeBackgroundSamples([sigMap[period][x] for x in channelBackground[channel]])
    if analysis in ['Hpp3l', 'Hpp4l']: plotter.initializeSignalSamples([sigMap[period][s]])
    plotter.initializeDataSamples([sigMap[period]['data']])
    intLumi = getIntLumiMap()[period]
    plotter.setIntLumi(intLumi)
    cutflowMap = defineCutFlowMap(analysis,finalStates,mass)
    cutflow = {}
    for c,l in zip(cutflowMap['cuts'],cutflowMap['labels_simple']):
        cutflow[l] = c
    cutflows = cutflowMap['labels_simple']
    allMC = [x for x in channelBackground[channel] if x not in ['Z','TT','T']] + ['datadriven'] if doDataDriven else channelBackground[channel]
    if analysis in ['Hpp3l', 'Hpp4l']: allMC += [s]

    if doCounts:
        print '{0} signal event counts'.format(analysis)
        tempCut = cut
        for chan in finalStates:
            num = plotter.getNumEntries('%s&channel=="%s"' %(tempCut,chan), sigMap[period][s], doUnweighted=True, do4l=do4l, bp=bp)
            print '%s: %10.2f' % (chan, num)
        print ''
        sys.stdout.flush()

    if doCategories:
        labels = {0:'F',1:'T'}
        cuts = {0:'{0}.PassTight{1}==0',1:'{0}.PassTight{1}==1'}
        yields = {}
        yieldsMC = {}
        yieldsWeightedData = {}
        yieldsWeightedMC = {}
        yieldsWeighted = {}
        nameMap = {
            0: 'z1.PassTight1',
            1: 'z1.PassTight2',
            2: 'w1.PassTight1',
        }
        effMap = {
            0: 'z1.LepEffTight1',
            1: 'z1.LepEffTight2',
            2: 'w1.LepEffTight1',
        }
        fakeMap = {
            0: 'z1.LepFake1',
            1: 'z1.LepFake2',
            2: 'w1.LepFake1',
        }

        if analysis in ['Hpp3l']:
            nameMap = {
                0: 'h1.PassTight1',
                1: 'h1.PassTight2',
                2: 'h2.PassTight1',
            }
            effMap = {
                0: 'h1.LepEffTight1',
                1: 'h1.LepEffTight2',
                2: 'h2.LepEffTight1',
            }
            fakeMap = {
                0: 'h1.LepFake1',
                1: 'h1.LepFake2',
                2: 'h2.LepFake1',
            }

        allCuts = 'finalstate.mass>100. && (z1.Pt1>20.&&z1.Pt2>10.) && z1.mass>60. && z1.mass<120. && w1.dR1_z1_1>0.1 && w1.dR1_z1_2>0.1 && w1.Pt1>20. && finalstate.met>30.'
        if analysis in ['Hpp3l']:
            allCuts = '1' # no special cuts at select.passTight
        myCut = cut
        for l in range(3):
            myCut = myCut.replace(nameMap[l],'1')
            myCut = myCut.replace('l{0}.PassTight'.format(l),'1')
        myCut = myCut.replace('select.passTight',allCuts)
        #isocuts = ' && ((l{0}.Iso<0.4 && l{0}Flv=="e") || (l{0}.Iso<1.0 && l{0}Flv=="m"))'
        #for l in range(3):
        #    myCut += isocuts.format(l+1)
        for chan in finalStates:
            yields[chan] = {}
            yieldsMC[chan] = {}
            yieldsWeighted[chan] = {}
            yieldsWeightedMC[chan] = {}
            yieldsWeightedData[chan] = {}
            for l1 in [0,1]:
                for l2 in [0,1]:
                    for l3 in [0,1]:
                        thisCut = '{0} && channel=="{1}" && '.format(myCut,chan) + ' && '.join([cuts[l1].format('z1','1'),cuts[l2].format('z1','2'),cuts[l3].format('w1','1')])
                        if analysis=='Hpp3l': thisCut = '{0} && channel=="{1}" && '.format(myCut,chan) + ' && '.join([cuts[l1].format('h1','1'),cuts[l2].format('h1','2'),cuts[l3].format('h2','1')])
                        #if not (l1 and l2 and l3): thisCut += ' && finalstate.bjetVeto30Medium==0'
                        thisName = ''.join([labels[l1],labels[l2],labels[l3]])
                        yields[chan][thisName]  = plotter.getDataEntries(thisCut)
                        num = '1' if thisName.count('F') in [1,3] else '-1'
                        for l in range(3):
                            if thisName[l] == 'F': num += '*{0}'.format(fakeMap[l])
                        scalefactor = num
                        if thisName=='TTT': scalefactor = '1'
                        yieldsWeightedData[chan][thisName] = plotter.getDataEntries(thisCut,customScale=scalefactor)
                        yieldsWeighted[chan][thisName] = yieldsWeightedData[chan][thisName]
                        if thisName=='TTT':
                            yieldsWeightedMC[chan][thisName] = 0.
                        else:
                            yieldsWeightedMC[chan][thisName] = 0.
                            for bg in allMC:
                                if bg in ['TT','Z','ZG','Sig']: continue
                                mcscale = plotter.getScaleFactor() + '*' + scalefactor
                                yieldsWeightedMC[chan][thisName] += plotter.getNumEntries(thisCut,sigMap[period][bg],customScale=mcscale)
                            yieldsWeighted[chan][thisName] -= yieldsWeightedMC[chan][thisName]
                        # yields in each channel
                        if doYields:
                            yieldsMC[chan][thisName] = {}
                            #print 'MC Yields (scaled to %i pb-1)' % intLumi
                            #theCut = '&'.join([cutflow[x] for x in cutflows])
                            #print 'Selection applied: %s' % thisCut
                            #print '%8s |    Channel |      Yield' % chan
                            for b in allMC:
                                val = plotter.getNumEntries('%s&channel=="%s"' %(thisCut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                                yieldsMC[chan][thisName][b] = val
                                #print '         | %10s | %10.2f' %(b,val)
                            val = plotter.getDataEntries('%s&channel=="%s"' %(thisCut,chan), do4l=do4l, bp=bp)
                            yieldsMC[chan][thisName]['data'] = val
                            #print '         | %10s | %10i' %('Data',val)
                            #print ''
                            sys.stdout.flush()

                            

        print 'Data events for each category'
        print 'CAT | {0}'.format(' | '.join(['{0:8}'.format(x) for x in finalStates + ['Total']]))
        for cat in ['TTT','TTF','TFT','FTT','TFF','FTF','FFT','FFF']:
            print '{0} | {1}'.format(cat,' | '.join(['{0:8d}'.format(int(yields[x][cat])) for x in finalStates])) + ' | {0:8d}'.format(sum([int(yields[x][cat]) for x in finalStates]))
        print ''

        print 'Weighted events: Data'
        print 'CAT | {0}'.format(' | '.join(['{0:8}'.format(x) for x in finalStates + ['Total']]))
        for cat in ['TTT','TTF','TFT','FTT','TFF','FTF','FFT','FFF']:
            print '{0} | {1}'.format(cat,' | '.join(['{0:8.2f}'.format(yieldsWeightedData[x][cat]) for x in finalStates])) + ' | {0:8.2f}'.format(sum([yieldsWeightedData[x][cat] for x in finalStates]))
        print ''
        
        print 'Weighted events: MC'
        print 'CAT | {0}'.format(' | '.join(['{0:8}'.format(x) for x in finalStates + ['Total']]))
        for cat in ['TTT','TTF','TFT','FTT','TFF','FTF','FFT','FFF']:
            print '{0} | {1}'.format(cat,' | '.join(['{0:8.2f}'.format(yieldsWeightedMC[x][cat]) for x in finalStates])) + ' | {0:8.2f}'.format(sum([yieldsWeightedMC[x][cat] for x in finalStates]))
        print ''

        print 'Weighted events: Total'
        print 'CAT | {0}'.format(' | '.join(['{0:8}'.format(x) for x in finalStates + ['Total']]))
        for cat in ['TTT','TTF','TFT','FTT','TFF','FTF','FFT','FFF']:
            print '{0} | {1}'.format(cat,' | '.join(['{0:8.2f}'.format(yieldsWeighted[x][cat]) for x in finalStates])) + ' | {0:8.2f}'.format(sum([yieldsWeighted[x][cat] for x in finalStates]))
        print ''
        sys.stdout.flush()

        if doYields:
            for chan in finalStates:
                print 'MC Yields: {0}'.format(chan)
                print 'CAT | {0}'.format(' | '.join(['{0:8}'.format(x) for x in allMC + ['Data']]))
                for cat in ['TTT','TTF','TFT','FTT','TFF','FTF','FFT','FFF']:
                    print '{0} | {1}'.format(cat,' | '.join(['{0:8.2f}'.format(yieldsMC[chan][cat][x]) for x in allMC + ['data']]))
                print ''


    # yields in each channel
    if doYields:
        print 'Yields (scaled to %i pb-1)' % intLumi
        #theCut = '&'.join([cutflow[x] for x in cutflows])
        print 'Selection applied: %s' % cut
        for chan in finalStates:
            print '%8s |    Channel |      Yield' % chan
            for b in allMC:
                val = plotter.getNumEntries('%s&channel=="%s"' %(cut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                print '         | %10s | %10.2f' %(b,val)
            val = plotter.getDataEntries('%s&channel=="%s"' %(cut,chan), do4l=do4l, bp=bp)
            print '         | %10s | %10i' %('Data',val)
            print ''
            sys.stdout.flush()

    # cut flow
    if doCutflow:
        print 'Cutflows'
        for chan in finalStates:
            print '%15s |        Sig |         BG |        S/B |       Data' % chan
            for c in range(len(cutflows)):
                sig = 0
                bg = 0
                theCut = '&'.join([cutflow[x] for x in cutflows[0:c+1]])
                for b in allMC:
                    val = plotter.getNumEntries('%s&channel=="%s"' %(theCut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                    if b==s: sig += val
                    elif 'data' in b: pass
                    else: bg += val
                data = plotter.getDataEntries('%s&channel=="%s"' %(theCut,chan), do4l=do4l, bp=bp)
                print '%15s | %10.2f | %10.2f | %10.2f | %10i' %(cutflows[c],sig,bg,sig/bg,data)
            print ''
            sys.stdout.flush()

    # electron charge id efficiency
    if doChargeId:
        print 'Charge ID efficiency'
        bgNum = 0.
        bgNum_err_2 = 0.
        bgDenom = 0.
        bgDenom_err_2 = 0.
        sigNum = 0.
        sigNum_err_2 = 0.
        sigDenom = 0.
        sigDenom_err_2 = 0.
        dataNum = 0.
        dataNum_err_2 = 0.
        dataDenom = 0.
        dataDenom_err_2 = 0.
        for b in allMC + ['data']:
            val01, err01 = plotter.getNumEntries('%s & l1Flv=="e" & l1.ChargeConsistent==0' %cut, sigMap[period][b], doError=True, do4l=do4l, bp=bp)
            val02, err02 = plotter.getNumEntries('%s & l2Flv=="e" & l2.ChargeConsistent==0' %cut, sigMap[period][b], doError=True, do4l=do4l, bp=bp)
            val03, err03 = plotter.getNumEntries('%s & l3Flv=="e" & l3.ChargeConsistent==0' %cut, sigMap[period][b], doError=True, do4l=do4l, bp=bp)
            val11, err11 = plotter.getNumEntries('%s & l1Flv=="e" & l1.ChargeConsistent==1' %cut, sigMap[period][b], doError=True, do4l=do4l, bp=bp)
            val12, err12 = plotter.getNumEntries('%s & l2Flv=="e" & l2.ChargeConsistent==1' %cut, sigMap[period][b], doError=True, do4l=do4l, bp=bp)
            val13, err13 = plotter.getNumEntries('%s & l3Flv=="e" & l3.ChargeConsistent==1' %cut, sigMap[period][b], doError=True, do4l=do4l, bp=bp)
            val0 = val01+val02+val03
            val1 = val11+val12+val13
            err0_2 = err01**2 + err02**2 + err03**2
            err1_2 = err11**2 + err12**2 + err13**2
            if b==s:
                sigNum += val1
                sigNum_err_2 += err1_2
                sigDenom += val0 + val1
                sigDenom_err_2 += err0_2 + err1_2
            elif 'data' in b:
                dataNum += val1
                dataNum_err_2 += err1_2
                dataDenom += val0 + val1
                dataDenom_err_2 += err0_2 + err1_2
            else:
                bgNum += val1
                bgNum_err_2 += err1_2
                bgDenom += val0 + val1
                bgDenom_err_2 += err0_2 + err1_2
        sigEff = sigNum/sigDenom
        sigErr = (sigNum_err_2/sigNum**2 + sigDenom_err_2/sigDenom**2)**0.5 * sigEff
        bgEff = bgNum/bgDenom
        bgErr = (bgNum_err_2/bgNum**2 + bgDenom_err_2/bgDenom**2)**0.5 * bgEff
        dataEff = float(dataNum)/dataDenom
        dataErr = (dataNum_err_2/dataNum**2 + dataDenom_err_2/dataDenom**2)**0.5 * dataEff
        print 'Sig Efficiency: %f +/- %f, BG Efficiency: %f +/- %f, Data Efficiency: %f +/- %f' % (sigEff, sigErr, bgEff, bgErr, dataEff, dataErr)
        ptBins = [20.,30.,40.,60.,100.,1000.]
        for p in range(len(ptBins)-1):
            bgNum = 0.
            bgNum_err_2 = 0.
            bgDenom = 0.
            bgDenom_err_2 = 0.
            sigNum = 0.
            sigNum_err_2 = 0.
            sigDenom = 0.
            sigDenom_err_2 = 0.
            dataNum = 0.
            dataNum_err_2 = 0.
            dataDenom = 0.
            dataDenom_err_2 = 0.
            ptlow = ptBins[p]
            pthigh = ptBins[p+1]
            for b in allMC + ['data']:
                val01, err01 = plotter.getNumEntries('%s && l1Flv=="e" && l1.ChargeConsistent==0 && l1.Pt > %f && l1.Pt < %f' %(cut,ptlow,pthigh), sigMap[period][b], doError=True, do4l=do4l, bp=bp)
                val02, err02 = plotter.getNumEntries('%s && l2Flv=="e" && l2.ChargeConsistent==0 && l2.Pt > %f && l2.Pt < %f' %(cut,ptlow,pthigh), sigMap[period][b], doError=True, do4l=do4l, bp=bp)
                val03, err03 = plotter.getNumEntries('%s && l3Flv=="e" && l3.ChargeConsistent==0 && l3.Pt > %f && l3.Pt < %f' %(cut,ptlow,pthigh), sigMap[period][b], doError=True, do4l=do4l, bp=bp)
                val11, err11 = plotter.getNumEntries('%s && l1Flv=="e" && l1.ChargeConsistent==1 && l1.Pt > %f && l1.Pt < %f' %(cut,ptlow,pthigh), sigMap[period][b], doError=True, do4l=do4l, bp=bp)
                val12, err12 = plotter.getNumEntries('%s && l2Flv=="e" && l2.ChargeConsistent==1 && l2.Pt > %f && l2.Pt < %f' %(cut,ptlow,pthigh), sigMap[period][b], doError=True, do4l=do4l, bp=bp)
                val13, err13 = plotter.getNumEntries('%s && l3Flv=="e" && l3.ChargeConsistent==1 && l3.Pt > %f && l3.Pt < %f' %(cut,ptlow,pthigh), sigMap[period][b], doError=True, do4l=do4l, bp=bp)
                val0 = val01+val02+val03
                val1 = val11+val12+val13
                err0_2 = err01**2 + err02**2 + err03**2
                err1_2 = err11**2 + err12**2 + err13**2
                if b==s:
                    sigNum += val1
                    sigNum_err_2 += err1_2
                    sigDenom += val0 + val1
                    sigDenom_err_2 += err0_2 + err1_2
                elif 'data' in b:
                    dataNum += val1
                    dataNum_err_2 += err1_2
                    dataDenom += val0 + val1
                    dataDenom_err_2 += err0_2 + err1_2
                else:
                    bgNum += val1
                    bgNum_err_2 += err1_2
                    bgDenom += val0 + val1
                    bgDenom_err_2 += err0_2 + err1_2
            sigEff = sigNum/sigDenom
            sigErr = (sigNum_err_2/sigNum**2 + sigDenom_err_2/sigDenom**2)**0.5 * sigEff
            bgEff = bgNum/bgDenom
            bgErr = (bgNum_err_2/bgNum**2 + bgDenom_err_2/bgDenom**2)**0.5 * bgEff
            dataEff = float(dataNum)/dataDenom
            dataErr = (dataNum_err_2/dataNum**2 + dataDenom_err_2/dataDenom**2)**0.5 * dataEff
            sf = dataEff/bgEff
            sfErr = ((dataErr/dataEff)**2 + (bgErr/bgEff)**2)**0.5 * sf
            print 'Pt range: [%i, %i], Sig Efficiency: %f +/- %f, BG Efficiency: %f +/- %f, Data Efficiency: %f +/- %f, Data/BG: %f +/- %f' % (int(ptlow),int(pthigh),sigEff,sigErr,bgEff,bgErr,dataEff,dataErr,sf,sfErr)
        print ''

    # efficiency of cuts
    if doEfficiency:
        print 'Cut efficiencies'
        sigPreCuts = {}
        bgPreCuts = {}
        sigFullCuts = {}
        bgFullCuts = {}
        allSigPre = 0
        allBgPre = 0
        allSigFull = 0
        allBgFull = 0
        for chan in finalStates:
            sigPre = 0
            bgPre = 0
            sigFull = 0
            bgFull = 0
            theFullCut = '&'.join([cutflow[x] for x in cutflows])
            for b in allMC:
                valPre = plotter.getNumEntries('%s&channel=="%s"' %(cut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                valFull = plotter.getNumEntries('%s&channel=="%s"' %(theFullCut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                if b==s:
                    sigPre += valPre
                    sigFull += valFull
                elif 'data' in b: pass
                else:
                    bgPre += valPre
                    bgFull += valFull
            sigPreCuts[chan] = sigPre
            bgPreCuts[chan] = bgPre
            sigFullCuts[chan] = sigFull
            bgFullCuts[chan] = bgFull
            allSigPre += sigPre
            allBgPre += bgPre
            allSigFull += sigFull
            allBgFull += bgFull
        sigPreCuts['all'] = allSigPre
        bgPreCuts['all'] = allBgPre
        sigFullCuts['all'] = allSigFull
        bgFullCuts['all'] = allBgFull


        for c in cutflows[1:]:
            print '%15s |  Sg Pre Eff |  BG Pre Eff | Sg Post Eff | BG Post Eff' % c
            allSigAllbut = 0
            allBgAllbut = 0
            allSigOnly = 0
            allBgOnly = 0
            for chan in finalStates:
                sigAllbut = 0
                bgAllbut = 0
                sigOnly = 0
                bgOnly = 0
                theCut = '&'.join([cutflow[x] for x in cutflows if x != c])
                theOnlyCut = '%s&%s' % (cut, cutflow[c])
                for b in allMC:
                    valAllbut = plotter.getNumEntries('%s&channel=="%s"' %(theCut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                    valOnly = plotter.getNumEntries('%s&channel=="%s"' %(theOnlyCut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                    if b==s:
                        sigAllbut += valAllbut
                        sigOnly += valOnly
                    elif 'data' in b: pass
                    else:
                        bgAllbut += valAllbut
                        bgOnly += valOnly
                allSigAllbut += sigAllbut
                allBgAllbut += bgAllbut
                allSigOnly += sigOnly
                allBgOnly += bgOnly
                sigEffPre = sigOnly/sigPreCuts[chan] if sigPreCuts[chan] else -1
                bgEffPre = bgOnly/bgPreCuts[chan] if bgPreCuts[chan] else -1
                sigEffPost = sigFullCuts[chan]/sigAllbut if sigAllbut else -1
                bgEffPost = bgFullCuts[chan]/bgAllbut if bgAllbut else -1
                print '%15s | %11.4f | %11.4f | %11.4f | %11.4f' % (chan, sigEffPre, bgEffPre, sigEffPost, bgEffPost)
            chan = 'all'
            sigEffPre = allSigOnly/sigPreCuts[chan] if sigPreCuts[chan] else -1
            bgEffPre = allBgOnly/bgPreCuts[chan] if bgPreCuts[chan] else -1
            sigEffPost = sigFullCuts[chan]/allSigAllbut if allSigAllbut else -1
            bgEffPost = bgFullCuts[chan]/allBgAllbut if allBgAllbut else -1
            print '%15s | %11.4f | %11.4f | %11.4f | %11.4f' % (chan, sigEffPre, bgEffPre, sigEffPost, bgEffPost)
            print ''

        print 'Cut efficiencies'
        allPre = {}
        allFull = {}
        for chan in finalStates:
            allPre[chan] = {}
            allFull[chan] = {}
            theFullCut = '&'.join([cutflow[x] for x in cutflows])
            for b in allMC:
                valPre = plotter.getNumEntries('%s&channel=="%s"' %(cut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                valFull = plotter.getNumEntries('%s&channel=="%s"' %(theFullCut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                allPre[chan][b] = valPre
                allFull[chan][b] = valFull
        allPre['all'] = {}
        allFull['all'] = {}
        for b in allMC:
            allPre['all'][b] = sum([allPre[x][b] for x in finalStates])
            allFull['all'][b] = sum([allFull[x][b] for x in finalStates])

        for c in cutflows[1:]:
            pre = ' | '.join(['{0:11}'.format(x+' Pre') for x in allMC])
            post = ' | '.join(['{0:11}'.format(x+' Post') for x in allMC])
            print '{0:15} | {1} | {2}'.format(c,pre,post)
            allBut = {}
            allOnly = {}
            allPreEff = {}
            allPostEff = {}
            for chan in finalStates:
                allBut[chan] = {}
                allOnly[chan] = {}
                allPreEff[chan] = {}
                allPostEff[chan] = {}
                theCut = '&'.join([cutflow[x] for x in cutflows if x != c])
                theOnlyCut = '%s&%s' % (cut, cutflow[c])
                for b in allMC:
                    valAllbut = plotter.getNumEntries('%s&channel=="%s"' %(theCut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                    valOnly = plotter.getNumEntries('%s&channel=="%s"' %(theOnlyCut,chan), sigMap[period][b], do4l=do4l, bp=bp)
                    allBut[chan][b] = valAllbut
                    allOnly[chan][b] = valOnly
                    allPreEff[chan][b] = allOnly[chan][b]/allPre[chan][b] if allPre[chan][b] else -1
                    allPostEff[chan][b] = allFull[chan][b]/allBut[chan][b] if allBut[chan][b] else -1
                pre = ' | '.join(['{0:11.4f}'.format(allPreEff[chan][x]) for x in allMC])
                post = ' | '.join(['{0:11.4f}'.format(allPostEff[chan][x]) for x in allMC])
                print '{0:15} | {1} | {2}'.format(chan, pre, post)
            allBut['all'] = {}
            allOnly['all'] = {}
            allPreEff['all'] = {}
            allPostEff['all'] = {}
            for b in allMC:
                allBut['all'][b] = sum([allBut[x][b] for x in finalStates])
                allOnly['all'][b] = sum([allOnly[x][b] for x in finalStates])
                allPreEff['all'][b] = allOnly['all'][b]/allPre['all'][b] if allPre['all'][b] else -1
                allPostEff['all'][b] = allFull['all'][b]/allBut['all'][b] if allBut['all'][b] else -1
            pre = ' | '.join(['{0:11.4f}'.format(allPreEff['all'][x]) for x in allMC])
            post = ' | '.join(['{0:11.4f}'.format(allPostEff['all'][x]) for x in allMC])
            print '{0:15} | {1} | {2}'.format('all', pre, post)
            print ''




    if doCorrelation:
        print 'Correlation matrix'
        for chan in finalStates:
            print ' | '.join(['%15s'%(chan+' denom.')] + ['%15s' % x for x in cutflows[1:]])
            for dc in cutflows[1:]:
                row = []
                for nc in cutflows[1:]:
                    numCut = '%s & %s & %s' % (cut, cutflow[dc], cutflow[nc])
                    denomCut = '%s & %s' % (cut, cutflow[dc])
                    numSig = plotter.getNumEntries('%s&channel=="%s"' %(numCut,chan), sigMap[period][s], do4l=do4l, bp=bp)
                    denomSig = plotter.getNumEntries('%s&channel=="%s"' %(denomCut,chan), sigMap[period][s], do4l=do4l, bp=bp)
                    eff = numSig/denomSig if denomSig else -1
                    row += [eff]
                print ' | '.join(['%15s'%dc] + ['%15.4f' % x for x in row])
            print ''


def parse_command_line(argv):
    parser = argparse.ArgumentParser(description="Plot a given channel and period")

    parser.add_argument('analysis', type=str, choices=['WZ','Hpp3l','Hpp4l'], help='Analysis to plot')
    parser.add_argument('channel', type=str, choices=['WZ','Hpp3l','Hpp4l','FakeRate'], help='Channel in analysis')
    parser.add_argument('period', type=int, choices=[7,8,13], help='Energy (TeV)')
    parser.add_argument('-m','--mass',nargs='?',type=int,const=500,default=500)
    parser.add_argument('-rt','--runTau',action='store_true',help='Run Tau finalStates (not implemented)')
    parser.add_argument('-ub','--unblind',action='store_false',help='Unblind signal channel')
    parser.add_argument('-dc','--doCounts',action='store_true',help='run sig counts')
    parser.add_argument('-dcat','--doCategories',action='store_true',help='categories for bg estimation')
    parser.add_argument('-dd','--doDataDriven',action='store_true',help='run datadriven')
    parser.add_argument('-dcf','--doCutflow',action='store_true',help='run cutflow')
    parser.add_argument('-de','--doEfficiency',action='store_true',help='run efficiencies')
    parser.add_argument('-dy','--doYields',action='store_true',help='run yields')
    parser.add_argument('-dco','--doCorrelation',action='store_true',help='run correlation')
    parser.add_argument('-df','--do4l',action='store_true',help='do PP part of Hpp3l')
    parser.add_argument('-dfs','--doFinalStates',action='store_true',help='do individual channels')
    parser.add_argument('-bp','--branchingPoint',nargs='?',type=str,const='',default='',choices=['ee100','em100','mm100','et100','mt100','tt100','BP1','BP2','BP3','BP4',''],help='Choose branching point for H++')
    parser.add_argument('-c','--cut',type=str,default='select.passTight',help='Cut to be applied to plots (default = "select.passTight").')
    args = parser.parse_args(argv)

    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_command_line(argv)

    if args.period == 7:
        print "7 TeV not implemented"

    # loose no iso
    # e: loose CBID
    # m: isLooseMuon

    # loose w/ iso
    # e: loose CBID and iso<0.2
    # m: isLooseMuon and iso<0.2

    # tight
    # e: medium CBID and iso<0.15
    # m: isTightMuon and iso<0.12

    sync(args.analysis,args.channel,args.period,
         runTau=args.runTau,
         blind=args.unblind,
         cut=args.cut,
         mass=args.mass,
         doEfficiency=args.doEfficiency,
         doCutflow=args.doCutflow,
         doYields=args.doYields,
         doCounts=args.doCounts,
         doCorrelation=args.doCorrelation,
         do4l=args.do4l,
         bp=args.branchingPoint,
         doFinalStates=args.doFinalStates,
         doCategories=args.doCategories,
         doDataDriven=args.doDataDriven
    )

    return 0


if __name__ == "__main__":
    main()
