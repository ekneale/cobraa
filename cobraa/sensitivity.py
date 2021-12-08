from decimal import *
from numpy import max

from ROOT import TH2D,TFile,TF1

from .load import *
from .coincidence import *
from .globals import *
from math import fabs
# This performs the final sensitivity calculations for the reactor analysis. 
# Author Liz Kneale (2021)
# Adapted from Watchmakers/sensitivity.py (Marc Bergevin)


setcontext(ExtendedContext)

def drange(start, stop, step):
    rii= start
    while rii<stop:
        yield rii
        rii+=step

def singlesrate():
    singlespersec = 0
    for _p in proc:
        for _loc in proc[_p]:
            for _element in d[_p][_loc]:
                if 'NA' in _p or 'RADIOGENIC' in _p:
                    singlespersec+=rates['%s_%s_%s'%(_element,_loc,_p)][0]
    print(singlespersec)
    return 0


def calculateSensitivity():

    hist = {}
    hacc = {}
    hrn = {}
    hfn = {}
    hibd = {}
    hibdBG = {}
    hli9 = {}
    hn17 = {}
    hgeo = {}

    print('Reading in root tree')
    # get the results of previous steps      
    if arguments['--core']:
        resultsstr = "core_root_files%s/coincidence_results.root"%(additionalString)
    else:
        resultsstr = "fred_root_files%s/coincidence_results.root"%(additionalString)
    resultsFile = TFile(resultsstr,"READ")
    print('reading in coincidence maps from %s'%(resultsstr))

    for delayed_nxcut,dTcut,maxEp,gcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT),drange(minEpmax,rangeEpmax,binwidthEpmax),drange(gmin,rangeGmax,binwidthG)):
    #for delayed_nxcut,dTcut,maxEp,gcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT),drange(gmin,rangeGmax,binwidthG)):

        hacc["hAcc%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)] = TH2D("hAcc%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10),'Accidental coincidence Rate',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hacc["hAcc%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetXTitle('distance from wall [m]')
        hacc["hAcc%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetYTitle('prompt %s cut'%(energyEstimator))
        hacc["hAcc%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetZTitle('accidental coincidence rate (Hz)')

        hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)] = TH2D("hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10),'IBD coincidence rate',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetXTitle('distance from wall [m]')
        hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetYTitle('prompt %s cut'%(energyEstimator))
        hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetZTitle('IBD rate (Hz)')

        hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)] = TH2D("hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10),'IBD BG coincidence rate',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetXTitle('distance from wall [m]')
        hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetYTitle('prompt %s cut'%(energyEstimator))
        hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetZTitle('IBD backgrounds rate (Hz)')

        hfn["hFN%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)] =  TH2D("hFN%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10),'Fast neutron coincidence rate ',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hfn["hFN%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetXTitle('distance from wall [m]')
        hfn["hFN%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetYTitle('prompt %s cut'%(energyEstimator))
        hfn["hFN%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetZTitle('Fast neutron rate (Hz)')

        # histograms for calculation of systematics
        hgeo["hGeo%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)] =  TH2D("hGeo%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10),'coincidence rate ',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hli9["hRNli%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)] =  TH2D("hRNli%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10),'coincidence rate ',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hn17["hRNn%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)] =  TH2D("hRNn%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10),'coincidence rate ',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        #hreac["hReac%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)] =  TH2D("hReac%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10),'coincidence rate ',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
      
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p][_loc]:
                    _tag = 'hist_%s_%s_%s_delayed%s_%d_%dus_maxEp%d_%d'%(_element,_loc,_p,energyEstimator,delayed_nxcut,dTcut,maxEp,gcut*10)
                    _tag = _tag.replace(" ","")
                    _str = "%s_%s_%s"%(_element,_loc,_p)
                    _str = _str.replace(" ","")
                    if 'pn_ibd' in _tag or 'A_Z' in _tag or 'FASTNEUTRONS' in _tag:
#                        print('correlated event, getting %s from %s\n'%(_tag,resultsstr))
                        hist[_tag] = resultsFile.Get(_tag)
                    elif 'singles' in _tag:
                        #print('uncorrelated event, getting %s from %s\n'%(_tag,resultsstr))
                        hist[_tag] = resultsFile.Get(_tag)
                    else:
                        continue

                    try:
                        #print(' entries:',hist[_tag].GetEntries(),' ... ', end = '')
                        if 'singles' in _tag: 
                            #print('identified as accidentals')
                            hacc["hAcc%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                        elif 'pn_ibd' in _tag:
                            #print('%s identified as IBD pair events\n\n'%(_tag))
                            if 'geo' in _tag:
                                hgeo["hGeo%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                            if arguments["--GSH"]:
                                if 'gravelines' in _tag or 'sizewell' in _tag or 'hinkley' in _tag:
                                    print('Adding %s to ibd signal'%(_tag))
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag:
                                    print('Adding %s to ibd bg'%(_tag))
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                            elif arguments["--SH"]:
                                if 'sizewell' in _tag or 'hinkley' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'gravelines' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                            elif arguments["--Gravelines"]:
                                if 'gravelines' in _tag:
                                    print('Adding %s to ibd signal'%(_tag))
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'sizewell' in _tag or 'hinkley' in _tag:
                                    print('Adding %s to ibd bg'%(_tag))
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                            elif arguments["--HinkleyC"]:
                                if 'hinkley' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'sizewell' in _tag or 'gravelines' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                            elif arguments["--Sizewell"]:
                                if 'sizewell' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'hinkley' in _tag or 'gravelines' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                            elif arguments["--Heysham"]:
                                if 'heysham_full' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'torness_full' in _tag or 'gravelines' in _tag or 'sizewell' in _tag or 'hinkley' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)

                            elif arguments["--HeyshamTorness"]:
                                if 'heysham_2' in _tag or 'torness_full' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'gravelines' in _tag or 'sizewell' in _tag or 'hinkley' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                    
                            elif arguments["--Heysham2"]:
                                if 'heysham_2' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'torness_full' in _tag or 'gravelines' in _tag or 'sizewell' in _tag or 'hinkley' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)

                            elif arguments["--Torness"]:
                                if 'torness_' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'heysham_2' in _tag or 'gravelines' in _tag or 'sizewell' in _tag or 'hinkley' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)

                            elif arguments['--Hartlepool1']:
                                if 'hartlepool_1' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'hartlepool_2' in _tag or 'heysham_full' in _tag or 'torness_full' in _tag or 'gravelines' in _tag or 'sizewell' in _tag or 'hinkley' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                            elif arguments['--Hartlepool2']:
                                if 'hartlepool_2' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'hartlepool_1' in _tag or 'heysham_full' in _tag or 'torness_full' in _tag or 'gravelines' in _tag or 'sizewell' in _tag or 'hinkley' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                            else:
                                if 'hartlepool_1' in _tag or 'hartlepool_2' in _tag:
                                    hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                                elif 'boulby_worldbg' in _tag or 'heysham_full' in _tag or 'torness_full' in _tag or 'gravelines' in _tag or 'sizewell' in _tag or 'hinkley' in _tag:
                                    hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)

                        elif 'FASTNEUTRONS' in _tag:
                            hfn["hFN%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                        elif 'A_Z' in _tag:
                            # scale by fraction surviving with 1 second after
                            # muon detected in fiducial (Study on the effect of varying veto
                            # thickness on the sensitivity of Gd-H2O filled tank, F. Sutanto),
                            # adjusted to 95% muon-detection efficiency with passive veto
                            # (Radionuclide Rates, E. Kneale 
                            # https://ait-neo.llnl.gov/confluence/pages/viewpage.action?pageId=10191240)
                            if 'li9' in _tag:
                                hist[_tag].Scale(0.069)
                                hli9["hRNli%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                            elif 'n17' in _tag:
                                hist[_tag].Scale(0.85)
                                hn17["hRNn%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Add(hist[_tag],1)
                        else:
                            print('not added to rates histogram')
                            continue

                    except:
                        print("Could not find ",_tag,". Skipping entry.")
    print("\nCompleted reading in of histogram, accidental and IBD identifitication. \n\n")

    optSignal,optBg,optSoB,optNXdelayed,optNXprompt,optDTW,optdT,optdR,optFidcut,optE,optG = -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
    optAcc,optAccErr,optFN,optFNerr,optRN,optRNerr,optIBDbg,optFnSysErr,optRnSysErr=-1,-1,-1,-1,-1,-1,-1,-1,-1
    optIBDbg, optSoB,optSigErr,optBgErr,optSoBErr,optTotSysErr,optt3sig = -1,-1,-1,-1,-1,-1,1000000.
    optReacibdBg, optGeoibdBg, optReacibdBgSysErr, optGeoibdBgSysErr = -1,-1,-1,-1

#    line,_line,_line2= ("",),"",""
    _histograms = {}

    for delayed_nxcut,dTcut,maxEp,gcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT),drange(minEpmax,rangeEpmax,binwidthEpmax),drange(gmin,rangeGmax,binwidthG)):
        
        binR = hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetNbinsX()
        binN = hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetNbinsY()

        _histograms["sOverB%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)]= hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Clone()
        _histograms["sOverB%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetZTitle('signal/sqrt(signal+background)')
        _histograms["sOverB%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetYTitle('%s cut on prompt'%(energyEstimator))
        _histograms["sOverB%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetTitle('hSoB - %s delayed %2d %d %dus %.1f'%(energyEstimator,delayed_nxcut,dTcut,maxEp,gcut))
        _histograms["sOverB%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetName('hSoB%d%d%d%d'%(delayed_nxcut,dTcut,maxEp,gcut*10))
        _histograms["sOverB%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Reset()

        _histograms["signal%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)]= hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Clone()
        _histograms["signal%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetZTitle('evts/day')
        _histograms["signal%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetYTitle('%s cut on prompt'%(energyEstimator))
        _histograms["signal%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetTitle('Signal - %s delayed %2d %d %dus %.1f'%(energyEstimator,delayed_nxcut,dTcut,maxEp,gcut))
        _histograms["signal%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetName('hSignal%d%d%d%d'%(delayed_nxcut,dTcut,maxEp,gcut*10))
        _histograms["signal%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Reset()
        _histograms["background%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)]= hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Clone()
        _histograms["background%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetZTitle('evts/day')
        _histograms["background%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetYTitle('%s cut on prompt'%(energyEstimator))
        _histograms["background%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetTitle('total backgrounds - %s delayed %2d %d %dus %1.f'%(energyEstimator,delayed_nxcut,dTcut,maxEp,gcut))
        _histograms["background%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetName('hBackground%d%d%d%d'%(delayed_nxcut,dTcut,maxEp,gcut*10))
        _histograms["background%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Reset()



        maxSoB,maxBg,maxSignal,maxBackground,maxNXd,maxNXp,maxFidcut,maxdT,maxdR,maxE,maxG = -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
        maxFnSysErr,maxIBDbgSysErr,maxRnSysErr,maxTotSysErr,maxt3sig,t3sig = -1,-1,-1,-1,1000000.,0

        for fidcut in drange(minFid,maxFid+binwidthFid,binwidthFid):
            for prompt_nxcut in drange(minNXprompt,maxNXprompt+binwidthNX,binwidthNX):
                # get the signal rates
                cutBin   = hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].FindBin(fidcut,prompt_nxcut)
                signal  = hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinContent(cutBin)
                signalError = hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinError(cutBin)

                # now get the background rates
                # first get the individual background contributions
                # at each cut value
                accRate  = hacc["hAcc%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinContent(cutBin)
                accError = hacc["hAcc%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinError(cutBin)
                fnRate   = hfn["hFN%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinContent(cutBin)
                fnError = hfn["hFN%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinError(cutBin)
                reacibdBgRate = hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinContent(cutBin)
                reacibdBgError = hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinError(cutBin)
                geoibdBgRate = hgeo["hGeo%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinContent(cutBin)
                geoibdBgError = hgeo["hGeo%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinError(cutBin)
                li9Rate = hli9["hRNli%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinContent(cutBin)
                li9Error = hli9["hRNli%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinError(cutBin)
                n17Rate = hn17["hRNn%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinContent(cutBin)
                n17Error = hn17["hRNn%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].GetBinError(cutBin)
                rnRate = li9Rate+n17Rate
                rnError = sqrt(li9Error**2+n17Error**2)

                # now we need to deal with systematic uncertainties
		# fractional systematic uncertainty on fast neutron rate (27%)
                fnSysError = fnRate*0.27
                # radionuclide systematics are somewhat more complex
                eMuSysError = (264/270.)**0.73*(log(264/270.)*0.1) #error on E=264GeV:0.73*264*sigmaEmu (not included); 0.1 error on 0.73 power: 264**0.73 * log(264)*0.1
                liSysError = li9Rate*sqrt(0.15e-8**2/4.09e-8+(eMuSysError**2)/((264/270.)**0.73))
                n17SysError = n17Rate*sqrt((0.15e-8**2)/4.09e-8+(eMuSysError**2)/((264/270.)**0.73))# error on muon flux and energy dependence
                rnSysError = sqrt(liSysError**2 + n17SysError**2)
                # fractional systematic uncertainties on IBD backgrounds (6% and 25%)
                reacibdBgSysError = reacibdBgRate*0.06
                geoibdBgSysError = geoibdBgRate*0.25

                # then get the total background rates (inc error rates)
                background = accRate + rnRate + reacibdBgRate +fnRate + geoibdBgRate # background rate
                bgError = sqrt(pow(accError,2)+pow(rnError,2)+pow(reacibdBgError,2)+pow(fnError,2)+pow(geoibdBgError,2))
                totSysError = sqrt(rnSysError*rnSysError+fnSysError*fnSysError+reacibdBgSysError*reacibdBgSysError+geoibdBgSysError*geoibdBgSysError)

		# calculate the signal over background and/or dwell time metric
                if arguments['--poissonpoisson']:
		    # poisson significance with poisson systematic uncertainty on backgrounds
                    sob = sqrt(2*((signal+background)*log((signal+background)*(background+ totSysError**2)/(background**2+(signal+background)*totSysError**2))-background**2/totSysError**2 * log(1 + (totSysError**2*signal/(background*(background + totSysError**2))))))
                elif arguments['--poisson']:
		    # poisson significance with gaussian systematic uncertainty on backgrounds
                    B_hat = 0.5*(background-pow(totSysError,2)+sqrt(pow(background,2)-2*totSysError*pow(totSysError,2)+4*(signal+background)*pow(totSysError,2)+pow(totSysError,4)))
                    sob = Z = sqrt(-2*((signal+background)*log(B_hat/(signal+background))-pow(background-B_hat,2)/(2*pow(totSysError,2))-B_hat+signal+background))
                #    t3sig = poissonTime(signal,background,totSysError)
                elif arguments['--knoll']:
		    # time to 3 sigma positive detection at 95% confidence level
                    if arguments['--2sigma']:
                        ton=((1.28*sqrt(2*background*(RonOff+1))+2*sqrt(signal*RonOff+2*background*(RonOff+1)))/(signal*RonOff))**2
                        toff=((1.28*sqrt(2*background*(1+1/RonOff))+2*sqrt(signal/RonOff+2*background*(1+1/RonOff)))/signal)**2
                    else:
                        ton=((1.645*sqrt(2*background*(RonOff+1))+3*sqrt(signal*RonOff+2*background*(RonOff+1)))/(signal*RonOff))**2
                        toff=((1.645*sqrt(2*background*(1+1/RonOff))+3*sqrt(signal/RonOff+2*background*(1+1/RonOff)))/signal)**2

                    t3sig=(2+RonOff)*toff
                    sob=-1
                else:
		    # gaussian significance
                    sob = signal/sqrt(background+totSysError**2)
                    t3sig = 9.*background/(signal**2-9.*totSysError**2)
                try:
		    # TODO get correct error for all metrics
                    sobError = sqrt(pow(signalError/signal,2)+pow(sqrt_s_plus_b_error/sqrt(signal+background),2))*sob
                except: 
                    sobError = 1

                _histograms["sOverB%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetBinContent(cutBin,sob)
                _histograms["signal%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetBinContent(cutBin,signal)
                _histograms["background%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetBinContent(cutBin,background)
                _histograms["sOverB%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetBinError(cutBin,sobError)
                _histograms["signal%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetBinError(cutBin,signalError)
                _histograms["background%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].SetBinError(cutBin,bgError)

                if arguments['--optimiseSoB']:
                    # find the maximum significance for each combination of NXd, maxEp, g cut
                    if sob >maxSoB:
                        maxt3sig = t3sig
                        maxTotSysErr = totSysError
                        maxSoB = sob
                        maxNXd = delayed_nxcut
                        maxNXp = prompt_nxcut
                        maxSignal = signal
                        maxBg = background
                        maxFidcut = fidcut
                        maxdT = dTcut
                        maxG = gcut
                        maxE = maxEp
                        line = ("beta/n: Wall-dist %4.1f m, %s cut (%d,%d) g cut %f, ibd rate : %4.2f per day;"\
                                %(maxFidcut,energyEstimator,maxNXp, maxNXd,maxG,signal),)
                        line2 =    (" acc. rate: %4.3f per day, rn rate: %4.3f per day, fn rate: %4.3f per day (all cuts)"\
                                %(accRate,rnRate,fnRate),)
                        # find the overall optimal significance
                    if sob>optSoB:
                        optt3sig = t3sig
                        optFnSysErr = fnSysError
                        optRnSysErr = rnSysError
                        optReacibdBgSysErr = reacibdBgSysError
                        optGeoibdBgSysErr = geoibdBgSysError
                        optSoB = sob
                        optSignal = signal
                        optBg     = background
                        optNXprompt  = prompt_nxcut
                        optNXdelayed  = delayed_nxcut
                        optFidcut = fidcut
                        optdT = dTcut
                        optG  = gcut
                        optE  = maxEp
                        optAcc = accRate
                        optRN = rnRate
                        optFN = fnRate
                        optReacibdBg = reacibdBgRate
                        optGeoibdBg = geoibdBgRate
                        optSigErr = signalError
                        optBgErr = bgError
                        optSoBErr = sobError
                        optTotSysErr = totSysError
                else:
                    # find the minimum dwell time for each combination of NXd, maxEp, g cut
                    if t3sig>0.5 and t3sig<maxt3sig:
                        maxt3sig = t3sig
                        maxTotSysErr = totSysError
                        maxSoB = sob
                        maxNXd = delayed_nxcut
                        maxNXp = prompt_nxcut
                        maxSignal = signal
                        maxBg = background
                        maxFidcut = fidcut
                        maxdT = dTcut
                        maxG = gcut
                        maxE = maxEp
                        line = ("beta/n: Wall-dist %4.1f m, %s cut (%d,%d) g cut %f, ibd rate : %4.2f per day;"\
                                %(maxFidcut,energyEstimator,maxNXp, maxNXd,maxG,signal),)
                        line2 =    (" acc. rate: %4.3f per day, rn rate: %4.3f per day, fn rate: %4.3f per day (all cuts)"\
                                %(accRate,rnRate,fnRate),)
                        # find the overall optimal s/sqrt(s+b)
                    if t3sig>0.5 and t3sig<optt3sig:
                        optt3sig = t3sig
                        optFnSysErr = fnSysError
                        optRnSysErr = rnSysError
                        optReacibdBgSysErr = reacibdBgSysError
                        optGeoibdBgSysErr = geoibdBgSysError
                        optSoB = sob
                        optSignal = signal
                        optBg     = background
                        optNXprompt  = prompt_nxcut
                        optNXdelayed  = delayed_nxcut
                        optFidcut = fidcut
                        optdT = dTcut
                        optG  = gcut
                        optE  = maxEp
                        optAcc = accRate
                        optRN = rnRate
                        optFN = fnRate
                        optReacibdBg = reacibdBgRate
                        optGeoibdBg = geoibdBgRate
                        optSigErr = signalError
                        optBgErr = bgError
                        optSoBErr = sobError
                        optTotSysErr = totSysError




        # print to screen the optimal values for each combination of delayed nx, goodness, maxEp and dT
        print('Dwell time:','{0:.1f}'.format(maxt3sig),'Delayed nx:',str(maxNXd).rjust(3,' '),'prompt nx:',str(maxNXp).rjust(3,' '),'dT:',str(maxdT).rjust(3,' '),'g:',str(maxG),',S/sqrt(B+sigma_b^2)','{0:.4f}'.format(maxSoB),',(S,B,dtw):','{0:.4f}'.format(maxSignal),'{0:.4f}'.format(maxBg),'{0:.1f}'.format(maxFidcut),'{0:.1f}'.format(maxE),")")
        #line += (_line + _line2,)


    print('\n\nMore info on the maximal sensitivity found:')
    # print line

    for _l in line:
        for i in range(len(_l)):
            print(_l[i], end=' ')
        print('')

    print("1-month significance: ",optSignal*30/sqrt(optBg*180+(optTotSysErr*180)**2))
    print("3-month significance: ",optSignal*90/sqrt(optBg*180+(optTotSysErr*180)**2))
    print("6-month significance: ",optSignal*180/sqrt(optBg*180+(optTotSysErr*180)**2))
    print("9-month significance: ",optSignal*270/sqrt(optBg*180+(optTotSysErr*180)**2))
    print("12-month significance: ",optSignal*360/sqrt(optBg*360+(optTotSysErr*360)**2))
    print("18-month significance: ",optSignal*540/sqrt(optBg*540+(optTotSysErr*540)**2))
    print("24-month significance: ",optSignal*720/sqrt(optBg*720+(optTotSysErr*720)**2))
    print("30-month significance: ",optSignal*900/sqrt(optBg*900+(optTotSysErr*900)**2))
    
    # calculate dwell time for 3 sigma significance, where sigma = s/sqrt(b)
    if optSignal >0:
        if arguments['--poisson']:
            T3SIGMA = poissonTime(optSignal,optBg,optTotSysErr)
        else:
            T3SIGMA = optt3sig #9.*optBg/(optSignal**2-9.*optTotSysErr**2)
    else:
        T3SIGMA = 1e999
    
    # write the optimal values to the results_coincidence* file
    result = "%s:\nfiducial:%.1f  nXprompt:%d  nXdelayed:%d  dT:%d maxE: %d  g:%.1f  \ns:%.3f+/-%1.2e(stat)  b:%.3f+/-%1.2e(stat) \nacc:%1.3e  rn:%1.2e +/- %1.2e(sys)  fn:%1.2e +/- %1.2e(sys)  \nreacIBDbg:%1.2e +/- %1.2e(sys)    geoIBDbg:%1.2e +/- %1.2e(sys) \ns/sqrt(b):%.3f+/-%1.2e(stat)  days to discovery:%.1f  total bg systematics: %1.2e" %(additionalString,optFidcut,optNXprompt,optNXdelayed,optdT,optE,optG,optSignal,optSigErr,optBg,optBgErr,optAcc,optRN,optRnSysErr,optFN,optFnSysErr,optReacibdBg,optReacibdBgSysErr,optGeoibdBg,optGeoibdBgSysErr,optSoB,optSoBErr,T3SIGMA,optTotSysErr)
    resultsfile = "results_coincidence%s.txt"%(additionalString)
    with open(resultsfile,'a') as file:
        file.write(result+'\n')
    # write the signal, background and s/sqrt(s+b) histograms to the analysis_results_coincidence* file
    if arguments['--core']:
        analysisfile = "core_root_files%s/sensitivity_results.root"%(additionalString)
    else:
        analysisfile = "fred_root_files%s/sensitivity_results.root"%(additionalString)
    print('\n\nWriting histograms to file',analysisfile)
    f_root = TFile(analysisfile,"recreate")
    for delayed_nxcut,dTcut,maxEp,gcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT),drange(minEpmax,rangeEpmax,binwidthEpmax),drange(gmin,rangeGmax,binwidthG)):
        _histograms["sOverB%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Write()
        _histograms["signal%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Write()
        _histograms["background%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Write()
        hacc["hAcc%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Write()
        hfn["hFN%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Write()
        hibd["hIBD%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Write()
        hibdBG["hIBDBG%d%d%d%d"%(delayed_nxcut,dTcut,maxEp,gcut*10)].Write()

    f_root.Close()

    print('\n\n')

def bhat(x,par):
    Bhat = 0.5*(par[2]*x[0]-pow(par[1]*x[0],2)+sqrt(pow(par[2]*x[0],2)-2*par[1]*x[0]*pow(par[1]*x[0],2)+4*(par[0]*x[0]+par[2]*x[0])*pow(par[1]*x[0],2)+pow(par[1]*x[0],4)))
    return Bhat

def poissonSig(x,par):
    B_hat = bhat(x,par)
    Z = sqrt(-2*((par[0]*x[0]+par[2]*x[0])*log(B_hat/(par[0]*x[0]+par[2]*x[0]))-pow(par[2]*x[0]-B_hat,2)/(2*pow(par[1]*x[0],2))-B_hat+par[0]*x[0]+par[2]*x[0]))
    return Z

def poissonTime(s,b,syserr):
    
    # with gaussian uncertainty on background
    Zfunc = TF1("Zfunc",poissonSig,0.1,100,3)
    Zfunc.SetParameter(0,s)
    Zfunc.SetParameter(1,syserr)
    Zfunc.SetParameter(2,b)
    Zfunc.SetParName(0,"sig")
    Zfunc.SetParName(1,"syserr")
    Zfunc.SetParName(2,"bg")
    t3sigma = Zfunc.GetX(3.)
    return t3sigma

def poissonpoissonTime(s,b,syserr):
        # with poisson uncertainty on background
    Z = TF1("Z","sqrt(2*(([0]*x+[2]*x)*log(([0]*x+[2]*x)*([2]*x + pow(x*[1],2))/(pow([2]*x,2)+([0]*x+[2]*x)*pow(x*[1],2)))-pow([2]*x,2)/pow([1]*x,2) * log(1 + (pow(x*[1],2)*x*[0]/([2]*x*([2]*x + pow(x*[1],2)))))))",0,50)
    Z.SetNpx(500)
    Z.SetParameter(0,s)
    Z.SetParameter(1,syserr)
    Z.SetParameter(2,b)
    t3sigma = Z.GetX(3.)
    return t3sigma

def gaussianTime(s,b,syserr):
    sigma = TF1("sigma","[0]*x/sqrt([2]*x+pow((x*[1]),2))",0,1000)
    sigma.SetNpx(10000)
    sigma.SetParameter(0,s)
    sigma.SetParameter(1,syserr)
    sigma.SetParameter(2,b)
    t3sigma = sigma.GetX(3.)
    return t3sigma

def ClTime(s,b,syserr):


    return t3sigma
