from decimal import *
from numpy import max

from ROOT import TH2D,TFile

from .load import *
from .coincidence import *
from .globals import *


# This performs the final sensitivity calculations for the reactor analysis. 
# Author Liz Kneale (2021)
# Adapted from Watchmakers/sensitivity.py (Marc Bergevin)


setcontext(ExtendedContext)

def drange(start, stop, step):
    rii= start
    while rii<stop:
        yield rii
        rii+=step


def calculateSensitivity():

    hist = {}
    hacc = {}
    hrn = {}
    hfn = {}
    hibd = {}
    hibdBG = {}

    print('Reading in root tree')
    ''' Following loop has two purposes. Read in the histogram from results file and scale them with appropriate rates'''
    # get the results of previous steps      
    resultsstr = "fred_root_files%s/coincidence_results.root"%(additionalString)
    resultsFile = TFile(resultsstr,"READ")
    print('reading in coincidence maps from %s'%(resultsstr))

    for delayed_nxcut,dTcut,dRcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT),drange(dRmin,rangedRmax,binwidthdR)):

        hacc["hAcc%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)] = TH2D("hAcc%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000),'Accidental coincidence Rate',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hacc["hAcc%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetXTitle('distance from wall [m]')
        hacc["hAcc%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetYTitle('prompt %s cut'%(energyEstimator))
        hacc["hAcc%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetZTitle('accidental coincidence rate (Hz)')

        hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)] = TH2D("hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000),'IBD coincidence rate',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetXTitle('distance from wall [m]')
        hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetYTitle('prompt %s cut'%(energyEstimator))
        hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetZTitle('IBD rate (Hz)')

        hibdBG["hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)] = TH2D("hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000),'IBD BG coincidence rate',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hibdBG["hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetXTitle('distance from wall [m]')
        hibdBG["hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetYTitle('prompt %s cut'%(energyEstimator))
        hibdBG["hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetZTitle('IBD backgrounds rate (Hz)')

        hrn["hRN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)] = TH2D("hRN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000),'Radionuclide coincidence rate ',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hrn["hRN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetXTitle('distance from wall [m]')
        hrn["hRN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetYTitle('prompt %s cut'%(energyEstimator))
        hrn["hRN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetZTitle('Radionuclide rate (Hz)')

        hfn["hFN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)] =  TH2D("hFN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000),'Fast neutron coincidence rate ',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hfn["hFN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetXTitle('distance from wall [m]')
        hfn["hFN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetYTitle('prompt %s cut'%(energyEstimator))
        hfn["hFN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetZTitle('Fast neutron rate (Hz)')
      
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p][_loc]:
                    _tag = 'hist_%s_%s_%s_delayed%scut%d_%dus_%dmm'%(_element,_loc,_p,energyEstimator,delayed_nxcut,dTcut,dRcut*1000)
                    _tag = _tag.replace(" ","")
                    _str = "%s_%s_%s"%(_element,_loc,_p)
                    _str = _str.replace(" ","")
                    if 'pn_ibd' in _tag or 'A_Z' in _tag or 'FASTNEUTRONS' in _tag:
                        print('correlated event, getting %s from %s\n'%(_tag,resultsstr))
                        hist[_tag] = resultsFile.Get(_tag)
                    elif 'SINGLES' in _tag:
                        print('uncorrelated event, getting %s from %s\n'%(_tag,resultsstr))
                        hist[_tag] = resultsFile.Get(_tag)
                    else:
                        continue

                    try:
                        print(' entries:',hist[_tag].GetEntries(),' ... ', end = '')
                        print(_tag)
                        if 'SINGLES' in _tag: 
                            print('identified as accidentals')
                            hacc["hAcc%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Add(hist[_tag],1)
                        if 'pn_ibd' in _tag:
                            print('%s identified as IBD pair events ,'%(_tag))
                            if arguments["--Heysham"]:
                                if 'heysham_signal' in _tag:
                                    hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Add(hist[_tag],1)
                                elif 'heysham_background' in _tag:
                                    hibdBG["hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Add(hist[_tag],1)
                            else:
                                if 'big_' in _tag or 'small' in _tag:
                                    hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Add(hist[_tag],1)
                                elif 'boulby' in _tag:
                                    hibdBG["hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Add(hist[_tag],1)
                        if 'FASTNEUTRONS' in _tag:
                            hfn["hFN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Add(hist[_tag],1)
                            print('identified as fast neutrons')
                        if 'A_Z' in _tag:
                            # if veto present, scale by fraction surviving with 1 second dead time
                            # of muon detected in veto (Study on the effect of varying veto
                            # thickness on the sensitivity of Gd-H2O ifilled tank, F. Sutanto)
                            # TODO currently removing for passive veto. Recalculate for active veto
                            #if 'li9' in _tag:
                            #    if 'default' in additionalString:
                            #        hist[_tag].Scale(0.03) currently removing for passive veto
                            #if 'n17' in _tag:
                            #    if 'default' in additionalString:
                            #        hist[_tag].Scale(0.85) currently removing for passive veto
                            print('identified as RN events, ')
                            hrn["hRN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Add(hist[_tag],1)
                        else:
                            print('not added to rates histogram')
                            continue

                    except:
                        print("Could not find ",_tag,". Skipping entry.")
    print("\nCompleted reading in of histogram, accidental and IBD identifitication. \n\n")

    optSignal,optBg,optSoB,optNXdelayed,optNXprompt,optDTW,optdT,optdR = -1,-1,-1,-1,-1,-1,-1,-1
    optAcc,optAccErr,optFN,optFNerr,optRN,optRNerr,optIBDbg=-1,-1,-1,-1,-1,-1,-1
    optIBDbg, optSoB,optSoB30days,optSoB30daysErr,optSigErr,optBgErr,optSoBErr = -1,-1,-1,-1,-1,-1,-1

    line,_line,_line2= ("",),"",""
    _histograms = {}

    for delayed_nxcut,dTcut,dRcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT),drange(dRmin,rangedRmax,binwidthdR)):
        
        binR = hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetNbinsX()
        binN = hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetNbinsY()

        _histograms["sOverB%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)]= hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Clone()
        _histograms["sOverB%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetZTitle('signal/sqrt(signal+background)')
        _histograms["sOverB%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetYTitle('%s cut on prompt'%(energyEstimator))
        _histograms["sOverB%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetTitle('hSoB - %s delayed %2d %dus %dmm'%(energyEstimator,delayed_nxcut,dTcut,dRcut*1000))
        _histograms["sOverB%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetName('hSoB%d%d%d'%(delayed_nxcut,dTcut,dRcut*1000))
        _histograms["sOverB%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Reset()

        _histograms["signal%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)]= hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Clone()
        _histograms["signal%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetZTitle('evts/day')
        _histograms["signal%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetYTitle('%s cut on prompt'%(energyEstimator))
        _histograms["signal%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetTitle('Signal - %s delayed %2d %dus %dmm'%(energyEstimator,delayed_nxcut,dTcut,dRcut*1000))
        _histograms["signal%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetName('hSignal%d%d%d'%(delayed_nxcut,dTcut,dRcut*1000))
        _histograms["signal%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Reset()
        _histograms["background%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)]= hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Clone()
        _histograms["background%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetZTitle('evts/day')
        _histograms["background%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetYTitle('%s cut on prompt'%(energyEstimator))
        _histograms["background%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetTitle('total backgrounds - %s delayed %2d %dus %dmm'%(energyEstimator,delayed_nxcut,dTcut,dRcut*1000))
        _histograms["background%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetName('hBackground%d%d%d'%(delayed_nxcut,dTcut,dRcut*1000))
        _histograms["background%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Reset()




        maxSoB,maxSignal,maxBackground,maxNXd,maxNXp,maxFidcut,maxdT,maxdR = -1,-1,-1,-1,-1,-1,-1,-1

        for fidcut in drange(minFid,maxFid+binwidthFid,binwidthFid):
            for prompt_nxcut in drange(minNXprompt,maxNXprompt+binwidthNX,binwidthNX):

                # get the signal rates
                cutBin   = hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].FindBin(fidcut,prompt_nxcut)
                ibdRate  = hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinContent(cutBin)
                signal   = ibdRate
                ibdError = hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinError(cutBin)
                signalError = ibdError

                # now get the background rates
                # first get the individual background contributions
                # at each cut value
                accRate  = hacc["hAcc%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinContent(cutBin)
                accError = hacc["hAcc%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinError(cutBin)
                fnRate   = hfn["hFN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinContent(cutBin)
                fnError = hfn["hFN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinError(cutBin)
                rnRate   = hrn["hRN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinContent(cutBin)
                rnError = hrn["hRN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinError(cutBin)
                ibdBgRate   = hibdBG["hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinContent(cutBin)
                ibdBgError = hibdBG["hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].GetBinError(cutBin)

                # then get the total background rates
                background = accRate + rnRate + ibdBgRate +fnRate # background rate
                bgError = sqrt(pow(accError,2)+pow(rnError,2)+pow(ibdBgError,2)+pow(fnError,2))

                # save the rates to the analysis results histograms
                sob = signal/sqrt(background)
                s_plus_b_error = sqrt(pow(signalError,2)+pow(bgError,2))
                sqrt_s_plus_b_error = 0.5*s_plus_b_error/(signal+background)*sqrt(signal+background)
                # error on s/sqrt(s+b)
                #sobError = sqrt(pow(signalError/signal,2)+pow(sqrt_s_plus_b_error/sqrt(signal+background),2))*sob
                # error on s/sqrt(b)
                sqrt_b_error = 0.5*bgError/background*sqrt(background)
                sobError = sqrt(pow(signalError/signal,2)+pow(sqrt_b_error/sqrt(background),2))
                _histograms["sOverB%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetBinContent(cutBin,sob)
                _histograms["signal%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetBinContent(cutBin,signal)
                _histograms["background%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetBinContent(cutBin,background)
                _histograms["sOverB%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetBinError(cutBin,sobError)
                _histograms["signal%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetBinError(cutBin,signalError)
                _histograms["background%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].SetBinError(cutBin,bgError)

                # find the maximum s/sqrt(b) for each value of delayed nxcut
                if sob >maxSoB:
                    maxSoB = sob
                    maxNXd = delayed_nxcut
                    maxNXp = prompt_nxcut
                    maxSignal = signal
                    maxBg = background
                    maxFidcut = fidcut
                    maxdT = dTcut
                    line = ("beta/n: Wall-dist %4.1f m, %s cut (%d,%d), ibd rate : %4.2f per day;"\
                            %(maxFidcut,energyEstimator,maxNXp, maxNXd,ibdRate),)
                    line2 =    (" acc. rate: %4.3f per day, rn rate: %4.3f per day, fn rate: %4.3f per day (all cuts)"\
                            %(accRate,rnRate,fnRate),)
                # find the overall optimal s/sqrt(s+b)
                if sob >optSoB:
                    optSoB = sob
                    optSignal = signal
                    optBg     = background
                    optNXprompt  = prompt_nxcut
                    optNXdelayed  = delayed_nxcut
                    optFidcut = fidcut
                    optdT = dTcut
                    optdR = dRcut
                    optAcc = accRate
                    optRN = rnRate
                    optFN = fnRate
                    optIBDbg = ibdBgRate
                    optSoB30daysErr = sobError*sqrt(30)
                    optSoB30days = sob*sqrt(30)
                    optSigErr = signalError
                    optBgErr = bgError
                    optSoBErr = sobError


        print('Delayed nx:',str(maxNXd).rjust(3,' '),'prompt nx:',str(maxNXp).rjust(3,' '),'dT:',str(maxdT).rjust(3,' '),',Found max S/sqrt(S+B)','{0:.4f}'.format(maxSoB),',(S,B,%s,dtw):('%(energyEstimator),'{0:.4f}'.format(maxSignal),'{0:.4f}'.format(maxBg),'{0:.1f}'.format(maxFidcut),')')
        line += (_line + _line2,)


    print('\n\nMore info on the maximal sensitivity found:')
    # print line

    for _l in line:
        for i in range(len(_l)):
            print(_l[i], end=' ')
        print('')

    # calculate dwell time for 3 sigma significance, where sigma = s/sqrt(b)
    if optSignal >0:
        T3SIGMA = 9/optSignal + 9*optBg/optSignal**2

    else:
        T3SIGMA = 1e999
    
    # write the optimal values to the results_coincidence* file
    result = "%s: \nOptimal cuts: %4.1f %3d %3d %dus %1.2fmm \nsig:%4.3f +/-%1.2e, bg:%4.3f +/-%1.2e, acc:%1.2e, rn:%1.2e, fn:%1.2e, ibdbg:%1.2e  \nsob:%4.3f +/-%1.2e sob30day:%4.3f +/-%1.2e dwell:%4.1f " %(additionalString,optFidcut,optNXprompt,optNXdelayed,optdT,optdR,optSignal,optSigErr,optBg,optBgErr,optAcc,optRN,optFN,optIBDbg,optSoB,optSoBErr,optSoB30days,optSoB30daysErr,T3SIGMA)
    resultsfile = "results_coincidence%s.txt"%(additionalString)
    with open(resultsfile,'a') as file:
        file.write(result+'\n')
    # write the signal, background and s/sqrt(s+b) histograms to the analysis_results_coincidence* file
    analysisfile = "fred_root_files%s/sensitivity_results.root"%(additionalString)
    print('\n\nWriting histograms to file',analysisfile)
    f_root = TFile(analysisfile,"recreate")
    for delayed_nxcut in drange(minNXdelayed,maxNXdelayed+binwidthNX,binwidthNX):
        _histograms["sOverB%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Write()
        _histograms["signal%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Write()
        _histograms["background%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Write()
        hacc["hAcc%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Write()
        hfn["hFN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Write()
        hrn["hRN%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Write()
        hibd["hIBD%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Write()
        hibdBG["hIBDBG%d%d%d"%(delayed_nxcut,dTcut,dRcut*1000)].Write()

    f_root.Close()

    print('\n\n')
