from decimal import *
from numpy import max

from ROOT import kOrange as kO,kBlue as kB,kGreen as kG
from ROOT import kMagenta as kM,kAzure as kA,kRed as kR
from ROOT import TCanvas,TH2D, gStyle,TFile,gROOT

from .load import *
from .coincidence import *
from .globals import *



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
    resultsstr = "core_root_files%s/coincidence_results.root"%(additionalString)
    resultsFile = TFile(resultsstr,"READ")
    print('reading in coincidence maps from %s'%(resultsstr))

    for delayed_nxcut,dTcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT)):

        hacc["hAcc%d%d"%(delayed_nxcut,dTcut)] = TH2D("hAcc%d%d"%(delayed_nxcut,dTcut),'Accidental coincidence Rate',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hacc["hAcc%d%d"%(delayed_nxcut,dTcut)].SetXTitle('distance from wall [m]')
        hacc["hAcc%d%d"%(delayed_nxcut,dTcut)].SetYTitle('prompt %s cut'%(energyEstimator))
        hacc["hAcc%d%d"%(delayed_nxcut,dTcut)].SetZTitle('accidental coincidence rate (Hz)')

        hibd["hIBD%d%d"%(delayed_nxcut,dTcut)] = TH2D("hIBD%d%d"%(delayed_nxcut,dTcut),'IBD coincidence rate',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].SetXTitle('distance from wall [m]')
        hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].SetYTitle('prompt %s cut'%(energyEstimator))
        hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].SetZTitle('IBD rate (Hz)')

        hibdBG["hIBDBG%d%d"%(delayed_nxcut,dTcut)] = TH2D("hIBDBG%d%d"%(delayed_nxcut,dTcut),'IBD BG coincidence rate',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hibdBG["hIBDBG%d%d"%(delayed_nxcut,dTcut)].SetXTitle('distance from wall [m]')
        hibdBG["hIBDBG%d%d"%(delayed_nxcut,dTcut)].SetYTitle('prompt %s cut'%(energyEstimator))
        hibdBG["hIBDBG%d%d"%(delayed_nxcut,dTcut)].SetZTitle('IBD backgrounds rate (Hz)')

        hrn["hRN%d%d"%(delayed_nxcut,dTcut)] = TH2D("hRN%d%d"%(delayed_nxcut,dTcut),'Radionuclide coincidence rate ',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hrn["hRN%d%d"%(delayed_nxcut,dTcut)].SetXTitle('distance from wall [m]')
        hrn["hRN%d%d"%(delayed_nxcut,dTcut)].SetYTitle('prompt %s cut'%(energyEstimator))
        hrn["hRN%d%d"%(delayed_nxcut,dTcut)].SetZTitle('Radionuclide rate (Hz)')

        hfn["hFN%d%d"%(delayed_nxcut,dTcut)] =  TH2D("hFN%d%d"%(delayed_nxcut,dTcut),'Fast neutron coincidence rate ',binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hfn["hFN%d%d"%(delayed_nxcut,dTcut)].SetXTitle('distance from wall [m]')
        hfn["hFN%d%d"%(delayed_nxcut,dTcut)].SetYTitle('prompt %s cut'%(energyEstimator))
        hfn["hFN%d%d"%(delayed_nxcut,dTcut)].SetZTitle('Fast neutron rate (Hz)')
      
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:
                    _tag = 'hist_%s_%s_%s_%sdelayed%d_%dus'%(_loc,_element,_p,energyEstimator,delayed_nxcut,dTcut)
                    _tag = _tag.replace(" ","")
                    _str = "%s_%s_%s"%(_element,_loc,_p)
                    _str = _str.replace(" ","")
                    print('Extracting coincidences from file ',_tag,'and scaling to daily coincidence rate\n')
                    if 'pn_ibd' in _tag or 'A_Z' in _tag or 'FASTNEUTRONS' in _tag:
                        print('correlated event, getting %s from %s'%(_tag,pairstr))
                        hist[_tag] = pairEffFile.Get(_tag)
                    elif 'SINGLES' in _tag:
                        print('uncorrelated event, getting %s from %s'%(_tag,pairstr))
                        hist[_tag] = pairEffFile.Get(_tag)
                    else:
                        continue

                    try:
                        print(' entries:',hist[_tag].GetEntries(),' ... ', end = '')
                        hist[_tag].Scale(_scale)
                        print(_tag)
                        if 'SINGLES' in _tag: 
                            print('identified as accidentals')
                            hacc["hAcc%d%d"%(delayed_nxcut,dTcut)].Add(hist[_tag],1)
                        if 'pn_ibd' in _tag:
                            print('%s identified as IBD pair events ,'%(_tag))
                            if arguments["--Heysham"]:
                                if 'heysham_signal' in _tag:
                                    hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].Add(hist[_tag],1)
                                elif 'heysham_background' in _tag:
                                    hibdBG["hIBDBG%d%d"%(delayed_nxcut,dTcut)].Add(hist[_tag],1)
                            else:
                                if 'big_' in _tag or 'small' in _tag:
                                    hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].Add(hist[_tag],1)
                                elif 'boulby' in _tag:
                                    hibdBG["hIBDBG%d%d"%(delayed_nxcut,dTcut)].Add(hist[_tag],1)
                        if 'FASTNEUTRONS' in _tag:
                            hfn["hFN%d%d"%(delayed_nxcut,dTcut)].Add(hist[_tag],1)
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
                            hrn["hRN%d%d"%(delayed_nxcut,dTcut)].Add(hist[_tag],1)
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

    for delayed_nxcut,dTcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT)):
        
        binR = hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].GetNbinsX()
        binN = hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].GetNbinsY()

        _histograms["sOverB%d%d"%(delayed_nxcut,dTcut)]= hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].Clone()
        _histograms["sOverB%d%d"%(delayed_nxcut,dTcut)].SetZTitle('signal/sqrt(signal+background)')
        _histograms["sOverB%d%d"%(delayed_nxcut,dTcut)].SetYTitle('%s cut on prompt'%(energyEstimator))
        _histograms["sOverB%d%d"%(delayed_nxcut,dTcut)].SetTitle('hSoB - %s delayed %2d %dus'%(energyEstimator,delayed_nxcut,dTcut))
        _histograms["sOverB%d%d"%(delayed_nxcut,dTcut)].SetName('hSoB%d%d'%(delayed_nxcut,dTcut))
        _histograms["sOverB%d%d"%(delayed_nxcut,dTcut)].Reset()

        _histograms["signal%d%d"%(delayed_nxcut,dTcut)]= hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].Clone()
        _histograms["signal%d%d"%(delayed_nxcut,dTcut)].SetZTitle('evts/day')
        _histograms["signal%d%d"%(delayed_nxcut,dTcut)].SetYTitle('%s cut on prompt'%(energyEstimator))
        _histograms["signal%d%d"%(delayed_nxcut,dTcut)].SetTitle('Signal - %s delayed %2d %dus'%(energyEstimator,delayed_nxcut,dTcut))
        _histograms["signal%d%d"%(delayed_nxcut,dTcut)].SetName('hSignal%d%d'%(delayed_nxcut,dTcut))
        _histograms["signal%d%d"%(delayed_nxcut,dTcut)].Reset()
        _histograms["background%d%d"%(delayed_nxcut,dTcut)]= hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].Clone()
        _histograms["background%d%d"%(delayed_nxcut,dTcut)].SetZTitle('evts/day')
        _histograms["background%d%d"%(delayed_nxcut,dTcut)].SetYTitle('%s cut on prompt'%(energyEstimator))
        _histograms["background%d%d"%(delayed_nxcut,dTcut)].SetTitle('total backgrounds - %s delayed %2d %dus'%(energyEstimator,delayed_nxcut,dTcut))
        _histograms["background%d%d"%(delayed_nxcut,dTcut)].SetName('hBackground%d%d'%(delayed_nxcut,dTcut))
        _histograms["background%d%d"%(delayed_nxcut,dTcut)].Reset()




        maxSoB,maxSignal,maxBackground,maxNXd,maxNXp,maxFidcut,maxdT,maxdR = -1,-1,-1,-1,-1,-1,-1,-1

        for fidcut in drange(minFid,maxFid+binwidthFid,binwidthFid):
            for prompt_nxcut in drange(minNXprompt,maxNXprompt+binwidthNX,binwidthNX):

                # get the signal rates
                cutBin   = hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].FindBin(fidcut,prompt_nxcut)
                ibdRate  = hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].GetBinContent(cutBin)
                signal   = ibdRate
                ibdError = hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].GetBinError(cutBin)
                signalError = ibdError

                # now get the background rates
                # first get the individual background contributions
                # at each cut value
                accRate  = hacc["hAcc%d%d"%(delayed_nxcut,dTcut)].GetBinContent(cutBin)
                accError = hacc["hAcc%d%d"%(delayed_nxcut,dTcut)].GetBinError(cutBin)
                fnRate   = hfn["hFN%d%d"%(delayed_nxcut,dTcut)].GetBinContent(cutBin)
                fnError = hfn["hFN%d%d"%(delayed_nxcut,dTcut)].GetBinError(cutBin)
                rnRate   = hrn["hRN%d%d"%(delayed_nxcut,dTcut)].GetBinContent(cutBin)
                rnError = hrn["hRN%d%d"%(delayed_nxcut,dTcut)].GetBinError(cutBin)
                ibdBgRate   = hibdBG["hIBDBG%d%d"%(delayed_nxcut,dTcut)].GetBinContent(cutBin)
                ibdBgError = hibdBG["hIBDBG%d%d"%(delayed_nxcut,dTcut)].GetBinError(cutBin)

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
                _histograms["sOverB%d%d"%(delayed_nxcut,dTcut)].SetBinContent(cutBin,sob)
                _histograms["signal%d%d"%(delayed_nxcut,dTcut)].SetBinContent(cutBin,signal)
                _histograms["background%d%d"%(delayed_nxcut,dTcut)].SetBinContent(cutBin,background)
                _histograms["sOverB%d%d"%(delayed_nxcut,dTcut)].SetBinError(cutBin,sobError)
                _histograms["signal%d%d"%(delayed_nxcut,dTcut)].SetBinError(cutBin,signalError)
                _histograms["background%d%d"%(delayed_nxcut,dTcut)].SetBinError(cutBin,bgError)

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
    result = "%s %4.1f %3d %3d %d %4.3f +/-%4.5f %4.3f +/-%4.3f %4.7f %4.3f %4.7f %4.7f  %4.3f +/-%4.5f %4.3f +/-%4.5f %4.1f " %(additionalString,optFidcut,optNXprompt,optNXdelayed,optdT,optSignal,optSigErr,optBg,optBgErr,optAcc,optRN,optFN,optIBDbg,optSoB,optSoBErr,optSoB30days,optSoB30daysErr,T3SIGMA)
    resultsfile = "results_coincidence%s.txt"%(additionalString)
    with open(resultsfile,'a') as file:
        file.write(result+'\n')
    # write the signal, background and s/sqrt(s+b) histograms to the analysis_results_coincidence* file
    analysisfile = "core_root_files%s/sensitivity_results.root"%(additionalString)
    print('\n\nWriting histograms to file',analysisfile)
    f_root = TFile(analysisfile,"recreate")
    for delayed_nxcut in drange(minNXdelayed,maxNXdelayed+binwidthNX,binwidthNX):
        _histograms["sOverB%d%d"%(delayed_nxcut,dTcut)].Write()
        _histograms["signal%d%d"%(delayed_nxcut,dTcut)].Write()
        _histograms["background%d%d"%(delayed_nxcut,dTcut)].Write()
        hacc["hAcc%d%d"%(delayed_nxcut,dTcut)].Write()
        hfn["hFN%d%d"%(delayed_nxcut,dTcut)].Write()
        hrn["hRN%d%d"%(delayed_nxcut,dTcut)].Write()
        hibd["hIBD%d%d"%(delayed_nxcut,dTcut)].Write()
        hibdBG["hIBDBG%d%d"%(delayed_nxcut,dTcut)].Write()

    f_root.Close()

    print('\n\n')
