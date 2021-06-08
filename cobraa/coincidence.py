from ROOT import TFile,TH2D 
import numpy as np
import time
from itertools import product
from numpy import multiply

from .load import *
from .globals import *

# This creates maps of signal and background coincidences in 5 dimensions. 
# It outputs 2D histograms of coincidence rates as a function of prompt-event
# energy threshold and fiducial cut for each value of delayed-event energy
# threshold and distance and time between interactions (coincidence cuts).
# Author Liz Kneale (2021)
# Adapted from Watchmakers/analysis.py (Marc Bergevin)

def drange(start, stop, step):
    rii= start
    while rii<stop:
        yield rii
        rii+=step

def coincidenceMap():

    # calls obtainCoincidences() to create and fill histograms of coincidences
    # as a function of nx_p and distance to nearest PMT over nx_d and dT ranges
    evtype = ""
    if arguments['--evtype']:
        evtype = arguments['--evtype']
        print('Generating histograms of coincidences for %s'%(evtype))
    else:
        print('Generating histograms of coincidences for all signal and background')

    _str = "fred_root_files%s/coincidence_results.root"%(additionalString)
    outfile = TFile(_str,"UPDATE")

    if arguments['--evtype']:
        for _p in proc:
            if _p==arguments['--evtype']:
                for _loc in proc[_p]:
                    for _element in d[_p][_loc]:
                        _tag = "%s_%s_%s"%(_element,_loc,_p)
                        _tag = _tag.replace(" ","")
                        _file = "fred_root_files%s/merged_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
                        _file = _file.replace(" ","")
                        print(_tag," from ",_file)
                        rate = rates[_tag][0]
                        if 'singles' in _tag:
                            obtainAccidentalCoincidences(_file,_tag,outfile,rate)
                        elif 'pn_ibd' in _tag or 'A_Z' in _tag or 'FAST' in _tag:
                            obtainCorrelatedCoincidences(_file,_tag,outfile,rate)
                        print('')
    else:
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p][_loc]:
                    _tag = "%s_%s_%s"%(_element,_loc,_p)
                    _tag = _tag.replace(" ","")
                    _file = "fred_root_files%s/merged_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
                    _file = _file.replace(" ","")
                    print(_tag," from ",_file)
                    rate = rates[_tag][0]
                    if 'singles' in _tag:
                        obtainAccidentalCoincidences(_file,_tag,outfile,rate)
                    elif 'pn_ibd' in _tag or 'A_Z' in _tag or 'FAST' in _tag:
                        obtainCorrelatedCoincidences(_file,_tag,outfile,rate)

                    print('')
    print('Saving outfile:',_str)
    outfile.Close()
    return 0

def obtainCorrelatedCoincidences(file,_tag,outfile,rate):

    # Performs the coincidence evaluation on the 
    # FRED root files
    # called by coincidenceMap()

    start_time = time.time()
    fredfile = TFile(file)
    print('Reading', file)
    hist = {}

    # check the root file is valid    
    try:
        runSummary = fredfile.Get('runSummary')
        runEntries = runSummary.GetEntries()
        
    except:
        print('File',file,'did not have run associated with it. Returning empty histogram.')
        return -1

    totalEvents  = 0
    for i in range(runEntries):
        runSummary.GetEntry(i)
        totalEvents += runSummary.nEvents
    print(totalEvents)

    # get only required branches from data tree for speed
    data         = fredfile.Get('data')
    dataEntries  = data.GetEntries()
    data.SetBranchStatus('*',0)
    data.SetBranchStatus('%s'%(energyEstimator),1)
    data.SetBranchStatus('%s_prev'%(energyEstimator),1)
    data.SetBranchStatus('closestPMT',1)
    data.SetBranchStatus('closestPMT_prev',1)
    data.SetBranchStatus('inner_hit',1)
    data.SetBranchStatus('inner_hit_prev',1)
    data.SetBranchStatus('veto_hit',1)
    data.SetBranchStatus('veto_hit_prev',1)
    data.SetBranchStatus('good_pos',1)
    data.SetBranchStatus('good_pos_prev',1)
    data.SetBranchStatus('dt_prev_us',1)
    data.SetBranchStatus('drPrevr',1)
    data.SetBranchStatus('timestamp',1)

    # now we can evaluate the event coincidence
    # and scale down to the day rate
    for delayed_nxcut,dTcut,dRcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT),drange(dRmin,rangedRmax,binwidthdR)):
        tag = _tag+'_%sdelayed%d_%dus_%dmm'%(energyEstimator,delayed_nxcut,dTcut,dRcut*1000)
        hist[tag] = TH2D('hist_%s'%(tag),'Coincidences -  %s '%(tag),binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hist[tag].SetXTitle('distance from wall [m]')
        hist[tag].SetYTitle('prompt %s cut'%(energyEstimator))
        hist[tag].SetZTitle('coincidences per day')

        for fidcut,prompt_nxcut in product(drange(minFid,rangeFidmax,binwidthFid),drange(minNXprompt,rangeNXpmax,binwidthNX)):

            # find the coincidence efficiency

            coincidences=0
            delayedtrigger  = "closestPMT/1000.>%f"%(fidcut)
            delayedtrigger  += "&& good_pos>%f " %(posGood)
            delayedtrigger  += "&& inner_hit > 4 &&  veto_hit < 4"
            delayedtrigger  += "&& %s > %f" %(energyEstimator,delayed_nxcut) 

            coincidencetrigger =  delayedtrigger
            coincidencetrigger += "&& closestPMT_prev/1000.>%f"%(fidcut)
            coincidencetrigger += "&& good_pos_prev>%f"%(posGood)
            coincidencetrigger += "&& inner_hit_prev > 4 && veto_hit_prev <4"
            coincidencetrigger += "&& %s_prev > %f"%(energyEstimator,prompt_nxcut)
            coincidencetrigger += "&& dt_prev_us > 0 && dt_prev_us < %f"%(dTcut) 
            coincidencetrigger += "&& drPrevr/1000. < %f"%(dRcut)

            # find all coincidences
            coincidences_tmp = data.Draw("timestamp:dt_prev_us",coincidencetrigger,"goff")
            t_delayed = data.GetV1()
            t_delayed = np.ndarray((coincidences_tmp),'d',t_delayed)
            dt_prev_us = data.GetV2() # time since previous event
            dt_prev_us = np.ndarray((coincidences_tmp),'d',dt_prev_us)
            t_prompt = t_delayed  - dt_prev_us 

            # do the multiplicity cut (for fast neutron multiplicity)
            # time between coincidence and next event
            dtnext = t_prompt[2:]-t_delayed[1:-1]
            # time between coincidence and previous event
            dtlast = t_prompt[1:-1]-t_delayed[:-2]
            # find all coincidences with no other event
            # immediately before and after pair
            coincidences = np.count_nonzero((dtlast>dTcut) & (dtnext>dTcut))
            
            # calculate statistical error and fill histogram
            coincidenceErr = 1/totalEvents*sqrt(coincidences*(1-coincidences/totalEvents))
            hist[tag].Fill(fidcut,prompt_nxcut,coincidences)
            errorbin = hist[tag].FindBin(fidcut,prompt_nxcut)
            hist[tag].SetBinError(errorbin, coincidenceErr)
            # end loop over prompt nx cuts
            # end loop over fiducial cuts
        # end loop over delayed nx cuts and dT time between triggers
        # scale to day rate
        ndays = float(totalEvents/rate/86400.)
        hist[tag].Scale(1/ndays)
        outfile.cd()
        hist[tag].Write()

    fredfile.Close()
    del data
    print("--- %s seconds ---" % (time.time() - start_time))
    return coincidences
   

def obtainAccidentalCoincidences(file,_tag,outfile,rate):

    start_time = time.time()
    fredfile = TFile(file)
    print('Reading', file)
    hist = {}

    # check the root file is valid
    try:
        runSummary = fredfile.Get('runSummary')
        Entries = runSummary.GetEntries()
        runSummary.GetEntry(Entries-1)
        events = 0
        eventsPerRun = runSummary.nEvents
    except:
        print('File',file,'did not have run associated with it.')
        return -1

    totalEvents  = Entries*eventsPerRun
    data         = fredfile.Get('data')
    dataEntries  = data.GetEntries()
    data.SetBranchStatus('*',0)
    data.SetBranchStatus('%s'%(energyEstimator),1)
    data.SetBranchStatus('%s_prev'%(energyEstimator),1)
    data.SetBranchStatus('closestPMT',1)
    data.SetBranchStatus('inner_hit',1)
    data.SetBranchStatus('veto_hit',1)
    data.SetBranchStatus('good_pos',1)
    data.SetBranchStatus('timestamp',1)
    data.SetBranchStatus('x',1)
    data.SetBranchStatus('y',1)
    data.SetBranchStatus('z',1)

    # now we can get the efficiencies of coincident singles events
    for delayed_nxcut,dTcut,dRcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT),drange(dRmin,rangedRmax,binwidthdR)):
        
        tag = _tag+'_%sdelayed%d_%dus_%dmm'%(energyEstimator,delayed_nxcut,dTcut,dRcut*1000)
        hist[tag] = TH2D('hist_%s'%(tag),'Coincidences -  %s '%(tag),binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hist[tag].SetXTitle('distance from wall [m]')
        hist[tag].SetYTitle('prompt %s cut'%(energyEstimator))
        hist[tag].SetZTitle('efficiency')

        for fidcut,prompt_nxcut in product(drange(minFid,rangeFidmax,binwidthFid),drange(minNXprompt,rangeNXpmax,binwidthNX)):
            # find which is the smaller out of the prompt and delayed cuts
            # (important for the negative scan)
            min_nxcut = min(prompt_nxcut,delayed_nxcut)
            max_nxcut = max(prompt_nxcut,delayed_nxcut)
            # now find the events which pass the minimal cuts plus
            # fiducial and smallest of the prompt and delayed nx cuts
            coincidences=0
            mintrigger  = "closestPMT/1000.>%f"%(fidcut)
            mintrigger  += "&& good_pos>%f " %(posGood)
            mintrigger  += "&& inner_hit > 4 &&  veto_hit < 4"
            mintrigger += "&& %s > %f" %(energyEstimator,min_nxcut) 

            # Save time and nx of all of the events which pass the min trigger
            evts = data.Draw("timestamp:%s"%(energyEstimator),mintrigger,"goff")
            t = data.GetV1()
            t = np.ndarray((evts),'d',t)
            dt = np.diff(t)
            nx = data.GetV2()
            nx = np.ndarray((evts),'d',nx)
            # Now save x, y and z of all of the events which pass the prompt trigger
            evts = data.Draw("x/1000.:y/1000.:z/1000.",mintrigger,"goff")
            x = data.GetV1()
            y = data.GetV2()
            z = data.GetV3()
            x = np.ndarray((evts),'d',x)
            y = np.ndarray((evts),'d',y)
            z = np.ndarray((evts),'d',z)
            # Find the distance between consecutive events
            dx = np.diff(x)
            dy = np.diff(y)
            dz = np.diff(z)
            dR2 = sum([multiply(dx,dx),multiply(dy,dy),multiply(dz,dz)])
            dR = sqrt(dR2)
            # find all of the subevents which pass the higher of the two nx cuts
            # and have a preceding event within dT and dR (for positiveScan)
            # OR 
            # find all of the sub events which pass the lower of the two nx cuts
            # and have a preciding event within dT and dR and passes the higher
            # of the two nxcuts (for negative scan)
            if min_nxcut==prompt_nxcut:
                # move the nx values of the delayed event left
                # to correspond with the dt and dR values
                nx = nx.pop(0)
            coincidences = np.count_nonzero((dt>0) & (dt<dTcut) & (dR<dRcut) & (nx>max_nxcut))
            

            # calculate statistical error and fill histogram
            coincidenceErr = 1/totalEvents*sqrt(coincidences*(1-coincidences/totalEvents))
            hist[tag].Fill(fidcut,prompt_nxcut,coincidences)
            errorbin = hist[tag].FindBin(fidcut,prompt_nxcut)
            hist[tag].SetBinError(errorbin, coincidenceErr)

            # end loop over prompt nx cuts
            # end loop over fiducial cuts
        # end loop over delayed nx cuts and dT time & dR distance between triggers
        # scale to day rate
        ndays = float(totalEvents/singlespersec/86400.)
        hist[tag].Scale(1/ndays)

        outfile.cd()
        hist[tag].Write()

