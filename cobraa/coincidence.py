from ROOT import gDirectory
from ROOT import TFile,TH2D
import numpy as np
import time
from itertools import product,combinations
from numpy import multiply
from math import fabs,log
import sys

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
   
    if arguments['--core']:
        _str = "core_root_files%s/coincidence_results.root"%(additionalString)
    else:
        _str = "fred_root_files%s/coincidence_results.root"%(additionalString)
#    if os.path.exists(_str):
#        proceed = intput("File %s already exists. Histograms will be added to the file. Do you want to proceed? Enter y or n:",_str)
#        if proceed ==n:
#            sys.exit("Exiting. Save coincidences results file before continuing.")
#        else:
#            print("Updating the coincidences results file.")

    outfile = TFile(_str,"UPDATE")

    if arguments['--evtype']:
        evtype = arguments['--evtype']
        if evtype=="li9":
            evtype="li 9"
        if evtype =="n17":
            evtype= "n 17"
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p][_loc]:
                    if _element==evtype:
                        _tag = "%s_%s_%s"%(_element,_loc,_p)
                        _tag = _tag.replace(" ","")
                        if arguments['--core']:
                            _file = "core_root_files%s/merged_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
                        else:
                            _file = "fred_root_files%s/merged_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
                        _file = _file.replace(" ","")
                        rate = rates[_tag][0]
                        if 'pn_ibd' in _tag or 'A_Z' in _tag or 'FAST' in _tag or 'singles' in _tag:
                            print(_tag," from ",_file)
                            obtainCoincidences(_file,_tag,outfile,rate)
                            print('')
    else:
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p][_loc]:
                    _tag = "%s_%s_%s"%(_element,_loc,_p)
                    _tag = _tag.replace(" ","")
                    if arguments['--core']:
                        _file = "core_root_files%s/merged_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
                    else:
                        _file = "fred_root_files%s/merged_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
                    _file = _file.replace(" ","")
                    if 'pn_ibd' in _tag or 'A_Z' in _tag or 'FAST' in _tag or 'singles' in _tag:
                        rate = rates[_tag][0]
                        print(_tag," from ",_file)
                        obtainCoincidences(_file,_tag,outfile,rate)
                        print('')

    print('Saving outfile:',_str)
    outfile.Close()
    return 0

def obtainCoincidences(file,_tag,outfile,rate):

    # Performs the coincidence evaluation on the 
    # FRED root files
    # called by coincidenceMap()
    start_time = time.time()
    recofile = TFile(file)
    print('Reading', file)
    hist = {}

    # check the root file is valid    
    try:
        runSummary = recofile.Get('runSummary')
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
    data         = recofile.Get('data')
    dataEntries  = data.GetEntries()
    data.SetBranchStatus('*',0)
    data.SetBranchStatus('%s'%(energyEstimator),1)
    data.SetBranchStatus('%s_prev'%(energyEstimator),1)
    data.SetBranchStatus('closestPMT',1)
    data.SetBranchStatus('inner_hit',1)
    data.SetBranchStatus('inner_hit_prev',1)
    data.SetBranchStatus('veto_hit',1)
    data.SetBranchStatus('veto_hit_prev',1)
    data.SetBranchStatus('good_pos',1)
    data.SetBranchStatus('good_pos_prev',1)
    data.SetBranchStatus('dt_prev_us',1)
    data.SetBranchStatus('timestamp',1)
    data.SetBranchStatus('gtid',1)
    if not arguments['--core']:
        data.SetBranchStatus('closestPMT_prev',1)

    # now we can evaluate the event coincidence
    # and scale down to the day rate
    for delayed_nxcut,dTcut,maxEp,gcut in product(drange(minNXdelayed,rangeNXdmax,binwidthNX),drange(dTmin,rangedTmax,binwidthdT),drange(minEpmax,rangeEpmax,binwidthEpmax),drange(gmin,rangeGmax,binwidthG)):
        tag = _tag+'_delayed%s_%d_%dus_maxEp%d_%d'%(energyEstimator,delayed_nxcut,dTcut,maxEp,gcut*10)
        histname = "hist_%s;1"%(tag)
        gDirectory.Delete("%s"%(histname));
        
        hist[tag] = TH2D('hist_%s'%(tag),'Coincidences -  %s '%(tag),binFid,rangeFidmin,rangeFidmax,binNX,rangeNXpmin,rangeNXpmax)
        hist[tag].SetXTitle('distance from wall [m]')
        hist[tag].SetYTitle('prompt %s cut'%(energyEstimator))
        hist[tag].SetZTitle('coincidences per day')
        hist[tag].Reset()
        
        for fidcut,prompt_nxcut in product(drange(minFid,rangeFidmax,binwidthFid),drange(minNXprompt,rangeNXpmax,binwidthNX)):
            if arguments['--positiveScan']:
                if prompt_nxcut>delayed_nxcut:
                    continue
            elif arguments['--negativeScan']:
                if prompt_nxcut<delayed_nxcut:
                    continue

            # define the prompt/delayed/coincidence cuts
            coincidences=0
            delayedtrigger  = "closestPMT/1000.>%f"%(fidcut)
            delayedtrigger  += "&& good_pos>%f " %(gcut)
            delayedtrigger  += "&& inner_hit > 10 &&  veto_hit < 4"
            delayedtrigger  += "&& %s > %f" %(energyEstimator,delayed_nxcut) 

            coincidencetrigger =  delayedtrigger
            coincidencetrigger += "&& good_pos_prev>%f"%(gcut)
            coincidencetrigger += "&& inner_hit_prev > 10 && veto_hit_prev <4"
            coincidencetrigger += "&& %s_prev > %f"%(energyEstimator,prompt_nxcut)
            coincidencetrigger += "&& dt_prev_us > 1 && dt_prev_us < %f"%(dTcut) 
            coincidencetrigger += "&& %s_prev < %f"%(energyEstimator,maxEp)
            if not arguments['--core']:
                # cut on the distance of the prompt event from PMTs if not using 
                # combined reconstruction
                coincidencetrigger += "&& closestPMT_prev/1000.>%f"%(fidcut)

            # find all coincidences
            coincidences_tmp = data.Draw("timestamp:dt_prev_us",coincidencetrigger,"goff")
            t_delayed = data.GetV1() # time of delayed-interaction trigger
            t_delayed = np.ndarray((coincidences_tmp),'d',t_delayed)
            dt_prev_us = data.GetV2() # time betweeen interactions in a pair
            dt_prev_us = np.ndarray((coincidences_tmp),'d',dt_prev_us)
            t_prompt = t_delayed  - dt_prev_us # time of prompt-interaction trigger
            dtnext = t_prompt[2:]-t_delayed[1:-1] # time to the first interaction in the next event pair
            dtlast = t_prompt[1:-1]-t_delayed[:-2] # time since the last interaction in the previous event pair
            # do the multiplicity cut (for fast neutron multiplicity)
            # check for coincidences immediately before and after pair
            coincidences = np.count_nonzero((dtnext>dTcut) & (dtlast>dTcut))
            # Poisson 95% UCL for zero if no coincidences are found
            if coincidences == 0:
                coincidences = 3.6889/totalEvents 
            # alternative since this is really a binomial distribution:
            # poisson approximation to binomial: -log(0.05)/totalEvents 
            # gives a lower value for the UCL so adopting the most cautious approach here
            # ref https://www.hilarispublisher.com/open-access/on-finding-the-upper-confidence-limit-for-a-binomial-proportion-when-zero-successes-are-observed-2155-6180-1000338.pdf
            coincidenceErr = 1/float(totalEvents)*sqrt(coincidences*(1-coincidences/float(totalEvents)))
            hist[tag].Fill(fidcut,prompt_nxcut,float(coincidences))
            errorbin = hist[tag].FindBin(fidcut,prompt_nxcut)
            hist[tag].SetBinError(errorbin, coincidenceErr)
            
            # end loop over prompt nx cuts
            # end loop over fiducial cuts
        # end loop over delayed nx cuts and dT time between triggers

        # scale no. of coincidences to day rate
        ndays = totalEvents/rate/86400.
        if 'singles' in _tag:
            ndays = totalEvents/float(singlespersec)/86400.
        hist[tag].Scale(1./float(ndays))
        
        outfile.cd()
        hist[tag].Write()
    
    # end loop over delayed nx, dT and dR cuts
    recofile.Close()
    del data
    print("--- %s seconds ---" % (time.time() - start_time))
    return
