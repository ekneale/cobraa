from ROOT import TFile,TH2D
import numpy as np
#import time
from itertools import product
from numpy import multiply

from .load import *
from .globals import *

# This creates maps of signal and background single-event efficiencies in 3 dimensions. 
# It outputs 2D histograms of coincidence rates as a function of prompt-event
# energy threshold and fiducial cut for each value of delayed-event energy
# threshold.
# Author Marc Bergevin

def efficiencyMapLassen():
    #calls obtainEventEfficiency() to create and fill histograms of efficiency as a function of nx and distance to nearest PMT

    print('Generating histograms of reconstruction efficiency for all background and signal processes')

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"

    for idx,_cover in enumerate(coverage):
        if (int(arguments["--QFIT"])==1):
            _strAdd = _cover+"_QFIT"
        else:
            _strAdd = _cover
        print("Saving with attachment",_strAdd)
        _str = "bonsai_root_files%s/results_%s.root"%(additionalMacStr,_strAdd)

        #outfile = TFile(_str,"UPDATE")

        ## This section needs to be revisited due to the partial obtion        

        if arguments["--procEle"] and arguments["--procType"]:
            _str = "bonsai_root_files%s/partial_%s_%s_%s_results_%s.root"%(additionalMacStr,arguments["--procType"],arguments["--procLoc"],arguments["--procEle"],_strAdd)
            outfile = TFile(_str,"RECREATE")
        elif  arguments["--procType"]:
            _str = "bonsai_root_files%s/partial_%s_%s_results_%s.root"%(additionalMacStr,arguments["--procType"],arguments["--procLoc"],_strAdd)
            outfile = TFile(_str,"RECREATE")
        else:
            outfile = TFile(_str,"RECREATE")

        if  arguments["--procType"]:
            proc = {}
            proc[arguments["--procType"]] = [arguments["--procLoc"]]
            print('Isolated proc:',proc)

        offset_binning = int(arguments["--offBin"])
        maxxoff = int(arguments["--scanMaxValue"])
        if arguments["--positiveScan"]:
             maxxoff_min = 0
             maxxoff_max = maxxoff+1
        elif arguments["--negativeScan"]:
             maxxoff_min = -maxxoff
             maxxoff_max = 0
        else:
             maxxoff_min = -maxxoff
             maxxoff_max = maxxoff+1     

        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p][_loc]:
                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _tag = _tag.replace(" ","")
                    _file = "bonsai_root_files%s/merged_%s_%s_%s_%s.root"%(additionalMacStr,_cover,_element,_loc,_p)
                    _file = _file.replace(" ","")
                    print(_tag," from ",_file)
                    if 'fast_neutron' in _tag or 'A_Z' in _tag or '_pn_' in _tag:
                        for _i in range(maxxoff_min,maxxoff_max,offset_binning):
                            print("Performing coincidence evaluation for offset: ",_i)
                            obtainDualEfficiency(_cover,_file,_tag,outfile,offset=_i)
                    else:
                        obtainEventEfficiency(_cover,_file,_tag,outfile)

                    print('')
    print('Saving outfile:',_str)
    outfile.Close()
    return 0


def obtainEventEfficiency(cover,file,_tag,outfile,_distance2pmt=1,_nx=10,_dist=30.0,\
_posGood=0.1,_dirGood=0.1,_pe=10,_nhit=10,_itr = 1.5):
    # covPCT  = coveragePCT[cover]

    _energyEstimator = arguments['--energyEst'] 
    para = testEnabledCondition(arguments)
    additionalString  = para[0]
    arbre = {}
    arbre["rfile"] = TFile(file)
    print('Reading', file)
    try:
        runSummary = arbre["rfile"].Get('runSummary')
        Entries = runSummary.GetEntries()
        runSummary.GetEntry(Entries-1)
        events = 0
        _eventPerRun = runSummary.nEvents
    except:
        print('File',file,'did not have run associated with it. Returning empty histogram.')
        binR,rangeRmin,rangeRmax = 31,0.45,3.55
        binwidthR = (rangeRmax-rangeRmin)/binR
        binN,rangeNmin,rangeNmax = 48,7.5,55.5
        binwidthN = (rangeNmax-rangeNmin)/binN
        h = TH2D('hist%s'%(_tag),'EMPTY - Rate of events -  %s '%(_tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        h.SetXTitle('distance from wall [m]')
        h.SetYTitle('%s cut'%(_energyEstimator))
        h.SetZTitle('efficiency')
        h.GetZaxis().SetTitleOffset(-.55);
        h.GetZaxis().SetTitleColor(1);
        h.GetZaxis().CenterTitle();
        #h.SaveAs("bonsai_root_files%s/%s/hist%s.C"%(additionalString,cover,_tag))
        return -1

    for i in range(10):
        events+= runSummary.subEventTally[i]
    totalEvents = float(Entries)*_eventPerRun

    arbre["data"]   = arbre["rfile"].Get('data')
    _someEntries = arbre["data"].GetEntries()

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 88,7.5,95.5
    binwidthN = (rangeNmax-rangeNmin)/binN
    

    h = TH2D('hist%s'%(_tag),'Rate of events -  %s '%(_tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
    #    h.SetName(_tag)
    h.SetXTitle('distance from wall [m]')
    h.SetYTitle('%s cut'%(_energyEstimator))
    h.SetZTitle('efficiency')
    h.GetZaxis().SetTitleOffset(-.55);
    h.GetZaxis().SetTitleColor(1);
    h.GetZaxis().CenterTitle();
    for _d in drange(rangeRmin+binwidthR/2.,rangeRmax,binwidthR):
        minAchieve = 0

        for _nx in range(int(rangeNmin+binwidthN/2.0),int(rangeNmax)):
            cond = "closestPMT/1000.>%f"%(_d)
            cond += "&& good_pos>%f " %(_posGood)
            cond += "&& inner_hit > 4 &&  veto_hit < 4"
            cond += "&& %s > %f" %(_energyEstimator,_nx)
            if int(arguments["--QFIT"])==1:
                cond = "closestPMTQFit/1000.>%f"%(_d)
                cond += "&& QFit ==1 "
                cond += "&& inner_hit > %f &&  veto_hit < 4" %(_nx)

            if _someEntries !=0:
                if minAchieve == 0:

                    evts = arbre["data"].Draw("",cond,"goff")
                    eff = evts/totalEvents

                    if eff == 0 and not arguments['--Heysham']:
                        eff = 1.0/totalEvents
                        minAchieve = 1
                    h.Fill(_d,_nx,eff)
                    effErr = 1./totalEvents*sqrt(evts*(1-evts/totalEvents))
                    errorbin = h.FindBin(_d,_nx)
                    h.SetBinError(errorbin, effErr)
                else:
                    h.Fill(_d,_nx,1./totalEvents)
                    eff = 1./totalEvents
            else:
                h.Fill(_d,_nx,1./totalEvents)
                eff =  1./totalEvents
    outfile.cd()
    ##Changes for multiple file access
    h.Write()
    #    h.SaveAs("bonsai_root_files%s/%s/hist%s.C"%(additionalString,cover,_tag))
    arbre["rfile"].Close()
    del arbre
    return eff


def obtainDualEfficiency(cover,file,_tag,outfile,_distance2pmt=1,_nx=10,_dist=30.0,\
		_posGood=0.1,_dirGood=0.1,_pe=10,_nhit=10,_itr = 1.5, offset=0,rate=1):

    _energyEstimator = arguments['--energyEst']
    para = testEnabledCondition(arguments)
    additionalString  = para[0]
    arbre = {}
    arbre["rfile"] = TFile(file)
    print('Reading', file)
    _offset = "%s"%(offset)
    _offset = _offset.replace("-","Minus")
    try:
        runSummary = arbre["rfile"].Get('runSummary')
        Entries = runSummary.GetEntries()
        runSummary.GetEntry(Entries-1)
        events = 0
        _eventPerRun = runSummary.nEvents
    except:
        print('File',file,'did not have run associated with it. Returning empty histogram.')
        binR,rangeRmin,rangeRmax = 31,0.45,3.55
        binwidthR = (rangeRmax-rangeRmin)/binR
        binN,rangeNmin,rangeNmax = 48,7.5,55.5
        binwidthN = (rangeNmax-rangeNmin)/binN
        h = TH2D('hist%s_%s'%(_tag,_offset),'EMPTY - Rate of events -  %s '%(_tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        h.SetXTitle('distance from wall [m]')
        h.SetYTitle('%s cut'%(_energyEstimator))
        h.SetZTitle('efficiency')
        h.GetZaxis().SetTitleOffset(-.55);
        h.GetZaxis().SetTitleColor(1);
        h.GetZaxis().CenterTitle();
        #h.SaveAs("bonsai_root_files%s/%s/hist%s.C"%(additionalString,cover,_tag))
        return -1

    for i in range(10):
        events+= runSummary.subEventTally[i]
    totalEvents = float(Entries)*_eventPerRun

    arbre["data"]   = arbre["rfile"].Get('data')
    _someEntries = arbre["data"].GetEntries()

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 88,7.5,95.5
    binwidthN = (rangeNmax-rangeNmin)/binN

    h = TH2D('hist%s_%s'%(_tag,_offset),'Rate of events -  %s_%s '%(_tag,_offset),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
    #    h.SetName(_tag)
    h.SetXTitle('distance from wall [m]')
    h.SetYTitle('%s cut'%(_energyEstimator))
    h.SetZTitle('efficiency')
    h.GetZaxis().SetTitleOffset(-.55);
    h.GetZaxis().SetTitleColor(1);
    h.GetZaxis().CenterTitle();
    for _d in drange(rangeRmin+binwidthR/2.,rangeRmax,binwidthR):
        minAchieve = 0

        for _nx in range(int(rangeNmin+binwidthN/2.0),int(rangeNmax)):
            cond = "closestPMT/1000.>%f && closestPMT_prev/1000.>%f "%(_d,_d)
            cond += "&& good_pos>%f " %(_posGood)
            cond += "&& inner_hit > 4 &&  veto_hit < 4"
            cond += "&& %s > %f && %s > %f" %(_energyEstimator,_nx,_energyEstimator,_nx+offset)
            cond += "&&  dt_prev_us < %s && drPrevr < %f" %(arguments["-t"],1000.*float(arguments["-t"]))
            if int(arguments["--QFIT"])==1:
                cond = "closestPMTQFit/1000.>%f"%(_d)
                cond += "&& QFit ==1 "
                cond += "&& inner_hit > %f &&  veto_hit < 4" %(_nx)

            if _someEntries !=0:
                if minAchieve == 0:
                    evts = arbre["data"].Draw("mcx",cond,"goff")
                    _multEvts = 0
                    _priorMCX = 0.
                    ### Simplified way of reject >double neutrons are to reject multiple events
                    _mcx = arbre["data"].GetV1()
                    _flag = 0
                    for __mcx in range(evts):
                        if _mcx[__mcx] == _priorMCX:
                            ###print('Found duplicate', _mcx[__mcx], _priorMCX, _multEvts)
                            if _flag == 1:
                                ## Must only remove the new pair from the total tally
                                _multEvts+=1
                            elif _flag == 0:
                                ## First event in burst must remove first and second pair
                                _multEvts+=2
                                _flag = 1
                            #print('Found duplicate', _mcx[__mcx], _priorMCX, _multEvts)
                        else:
                            #print('Not   duplicate', _mcx[__mcx], _priorMCX, _multEvts)
                            _flag = 0
                        _priorMCX = _mcx[__mcx]
                    eff = (evts-_multEvts)/totalEvents
                    #print('Breakdown. Offset:',offset,'Events (tot,mult,subtracted..)',evts,_multEvts,evts-_multEvts,totalEvents,(evts-_multEvts)/totalEvents)
                    if eff == 0 and not arguments['--Heysham']:
                        eff = 1.0/totalEvents
                        minAchieve = 1
                    h.Fill(_d,_nx,eff)
                    effErr = 1./totalEvents*sqrt(evts*(1-evts/totalEvents))
                    errorbin = h.FindBin(_d,_nx)
                    h.SetBinError(errorbin, effErr)
                else:
                    h.Fill(_d,_nx,1./totalEvents)
                    eff = 1./totalEvents
            else:
                h.Fill(_d,_nx,1./totalEvents)
                eff =  1./totalEvents
    #hist = h.Clone()
    outfile.cd()
    #Changes for multiple file access
    h.Write()
    #h.SaveAs("bonsai_root_files%s/hist_%s_%s_%d_.C"%(additionalString,cover,_tag,offset))
    arbre["rfile"].Close()
    del arbre
    return h


def EfficiencyMapInPMTVol():
    '''For the given configuration on initializing watchmakers,
    Generates PMT Volume Efficiency histograms for all merged
    files.'''
    nx = int(arguments['--minN9'])
    good_pos = float(arguments['-g'])
    good_dir = float(arguments['-G'])
    print("Evaluating sensitivity in PMT Volume for all WaterVolume types. Using given minimum parameter requirements")
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"
    #d,proc,coverage = loadSimulationParametersNew()
    for _p in proc:
        for _loc in proc[_p]:
            if _loc != "WaterVolume":
                continue
            for idx,_cover in enumerate(coverage):
                for _element in d[_p]:
                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _file = "bonsai_root_files%s/%s/merged_%s_%s_%s.root"%(additionalMacStr,_cover,_loc,_element,_p)
                    print(_tag)
                    obtainEfficiencyInPMTVol(_cover,_file,_tag, _nx=nx, _posGood=good_pos, _dirGood=good_dir)
                    print('')


def obtainEfficiencyInPMTVol(cover,file,_tag,_nx=8,\
_posGood=0.1,_dirGood=0.1):
    '''For the given merged bonsai file, will generate a histogram giving the
    efficiency inside of the defined fiducial volume'''
    # covPCT  = coveragePCT[cover]
    _energyEstimator = arguments['--energyEst']
    para = testEnabledCondition(arguments)
    additionalString  = para[0]
    arbre = {}
    arbre["rfile"] = TFile(file)
    print('Reading', file)
    try:
        runSummary = arbre["rfile"].Get('runSummary')
        Entries = runSummary.GetEntries()
        runSummary.GetEntry(Entries-1)
        events = 0
        _eventPerRun = runSummary.nEvents
    except:
        print('File',file,'did not have run associated with it. Returning empty histogram.')
        binR,rangeRmin,rangeRmax = 31,0.0, 1.0
        binwidthR = (rangeRmax-rangeRmin)/binR
        binN,rangeNmin,rangeNmax = 48,-1.0*pmtHeight,pmtHeight
        binwidthN = (rangeNmax-rangeNmin)/binN
        h = TH2D('hist%s'%(_tag),'EMPTY - Acceptance in PMT Volume -  %s '%(_tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        h.SetXTitle(r'($\rho$/$\rho_{tank}$)$^{2}$')
        h.SetYTitle('Z (m)')
        h.SetZTitle('Acceptance Fraction')
        h.GetZaxis().SetTitleOffset(-.55);
        h.GetZaxis().SetTitleColor(1);
        h.GetZaxis().CenterTitle();
        h.SaveAs("bonsai_root_files%s/%s/PMTVolEff%s.C"%(additionalString,cover,_tag))
        return -1

    for i in range(10):
        events+= runSummary.subEventTally[i]
    totalEvents = float(Entries)*_eventPerRun

    arbre["data"]   = arbre["rfile"].Get('data')
    _someEntries = arbre["data"].GetEntries()

    binR,rangeRmin,rangeRmax = 10,0.0,1.0
    binN,rangeNmin,rangeNmax = 10,-1.*pmtHeight,pmtHeight
    h = TH2D('delEff%s'%(_tag),'Acceptance in PMT Volume -  %s '%(_tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
    h.SetXTitle(r'($\rho_{true}$/$\rho_{PMT}$)$^{2}$')
    h.SetYTitle('True Z (m)')
    h.SetZTitle('Acceptance Fraction')
    h.GetZaxis().SetTitleOffset(-.55)
    h.GetZaxis().SetTitleColor(1)
    h.GetZaxis().CenterTitle()

    #Now, we need to go entry by entry and fill our total histogram (denom)
    #And also fill our numerator if events pass the condition
    mcr = "sqrt(mcx**2 + mcy**2)"
    r = "sqrt(x**2 + y**2)"
    rho2eqn ="(%s/%f)**2"%(mcr,pmtRadius)
    thedraw = "mcz:%s"%(rho2eqn)

    MCPMTVolCond = "mcz>%f"%(rangeNmin)
    MCPMTVolCond += "&& mcz<%f"%(rangeNmax)
    MCPMTVolCond += "&& %s<%f"%(mcr,pmtRadius)

    PMTVolCond = "z>%f"%(rangeNmin)
    PMTVolCond += "&& z<%f"%(rangeNmax)
    PMTVolCond += "&& %s<%f"%(r,pmtRadius)

    effcond = "good_pos>%f " %(_posGood)
    effcond += "&& good_dir>%f " %(_dirGood)
    effcond += "&& %s > %i" %(_energyEstimator,_nx)

    #First, draw that sweet, sweet total events in FV
    arbre["data"].Draw("%s>>h_effdenominator(%i,%f,%f,%i,%f,%f)"%(thedraw,\
          binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax),\
          MCPMTVolCond,"goff")
    hdenom = gDirectory.Get("h_effdenominator")
    arbre["data"].Draw("%s>>h_effnumerator(%i,%f,%f,%i,%f,%f)"%(thedraw,\
          binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax) ,"%s && %s && %s"%(MCPMTVolCond,PMTVolCond,effcond),"goff")
    hnum = gDirectory.Get("h_effnumerator")
    #Now, we fill the actual efficiency histogram with the division of the two
    h.Divide(hnum,hdenom,1.,1.,"b")
    h.SaveAs("bonsai_root_files%s/%s/PMTVolEff%s_%s_%i_goodpos_%f_gooddir_%f.C"%(additionalString,\
            cover,_tag,_energyEstimator,_nx,_posGood,_dirGood))
    arbre["rfile"].Close()
    del arbre

