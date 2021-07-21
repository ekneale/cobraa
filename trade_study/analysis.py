from .load import *
from .io_operations import testEnabledCondition
from ROOT import gDirectory
from ROOT import TFile,TH2D,TH3D,Double
import warnings
import glob
from math import fabs

d,proc,coverage,rates = loadSimulationParametersNew()
_energyEstimator = arguments['--energyEst']

fidRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])-float(arguments['--fidThick'])
fidHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])-float(arguments['--fidThick'])

pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])
pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])

detectorRadius  = float(arguments['--tankRadius'])-float(arguments['--steelThick'])
detectorHeight  = float(arguments['--halfHeight'])-float(arguments['--steelThick'])

if arguments["--letterboxSize"]:
    if arguments["--letterboxSize"]=='80':
        detectorRadius = 4000.
        detectorHeight = 40000.
    else: 
        detectorRadius = 4000.
        detectorHeight = 25000.

if arguments["--cylinderSize"]:
    if arguments["--cylinderSize"]=='12':
        detectorRadius = 6000.
        detectorHeight = 6000.
    else:
        detectorRadius = 5000.
        detectorHeight = 5000.

_posGood = float(arguments['-g'])
_dirGood = float(arguments['-G'])
minNX = float(arguments['--minNX'])
maxNX = float(arguments['--maxNX'])

warnings.simplefilter("ignore")

def drange(start, stop, step):
    rii= start
    while rii<stop:
        yield rii
        rii+=step

def logx_logy_array(nbins = 500,xmin = 1e-2,xmax = 30.,ymin = 1e-9,ymax = 1e3):
    #x-axis
    logxmin = log10(xmin)
    logxmax = log10(xmax)
    xbinwidth= (logxmax-logxmin)/float(nbins)
    xbins	= zeros(nbins,dtype=float)
    xbins[0] = xmin
    #y-axis
    logymin = log10(ymin)
    logymax = log10(ymax)
    ybinwidth= (logymax-logymin)/float(nbins)
    ybins	= zeros(nbins,dtype=float)
    ybins[0] = ymin
    for i in range(nbins):
        xbins[i] = xmin + pow(10,logxmin+i*xbinwidth)
        ybins[i] = ymin + pow(10,logymin+i*ybinwidth)
    return nbins,xbins,ybins


def obtainEventEfficiency(cover,file,_tag,outfile,_distance2pmt=1,_nx=8,_dist=30.0,\
_posGood=0.1,_dirGood=0.1,_pe=8,_nhit=8,_itr = 1.5):
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
        h.SaveAs("bonsai_root_files%s/%s/hist%s.C"%(additionalString,cover,_tag))
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
                        eff = 1/totalEvents
                        minAchieve = 1
                    h.Fill(_d,_nx,eff)
                    effErr = 1./totalEvents*sqrt(evts*(1-evts/totalEvents))
                    errorbin = h.FindBin(_d,_nx)
                    h.SetBinError(errorbin, effErr)
                else:
                    h.Fill(_d,_nx,1./totalEvents)
                    eff = 1./totalEvents
                    effErr = 1./totalEvents
                    h.SetBinError(errorbin,effErr)
            else:
                h.Fill(_d,_nx,1./totalEvents)
                eff =  1./totalEvents
                effErr = 1./totalEvents
                h.SetBinError(errorbin,effErr)
    outfile.cd()
    h.Write()
    #    h.SaveAs("bonsai_root_files%s/%s/hist%s.C"%(additionalString,cover,_tag))
    arbre["rfile"].Close()
    del arbre
    return eff

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


def obtainPairEfficiency(cover,file,_tag,outfile,dtcut,dscut):
    print('cut on goodness %s' %(_posGood))
    para = testEnabledCondition(arguments)
    additionalString  = para[0]
    arbre = {}
    arbre["rfile"] = TFile(file)
    print('Reading', file)
    hist = {}
    try:
        runSummary = arbre["rfile"].Get('runSummary')
        Entries = runSummary.GetEntries()
        runSummary.GetEntry(Entries-1)
        events = 0
        _eventPerRun = runSummary.nEvents
    except:
        print('File',file,'did not have run associated with it. Returning no histogram.')
        return -1

    totalEvents = Entries*_eventPerRun

    arbre["data"]    = arbre["rfile"].Get('data')

    _someEntries = arbre["data"].GetEntries()

    binwidthR = float(arguments["--binwidthR"])#(rangeRmax-rangeRmin)/binR
    rangeRmin,rangeRmax = 0.5-binwidthR/2.,3.0+binwidthR/2.
    binR = int((rangeRmax-rangeRmin)/binwidthR)
    binwidthN = float(arguments["--binwidthN"])#(rangeNmax-rangeNmin)/binN
    rangeNmin,rangeNmax = minNX-binwidthN/2.,maxNX+binwidthN/2.
    binN = int((rangeNmax-rangeNmin)/binwidthN)
    binwidthT = 10
    binwidthS = 0.1 
    scanStep = int(arguments['--scanStep'])
    scanMax = int(arguments["--scanMaxValue"])
    scanMin = int(arguments["--scanMinValue"])
    if  arguments["--positiveScan"]:
        offsets_nx = np.arange(scanMin,scanMax+1,scanStep,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        if scanMin==0:
            offsets_nx = np.arange(scanMax,scanMin+1,scanStep,dtype=int)## Correction to add last value
        else:
            offsets_nx = np.arange(-scanMax,-scanMin+1,scanStep,dtype=int)## Correction to add last value
    else:
        offsets_nx = np.arange(-scanMax,scanMax+1,scanStep,dtype=int)## Correction to add last value


    # now we can get the efficiencies of coincident correlated events
    for offset in offsets_nx:
        _offset = str(offset).replace("-","Minus")
        
        tag = _tag+'_%soffset%s'%(_energyEstimator,_offset)
        hist[tag] = TH2D('hist%s'%(tag),'Coincidence Efficiency -  %s '%(tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        hist[tag].SetXTitle('distance from wall [m]')
        hist[tag].SetYTitle('prompt %s cut'%(_energyEstimator))
        hist[tag].SetZTitle('efficiency')
        hist[tag].GetZaxis().SetTitleOffset(-.55);
        hist[tag].GetZaxis().SetTitleColor(1);
        hist[tag].GetZaxis().CenterTitle();

        for _d in drange(rangeRmin+binwidthR/2.,rangeRmax,binwidthR):
            minAchieve = 0

            for _nx in drange(int(rangeNmin+binwidthN/2.0),int(rangeNmax),binwidthN):
                coincidences=0
                delayedtrigger  = "closestPMT/1000.>%f"%(_d)
                delayedtrigger  += "&& good_pos>%f " %(_posGood)
                delayedtrigger  += "&& inner_hit > 4 &&  veto_hit < 4"
                delayedtrigger  += "&& %s > %f" %(_energyEstimator,_nx+offset) # all of the events pass the delayed nx cut
                coincidences=0
                eff = 0
                # find the coincidence efficiency
                if _someEntries !=0:
                    if minAchieve ==0:
                        coincidencetrigger =  delayedtrigger
                        coincidencetrigger += "&& good_pos_prev>%f"%(_posGood)
                        coincidencetrigger += "&& inner_hit_prev > 4 && veto_hit_prev <4"
                        coincidencetrigger += "&& %s_prev > %f"%(_energyEstimator,_nx) # preceding event passes the prompt nx cut
                        coincidencetrigger += "&& dt_prev_us > 0 && dt_prev_us < %f"%(dtcut) # time coincidence
                                 
                        coincidences_tmp = arbre["data"].Draw("timestamp",coincidencetrigger,"goff")
                        t = arbre["data"].GetV1()
                        # do the multiplicity cut (for fast neutrons)
                        for coincidence in range(coincidences_tmp-2):
                            dt_next = t[coincidence+1]-t[coincidence]
                            if coincidence<2:
                                if dt_next>dtcut and dt_next>0:
                                    coincidences+=1
                            else:
                                dt_prev_prev = t[coincidence-1]-t[coincidence-2]
                                if dt_next>dtcut and dt_next>0 and dt_prev_prev>dtcut:
                                    coincidences+=1
                        minAchive = 1
                        if coincidences==0:
                            eff = 1.e-6/totalEvents
                            effErr = 1.e-6/sqrt(totalEvents)
                        else:
                            eff = coincidences/totalEvents
                            effErr = 1/totalEvents*sqrt(coincidences*(1-coincidences/totalEvents))
                        hist[tag].Fill(_d,_nx,eff)
                        errorbin = hist[tag].FindBin(_d,_nx)
                        hist[tag].SetBinError(errorbin, effErr)
                    else:
                        eff = 1.e-6/totalEvents
                        effErr = 1/totalEvents*sqrt(1-1e-6*(1-1e-6/totalEvents))
                        hist[tag].Fill(_d,_nx,eff)
                        errorbin = hist[tag].FindBin(_d,_nx)
                        hist[tag].SetBinError(errorbin, effErr)
                else:
                    eff = 1.e-6/totalEvents
                    effErr = 1/totalEvents*sqrt(1-1e-6*(1-1e-6/totalEvents))
                    hist[tag].Fill(_d,_nx,eff)
                    errorbin = hist[tag].FindBin(_d,_nx)
                    hist[tag].SetBinError(errorbin, effErr)
        # offset done

        outfile.cd()
        hist[tag].Write()

    arbre["rfile"].Close()
    del arbre
    return eff


def obtainCorrelatedEventEfficiency(cover,file,_tag,outfile,dtcut,dscut):
    
    para = testEnabledCondition(arguments)
    additionalString  = para[0]
    arbre = {}
    arbre["rfile"] = TFile(file)
    print('Reading', file)
    hist = {}
    try:
        runSummary = arbre["rfile"].Get('runSummary')
        Entries = runSummary.GetEntries()
        runSummary.GetEntry(Entries-1)
        events = 0
        _eventPerRun = runSummary.nEvents
    except:
        print('File',file,'did not have run associated with it. Returning empty histogram.')
        rangeRmin,rangeRmax = 0.45,3.05
        binwidthR = float(arguments["--binwidthR"])
        binR = int((rangeRmax-rangeRmin)/binwidthR)
        rangeNmin,rangeNmax = 7.5,75.5
        binwidthN = float(arguments["--binwidthN"])
        binN = int((rangeNmax-rangeNmin)/binwidthN)
        hist[_tag] = TH2D('hist%s'%(_tag),'EMPTY - Coincidence efficiency -  %s '%(_tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        hist[_tag].SetXTitle('distance from wall [m]')
        hist[_tag].SetYTitle('%s cut'%(_energyEstimator))
        hist[_tag].SetZTitle('efficiency *1e6')
        hist[_tag].GetZaxis().SetTitleOffset(-.55);
        hist[_tag].GetZaxis().SetTitleColor(1);
        hist[_tag].GetZaxis().CenterTitle();
        hist[_tag].SaveAs("bonsai_root_files%s/%s/hist%s.C"%(additionalString,cover,_tag))
        return -1

    totalEvents = Entries*_eventPerRun

    arbre["data"]    = arbre["rfile"].Get('data')

    _someEntries = arbre["data"].GetEntries()

    binwidthR = float(arguments["--binwidthR"])
    rangeRmin,rangeRmax = 0.5-binwidthR/2.,3.0+binwidthR/2.
    binR = int((rangeRmax-rangeRmin)/binwidthR)
    binwidthN = float(arguments["--binwidthN"])
    rangeNmin,rangeNmax = minNX-binwidthN/2.,maxNX+binwidthN/2.
    binN = int((rangeNmax-rangeNmin)/binwidthN)
    binwidthT = 10
    binwidthS = 0.1 
    scanStep = int(arguments['--scanStep'])
    scanMax = int(arguments["--scanMaxValue"])
    scanMin = int(arguments["--scanMinValue"])
    if  arguments["--positiveScan"]:
        offsets_nx = np.arange(scanMin,scanMax+1,scanStep,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        offsets_nx = np.arange(-scanMax,-scanMin+1,scanStep,dtype=int)## Correction to add last value
    else:
        offsets_nx = np.arange(-scanMax,scanMax+1,scanStep,dtype=int)## Correction to add last value


    # now we can get the efficiencies of coincident correlated events
    for offset in offsets_nx:
        _offset = str(offset).replace("-","Minus")
        tag = _tag+'_%soffset%s'%(_energyEstimator,_offset)
        hist[tag] = TH2D('hist%s'%(tag),'Coincidence Efficiency -  %s '%(tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        hist[tag].SetXTitle('distance from wall [m]')
        hist[tag].SetYTitle('prompt %s cut'%(_energyEstimator))
        hist[tag].SetZTitle('efficiency')
        hist[tag].GetZaxis().SetTitleOffset(-.55);
        hist[tag].GetZaxis().SetTitleColor(1);
        hist[tag].GetZaxis().CenterTitle();

        for _d in drange(rangeRmin+binwidthR/2.,rangeRmax,binwidthR):
            minAchieve = 0

            for _nx in drange(int(rangeNmin+binwidthN/2.0),int(rangeNmax),binwidthN):
                coincidences=0
                trigger  = "closestPMT/1000.>%f"%(_d)
                trigger  += "&& good_pos>%f " %(_posGood)
                trigger  += "&& inner_hit > 4 &&  veto_hit < 4"
                delayedtrigger = trigger
                delayedtrigger  += "&& %s > %f" %(_energyEstimator,_nx+offset) # all of the events pass the delayed nx cut
                if int(arguments["--QFIT"])==1:
                    delayedtrigger = "closestPMTQFit/1000.>%f"%(_d)
                    delayedtrigger += "&& QFit ==1 "
                    delayedtrigger += "&& inner_hit > %f &&  veto_hit < 4" %(_nx+offset)

                # find the coincidence efficiency
                if _someEntries !=0:
                    if minAchieve ==0:
                        coincidencetrigger =  delayedtrigger
                        coincidencetrigger += "&& closestPMT_prev/1000.>%f"%(_d)
                        coincidencetrigger += "&& good_pos_prev>%f"%(_posGood)
                        coincidencetrigger += "&& inner_hit_prev > 4 && veto_hit_prev <4"
                        coincidencetrigger += "&& %s_prev > %f"%(_energyEstimator,_nx) # preceding event passes the prompt nx cut
                        coincidencetrigger += "&& dt_prev_us > 0 && dt_prev_us < %f"%(dtcut) # time coincidence
                        coincidencetrigger += "&& drPrevr/1000. < %f"%(dscut) # spatial coincidence
                        if 'wbls_l' in arguments["--detectMedia"]:
                            coincidencetrigger += "&& best_like-average_like_05m > 20"
                        if int(arguments["--QFIT"])==1:
                            coincidencetrigger = delayedtrigger
                            coincidencetrigger += "&& closestPMTQFit_prev/1000.>%f"%(_d)
                            coincidencetrigger += "&& inner_hit_prev > %f &&  veto_hit < 4" %(_nx)
                            coincidencetrigger += "&& dt_prev_us > 0 && dt_prev_us < %f"%(dtcut) # time coincidence
                            coincidencetrigger += "&& drPrevrQFit/1000. < %f"%(dscut) # spatial coincidence
                        coincidences_tmp = arbre["data"].Draw("timestamp",coincidencetrigger,"goff")
                        t = arbre["data"].GetV1()
                        # do the multiplicity cut (for fast neutrons)
                        for coincidence in range(coincidences_tmp-2):
                            dt_next = t[coincidence+1]-t[coincidence]
                            if coincidence<2:
                                if dt_next>dtcut and dt_next>0:
                                    coincidences+=1
                            else:
                                dt_prev_prev = t[coincidence-1]-t[coincidence-2]
                                if dt_next>dtcut and dt_next>0 and dt_prev_prev>dtcut:
                                    coincidences+=1
                        minAchive = 1
                        if coincidences==0:
                            coincidences = 1.e-6
                        eff = coincidences/totalEvents
                        effErr = 1/totalEvents*sqrt(coincidences*(1-coincidences/totalEvents))
                        hist[tag].Fill(_d,_nx,eff)
                        errorbin = hist[tag].FindBin(_d,_nx)
                        hist[tag].SetBinError(errorbin, effErr)
                    else:
                        coincidences = 1.e-6
                        eff = coincidences/totalEvents
                        effErr = 1/totalEvents*sqrt(coincidences*(1-coincidences/totalEvents))
                        hist[tag].Fill(_d,_nx,eff)
                        errorbin = hist[tag].FindBin(_d,_nx)
                        hist[tag].SetBinError(errorbin, effErr)
                else:
                    coincidences = 1.e-6
                    eff = coincidences/totalEvents
                    hist[tag].Fill(_d,_nx,eff)
                    effErr = 1/totalEvents*sqrt(coincidences*(1-coincidences/totalEvents))
                    errorbin = hist[tag].FindBin(_d,_nx)
                    hist[tag].SetBinError(errorbin, effErr)
        # offset done

        outfile.cd()
        hist[tag].Write()

    arbre["rfile"].Close()
    del arbre
    return eff


def obtainAccidentalCoincidenceEfficiency(cover,file,_tag,pdffile,outfile,dtcut,dscut):

    para = testEnabledCondition(arguments)
    additionalString  = para[0]
    arbre = {}
    arbre["rfile"] = TFile(file)
    print('Reading', file)
    hist = {}
    try:
        runSummary = arbre["rfile"].Get('runSummary')
        Entries = runSummary.GetEntries()
        runSummary.GetEntry(Entries-1)
        events = 0
        _eventPerRun = runSummary.nEvents
    except:
        print('File',file,'did not have run associated with it. Returning empty histogram.')
        rangeRmin,rangeRmax = 0.45,3.05
        binwidthR = float(arguments["--binwidthR"])#(rangeRmax-rangeRmin)/binR
        binR = int((rangeRmax-rangeRmin)/binwidthR)
        rangeNmin,rangeNmax = 7.5,75.5
        binwidthN = float(arguments["--binwidthN"])#(rangeNmax-rangeNmin)/binN
        binN = int((rangeNmax-rangeNmin)/binwidthN)
        hist[_tag] = TH2D('hist%s'%(_tag),'EMPTY - Coincidence efficiency -  %s '%(_tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        hist[_tag].SetXTitle('distance from wall [m]')
        hist[_tag].SetYTitle('%s cut'%(_energyEstimator))
        hist[_tag].SetZTitle('efficiency *1e6')
        hist[_tag].GetZaxis().SetTitleOffset(-.55);
        hist[_tag].GetZaxis().SetTitleColor(1);
        hist[_tag].GetZaxis().CenterTitle();
        hist[_tag].SaveAs("bonsai_root_files%s/%s/hist%s.C"%(additionalString,cover,_tag))
        return -1

    totalEvents = float(Entries*_eventPerRun)

    arbre["data"]    = arbre["rfile"].Get('data')

    _someEntries = arbre["data"].GetEntries()

    binwidthR = float(arguments["--binwidthR"])
    rangeRmin,rangeRmax = 0.5-binwidthR/2.,3.0+binwidthR/2.
    binR = int((rangeRmax-rangeRmin)/binwidthR)
    binwidthN = float(arguments["--binwidthN"])
    rangeNmin,rangeNmax = minNX-binwidthN/2.,maxNX+binwidthN/2.
    binN = int((rangeNmax-rangeNmin)/binwidthN)
    binwidthT = 10
    binwidthS = 0.1 
    scanStep = int(arguments['--scanStep'])
    scanMax = int(arguments["--scanMaxValue"])
    scanMin = int(arguments["--scanMinValue"])
    if  arguments["--positiveScan"]:
        offsets_nx = np.arange(scanMin,scanMax+1,scanStep,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        offsets_nx = np.arange(-scanMax,-scanMin+1,scanStep,dtype=int)## Correction to add last value
    else:
        offsets_nx = np.arange(-scanMax,scanMax+1,scanStep,dtype=int)## Correction to add last value

    # now we can get the efficiency of single events
    for offset in offsets_nx:
        _offset = str(offset).replace("-","Minus")
        tag = _tag+'_%soffset%s'%(_energyEstimator,_offset)
        hist[tag] = TH2D('hist%s'%(tag),'Singles Efficiency -  %s '%(tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        hist[tag].SetXTitle('distance from wall [m]')
        hist[tag].SetYTitle('prompt %s cut'%(_energyEstimator))
        hist[tag].SetZTitle('efficiency')
        hist[tag].GetZaxis().SetTitleOffset(-.55);
        hist[tag].GetZaxis().SetTitleColor(1);
        hist[tag].GetZaxis().CenterTitle();
    

        for _d in drange(rangeRmin+binwidthR/2.,rangeRmax,binwidthR):
            minAchieve = 0

            for _nx in drange(int(rangeNmin+binwidthN/2.0),int(rangeNmax),binwidthN):
                trigger  = "closestPMT/1000.>%f"%(_d)
                trigger  += "&& good_pos>%f " %(_posGood)
                trigger  += "&& inner_hit > 4 &&  veto_hit < 4"
                trigger += "&& %s > %f" %(_energyEstimator,_nx) # all of the events pass the prompt nx cut
                if int(arguments["--QFIT"])==1:
                    trigger = "closestPMTQFit/1000.>%f"%(_d)
                    trigger += "&& QFit ==1 "
                    trigger += "&& inner_hit > %f &&  veto_hit < 4" %(_nx+offset)

                # find the singles efficiency
                if _someEntries !=0:
                    coincidences=0
                    if minAchieve ==0:


                        if int(arguments["--QFIT"])==1:
                            evts             = arbre["data"].Draw("xQFit/1000.:yQFit/1000.:zQFit/1000.",delayedtrigger,"goff")
                        else:
                            evts             = arbre["data"].Draw("x/1000.:y/1000.:z/1000.",trigger,"goff")
                        x             = arbre["data"].GetV1()
                        y             = arbre["data"].GetV2()
                        z             = arbre["data"].GetV3()
                        if evts!=0:
                            coincidences = coincidenceCalculation(pdffile,evts,x,y,z,_nx+offset,_d,dscut)

                        if coincidences == 0:
                            coincidences = 0.01
                            minAchieve = 1
                        effErr = 1./totalEvents*sqrt(coincidences*(1-coincidences/totalEvents))
                        eff = coincidences/totalEvents
                        hist[tag].Fill(_d,_nx,eff)
                        errorbin = hist[tag].FindBin(_d,_nx)
                        hist[tag].SetBinError(errorbin, effErr)

                    else:
                        coincidences = 0.01
                        eff = coincidences/totalEvents
                        hist[tag].Fill(_d,_nx,eff)
                        effErr = 1./totalEvents*sqrt(coincidences*(1-coincidences/totalEvents))
                        errorbin = hist[tag].FindBin(_d,_nx)
                        hist[tag].SetBinError(errorbin, effErr)
                else:
                    coincidences = 0.01
                    eff = coincidences/totalEvents
                    hist[tag].Fill(_d,_nx,eff)
                    effErr = 1./totalEvents*sqrt(coincidences*(1-coincidences/totalEvents))
                    errorbin = hist[tag].FindBin(_d,_nx)
                    hist[tag].SetBinError(errorbin, effErr)
        outfile.cd()
        hist[tag].Write()

    arbre["rfile"].Close()
    del arbre
    return eff

def obtainEventPDFSingles(singlesfile,outfile):
   
    data       = singlesfile.Get('data')

    minNX = float(arguments['--minNX'])
    maxNX = float(arguments['--maxNX'])

    scanStep = int(arguments['--scanStep'])
    scanMax = int(arguments["--scanMaxValue"])
    scanMin = int(arguments["--scanMinValue"])
    if  arguments["--positiveScan"]:
        nxRange = np.arange(minNX+scanMin,maxNX+scanMax+1,scanStep,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        nxRange = np.arange(minNX-scanMax,maxNX-scanMin+1,scanStep,dtype=int)## Correction to add last value
    else:
        offsets_nx = np.arange(minNX-scanMax,maxNX+scanMax+1,scanStep,dtype=int)## Correction to add last value
    

    # create the bins in x, y and z
    # variable bins to allow for variable detector size
    xMin,xMax = -detectorRadius/1000., detectorRadius/1000.
    yMin,yMax = -detectorRadius/1000., detectorRadius/1000.
    zMin,zMax = -detectorHeight/1000., detectorHeight/1000.
    binwidthxyz = 0.1
    binx = int((xMax-xMin)/binwidthxyz)
    biny = int((yMax-yMin)/binwidthxyz)
    binz = int((zMax-zMin)/binwidthxyz)
    hist={}
    for nxval in nxRange: 
   
        tag = '%s_%d'%(_energyEstimator,nxval)
        hist[tag] = TH3D('hist%s'%(tag), 'frequency distribution of singles events vertices-  %s '%(tag),binx,xMin,xMax,biny,yMin,yMax,binz,zMin,zMax)
        hist[tag].SetXTitle('x [m]')
        hist[tag].SetYTitle('y [m]')
        hist[tag].SetZTitle('z [m]')
        hist[tag].GetZaxis().SetTitleOffset(-.55);
        hist[tag].GetZaxis().SetTitleColor(1);
        hist[tag].GetZaxis().CenterTitle();

        # get the data for the vertex frequency distribution
        trigger  = "good_pos>%f " %(_posGood)
        trigger += "&& inner_hit > 4 &&  veto_hit < 4"
        trigger += "&& x>-9999"
        trigger += "&& %s > %d"%(_energyEstimator,nxval)
        # get the reconstructed rho^2, timestamp nx from the bonsai file
        evts     = data.Draw("x/1000.:y/1000.:z/1000.",trigger,"goff")
        x        = data.GetV1()
        y        = data.GetV2()
        z        = data.GetV3()
        for ev in range(evts):
            hist[tag].Fill(x[ev],y[ev],z[ev],1)
        outfile.cd()
        hist[tag].Write()

    # end loop over events which pass the cut
    del data
    outfile.Write()
    
    return 0

def obtainEventPDF(files,outfile):

    print('detector radius: ', detectorRadius,', detector height: ',detectorHeight)
    # makes pdf of all accidental background events in rho2,z2 and nx
    # also makes pdf of IBD positrons for signal evaluation using IBD singles if required

    minNX = float(arguments['--minNX'])
    maxNX = float(arguments['--maxNX'])

    scanStep = int(arguments["--scanStep"])
    scanMax = int(arguments["--scanMaxValue"])
    scanMin = int(arguments["--scanMinValue"])
  
    if  arguments["--positiveScan"]:
        nxRange = np.arange(minNX+scanMin,maxNX+scanMax+1,scanStep,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        nxRange = np.arange(minNX-scanMax,minNX-scanMin+1,scanStep,dtype=int)## Correction to add last value
    else:
        nxRange = np.arange(minNX-scanMax,maxNX+scanMax+1,scanStep,dtype=int)## Correction to add last value
 

    # create the bins in x, y and z
    # variable bins to allow for variable detector size
    xMin,xMax = -detectorRadius/1000., detectorRadius/1000.
    yMin,yMax = -detectorRadius/1000., detectorRadius/1000.
    zMin,zMax = -detectorHeight/1000., detectorHeight/1000.
    binwidthxyz = 0.1
    binx = int((xMax-xMin)/binwidthxyz)
    biny = int((yMax-yMin)/binwidthxyz)
    binz = int((zMax-zMin)/binwidthxyz)
    
    # find the minimum time simulated and set as the maximum time to look at
    maxRate=0.
    if arguments["--lightSim"]:
        newRates = {k:v for (k,v) in rates.items() if 'PMT' in k} # assumes that the PMT rates are highest
        maxRate = max(newRates.values())[0]
        print(maxRate)
    else:
        maxRate = max(rates.values())[0]
        print(maxRate)
    firstgo=1

#    hist = TH3D('hist', 'frequency distribution of singles events vertices'   ,binx,xMin,xMax,biny,yMin,yMax,binz,zMin,zMax)
    hist={}
    for nxval in nxRange: 
        print(nxval)   
        tag = '%s_%d'%(_energyEstimator,nxval)
        hist[tag] = TH3D('hist%s'%(tag), 'frequency distribution of singles events vertices-  %s '%(tag),binx,xMin,xMax,biny,yMin,yMax,binz,zMin,zMax)
        hist[tag].SetXTitle('x [m]')
        hist[tag].SetYTitle('y [m]')
        hist[tag].SetZTitle('z [m]')
        hist[tag].GetZaxis().SetTitleOffset(-.55);
        hist[tag].GetZaxis().SetTitleColor(1);
        hist[tag].GetZaxis().CenterTitle();
        for filename in glob.glob(files):
            bonsaifile = TFile(filename)
            try:
                runSummary = bonsaifile.Get('runSummary')
                Entries = runSummary.GetEntries()
                runSummary.GetEntry(Entries-1)
                events = 0
                _eventPerRun = runSummary.nEvents
            except:
                print('File',filename,'did not have run associated with it. No events added to histogram.')
                continue
            if firstgo:
                maxTime = float(Entries)*float(_eventPerRun)/maxRate*1e6 #simulation time in us at fastest simulated rate
                firstgo=0

            if 'ibd' in filename or 'fast' in filename or 'A_Z' in filename:
                continue
            print('Adding events from ', filename)

            data       = bonsaifile.Get('data')

            # get the data for the vertex frequency distribution
            trigger  = "good_pos>%f " %(_posGood)
            trigger += "&& inner_hit > 4 &&  veto_hit < 4"
            trigger += "&& timestamp<%f"%(maxTime)
            trigger += "&& %s>%d"%(_energyEstimator,nxval)
            # get the reconstructed rho^2,  nx from the bonsai file
            if int(arguments["--QFIT"])==1:
                trigger = "inner_hit > 4 &&  veto_hit < 4"
                trigger += "&& timestamp<%f"%(maxTime)
                trigger += "&& inner_hit>%d"%(nxval)
                evts     = data.Draw("xQFit/1000.:yQFit/1000.:zQFit/1000.",trigger,"goff")
                print("using QFit results")
            else:
                evts     = data.Draw("x/1000.:y/1000.:z/1000.",trigger,"goff")
            x        = data.GetV1()
            y        = data.GetV2()
            z        = data.GetV3()
            for ev in range(evts):
                hist[tag].Fill(x[ev],y[ev],z[ev],1)
            # end loop over events which pass the cut
            bonsaifile.Close()
            del data
        outfile.cd()
        hist[tag].Write()

    
    #outfile.Write()
    return 0


def obtainCoincidenceMap(singlesfile,outfile,dscut,dtcut):
    minNX = float(arguments['--minNX'])
    maxNX = float(arguments['--maxNX'])
    scanStep = int(arguments['--scanStep'])

    print('Obtaining coincidences in range around %d us dt cut'%(dtcut)) 
    hist = {}
    binwidthR = float(arguments["--binwidthR"])
    rangeRmin,rangeRmax = 0.5-binwidthR/2.,2.5+binwidthR/2.
    binR = int((rangeRmax-rangeRmin)/binwidthR)
    binwidthN = float(arguments["--binwidthN"])
    rangeNmin,rangeNmax = minNX-binwidthN/2.,maxNX+binwidthN/2.
    binN = int((rangeNmax-rangeNmin)/binwidthN)
    binwidthS = 0.5

    scanStep = int(arguments['--scanStep'])
    scanMax = int(arguments["--scanMaxValue"])
    scanMin = int(arguments["--scanMinValue"])
    if  arguments["--positiveScan"]:
        offsets_nx = np.arange(scanMin,scanMax+1,scanStep,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        offsets_nx = np.arange(-scanMax,-scanMin+1,scanStep,dtype=int)## Correction to add last value
    else:
        offsets_nx = np.arange(-scanMax,scanMax+1,scanStep,dtype=int)## Correction to add last value

    singlesfile = TFile(singlesfile)
    singlesdata    = singlesfile.Get('data')
    runSummary = singlesfile.Get('runSummary')
    entries = runSummary.GetEntries()
    runSummary.GetEntry(entries-1)
    eventsPerRun = runSummary.nEvents
    totalEvents = entries*eventsPerRun
    
    # get the probability of a prompt event
    # having a delayed event in coincidence
    for offset in offsets_nx:
        _offset = str(offset).replace("-","Minus")
        tag = '%soffset%s_%d_%d'%(_energyEstimator,_offset,dtcut,dscut*1000)
        hist[tag] = TH2D('hist%s'%(tag),'Coincidence Probability -  %s '%(tag),binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
        hist[tag].SetXTitle('distance from wall [m]')
        hist[tag].SetYTitle('prompt %s cut'%(_energyEstimator))
        hist[tag].SetZTitle('probability')
        hist[tag].GetZaxis().SetTitleOffset(-.55);
        hist[tag].GetZaxis().SetTitleColor(1);
        hist[tag].GetZaxis().CenterTitle();
        for fidcut in drange(rangeRmin+binwidthR/2.,rangeRmax,binwidthR):
            minAchieve = 0

            for nxcut in drange(int(rangeNmin+binwidthN/2.0),int(rangeNmax),binwidthN):

                nCoincidences=0

                prompttrigger  = "closestPMT/1000.>%f"%(fidcut)
                prompttrigger  += "&& good_pos>%f " %(_posGood)
                prompttrigger  += "&& inner_hit > 4 &&  veto_hit < 4"
                prompttrigger  += "&& %s > %f" %(_energyEstimator,nxcut)
                if 'wbls_l' in arguments["--detectMedia"]:
                    print('adding condition from Kat to prompt trigger')
                    prompttrigger += best_like-average_like_05m > 20
                if int(arguments["--QFIT"])==1:
                    evts        = singlesdata.Draw("inner_hit:timestamp",prompttrigger,"goff")
                else:
                    evts        = singlesdata.Draw("%s:timestamp"%(_energyEstimator),prompttrigger,"goff")

                nx              = singlesdata.GetV1()
                time            = singlesdata.GetV2()

                if evts!=0:
                    for ev in range(0,evts-1):
                        for ev_d in range(ev+1,evts-1):
                            dt = time[ev_d]-time[ev]
                            if dt>dtcut:
                                continue
                            if dt > 0 and dt < dtcut:# and nx[ev_d]>(nxcut+offset): 
                               nCoincidences+=1
                                   
                if nCoincidences!=0:
                    # find probability of a coincidence
                    coincidenceEff = nCoincidences/float(totalEvents)
                    error = 1./totalEvents*sqrt(nCoincidences*(1-nCoincidences/float(totalEvents)))
                    hist[tag].Fill(fidcut,nxcut,coincidenceEff)
                    errorbin = hist[tag].FindBin(fidcut,nxcut)
                    hist[tag].SetBinError(errorbin,error)
                else:
                    coincidenceEff = 0.1/totalEvents
                    error = 1./totalEvents*sqrt(0.1*(1-0.1/float(totalEvents)))
                    hist[tag].Fill(fidcut,nxcut,coincidenceEff)
                    errorbin = hist[tag].FindBin(fidcut,nxcut)
                    hist[tag].SetBinError(errorbin,error)
            else:
                coincidenceEff = 0.1/totalEvents
                error = 1./totalEvents*sqrt(0.1*(1-0.1/float(totalEvents)))
                hist[tag].Fill(fidcut,nxcut,coincidenceEff)
                errorbin = hist[tag].FindBin(fidcut,nxcut)
                hist[tag].SetBinError(errorbin,error)
        outfile.cd()
        hist[tag].Write()
    # end loops to create histogram of coincidence probabilities

    return 0

def coincidenceCalculation(pdffile,evts,x,y,z,nxcut,fidcut,dscut):
    
    _energyEstimator = arguments['--energyEst'] 
    h = pdffile.Get('hist%s_%d'%(_energyEstimator,nxcut))
    fidRadius = pmtRadius/1000.-fidcut
    fidHeight = pmtHeight/1000.-fidcut
    
    # only look at events which have passed the prompt trigger
    xMin,xMax           = -fidRadius/1000.,fidRadius/1000.
    yMin,yMax           = -fidRadius/1000.,fidRadius/1000.
    zMin,zMax           = -fidHeight/1000.,fidHeight/1000. # PSUP base to top
    nxMin               = nxcut
    dPMTbinMin          = dscut
    xBinMin,xBinMax     = h.GetXaxis().FindBin(xMin), h.GetXaxis().FindBin(xMax)
    yBinMin,yBinMax     = h.GetYaxis().FindBin(yMin), h.GetYaxis().FindBin(yMax)
    zBinMin, zBinMax      = h.GetZaxis().FindBin(zMin), h.GetZaxis().FindBin(zMax)
    h.GetXaxis().SetRange(xBinMin,xBinMax)
    h.GetYaxis().SetRange(yBinMin, yBinMax)
    h.GetZaxis().SetRange(zBinMin, zBinMax)
    ncoincidences=0
    for ev in range(evts):
        count=0
        if (ev % 20000) == 0 and ev>0:
            print('getting coincidences for ev %d, %s %d, closestPMT %f'%(ev,_energyEstimator,nxcut,fidcut))
        # for each event, randomly generate an event from the pdf 
        # assign as a coincidence if it passes coincidence cuts
        xRandom=Double(0.)
        yRandom=Double(0.)
        zRandom=Double(0.)
        h.GetRandom3(xRandom,yRandom,zRandom)
        ds = sqrt(pow(xRandom-x[ev],2)+pow(yRandom-y[ev],2)+pow(zRandom-z[ev],2))
        # check whether the random event is in spatial coincidence
        if ds<dscut:
            ncoincidences+=1

    return ncoincidences


def integralCoincidence(R,lowerBound,upperBound):
    low = -exp(-lowerBound) * (1+lowerBound )
    up  = -exp(-upperBound) * (1+upperBound)
    return up - low


def histIntegral(s,f,cut):
    # H = {}
    Hs = f.Get(s)

    if Hs!= None:
        a = Hs.Integral(0,int(cut*10))
        N = Hs.Integral(0,5000)
    else:
        a = 1.0
        N = 1.0

    # print 'histInt',cut,s,Hs,a,N
    if N !=0:
        return (1.0 - a/N )
    else:
        return 0
