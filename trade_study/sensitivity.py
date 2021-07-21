#import watchmakers as PR

from watchmakers.load import *
from watchmakers.analysis import *

from ROOT import kOrange as kO,kBlue as kB,kGreen as kG
from ROOT import kMagenta as kM,kAzure as kA,kRed as kR
from ROOT import TCanvas,TH2D, gStyle,TFile,gROOT
#from .NeutrinoOscillation import nuOsc


from .io_operations import testEnabledCondition
   
import watchmakers.NeutrinoOscillation as nuOsc

from decimal import *
setcontext(ExtendedContext)

from numpy import max

t = arguments['--timeScale']

fidRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])-float(arguments['--fidThick'])
fidHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])-float(arguments['--fidThick'])

pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])
pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])

detectorRadius  = float(arguments['--tankRadius'])-float(arguments['--steelThick'])
detectorHeight  = float(arguments['--halfHeight'])-float(arguments['--steelThick'])

_energyEstimator = arguments['--energyEst']
dtcut = float(arguments['-t'])
dscut = float(arguments['-d'])
minNX = float(arguments['--minNX'])
maxNX = float(arguments['--maxNX'])

d,proc,coverage,rates = loadSimulationParametersNew()

def drange(start, stop, step):
    rii= start
    while rii<stop:
        yield rii
        rii+=step

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
        outfile = TFile(_str,"RECREATE")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:

                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _tag = _tag.replace(" ","")
                    _file = "bonsai_root_files%s/merged_%s_%s_%s_%s.root"%(additionalMacStr,_cover,_element,_loc,_p)
                    _file = _file.replace(" ","")
                    print(_tag," from ",_file)
                    obtainEventEfficiency(_cover,_file,_tag,outfile)

                    print('')
    print('Saving outfile:',_str)
    outfile.Close()
    return 0

def pdfLassenSingles():

    print('Generating pdfs for accidentals')#IBDs, fast neutrons and radionuclides will be added later

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"

    for idx,_cover in enumerate(coverage):
        _str = "bonsai_root_files%s/pdfs_%s.root"%(additionalMacStr,_cover)
        outfile = TFile(_str,"RECREATE")
        singlestr = "bonsai_root_files%s/merged_%s_singles_SINGLES_singles.root"%(additionalMacStr,_cover)
        singlestr = singlestr.replace(" ","")
        _file = TFile(singlestr)
        print(_file)
        obtainEventPDFSingles(_file,outfile)

    print('Saving outfile:',_str)
    outfile.Close()
    return 0


def pdfLassen():
    #calls obtainPDF() to create and fill pdf as a function of nx, rho2 and z

    print('Generating pdfs for accidentals')#IBDs, fast neutrons and radionuclides will be added later

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"

    for idx,_cover in enumerate(coverage):
        _str = "bonsai_root_files%s/pdfs_%s.root"%(additionalMacStr,_cover)
        outfile = TFile(_str,"RECREATE")
        _files = "bonsai_root_files%s/merged_%s_*.root"%(additionalMacStr,_cover)
        _files = _files.replace(" ","")
        print(_files)
        obtainEventPDF(_files,outfile)

    print('Saving outfile:',_str)
    outfile.Close()
    return 0


def coincidenceMapLassen():
    # calls obtainCoincidenceEfficiency() to create coincidence efficiency map 
    # as a function of nx and closestPMT

    print('Generating coincidence efficiency map')

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"
    for idx,_cover in enumerate(coverage):
        file = "bonsai_root_files%s/merged_%s_singles_SINGLES_singles.root"%(additionalMacStr,_cover)
        _str = "bonsai_root_files%s/results_%s_coincidenceMap_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        outfile = TFile(_str,"RECREATE")
        obtainCoincidenceMap(file,outfile,dscut,dtcut)
    print('Saving outfile:',_str)
    outfile.Close()
    return 0


def pairEfficiencyMapLassen():
    #calls obtainCorrelationEfficiency() to create and fill histograms of coincidence 
    #efficiency as a function of nx and distance to nearest PMT for correlated events

    print('Generating histograms of reconstruction efficiency for all background and signal processes')

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"

    for idx,_cover in enumerate(coverage):
        _str = "bonsai_pair_root_files%s/results_%s_pair_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        outfile = TFile(_str,"UPDATE")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:

                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _tag = _tag.replace(" ","")
                    if 'fast' in _tag or 'singles' in _tag or 'A_Z' in _tag or 'pn_ibd' in _tag:
                        _file = "bonsai_pair_root_files%s/merged_%s_%s_%s_%s.root"%(additionalMacStr,_cover,_element,_loc,_p)
                        _file = _file.replace(" ","")
                        print(_tag," from ",_file)
                        obtainPairEfficiency(_cover,_file,_tag,outfile,dtcut,dscut)
                        print('')
    print('Saving outfile:',_str)
    outfile.Close()
    return 0

def ibdEventEfficiencyMapLassen():
    #calls obtainCorrelationEfficiency() to create and fill histograms of coincidence 
    #efficiency as a function of nx and distance to nearest PMT for correlated events

    print('Generating histograms of reconstruction efficiency for all background and signal processes')

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"

    for idx,_cover in enumerate(coverage):
        _str = "bonsai_root_files%s/results_%s_coincidence_correlated_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        outfile = TFile(_str,"UPDATE")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:

                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _tag = _tag.replace(" ","")
                    if 'pn_ibd' in _tag: # not including fast neutrons
                        if arguments['--Heysham']:
                            if 'hartlepool' not in _tag:
                                _file = "bonsai_root_files%s/merged_%s_%s_%s_%s.root"%(additionalMacStr,_cover,_element,_loc,_p)
                                _file = _file.replace(" ","")
                                print(_tag," from ",_file)
                                obtainCorrelatedEventEfficiency(_cover,_file,_tag,outfile,dtcut,dscut)
                                print('')
                        else:
                            if 'heysham' not in _tag:
                                _file = "bonsai_root_files%s/merged_%s_%s_%s_%s.root"%(additionalMacStr,_cover,_element,_loc,_p)
                                _file = _file.replace(" ","")
                                print(_tag," from ",_file)
                                obtainCorrelatedEventEfficiency(_cover,_file,_tag,outfile,dtcut,dscut)
                                print('')
    print('Saving outfile:',_str)
    outfile.Close()
    return 0

def radionuclideEventEfficiencyMapLassen():
    #calls obtainCorrelationEfficiency() to create and fill histograms of coincidence 
    #efficiency as a function of nx and distance to nearest PMT for correlated events

    print('Generating histograms of reconstruction efficiency for all background and signal processes')

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"

    for idx,_cover in enumerate(coverage):
        _str = "bonsai_root_files%s/results_%s_coincidence_correlated_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        outfile = TFile(_str,"UPDATE")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:

                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _tag = _tag.replace(" ","")
                    if 'A_Z' in _tag:
                        _file = "bonsai_root_files%s/merged_%s_%s_%s_%s.root"%(additionalMacStr,_cover,_element,_loc,_p)
                        _file = _file.replace(" ","")
                        print(_tag," from ",_file)
                        obtainCorrelatedEventEfficiency(_cover,_file,_tag,outfile,dtcut,dscut)
                        print('')
    print('Saving outfile:',_str)
    outfile.Close()
    return 0

def fastneutronEventEfficiencyMapLassen():
    #calls obtainCorrelationEfficiency() to create and fill histograms of coincidence 
    #efficiency as a function of nx and distance to nearest PMT for correlated events

    print('Generating histograms of reconstruction efficiency for all background and signal processes')

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"

    for idx,_cover in enumerate(coverage):
        _str = "bonsai_root_files%s/results_%s_coincidence_correlated_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        outfile = TFile(_str,"UPDATE")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:

                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _tag = _tag.replace(" ","")
                    if 'FASTNEUTRONS' in _tag:
                        _file = "bonsai_root_files%s/merged_%s_%s_%s_%s.root"%(additionalMacStr,_cover,_element,_loc,_p)
                        _file = _file.replace(" ","")
                        print(_tag," from ",_file)
                        obtainCorrelatedEventEfficiency(_cover,_file,_tag,outfile,dtcut,dscut)
                        print('')
    print('Saving outfile:',_str)
    outfile.Close()
    return 0



def accidentalCoincidenceEfficiencyMapLassen():
    #calls obtainAccidentalCoincidenceEfficiency() to create and fill histograms of coincidence 
    #efficiency as a function of nx and distance to nearest PMT for uncorrelated backgrounds

    print('Generating histograms of reconstruction efficiency for all background and signal processes')

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"

    for idx,_cover in enumerate(coverage):
        _str = "bonsai_root_files%s/results_%s_coincidence_accidentals_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        outfile = TFile(_str,"RECREATE")
#        _singlestr = "bonsai_root_files%s/merged_%s_singles_SINGLES_singles.root"%(additionalMacStr,_cover)
#        singlesfile = TFile(_singlestr)
        _pdfstr = "bonsai_root_files%s/pdfs_%s.root"%(additionalMacStr,_cover)
        pdffile = TFile(_pdfstr)
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:

                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _tag = _tag.replace(" ","")
                    if 'ibd' not in _tag and 'A_Z' not in _tag and 'singles' not in _tag and 'mono' not in _tag and 'FASTNEUTRONS' not in _tag:
                        _file = "bonsai_root_files%s/merged_%s_%s_%s_%s.root"%(additionalMacStr,_cover,_element,_loc,_p)
                        _file = _file.replace(" ","")
                        print(_tag," from ",_file)
                        obtainAccidentalCoincidenceEfficiency(_cover,_file,_tag,pdffile,outfile,dtcut,dscut)
                        print('')
    print('Saving outfile:',_str)
    outfile.Close()
    return 0



def sensitivityMapPass2New():

    site = arguments["--site"]

    # Need to fix this for future running

    OnOffRatio = float(arguments["--OnOff"])
    print(site,'with on-off ratio of ',OnOffRatio)

    cores = int(arguments["--cores"])

    if arguments["--RNRedux"]:
        rnRedux = float(arguments["--RNRedux"])
        if rnRedux>1:
            print("Value of reduction of radionuclide greater than 1, setting to 0")
            rnRedux = 0.0
    else:
        rnRedux = 0.0

    if t == 'sec':
        timeAdjustment = 24*3600.
    if t == 'day':
        timeAdjustment = 1.0
    if t == 'month':
        timeAdjustment = 1./31.
    if t == 'year':
        timeAdjustment = 1./365.
    maxTime = 14400.*timeAdjustment

    print('\nEvaluation based on geoneutrinos.org')
    #parameters  = loadAnalysisParametersNew(t)
    rates       = {}#parameters[11]
    rates["boulby_S"]=1.0
    rates["imb_S"]=1.0
    print('Wrong rates for now')
    sizeDetc    = 2.*pi*pow(fidRadius/1000.,2)*fidHeight/1000./1000.
    sizeTank    = 2.*pi*pow(detectorRadius/1000.,2)*detectorHeight/1000./1000.
    FVkTonRatio = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
    print('Double check:',fidRadius,fidHeight,detectorRadius,detectorHeight)
    boulbyRate,imbRate = rates["boulby_S"]*FVkTonRatio,rates["imb_S"]*FVkTonRatio
    print(' boulby rates: %4.2e per %s per %4.2f kton; [r: %4.2f m; z: %4.2f m]'\
    %(boulbyRate,t,sizeDetc,fidRadius/1000.,fidHeight/1000.))
    #fast neutrons
    # print 'Debug',rates["boulby_S"],FVkTonRatio
    print('\nEvaluation not based on geoneutrinos.org')
    detectorMedium,detectorMass,reactorPower,reactorStandoff = 1,sizeDetc*1000.,1.5,24.98
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff,visual=1)
    preOsc,afterOsc = experiment.FindRate()
    print(' Neutrino rate pre osc: %4.2f; neutrino rate post osc: %4.2f at %4.2f GWth, at %4.2f km, for %4.2f kton' %(preOsc,afterOsc,reactorPower,reactorStandoff,detectorMass/1000.))
    detectorMedium,detectorMass,reactorPower,reactorStandoff = 1,sizeDetc*1000.,1.575,24.98
    experiment = nuOsc.NeutrinoOscillation(detectorMedium,detectorMass,reactorPower,reactorStandoff)
    preOsc,afterOsc = experiment.FindRate()
    print(' Neutrino rate pre osc: %4.2f; neutrino rate post osc: %4.2f at %4.2f GWth, at %4.2f km, for %4.2f kton' %(preOsc,afterOsc,reactorPower,reactorStandoff,detectorMass/1000.))
    print('')

    proc,loca,type,color,lineS,acc,scale   = [],[],[],[],[],[],[]

    proc        += ['QGSP_BERT_EMV','QGSP_BERT_EMX','QGSP_BERT','QGSP_BIC',\
    'QBBC','QBBC_EMZ','FTFP_BERT','QGSP_FTFP_BERT']
    _t          = 'FN%s' % (site)
    loca        += [_t,_t, _t,_t,_t,_t,_t,_t]
    type        += ['si','si','si','si','si','si','si','si']
    acc         += ['corr','corr','corr','corr','corr','corr','corr', 'corr']
    color       += [kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6,kO+6]
    lineS       += [2,2,2,2,2,2,2,2]
    scale       += [1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.,1./8.]

    #radionuclides
    proc        += ['9003','11003']
    _t          =  'RN%s' % (site)
    loca        += [_t,_t]
    type        += ['ei','ei']
    acc         += ['di','di']
    color       += [kG+3,kG+3]
    lineS       += [1,1]
    scale       += [1.,1.]
    #ibds
    if site == 'boulby':
        proc    += ['boulby','boulby','neutron']
    else:
        proc    += ['imb','imb','neutron']
    loca        += ['S','S','N%s'%(site)]
    type        += ['ei','ei','ei']
    acc         += ['di','corr','corr']
    color       += [kA+0,kA-0,kA-0]
    lineS       += [1,2,2]
    scale       += [-1.0,0.0,0.0]

    c1 = TCanvas('c1','c1',1618,1000)
    c1.SetRightMargin(0.53)
    c1.Divide(2,2)

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    if additionalString == "":
        additionalString = "_default"
    location = 'Boulby '

    hist = TH2D('hist','3#sigma discovery phase space -  %s '%(location),31,9.5,40.5,30,9.5,39.5)
    hist.SetXTitle('photocoverage [%]')
    hist.SetYTitle('photoelectron threhsold cut [p.e.]')
    hist.SetZTitle('off-time [%s] for 3 #sigma discovery'%(t))
    hist.GetZaxis().SetTitleOffset(-.55);
    hist.GetZaxis().SetTitleColor(1);
    hist.GetZaxis().CenterTitle();
    gStyle.SetOptStat(0)
    gStyle.SetPalette(55)

    h = {}

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 88,7.5,95.5
    binwidthN = (rangeNmax-rangeNmin)/binN

    _cov = arguments['-C']

    #d,proc,coverage = loadSimulationParametersNew()

    minAchieve = 0
    
    for idx,_cover in enumerate(coverage):
        if (arguments["--QFIT"]==1):
            _strAdd = _cover+"_QFIT"
        else:
            _strAdd = _cover
        _str = "bonsai_root_files%s/%s/results.root"%(additionalMacStr,_strAdd)
        outfile = TFile(_str,"RECREATE")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:
                    _tag = "%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _file = "bonsai_root_files%s/%s/merged_%s_%s_%s.root"%(additionalMacStr,_cover,_loc,_element,_p)
                    print(_tag)
                    obtainEventEfficiency(_cover,_file,_tag,outfile)

                    print('')
        print('Saving outfile:',_str)
        outfile.Close()
    return 0

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


def readPairEfficiencyHistogramLassen():

    hist = {}
    hacc = {}
    hrn = {}
    hfn = {}
    hibd = {}
    hibdBG = {}

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    _energyEstimator = arguments['--energyEst']
    dtcut = float(arguments['-t'])
    dscut = float(arguments['-d'])
    minNX = float(arguments['--minNX'])
    maxNX = float(arguments['--maxNX'])
    scanMax = int(arguments["--scanMaxValue"])
    scanMin = int(arguments["--scanMinValue"])
    if  arguments["--positiveScan"]:
        offsets_nx = np.arange(scanMin,scanMax+1,scanStep,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        offsets_nx = np.arange(-scanMax,-scanMin+1,scanStep,dtype=int)## Correction to add last value
    else:
        offsets_nx = np.arange(-scanMax,scanMax+1,scanStep,dtype=int)## Correction to add last value


    print('Reading in root tree')
    ''' Following loop has two purposes. Read in the histogram from results file and scale them with appropriate rates'''
    for idx,_cover in enumerate(coverage):
        # get the results of previous steps      
        _pairstr = "bonsai_pair_root_files%s/results_%s_pair_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        pairEffFile = TFile(_pairstr,"READ")

        for _offset in offsets_nx:
            offset = str(_offset).replace("-","Minus")
            firstGo  = 1 
            for _p in proc:
                for _loc in proc[_p]:
                    for _element in d[_p]:
                        _tag = "hist%s_%s_%s_%s_%soffset%s"%(_cover,_loc,_element,_p,_energyEstimator,offset)
                        _tag = _tag.replace(" ","")
                        if 'pn_ibd' in _tag or 'A_Z' in _tag or 'singles' in _tag or 'FASTNEUTRONS' in _tag:
                            hist[_tag] = pairEffFile.Get(_tag)
                            _str = "%s_%s_%s"%(_element,_loc,_p)
                            _str = _str.replace(" ","")
                            _scale = rates[_str][0]
                            print('Extracting efficiency histogram from file ',_tag,'and applying rate: ',_scale,'evts/s ... ',_offset, 'offset',end = '\n')

                            try:
                                print(' entries:',hist[_tag].GetEntries(),' ... ', end = '')
                                if firstGo==1:
                                    hacc["hAcc%s"%(offset)] = hist[_tag].Clone()
                                    hacc["hAcc%s"%(offset)].SetZTitle('accidentals rate (Hz)')
                                    hacc["hAcc%s"%(offset)].SetTitle('Accidentals rate')
                                    hacc["hAcc%s"%(offset)].SetName('hAccidentalsRate%s'%(_offset))
                                    hacc["hAcc%s"%(offset)].Reset()
                                    hibd["hIBD%s"%(offset)] = hist[_tag].Clone()
                                    hibd["hIBD%s"%(offset)].SetZTitle('IBD rate (Hz)')
                                    hibd["hIBD%s"%(offset)].SetTitle('rate')
                                    hibd["hIBD%s"%(offset)].SetName('hIBDRate%s'%(_offset))
                                    hibd["hIBD%s"%(offset)].Reset()
                                    hibdBG["hIBDBG%s"%(offset)] = hist[_tag].Clone()
                                    hibdBG["hIBDBG%s"%(offset)].SetZTitle('IBD backgrounds rate (Hz)')
                                    hibdBG["hIBDBG%s"%(offset)].SetTitle('rate')
                                    hibdBG["hIBDBG%s"%(offset)].SetName('hIBDBGRate%s'%(_offset))
                                    hibdBG["hIBDBG%s"%(offset)].Reset()
                                    hrn["hRN%s"%(offset)] = hist[_tag].Clone()
                                    hrn["hRN%s"%(offset)].SetZTitle('Radionuclide rate (Hz)')
                                    hrn["hRN%s"%(offset)].SetTitle('rate')
                                    hrn["hRN%s"%(offset)].SetName('hRNRate%s'%(_offset))
                                    hrn["hRN%s"%(offset)].Reset()
                                    hfn["hFN%s"%(offset)] = hist[_tag].Clone()
                                    hfn["hFN%s"%(offset)].SetZTitle('Fast neutron rate (Hz)')
                                    hfn["hFN%s"%(offset)].SetTitle('rate')
                                    hfn["hFN%s"%(offset)].SetName('hFNRate%s'%(_offset))
                                    hfn["hFN%s"%(offset)].Reset()
                                    firstGo =0
                                hist[_tag].Scale(_scale)
                                if 'pn_ibd' in _tag:
                                    print('%s identified as IBD pair events ,'%(_tag))
                                    if arguments["--Heysham"]:
                                        if 'heysham_signal' in _tag:
#                                           hist[_tag].Smooth()
                                            hibd["hIBD%s"%(offset)].Add(hist[_tag],1)
                                        elif 'heysham_background' in _tag: #or 'boulby_geo' in _tag:
#                                           hist[_tag].Smooth()
                                            hibdBG["hIBDBG%s"%(offset)].Add(hist[_tag],1)
                                    else:
                                        if 'big_' in _tag or 'small' in _tag:
#                                           hist[_tag].Smooth()
                                            hibd["hIBD%s"%(offset)].Add(hist[_tag],1)
                                        elif 'boulby' in tag:
                                            print('identified as boulby background ',_tag)
#                                           hist[_tag].Smooth()
                                            hibdBG["hIBDBG%s"%(offset)].Add(hist[_tag],1)    
                                if 'A_Z' in _tag:
                                    # scale by fraction surviving with 1 second dead time
                                    # of muon detected in veto (Study on the effect of varying veto 
                                    # thickness on the sensitivity of Gd-H2O filled tank, F. Sutanto)
                                    if 'li9' in _tag:
                                        hist[_tag].Scale(0.03)
                                    if 'n17' in _tag:
                                        hist[_tag].Scale(0.85)
                                    print('identified as RN events, ')
#                                   hist[_tag].Smooth()
                                    hrn["hRN%s"%(offset)].Add(hist[_tag],1)
                                if 'FASTNEUTRONS' in _tag:
                                    print('identified as fast neutrons')
                                    hfn["hFN%s"%(offset)].Add(hist[_tag],1)
                                if 'singles' in _tag:
                                    print('identified as accidentals')
#                                   hist[_tag].Smooth()
                                    hacc["hAcc%s"%(offset)].Add(hist[_tag],1) 
                                else:
                                    continue
     
                            except:
                                print("Could not find ",_tag,". Skipping entry.")


    print("\nCompleted reading in of histogram, accidental and IBD identifitication. \n\n")
    
    signalEff = float(arguments["--se"])
    print("Default signal efficiency after proximity cuts (--se) option",signalEff,"\n")

    binwidthR = float(arguments["--binwidthR"])#(rangeRmax-rangeRmin)/binR
    rangeRmin,rangeRmax = 0.5-binwidthR/2.,2.5+binwidthR/2.
    binR = int((rangeRmax-rangeRmin)/binwidthR)
    binwidthN = float(arguments["--binwidthN"])#(rangeNmax-rangeNmin)/float(binN)
    rangeNmin,rangeNmax = minNX-binwidthN/2.,maxNX+binwidthN/2.
    binN = int((rangeNmax-rangeNmin)/binwidthN)
    _cov = arguments['-C']

    _maxSignal,_maxBkgd,_maxSoverB,_maxOffnx,_maxOffdtw = -1,-1,-1,-1,-1
    _maxSignal2,_maxBkgd2,_maxSoverB2,_maxOffnx2,_maxOffdtw2,_max2 = -1,-1,-1,-1,-1,-1
    _maxAcc,_maxAccErr,_maxFN,_maxFNerr,_maxRN,_maxRNerr,_maxIBDbg=-1,-1,-1,-1,-1,-1,-1
    _maxIBDbg, _maxSOB,_maxSOB30days,_maxSOB30daysErr,_maxSigErr,_maxBgErr,_maxSoBErr = -1,-1,-1,-1,-1,-1,-1
    _maxOffset,_maxOffset2=-1,-1

    line,_line,_line2= ("",),"",""
    _histograms = {}

    for offset in offsets_nx:

        _offset = str(offset)
        _offset = _offset.replace("-","Minus")

        _histograms["sOverB%s"%(_offset)]= hibd["hIBD%s"%(_offset)].Clone()
        _histograms["sOverB%s"%(_offset)].SetZTitle('signal/sqrt(signal+background)')
        _histograms["sOverB%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
        _histograms["sOverB%s"%(_offset)].SetTitle('hSoB - %s offset %2d'%(_energyEstimator,offset))
        _histograms["sOverB%s"%(_offset)].SetName('hSoB%s'%(_offset))
        _histograms["sOverB%s"%(_offset)].Reset()

        _histograms["signal%s"%(_offset)]= hibd["hIBD%s"%(_offset)].Clone()
        _histograms["signal%s"%(_offset)].SetZTitle('evts/day')
        _histograms["signal%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
        _histograms["signal%s"%(_offset)].SetTitle('Signal - %s offset %2d'%(_energyEstimator,offset))
        _histograms["signal%s"%(_offset)].SetName('hSignal%s'%(_offset))
        _histograms["signal%s"%(_offset)].Reset()

        _histograms["background%s"%(_offset)]= hibd["hIBD%s"%(_offset)].Clone()
        _histograms["background%s"%(_offset)].SetZTitle('evts/day')
        _histograms["background%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
        _histograms["background%s"%(_offset)].SetTitle('total backgrounds - %s offset %2d'%(_energyEstimator,offset))
        _histograms["background%s"%(_offset)].SetName('hBackground%s'%(_offset))
        _histograms["background%s"%(_offset)].Reset() 



        _proc = '_%s'%(_cov)

        _maxSignal,_maxBkgd,_maxSoverB,_maxnx,_max_dtw = -1,-1,-1,-1,-1

        for _d in range(binR-1):
            for _n in range(binN-offset-1):
                _db,_nb,offset=int(_d+1),int(_n+1),int(offset)
                # get the signal rates 
                _p_d  = hibd["hIBD%s"%(_offset)].GetXaxis().GetBinCenter(_db)
                _p_nx = hibd["hIBD%s"%(_offset)].GetYaxis().GetBinCenter(_nb)
                _n_d  = hibd["hIBD%s"%(_offset)].GetXaxis().GetBinCenter(_db)
                _n_nx = hibd["hIBD%s"%(_offset)].GetYaxis().GetBinCenter(_nb+offset)
                _ibd_rate     = hibd["hIBD%s"%(_offset)].GetBinContent(_db,_nb) #rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                _signal = _ibd_rate*86400.*signalEff
                _ibd_error = hibd["hIBD%s"%(_offset)].GetBinError(_db,_nb)
                _signal_error = _ibd_error*86400.*signalEff

                # now get the background rates 
                _p_d_1  = hacc["hAcc%s"%(_offset)].GetXaxis().GetBinCenter(_db)
                _p_nx_1 = hacc["hAcc%s"%(_offset)].GetYaxis().GetBinCenter(_nb)
                _n_d_1  = hacc["hAcc%s"%(_offset)].GetXaxis().GetBinCenter(_db)
                _n_nx_1 = hacc["hAcc%s"%(_offset)].GetYaxis().GetBinCenter(_nb+offset)
                # first get the individual background contributions
                _acc_rate  = hacc["hAcc%s"%(_offset)].GetBinContent(_db,_nb)*86400*signalEff  #accidentals rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                _acc_error = hacc["hAcc%s"%(_offset)].GetBinError(_db,_nb)*86400
                _rn_rate   = hrn["hRN%s"%(_offset)].GetBinContent(_db,_nb)*86400.   #radionuclides rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                _rn_error = hrn["hRN%s"%(_offset)].GetBinError(_db,_nb)*86400.
                _fn_rate   = hfn["hFN%s"%(_offset)].GetBinContent(_db,_nb)*86400.   #radionuclides rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                _fn_error = hfn["hFN%s"%(_offset)].GetBinError(_db,_nb)*86400.
                _ibd_bg_rate   = hibdBG["hIBDBG%s"%(_offset)].GetBinContent(_db,_nb)*86400. #radionuclides rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                _ibd_bg_error = hibdBG["hIBDBG%s"%(_offset)].GetBinError(_db,_nb)*86400.

                # then get the total background rates
                _background = _acc_rate + _rn_rate + _ibd_bg_rate +_fn_rate# background rate (not including fn)
                _background_error = sqrt(pow(_acc_error,2)+pow(_rn_error,2)+pow(_ibd_bg_error,2)+pow(_fn_error,2))

                # save the rates to the analysis results histograms
                sob = _signal/sqrt(_signal+_background)
                sob30 = _signal*30/sqrt(_signal*30+_background*30)
                _s_plus_b_error = sqrt(pow(_signal_error,2)+pow(_background_error,2))
                sqrt_s_plus_b_error = 0.5*_s_plus_b_error/(_signal+_background)*sqrt(_signal+_background)
                sob_error = sqrt(pow(_signal_error/_signal,2)+pow(sqrt_s_plus_b_error/sqrt(_signal+_background),2))*sob
                _histograms["sOverB%s"%(_offset)].SetBinContent(_db,_nb,sob)
                _histograms["signal%s"%(_offset)].SetBinContent(_db,_nb,_signal)
                _histograms["background%s"%(_offset)].SetBinContent(_db,_nb,_background)
                _histograms["sOverB%s"%(_offset)].SetBinError(_db,_nb,sob_error)
                _histograms["signal%s"%(_offset)].SetBinError(_db,_nb,_signal_error)
                _histograms["background%s"%(_offset)].SetBinError(_db,_nb,_background_error)

                # find the optimal significance as a function of accidentals
                if sob >_maxSoverB:
                    _maxSoverB = sob
                    _maxSignal = _signal
                    _maxBkgd   = _background
                    _maxOffnx   = _n_nx
                    _maxOff_dtw = _n_d
                    _line = ("Offset:%3d, beta/n: Wall-dist (%4.1f,%4.1f) m, %s cut (%d,%d), ibd rate : %4.2f per day;"\
                    %(offset,_p_d,_n_d,_energyEstimator\
                    ,_p_nx,_n_nx\
                    ,_ibd_rate*86400.),)
                    _line2 =    (" acc. rate: %4.3f per day, rn rate: %4.3f per day (all cuts)"%(_acc_rate,_rn_rate),)
                if sob >_maxSoverB2:
                    _maxSoverB2 = sob
                    _maxSignal2 = _signal
                    _maxBkgd2   = _background
                    _maxOffnx2   = _n_nx
                    _maxOff_dtw2 = _n_d
                    _maxOffset2 = offset
                    _maxAcc = _acc_rate
                    _maxFN = _fn_rate
                    _maxFNerr = _fn_error 
                    _maxRN = _rn_rate
                    _maxRNerr = _rn_error
                    _maxIBDbg = _ibd_bg_rate
                    _maxSOB30daysErr = sob_error*sqrt(30)
                    _maxSOB30days = _signal*30/sqrt(_signal*30+_background*30)
                    _maxSigErr = _signal_error
                    _maxBgErr = _background_error
                    _maxSoBErr = sob_error


        print('Offset:',str(offset).rjust(3,' '),',Found max S/sqrt(S+B)','{0:.4f}'.format(_maxSoverB),',(S,B,%s,dtw):('%(_energyEstimator),'{0:.4f}'.format(_maxSignal),'{0:.4f}'.format(_maxBkgd),_maxOffnx,'{0:.1f}'.format(_maxOff_dtw),')')
        line += (_line + _line2,)


    print('\n\nMore info on the maximal sensitivity found:')
    # print line

    for _l in line:
        for i in range(len(_l)):
            print(_l[i], end=' ')
        print('')

    if _maxSignal2 >0:
        T3SIGMA = 9/_maxSignal2 + 9*_maxBkgd2/_maxSignal2**2
    else:
        T3SIGMA = 1e999

    _res = "%s %4.1f %3d %3d %4.3f +/-%4.5f %4.3f +/-%4.7f %4.10f %4.3f %4.3f %4.7f %4.3f +/-%4.5f %4.3f +/-%4.5f %4.1f" %(_cover,_maxOff_dtw2,_maxOffnx2,_maxOffnx2-_maxOffset2,_maxSignal2,_maxSigErr,_maxBkgd2,_maxBgErr,_maxAcc,_maxRN,_maxFN,_maxIBDbg,_maxSOB,_maxSoBErr,_maxSOB30days,_maxSOB30daysErr,T3SIGMA)
    _strRes = "results_DTW_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.txt"%(float(arguments['--vetoThickR']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))
    # _strRes = 'res.txt'
    with open(_strRes,'a') as file:
        file.write(_res+'\n')

    _str = "bonsai_pair_root_files%s/analysis_results_%s_pair_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
    print('\n\nWriting histograms to file',_str)
    f_root = TFile(_str,"recreate")
    for offset in offsets_nx:
        _offset = str(offset).replace("-","Minus")
        _histograms["sOverB%s"%(_offset)].Write()
        _histograms["signal%s"%(_offset)].Write()
        _histograms["background%s"%(_offset)].Write()
        hacc["hAcc%s"%(_offset)].Write()
        hrn["hRN%s"%(_offset)].Write()
        hfn["hFN%s"%(_offset)].Write()
        hibd["hIBD%s"%(_offset)].Write()
        hibdBG["hIBDBG%s"%(_offset)].Write()

    f_root.Close()

    print('\n\n')


def readCoincidenceEfficiencyHistogramLassen():
    timeWindow = float(arguments["-t"])/1e6 # time in microsecond
    timeAcc = timeWindow*86400.

    hist = {}
    hacc = {}
    hrn = {}
    hfn = {}
    hibd = {}
    hibdBG = {}

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    scanMax = int(arguments["--scanMaxValue"])
    scanMin = int(arguments["--scanMinValue"])
    scanStep = int(arguments["--scanStep"])
    if  arguments["--positiveScan"]:
        offsets_nx = np.arange(scanMin,scanMax+1,scanStep,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        offsets_nx = np.arange(-scanMax,-scanMin+1,scanStep,dtype=int)## Correction to add last value
    else:
        offsets_nx = np.arange(-scanMax,scanMax+1,scanStep,dtype=int)## Correction to add last value

    binwidthR = float(arguments["--binwidthR"])# (rangeRmax-rangeRmin)/binR
    rangeRmin,rangeRmax = 0.5-binwidthR/2.,2.5+binwidthR/2.
    binR=int((rangeRmax-rangeRmin)/binwidthR)
    binwidthN = float(arguments["--binwidthN"])#(rangeNmax-rangeNmin)/float(binN)
    rangeNmin,rangeNmax = minNX-binwidthN/2.,maxNX+binwidthN/2.
    binN = int((rangeNmax-rangeNmin)/binwidthN)


    print('Reading in root tree')
    ''' Following loop has two purposes. Read in the histogram from results file and scale them with appropriate rates'''
    for idx,_cover in enumerate(coverage):
        # get the results of previous steps      
        _corrstr = "bonsai_root_files%s/results_%s_coincidence_correlated_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        corrEffFile = TFile(_corrstr,"READ")
        _accstr = "bonsai_root_files%s/results_%s_coincidence_accidentals_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        accEffFile = TFile(_accstr,"READ")
        _mapstr = "bonsai_root_files%s/results_%s_coincidenceMap_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
        print('reading in coincidence map %s, correlated efficiencies %s and uncorrelated efficiencies %s'%(_mapstr,_corrstr,_accstr))
        mapFile = TFile(_mapstr,"READ")

        for _offset in offsets_nx:
            offset = str(_offset).replace("-","Minus")
            firstGo  = 1
            hacc["hAcc%s"%(offset)] = TH2D('hist','Coincidence Probability ',binR,rangeRmin,rangeRmax,binN,rangeNmin,rangeNmax)
            hacc["hAcc%s"%(offset)].SetZTitle('accidentals rate (Hz)')
            hacc["hAcc%s"%(offset)].SetTitle('Accidentals rate')
            hacc["hAcc%s"%(offset)].SetName('hAccidentalsRate%s'%(_offset))
            hacc["hAcc%s"%(offset)].Reset()
            for _p in proc:
                for _loc in proc[_p]:
                    for _element in d[_p]:
                        _tag = 'hist%s_%s_%s_%s_%soffset%s'%(_cover,_loc,_element,_p,_energyEstimator,offset)
                        _tag = _tag.replace(" ","")
                        if 'singles' not in _tag:
                            if 'pn_ibd' in _tag or 'A_Z' in _tag or 'FASTNEUTRONS' in _tag:
                                print('correlated event, getting %s from %s'%(_tag,_corrstr))
                                hist[_tag] = corrEffFile.Get(_tag)
                            else:
                                print('uncorrelated event, getting %s from %s'%(_tag,_accstr))
                                hist[_tag] = accEffFile.Get(_tag)
                            _str = "%s_%s_%s"%(_element,_loc,_p)
                            _str = _str.replace(" ","")
                            _scale = rates[_str][0]
                            print('Extracting efficiency histogram from file ',_tag,'and applying rate: ',_scale,'evts/s ... ',_offset, 'offset',end = '\n')

                            try:
                                print(' entries:',hist[_tag].GetEntries(),' ... ', end = '')
                                if firstGo:
#                                    hacc["hAcc%s"%(offset)] = hist[_tag].Clone()
#                                    hacc["hAcc%s"%(offset)].SetZTitle('accidentals rate (Hz)')
#                                    hacc["hAcc%s"%(offset)].SetTitle('Accidentals rate')
#                                    hacc["hAcc%s"%(offset)].SetName('hAccidentalsRate%s'%(_offset))
#                                    hacc["hAcc%s"%(offset)].Reset()
                                    hibd["hIBD%s"%(offset)] = hist[_tag].Clone()
                                    hibd["hIBD%s"%(offset)].SetZTitle('IBD rate (Hz)')
                                    hibd["hIBD%s"%(offset)].SetTitle('rate')
                                    hibd["hIBD%s"%(offset)].SetName('hIBDRate%s'%(_offset))
                                    hibd["hIBD%s"%(offset)].Reset()
                                    hibdBG["hIBDBG%s"%(offset)] = hist[_tag].Clone()
                                    hibdBG["hIBDBG%s"%(offset)].SetZTitle('IBD backgrounds rate (Hz)')
                                    hibdBG["hIBDBG%s"%(offset)].SetTitle('rate')
                                    hibdBG["hIBDBG%s"%(offset)].SetName('hIBDBGRate%s'%(_offset))
                                    hibdBG["hIBDBG%s"%(offset)].Reset()
                                    hrn["hRN%s"%(offset)] = hist[_tag].Clone()
                                    hrn["hRN%s"%(offset)].SetZTitle('Radionuclide rate (Hz)')
                                    hrn["hRN%s"%(offset)].SetTitle('rate')
                                    hrn["hRN%s"%(offset)].SetName('hRNRate%s'%(_offset))
                                    hrn["hRN%s"%(offset)].Reset()
                                    hfn["hFN%s"%(offset)] = hist[_tag].Clone()
                                    hfn["hFN%s"%(offset)].SetZTitle('Fast neutron rate (Hz)')
                                    hfn["hFN%s"%(offset)].SetTitle('rate')
                                    hfn["hFN%s"%(offset)].SetName('hFNRate%s'%(_offset))
                                    hfn["hFN%s"%(offset)].Reset()
                                    firstGo =0
                                hist[_tag].Scale(_scale)
                                if 'pn_ibd' in _tag:
                                    print('%s identified as IBD pair events ,'%(_tag))
                                    if arguments["--Heysham"]:
                                        if 'heysham_signal' in _tag:
                                            hibd["hIBD%s"%(offset)].Add(hist[_tag],1)
                                        elif 'heysham_background' in _tag: #or 'boulby_geo' in _tag:
                                            hibdBG["hIBDBG%s"%(offset)].Add(hist[_tag],1)
                                    else:
                                        if 'big_' in _tag or 'small' in _tag:
                                            hibd["hIBD%s"%(offset)].Add(hist[_tag],1)
                                        elif 'boulby' in _tag:
                                            hibdBG["hIBDBG%s"%(offset)].Add(hist[_tag],1)    
                                if 'FASTNEUTRONS' in _tag:
                                    hfn["hFN%s"%(offset)].Add(hist[_tag],1)
                                if 'A_Z' in _tag:
                                    # scale by fraction surviving with 1 second dead time
                                    # of muon detected in veto (Study on the effect of varying veto 
                                    # thickness on the sensitivity of Gd-H2O filled tank, F. Sutanto)
                                    if 'li9' in _tag:
                                        hist[_tag].Scale(0.03)
                                    if 'n17' in _tag:
                                        hist[_tag].Scale(0.85)
                                    print('identified as RN events, ')
                                    hrn["hRN%s"%(offset)].Add(hist[_tag],1)
                                if 'A_Z' not in _tag and 'singles' not in _tag and 'FASTNEUTRONS' not in _tag and 'ibd_p' not in _tag and 'ibd_n' not in _tag and 'mono' not in _tag:
                                    print('identified and adding to accidentals')
                                    hacc["hAcc%s"%(offset)].Add(hist[_tag],1) 
                                else:
                                    continue
     
                            except:
                                print("Could not find ",_tag,". Skipping entry.")


    print("\nCompleted reading in of histogram, accidental and IBD identifitication. \n\n")
    
    signalEff = float(arguments["--se"])
    print("Default signal efficiency after proximity cuts (--se) option",signalEff,"\n")

    _cov = arguments['-C']

    _maxSignal,_maxBkgd,_maxSoverB,_maxOffnx,_maxOffdtw = -1,-1,-1,-1,-1
    _maxSignal2,_maxBkgd2,_maxSoverB2,_maxOffnx2,_maxOffdtw2,_max2 = -1,-1,-1,-1,-1,-1
    _maxAcc,_maxAccErr,_maxFN,_maxFNerr,_maxRN,_maxRNerr,_maxIBDbg=-1,-1,-1,-1,-1,-1,-1
    _maxIBDbg, _maxSOB,_maxSOB30days,_maxSOB30daysErr,_maxSigErr,_maxBgErr,_maxSoBErr = -1,-1,-1,-1,-1,-1,-1
    _maxOffset,_maxOffset2=-1,-1

    line,_line,_line2= ("",),"",""
    _histograms = {}

    for offset in offsets_nx:

        binR = hibd["hIBD%s"%(_offset)].GetNbinsX()
        binN = hibd["hIBD%s"%(_offset)].GetNbinsY()
        _offset = str(offset)
        _offset = _offset.replace("-","Minus")

        _histograms["sOverB%s"%(_offset)]= hibd["hIBD%s"%(_offset)].Clone()
        _histograms["sOverB%s"%(_offset)].SetZTitle('signal/sqrt(signal+background)')
        _histograms["sOverB%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
        _histograms["sOverB%s"%(_offset)].SetTitle('hSoB - %s offset %2d'%(_energyEstimator,offset))
        _histograms["sOverB%s"%(_offset)].SetName('hSoB%s'%(_offset))
        _histograms["sOverB%s"%(_offset)].Reset()

        _histograms["signal%s"%(_offset)]= hibd["hIBD%s"%(_offset)].Clone()
        _histograms["signal%s"%(_offset)].SetZTitle('evts/day')
        _histograms["signal%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
        _histograms["signal%s"%(_offset)].SetTitle('Signal - %s offset %2d'%(_energyEstimator,offset))
        _histograms["signal%s"%(_offset)].SetName('hSignal%s'%(_offset))
        _histograms["signal%s"%(_offset)].Reset()

        _histograms["background%s"%(_offset)]= hibd["hIBD%s"%(_offset)].Clone()
        _histograms["background%s"%(_offset)].SetZTitle('evts/day')
        _histograms["background%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
        _histograms["background%s"%(_offset)].SetTitle('total backgrounds - %s offset %2d'%(_energyEstimator,offset))
        _histograms["background%s"%(_offset)].SetName('hBackground%s'%(_offset))
        _histograms["background%s"%(_offset)].Reset() 



        _proc = '_%s'%(_cov)

        _maxSignal,_maxBkgd,_maxSoverB,_maxnx,_max_dtw = -1,-1,-1,-1,-1

        for _d in range(binR-1):
            for _n in range(binN-1):
                _db,_nb,offset=int(_d+1),int(_n+1),int(offset)
                # get the signal rates 
                _p_d  = hibd["hIBD%s"%(_offset)].GetXaxis().GetBinCenter(_db)
                _p_nx = hibd["hIBD%s"%(_offset)].GetYaxis().GetBinCenter(_nb)
                _n_d  = _p_d 
                _n_nx = _p_nx+offset
                _ibd_bin = hibd["hIBD%s"%(_offset)].FindBin(_p_d,_p_nx)
                _ibd_rate     = hibd["hIBD%s"%(_offset)].GetBinContent(_ibd_bin) #rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                _signal = _ibd_rate*86400.*signalEff
                _ibd_error = hibd["hIBD%s"%(_offset)].GetBinError(_db,_nb)
                _signal_error = _ibd_error*86400.*signalEff

                # now get the background rates 
                _p_d_1  = hacc["hAcc%s"%(_offset)].GetXaxis().GetBinCenter(_db)
                _p_nx_1 = hacc["hAcc%s"%(_offset)].GetYaxis().GetBinCenter(_nb)
                _n_d_1  = _p_d_1
                _n_nx_1 = _p_nx_1
                # first get the individual background contributions
                accmap = mapFile.Get("hist%soffset%s_%d_%d"%(_energyEstimator,_offset,dtcut,dscut*1000))
                accEff = accmap.GetBinContent(_db,_nb)
                _acc_rate  = hacc["hAcc%s"%(_offset)].GetBinContent(_db,_nb)*86400*accEff  #accidentals rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                try:
                    _acc_error = sqrt(pow(hacc["hAcc%s"%(_offset)].GetBinError(_db,_nb)/hacc["hAcc%s"%(_offset)].GetBinContent(_db,_nb),2)+pow(accmap.GetBinError(_db,_nb)/accmap.GetBinContent(_db,_nb),2))*_acc_rate
                except:
                    _acc_error = 10.
                    print(_p_d,", ",_p_nx)
                _fn_rate   = hfn["hFN%s"%(_offset)].GetBinContent(_db,_nb)*86400.*signalEff   #fast neutrons rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                _fn_error = hfn["hFN%s"%(_offset)].GetBinError(_db,_nb)*86400.*signalEff
                _rn_rate   = hrn["hRN%s"%(_offset)].GetBinContent(_db,_nb)*86400.*signalEff   #radionuclides rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                _rn_error = hrn["hRN%s"%(_offset)].GetBinError(_db,_nb)*86400.*signalEff
                _ibd_bg_rate   = hibdBG["hIBDBG%s"%(_offset)].GetBinContent(_db,_nb)*86400.*signalEff  #radionuclides rate at fidcut=_db, prompt nxcut = _nb, delayed nxcut = _nb+offset
                _ibd_bg_error = hibdBG["hIBDBG%s"%(_offset)].GetBinError(_db,_nb)*86400.*signalEff
               
                # then get the total background rates
                _background = _acc_rate + _rn_rate + _ibd_bg_rate +_fn_rate # background rate
                _background_error = sqrt(pow(_acc_error,2)+pow(_rn_error,2)+pow(_ibd_bg_error,2)+pow(_fn_error,2))

                # save the rates to the analysis results histograms
                sob = _signal/sqrt(_signal+_background)
                _s_plus_b_error = sqrt(pow(_signal_error,2)+pow(_background_error,2))
                sqrt_s_plus_b_error = 0.5*_s_plus_b_error/(_signal+_background)*sqrt(_signal+_background)
                sob_error = sqrt(pow(_signal_error/_signal,2)+pow(sqrt_s_plus_b_error/sqrt(_signal+_background),2))*sob
                _histograms["sOverB%s"%(_offset)].SetBinContent(_db,_nb,sob)
                _histograms["signal%s"%(_offset)].SetBinContent(_db,_nb,_signal)
                _histograms["background%s"%(_offset)].SetBinContent(_db,_nb,_background)
                _histograms["sOverB%s"%(_offset)].SetBinError(_db,_nb,sob_error)
                _histograms["signal%s"%(_offset)].SetBinError(_db,_nb,_signal_error)
                _histograms["background%s"%(_offset)].SetBinError(_db,_nb,_background_error)

                # find the optimal significance as a function of accidentals
                if sob >_maxSoverB:
                    _maxSoverB = sob
                    _maxSignal = _signal
                    _maxBkgd   = _background
                    _maxOffnx   = _n_nx
                    _maxOff_dtw = _n_d
                    _line = ("Offset:%3d, beta/n: Wall-dist (%4.1f,%4.1f) m, %s cut (%d,%d), ibd rate : %4.2f per day;"\
                    %(offset,_p_d,_n_d,_energyEstimator\
                    ,_p_nx,_n_nx\
                    ,_ibd_rate*86400.),)
                    _line2 =    (" acc. rate: %4.3f per day, rn rate: %4.3f per day (all cuts)"%(_acc_rate,_rn_rate),)
                if sob >_maxSoverB2:
                    _maxSoverB2 = sob
                    _maxSignal2 = _signal
                    _maxBkgd2   = _background
                    _maxOffnx2   = _n_nx
                    _maxOff_dtw2 = _n_d
                    _maxOffset2 = offset
                    _maxAcc = _acc_rate
                    _maxRN = _rn_rate
                    _maxFN = _fn_rate
                    _maxIBDbg = _ibd_bg_rate
                    _maxSOB30daysErr = sob_error*sqrt(30)
                    _maxSOB30days = _signal*30/sqrt(_signal*30+_background*30)
                    _maxSigErr = _signal_error
                    _maxBgErr = _background_error
                    _maxSoBErr = sob_error


        print('Offset:',str(offset).rjust(3,' '),',Found max S/sqrt(S+B)','{0:.4f}'.format(_maxSoverB),',(S,B,%s,dtw):('%(_energyEstimator),'{0:.4f}'.format(_maxSignal),'{0:.4f}'.format(_maxBkgd),_maxOffnx,'{0:.1f}'.format(_maxOff_dtw),')')
        line += (_line + _line2,)


    print('\n\nMore info on the maximal sensitivity found:')
    # print line

    for _l in line:
        for i in range(len(_l)):
            print(_l[i], end=' ')
        print('')

    if _maxSignal2 >0:
        T3SIGMA = 9/_maxSignal2 + 9*_maxBkgd2/_maxSignal2**2

    else:
        T3SIGMA = 1e999

    _res = "%s %4.1f %3d %3d %4.3f +/-%4.5f %4.3f +/-%4.3f %4.7f %4.3f %4.7f %4.7f  %4.3f +/-%4.5f %4.3f +/-%4.5f %4.1f " %(_cover,_maxOff_dtw2,_maxOffnx2,_maxOffnx2-_maxOffset2,_maxSignal2,_maxSigErr,_maxBkgd2,_maxBgErr,_maxAcc,_maxRN,_maxFN,_maxIBDbg,_maxSoverB2,_maxSoBErr,_maxSOB30days,_maxSOB30daysErr,T3SIGMA)
    _strRes = "results_coincidence_DTW_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.txt"%(float(arguments['--vetoThickR']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))
    # _strRes = 'res.txt'
    with open(_strRes,'a') as file:
        file.write(_res+'\n')

    _str = "bonsai_root_files%s/analysis_results_%s_coincidence_%dus_%dmm.root"%(additionalMacStr,_cover,dtcut,dscut*1000)
    print('\n\nWriting histograms to file',_str)
    f_root = TFile(_str,"recreate")
    for offset in offsets_nx:
        _offset = str(offset).replace("-","Minus")
        _histograms["sOverB%s"%(_offset)].Write()
        _histograms["signal%s"%(_offset)].Write()
        _histograms["background%s"%(_offset)].Write()
        hacc["hAcc%s"%(_offset)].Write()
        hfn["hFN%s"%(_offset)].Write()
        hrn["hRN%s"%(_offset)].Write()
        hibd["hIBD%s"%(_offset)].Write()
        hibdBG["hIBDBG%s"%(_offset)].Write()

    f_root.Close()

    print('\n\n')


def readEfficiencyHistogramLassen():

    hist = {}

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    firstGo  = 1 
    print('Reading in root tree')
    ''' Following loop has two purposes. Read in the histogram from results file and scale them with appropriate rates'''
    for idx,_cover in enumerate(coverage):
        if (int(arguments["--QFIT"])==1):
            _strAdd = _cover+"_QFIT"
        else:
            _strAdd = _cover
        _str = "bonsai_root_files%s/results_%s.root"%(additionalMacStr,_strAdd)
        outfile = TFile(_str,"READ")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:
                    _tag = "hist%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    _tag = _tag.replace(" ","")
                    if 'singles' not in _tag:
                      hist[_tag] =  outfile.Get(_tag)
                      #print('Extraction efficiency histogram from file ',_tag,'and applying rate: ',,'evts/s')
                      _str = "%s_%s_%s"%(_element,_loc,_p)
                      _str = _str.replace(" ","")
                      _scale = rates[_str][0]
                      print('Extraction efficiency histogram from file ',_tag,'and applying rate: ',_scale,'evts/s ... ', end = '')

                      try:
                        print(' entries:',hist[_tag].GetEntries(),' ... ', end = '')
                        #hist[_tag].Scale(_scale)
                        if firstGo:
                            h = hist[_tag].Clone()
                            h.SetZTitle('singles rate (Hz)')
                            h.SetTitle('Singles rate')
                            h.SetName('hSinglesRate')
                            h.Sumw2()
                            h.Reset()
                            hn = hist[_tag].Clone()
                            hn.SetZTitle('neutron rate (Hz)')
                            hn.SetTitle('rate')
                            hn.SetName('hNeutronRate')
                            hn.Reset()
                            hee = hist[_tag].Clone()
                            hee.SetZTitle('positron efficiency')
                            hee.SetTitle('efficiency')
                            hee.SetName('hPositronEfficiency')
                            hee.Reset()
                            hee_hs = hist[_tag].Clone()
                            hee_hs.SetZTitle('positron efficiency')
                            hee_hs.SetTitle('efficiency')
                            hee_hs.SetName('hPositronHeyshamSigEfficiency')
                            hee_hs.Reset()
                            hee_hb = hist[_tag].Clone()
                            hee_hb.SetZTitle('positron efficiency')
                            hee_hb.SetTitle('efficiency')
                            hee_hb.SetName('hPositronHeyshamBkgEfficiency')
                            hee_hb.Reset()
                            hne = hist[_tag].Clone()
                            hne.SetZTitle('neutron efficiency')
                            hne.SetTitle('efficiency')
                            hne.SetName('hNeutronEfficiency')
                            hne.Reset()
                            firstGo =0
                        if 'Neutron' in _tag:
                            print('identified as neutron events ...')
                            hne.Add(hist[_tag],1)
                            hist[_tag].Scale(_scale)
                            hn.Add(hist[_tag],1)
                        elif 'Positron' in _tag:
                            print('identified as positron events  ,')
                            if 'HeyshamSig' in _tag:
                                hee_hs.Add(hist[_tag],1)
                            elif 'HeyshamBkg' in _tag:
                                hee_hb.Add(hist[_tag],1)
                            else:
                                hee.Add(hist[_tag],1)
                        elif '238U' in _tag or '232Th' in _tag or '40K' in _tag or '222Rn' in _tag:
                            print('identified and adding to accidentals')
                            hist[_tag].Scale(_scale)
                            h.Add(hist[_tag],1) 
                        else:
                            print('not identified and not used in analsysis.')
			    
                      except:
                        print("Could not find ",_tag,". Skipping entry.")
    '''print(h)
    h.SaveAs('h.C')
    hn.SaveAs('hn.C')
    hne.SaveAs('hne.C')
    hee.SaveAs('hee.C')'''


    print("\nCompleted reading in of histogram, accidental and IBD identifitication. \n\n")
    
    timeWindow = float(arguments["-t"])/1e6 # time in microsecond
    timeAcc = timeWindow*86400.
    accidentalContamination = float(arguments["--acc"])
    signalEff = float(arguments["--se"])
    print("\nPerforming analyis with time window",timeWindow, "seconds (-t option); daily coincidence factor : ", timeAcc, "evts/days.")
    print("Default accidental proximity contamination (--acc option): ",accidentalContamination)
    print("Default signal efficiency after proximity cuts (--se) option",signalEff,"\n")
    scanMax = int(arguments["--scanMaxValue"])
    if  arguments["--positiveScan"]:
        offsets_nx = np.arange(0,scanMax+1,scanStep,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        offsets_nx = np.arange(-scanMax,0+1,dtype=int)## Correction to add last value
    else:
        offsets_nx = np.arange(-scanMax,scanMax+1,scanStep,dtype=int)## Correction to add last value
        #offsets_nx = [-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]  ## bin numbers
    offsets_dtw = [0]       ## bin numbers

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 88,7.5,95.5
    binwidthN = (rangeNmax-rangeNmin)/binN

    _cov = arguments['-C']


    pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])
    pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])
    detectorRadius  = float(arguments['--tankRadius'])-float(arguments['--steelThick'])
    detectorHeight  = float(arguments['--halfHeight'])-float(arguments['--steelThick'])

    _maxSignal,_maxBkgd,_maxSoverB,_maxOffnx,_maxOff_dtw = -1,-1,-1,-1,-1
    _maxSignal2,_maxBkgd2,_maxSoverB2,_maxOffnx2,_maxOff_dtw2,_maxOffset2 = -1,-1,-1,-1,-1,-1

    heysham_signal_rate_ratio = rates['IBDPositronHeyshamSig_LIQUID_ibd_p_hs'][0]/rates['IBDPositron_LIQUID_ibd_p'][0]
    heysham_background_rate_ratio = rates['IBDPositronHeyshamBkg_LIQUID_ibd_p_hb'][0]/rates['IBDPositron_LIQUID_ibd_p'][0]


    line,_line,_line2= ("",),"",""
    _histograms = {}
    for offset in offsets_nx:
        for fv_offset in offsets_dtw:
            _offset = str(offset)
            _offset = _offset.replace("-","Minus")
            _histograms["sOverB_%d"%(offset)]= h.Clone()
            _histograms["sOverB_%d"%(offset)].SetZTitle('signal/sqrt(signal+background)')
            _histograms["sOverB_%d"%(offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
            _histograms["sOverB_%d"%(offset)].SetTitle('hSoB - offset %2d'%(offset))
            _histograms["sOverB_%d"%(offset)].SetName('hSoB%s'%(_offset))
            _histograms["sOverB_%d"%(offset)].Reset()

            _histograms["signal_%d"%(offset)]= h.Clone()
            _histograms["signal_%d"%(offset)].SetZTitle('evts/day')
            _histograms["signal_%d"%(offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
            _histograms["signal_%d"%(offset)].SetTitle('Signal - offset %2d'%(offset))
            _histograms["signal_%d"%(offset)].SetName('hSignal%s'%(_offset))
            _histograms["signal_%d"%(offset)].Reset()

            _histograms["background_%d"%(offset)]= h.Clone()
            _histograms["background_%d"%(offset)].SetZTitle('evts/day')
            _histograms["background_%d"%(offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
            _histograms["background_%d"%(offset)].SetTitle('accidentals - offset %2d'%(offset))
            _histograms["background_%d"%(offset)].SetName('hBackground%s'%(_offset))
            _histograms["background_%d"%(offset)].Reset() 



            _proc = '_%d_%d_%s'%(offset,fv_offset,_cov)

            _maxSignal,_maxBkgd,_maxSoverB,_maxOffnx,_maxOff_dtw = -1,-1,-1,-1,-1
            for _d in range(binR-fv_offset-1):
                for _n in range(binN-offset-1):
                    _db,_nb,offset=int(_d+1),int(_n+1),int(offset)
                    _p_d  = hee.GetXaxis().GetBinCenter(_db)
                    _p_nx = hee.GetYaxis().GetBinCenter(_nb)
                    _n_d  = hne.GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_nx = hne.GetYaxis().GetBinCenter(_nb+offset)
                    fidRadius,fidHeight = pmtRadius-_p_d*1e3,pmtHeight-_p_d*1e3
                    ratioScaling = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
                    _p_v        = hee.GetBinContent(_db,_nb)/ratioScaling
                    _p_v_err    = hee.GetBinError(_db,_nb)/ratioScaling
                    _n_v        = hne.GetBinContent(_db+fv_offset,_nb+offset)/ratioScaling
                    _n_v_error  = hne.GetBinError(_db+fv_offset,_nb+offset)/ratioScaling
                    _rate_v     = hn.GetBinContent(_db+fv_offset,_nb+offset)
                    _rate_v_err = hn.GetBinError(_db+fv_offset,_nb+offset)

                    if arguments['--Heysham']:
                        _p_v = hee_hs.GetBinContent(_db,_nb)/ratioScaling
                        _p_v_err = hee_hs.GetBinError(_db,_nb)/ratioScaling
                        _rate_v = _rate_v * heysham_signal_rate_ratio
                        _rate_v_err = _rate_v_err * heysham_signal_rate_ratio

                    _signal = _rate_v*_p_v*86400.*signalEff

                    ###print(_rate_v, _rate_v_err, _p_v, _p_v_err)
                    if _rate_v == 0:
                        _rate_v_err = 1
                        _rate_v = 1
                    if _p_v == 0:
                        _p_v_err = 1
                        _p_v  = 1

                    _signal_err = sqrt(_rate_v_err*_rate_v_err/(_rate_v*_rate_v) + _p_v_err*_p_v_err/(_p_v*_p_v))*_signal 

                    _p_d_1  = h.GetXaxis().GetBinCenter(_db)
                    _p_nx_1 = h.GetYaxis().GetBinCenter(_nb)
                    _n_d_1  = h.GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_nx_1 = h.GetYaxis().GetBinCenter(_nb+offset)
                    _p_v_1  = h.GetBinContent(_db,_nb)
                    _p_v_1_err = h.GetBinError(_db,_nb)
                    _n_v_1  = h.GetBinContent(_db+fv_offset,_nb+offset)
                    _n_v_1_err  = h.GetBinError(_db+fv_offset,_nb+offset)

                    _background = _p_v_1*_n_v_1*timeAcc*accidentalContamination + _signal*0.30
                    if _p_v_1 == 0:
                        _p_v_1 = 1
                        _p_v_1_err = 1
                    if _n_v_1 == 0:
                        _n_v_1 = 1
                        _n_v_1_err = 1
                    _background_err = sqrt( _p_v_1_err*_p_v_1_err/(_p_v_1*_p_v_1) +  _n_v_1_err*_n_v_1_err/(_n_v_1*_n_v_1))*_background

                    if arguments['--Heysham']:
                        _p_v_bkg = hee_hb.GetBinContent(_db,_nb)/ratioScaling
                        _rate_v_bkg = hn.GetBinContent(_db+fv_offset,_nb+offset)*heysham_background_rate_ratio
                        _background += _rate_v_bkg*_p_v_bkg*86400.*signalEff*0.9



                    sob = _signal/sqrt(_signal+_background)
                    sob_err = sqrt( pow(_signal_err*(2.*_background+_signal)/2.0/pow(_signal+_background,3./2.),2) + pow(_background_err*_signal/2.0/pow(_signal+_background,3./2.),2)  ) * sob

                    _histograms["sOverB_%d"%(offset)].SetBinContent(_db,_nb,sob)
                    _histograms["signal_%d"%(offset)].SetBinContent(_db,_nb,_signal)
                    _histograms["background_%d"%(offset)].SetBinContent(_db,_nb,_background)

                    if sob >_maxSoverB:
                        _maxSoverB = _signal/sqrt(_signal+_background)
                        _maxSignal = _signal
                        _maxSignalErr = _signal_err
                        _maxBkgd   = _background
                        _maxBkgdErr   = _background_err
                        _maxOffnx   = _n_nx
                        _maxOff_dtw = _n_d
                        _line = ("Offset:%3d, beta/n: Wall-dist (%4.1f,%4.1f) m, %s cut (%d,%d), rel. (e+,n) eff (%4.2f,%4.2f), neutron rate :(%4.2f per day), combined rate : %4.2f per day;"\
                        %(offset,_p_d,_n_d\
                        ,_energyEstimator\
                        ,_p_nx,_n_nx\
                        ,_p_v,_n_v\
                        ,_rate_v*86400.\
                        ,_rate_v*_p_v*86400.),)
                        _line2 =    ("acc. rel. (e+,n) rate (%5.3f,%5.3f): acc. combined rate: %4.3f per day (all cuts)"\
                         %(_p_v_1,_n_v_1,_p_v_1*_n_v_1*timeAcc*accidentalContamination),)
                        if arguments['--Heysham']:
                          _line2 =    ("\nacc. rel. (e+,n) rate (%5.3f,%5.3f): acc. combined rate: %4.3f per day (all cuts)"\
                            %(_p_v_1,_n_v_1,_p_v_1*_n_v_1*timeAcc*accidentalContamination),) +\
                            ("\nIBD background (e+, n) eff (%4.2f,%4.2f), neutron rate: %4.2f per day, combined rate: %4.2f per day"\
                                      %(_p_v_bkg, _n_v, _rate_v_bkg*86400*hb_ratio, _rate_v_bkg*_p_v_bkg*86400.),)

                    if sob >_maxSoverB2:
                        _maxSoverB2    = _signal/sqrt(_signal+_background)
                        _maxSoverBErr2 = sob_err
                        _maxSignal2    = _signal
                        _maxSignalErr2 = _signal_err
                        _maxBkgd2      = _background
                        _maxBkgdErr2   = _background_err
                        _maxOffnx2     = _n_nx
                        _maxOff_dtw2   = _n_d
                        _maxOffset2    = offset


            print('Offset:',str(offset).rjust(3,' '),',Found max S/sqrt(S+B)','{0:.4f}'.format(_maxSoverB),',(S,B,%s,dtw):('%(_energyEstimator),'{0:.4f}'.format(_maxSignal),'{0:.4f}'.format(_maxBkgd),_maxOffnx,'{0:.1f}'.format(_maxOff_dtw),')')
            line += (_line + _line2,)


    print('\n\nMore info on the maximal sensitivity found:')
    # print line

    for _l in line:
        for i in range(len(_l)):
            print(_l[i], end=' ')
        print('')

    cut_signal = 0.9
    sigma = 4.65 # 3-sigma 95% of the time, according to Owen report
    _S = _maxSignal2*cut_signal
    _B = _maxSignal2*cut_signal*1.15 + _maxBkgd2 # Includes a 15% other reactor component
    OnOffRatio = 1.5

    if _S >0:
        T3SIGMA = sigma**2*(_B +(_S+_B)/OnOffRatio)/_S/_S

    else:
        T3SIGMA = 1e999
    metric = T3SIGMA*OnOffRatio + T3SIGMA

    if arguments['--Heysham']:
        _B = _maxBkgd2
        if _S > 0:
            T3SIGMA = _B/pow(_S/3,2)
        else:
            T3SIGMA = 1e999
        # Dwell time for median 3 sigma detection given no systematic background uncertainty
        metric = T3SIGMA
    home = os.environ['HOME']
    _res = "%s %4.1f %3d %3d %4.3f %4.3f %4.3f %4.3f %4.1f %4.1f %4.2f %4.2f" %(_cover,_maxOff_dtw2,_maxOffnx2,_maxOffnx2-_maxOffset2,_maxSignal2,_maxSignalErr2,_maxBkgd2,_maxBkgdErr2,T3SIGMA,metric,_maxSoverB2*5.47722,_maxSoverBErr2*5.47722)
    _strRes = "results_DTW_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.txt"%(float(arguments['--vetoThickR']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))
    _strRes2 = "%s/results_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.txt"%(home,float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))

    # _strRes = 'res.txt'
    with open(_strRes,'a') as file:
        file.write(_res+'\n')

    with open(_strRes2,'a') as file:
        file.write(_res+'\n')

    _str = "bonsai_root_files%s/analysis_results_%s.root"%(additionalMacStr,_cover)
    print('\n\nWriting histograms to file',_str)
    f_root = TFile(_str,"recreate")
    h.Write()
    hn.Write()
    hne.Write()
    hee.Write()
    hee_hs.Write()
    hee_hb.Write()
    for offset in offsets_nx:
        _histograms["sOverB_%d"%(offset)].Write()
        _histograms["signal_%d"%(offset)].Write()
        _histograms["background_%d"%(offset)].Write()

    f_root.Close()


    print('\n\n')


def readEfficiencyHistogram():

    hist = {}

    #d,proc,coverage = loadSimulationParametersNew()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    print('Reading in root tree')
    for idx,_cover in enumerate(coverage):
        if (arguments["--QFIT"]==1):
            _strAdd = _cover+"_QFIT"
        else:
            _strAdd = _cover
        _str = "bonsai_root_files%s/results_%s.root"%(additionalMacStr,_strAdd)
        outfile = TFile(_str,"READ")
        for _p in proc:
            for _loc in proc[_p]:
                for _element in d[_p]:
                    _tag = "hist%s_%s_%s_%s"%(_cover,_loc,_element,_p)
                    hist[_tag] =  outfile.Get(_tag)
                    print('Reading file ',_tag)
                    try:
                        print(hist[_tag].GetEntries())
                    except:
                        print("Could not find ",_tag,". Skipping entry.")
                        continue

    print('done processing result file. What is in the directory:')
    gROOT.ProcessLine('.ls')

    print('\nLoading PMT activity:')
    mPMTs,mPMTsU238,mPMTsTh232,mPMTsK40 = loadPMTActivity()
    print('done.')

    print('\nLoading Veto activity:')
    mVETOs,mVETOsU238,mVETOsTh232,mVETOsK40 = loadVETOActivity()
    print('done.')

    print('\nLoading Gd activity:')
    GdU238,GdTh232,GdU235,GdU238_l,GdTh232_l,GdU235_l = loadGdActivity()
    print('done.')

    tankRadius  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])
    tankHeight  = float(arguments["--halfHeight"])-float(arguments['--steelThick'])
    nKiloTons = pi*pow(tankRadius/1000.,2)*(2.*tankHeight/1000.)
    rRn222 = float(arguments["--Rn222"])*nKiloTons
    print('\nLoaded Rn-222 activity of ',rRn222,'Bq per water volume, assumed a rate of %4.3e Bq/m^3'%(float(arguments["--Rn222"])))


    FreeProtons = 0.668559
    TNU         = FreeProtons* nKiloTons /1000./365./24./3600.
    boulbyIBDRate   = 721.173964*TNU
    print('\nLoaded an IBD rate of ',boulbyIBDRate,' events per water volume per second, assumed a rate of %4.3e TNU'%(boulbyIBDRate/TNU))

    innerRad = 12.5 #meters
    outerRad = 13.5 #meters
    rockMass = (pi*pow(outerRad,2)*(2.*outerRad)-pi*pow(innerRad,2)*(2.*innerRad))*power(100.,3)*2.39
    # volumeR         = power(22.,3)-power(20.,3)# Rock cavern e.g. (22m x 22m x 22m) - (20m x 20m x 20m)
    # density         = 2.39 #from McGrath
    # rockMass        = volumeR*power(100.,3)*density
    #Mass of rock evalyated
    avgMuon         = npa([180.,264.])
    avgMuonNC       = power(avgMuon,0.849)
    avgNFluxMag     = 1e-6
    muonRate        = npa([7.06e-7,4.09e-8]) # mu/cm2/s
    tenMeVRatio     = npa([7.51/34.1,1.11/4.86])
    fastNeutrons    = rockMass*avgMuonNC*avgNFluxMag*muonRate*tenMeVRatio
    FN_boulby = fastNeutrons[1]

    avgRNYieldRC    = power(avgMuon,0.73)
    skRNRate        = 0.5e-7 # 1/mu/g cm2
    avgMuonSK       = power(219.,0.73)
    skMuFlux        = 1.58e-7 #mu/cm2/sec
    radionuclideRate= (skRNRate*avgRNYieldRC/avgMuonSK)*muonRate*nKiloTons*1e9
    RN_boulby        = radionuclideRate[1]
    print('\nLoaded mass of rock %e g. Fast Neutron Yield %e per sec; radionuclide yield %e per sec'%(rockMass,FN_boulby,RN_boulby))


    print('\n What are the maximum efficiency/rate found in each histogram:')
    _sing = 0.0
    lineU238PMT,lineTh232PMT,lineKPMT = '','',''
    lineU238VETO,lineTh232VETO,lineKVETO = '','',''
    lineFNROCK = ''
    lineU238GUN,lineTh232GUN,lineKGUN = '','',''
    lineU238ROCK,lineTh232ROCK,lineKROCK = '','',''
    lineU238CONC,lineTh232CONC,lineKCONC = '','',''
    lineRn222WaterVolume = ''
    lineU238GD,lineU235GD,lineTH232GD = '','',''
    linePromptWaterVolume,lineDelayedWaterVolume,linePromptDelayedWaterVolume = '','',''
    lineELSE = ''

    firstGo = 1
    for _t in hist:
        if firstGo:
            h = hist[_t].Clone()
            h.SetZTitle('singles rate (Hz)')
            h.SetTitle('Singles rate')
            h.SetName('hSinglesRate')
            h.Sumw2()
            h.Reset()
            hn = hist[_t].Clone()
            hn.SetZTitle('neutron rate (Hz)')
            hn.SetTitle('rate')
            hn.SetName('hNeutronRate')
            hn.Reset()
            hee = hist[_t].Clone()
            hee.SetZTitle('positron efficiency')
            hee.SetTitle('efficiency')
            hee.SetName('hPositronEfficiency')
            hee.Reset()
            hne = hist[_t].Clone()
            hne.SetZTitle('neutron efficiency')
            hne.SetTitle('efficiency')
            hne.SetName('hNeutronEfficiency')
            hne.Reset()

            firstGo =0
        if 'PMT' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                try:
                    _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                    h.Add(hist[_t],mPMTsU238[0]*0.002)
                    lineU238PMT += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
                except:
                    print("Could not read,",_t)
            else:
                try:
                    _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                    h.Add(hist[_t],mPMTsU238[0])
                    lineU238PMT+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
                except:
                    print("Could not read,",_t)
        elif 'PMT' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            h.Add(hist[_t],mPMTsTh232[0])
            lineTh232PMT += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'PMT' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            h.Add(hist[_t],mPMTsK40[0])
            lineKPMT += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])

        elif 'VETO' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mVETOsU238*0.002
                h.Add(hist[_t],mVETOsU238*0.002)
                lineU238VETO += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mVETOsU238*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mVETOsU238
                h.Add(hist[_t],mVETOsU238)
                lineU238VETO+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mVETOsU238)
        elif 'VETO' in _t and 'CHAIN_232Th_NA' in _t:
            print(mVETOsTh232)
            _sing+=hist[_t].GetMaximum()*mVETOsTh232
            h.Add(hist[_t],mVETOsTh232)
            lineTh232VETO += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mVETOsTh232)
        elif 'VETO' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mVETOsK40
            h.Add(hist[_t],mVETOsK40)
            lineKVETO += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mVETOsK40)

        elif 'GUNITE' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                lineU238GUN += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                lineU238GUN+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'GUNITE' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            lineTh232GUN += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'GUNITE' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            lineKGUN += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])

        elif 'ROCK' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                lineU238ROCK += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                lineU238ROCK+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'ROCK' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            lineTh232ROCK += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'ROCK' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            lineKROCK += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])
        elif 'ROCK' in _t and '_FN' in _t:
            lineFNROCK += "%50s %e %15.10f (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*FN_boulby,hist[_t].GetMaximum()*FN_boulby*3600.*24.)
        elif 'CONC' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]*0.002
                lineU238CONC += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0]*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*mPMTsU238[0]
                lineU238CONC+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsU238[0])
        elif 'CONC' in _t and 'CHAIN_232Th_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsTh232[0]
            lineTh232CONC += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsTh232[0])
        elif 'CONC' in _t and '40K_NA' in _t:
            _sing+=hist[_t].GetMaximum()*mPMTsK40[0]
            lineKCONC += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*mPMTsK40[0])

        elif 'GD' in _t and 'CHAIN_238U_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*GdU238_l*0.002
                h.Add(hist[_t],GdU238_l*0.002)
                lineU238GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU238_l*0.002)
            elif '234Pa' in _t:
                _sing+=hist[_t].GetMaximum()*GdU238
                h.Add(hist[_t],GdU238)
                lineU238GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU238)
            else:
                _sing+=hist[_t].GetMaximum()*GdU238_l
                h.Add(hist[_t],GdU238_l)
                lineU238GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU238_l)
        elif 'GD' in _t and 'CHAIN_232Th_NA' in _t:
            if '228Ac' in _t:
                _sing+=hist[_t].GetMaximum()*GdTh232
                h.Add(hist[_t],GdTh232)
                lineTH232GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdTh232)
            else:
                _sing+=hist[_t].GetMaximum()*GdTh232_l
                h.Add(hist[_t],GdTh232_l)
                lineTH232GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdTh232_l)
        elif 'GD' in _t and 'CHAIN_235U_NA' in _t:
            if '231Th' in _t:
                _sing+=hist[_t].GetMaximum()*GdU235
                h.Add(hist[_t],GdU235)
                lineU235GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU235)
            else:
                _sing+=hist[_t].GetMaximum()*GdU235_l
                h.Add(hist[_t],GdU235_l)
                lineU235GD += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*GdU235_l)


        elif 'WaterVolume' in _t and 'CHAIN_222Rn_NA' in _t:
            if '210Tl' in _t:
                _sing+=hist[_t].GetMaximum()*rRn222*0.002
                h.Add(hist[_t],rRn222*0.002)
                lineRn222WaterVolume += "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*rRn222*0.002)
            else:
                _sing+=hist[_t].GetMaximum()*rRn222
                h.Add(hist[_t],rRn222)
                lineRn222WaterVolume+= "%50s %e %15.10f\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*rRn222)

        elif 'WaterVolume' in _t and 'promptPositron' in _t:
            _day=hist[_t].GetMaximum()*boulbyIBDRate*3600.*24
            hee.Add(hist[_t],1.)#Only add efficiency for positron
            linePromptWaterVolume += "%50s %e %15.10f per sec (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*boulbyIBDRate,_day)

        elif 'WaterVolume' in _t and 'delayedNeutron' in _t:
            _day=hist[_t].GetMaximum()*boulbyIBDRate*3600.*24
            hn.Add(hist[_t],boulbyIBDRate)
            hne.Add(hist[_t],1.)
            linePromptWaterVolume += "%50s %e %15.10f per sec (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*boulbyIBDRate,_day)

        elif 'WaterVolume' in _t and 'promptDelayedPair' in _t:
            _day=hist[_t].GetMaximum()*boulbyIBDRate*3600.*24
            linePromptWaterVolume += "%50s %e %15.10f per sec (%15.10f per day)\n"%(_t,hist[_t].GetMaximum(),hist[_t].GetMaximum()*boulbyIBDRate,_day)

        else:
            try:
                lineELSE += "%50s %e\n"%(_t,hist[_t].GetMaximum())
            except:
                print('Could not find:',_t)

                    # hist =
                    # print _hist.GetMaximum()
    print('')
    print('PMT U-238 \n', lineU238PMT)
    print('PMT Th-232\n', lineTh232PMT)
    print('PMT K\n', lineKPMT,'\n')
    print('VETO U-238 \n', lineU238VETO)
    print('VETO Th-232 \n', lineTh232VETO)
    print('VETO K \n', lineKVETO,'\n')
    print('Water Rn-222\n',lineRn222WaterVolume,'\n')
    print('Gunite U-238 \n', lineU238GUN)
    print('Gunite Th-232\n', lineTh232GUN)
    print('Gunite K\n', lineKGUN,'\n')
    print('Concrete U-238 \n', lineU238CONC)
    print('Concrete Th-232\n', lineTh232CONC)
    print('Concrete K\n', lineKCONC,'\n')
    print('ROCK U-238 \n', lineU238ROCK)
    print('ROCK Th-232\n', lineTh232ROCK)
    print('ROCK K\n', lineKROCK,'\n')
    print('ROCK Fast Neutron\n',lineFNROCK,'\n')
    print('Else  \n', lineELSE,'\n')
    print('Total singles rate:\t\t\t',_sing,'events per sec at minimum buffer distance of 0.5 m\n')
    print('Signal information')
    print('Prompt positron Water volume \n', linePromptWaterVolume)
    signal = ['WaterVolume_delayedNeutron_ibd_n','WaterVolume_promptPositron_ibd_p']
    _str =  "bonsai_root_files%s/%s/histograms_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.root"%(additionalMacStr,_cover,float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))


    timeAcc = 0.0001*86400.

    offsets_nx = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]  ## bin numbers
    offsets_dtw = [0]       ## bin numbers

    binR,rangeRmin,rangeRmax = 31,0.45,3.55
    binwidthR = (rangeRmax-rangeRmin)/binR
    binN,rangeNmin,rangeNmax = 88,7.5,95.5
    binwidthN = (rangeNmax-rangeNmin)/binN

    _cov = arguments['-C']


    pmtRadius = float(arguments['--tankRadius'])-float(arguments['--steelThick'])-float(arguments['--vetoThickR'])
    pmtHeight = float(arguments['--halfHeight'])-float(arguments['--steelThick'])-float(arguments['--vetoThickZ'])
    detectorRadius  = float(arguments['--tankRadius'])-float(arguments['--steelThick'])
    detectorHeight  = float(arguments['--halfHeight'])-float(arguments['--steelThick'])

    _maxSignal,_maxBkgd,_maxSoverB,_maxOffnx,_maxOff_dtw = -1,-1,-1,-1,-1
    _maxSignal2,_maxBkgd2,_maxSoverB2,_maxOffnx2,_maxOff_dtw2,_maxOffset2 = -1,-1,-1,-1,-1,-1


    line,_line,_line2= ("",),"",""
    _histograms = {}
    for offset in offsets_nx:
        for fv_offset in offsets_dtw:
            _histograms["sOverB_%d"%(offset)]= h.Clone()
            _histograms["sOverB_%d"%(offset)].SetZTitle('signal/sqrt(signal+background)')
            _histograms["sOverB_%d"%(offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
            _histograms["sOverB_%d"%(offset)].SetTitle('hSoB - offset %2d'%(offset))
            _histograms["sOverB_%d"%(offset)].SetName('hSoB%d'%(offset))
            _histograms["sOverB_%d"%(offset)].Reset()

            _histograms["signal_%d"%(offset)]= h.Clone()
            _histograms["signal_%d"%(offset)].SetZTitle('evts/day')
            _histograms["signal_%d"%(offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
            _histograms["signal_%d"%(offset)].SetTitle('Signal - offset %2d'%(offset))
            _histograms["signal_%d"%(offset)].SetName('hSignal%d'%(offset))
            _histograms["signal_%d"%(offset)].Reset()

            _histograms["background_%d"%(offset)]= h.Clone()
            _histograms["background_%d"%(offset)].SetZTitle('evts/day')
            _histograms["background_%d"%(offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
            _histograms["background_%d"%(offset)].SetTitle('accidentals - offset %2d'%(offset))
            _histograms["background_%d"%(offset)].SetName('hBackground%d'%(offset))
            _histograms["background_%d"%(offset)].Reset()



            _proc = '_%d_%d_%s'%(offset,fv_offset,_cov)

            _maxSignal,_maxBkgd,_maxSoverB,_maxOffnx,_maxOff_dtw = -1,-1,-1,-1,-1
            for _d in range(binR-fv_offset-1):
                for _n in range(binN-offset-1):
                    _db,_nb=_d+1,_n+1
                    _p_d  = hee.GetXaxis().GetBinCenter(_db)
                    _p_nx = hee.GetYaxis().GetBinCenter(_nb)
                    _n_d  = hne.GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_nx = hne.GetYaxis().GetBinCenter(_nb+offset)
                    fidRadius,fidHeight = pmtRadius-_p_d*1e3,pmtHeight-_p_d*1e3
                    ratioScaling = (pow(fidRadius,2)*fidHeight)/(pow(detectorRadius,2)*detectorHeight)
                    _p_v  = hee.GetBinContent(_db,_nb)/ratioScaling
                    _n_v  = hne.GetBinContent(_db+fv_offset,_nb+offset)/ratioScaling
                    _rate_v  = hn.GetBinContent(_db+fv_offset,_nb+offset)

                    _signal = _rate_v*_p_v*86400.

                    _p_d_1  = h.GetXaxis().GetBinCenter(_db)
                    _p_nx_1 = h.GetYaxis().GetBinCenter(_nb)
                    _n_d_1  = h.GetXaxis().GetBinCenter(_db+fv_offset)
                    _n_nx_1 = h.GetYaxis().GetBinCenter(_nb+offset)
                    _p_v_1  = h.GetBinContent(_db,_nb)
                    _n_v_1  = h.GetBinContent(_db+fv_offset,_nb+offset)

                    _background = _p_v_1*_n_v_1*timeAcc*0.05
                    sob = _signal/sqrt(_signal+_background)

                    _histograms["sOverB_%d"%(offset)].SetBinContent(_db,_nb,sob)
                    _histograms["signal_%d"%(offset)].SetBinContent(_db,_nb,_signal)
                    _histograms["background_%d"%(offset)].SetBinContent(_db,_nb,_background)

                    if sob >_maxSoverB:
                        _maxSoverB = _signal/sqrt(_signal+_background)
                        _maxSignal = _signal
                        _maxBkgd   = _background
                        _maxOffnx   = _n_nx
                        _maxOff_dtw = _n_d
                        _line = ("Offset:%3d, beta/n: Wall-dist (%4.1f,%4.1f) m, %s cut (%d,%d), rel. efficiency (%4.2f,%4.2f), neutron rate :(%4.2f per day), combined eff/rate : %4.2f per day;"\
                            %(offset,_p_d,_n_d\
                            ,_energyEstimator\
                            ,_p_nx,_n_nx\
                            ,_p_v,_n_v\
                            ,_rate_v*86400.\
                            ,_rate_v*_p_v*86400.),)
                        _line2 =    ("acc. rate (%5.3f,%5.3f): acc. combined rate: %4.3f per day (pre-prox)"\
                         %(_p_v_1,_n_v_1,_p_v_1*_n_v_1*timeAcc),)

                    if sob >_maxSoverB2:
                        _maxSoverB2 = _signal/sqrt(_signal+_background)
                        _maxSignal2 = _signal
                        _maxBkgd2   = _background
                        _maxOffnx2   = _n_nx
                        _maxOff_dtw2 = _n_d
                        _maxOffset2 = offset


            print('Offset:',str(offset).rjust(3,' '),',Found max S/sqrt(S+B)',_maxSoverB,',(S,B,%s,dtw):('%(_energyEstimator),_maxSignal,_maxBkgd,_maxOffnx,_maxOff_dtw,')')
            line += (_line + _line2,)

    print('\n\nMore info on the maximal sensitivity found:')
    # print line

    for _l in line:
        for i in range(len(_l)):
            print(_l[i], end=' ')
        print('')

    cut_signal = 0.9
    sigma = 4.65 # 3-sigma 95% of the time, according to Owen report
    _S = _maxSignal2*cut_signal
    _B = _maxSignal2*cut_signal*1.15 + _maxBkgd2 # Includes a 15% other reactor component
    OnOffRatio = 1.5

    if _S >0:
        T3SIGMA = sigma**2*(_B +(_S+_B)/OnOffRatio)/_S/_S

    else:
        T3SIGMA = 1e999
    metric = T3SIGMA*OnOffRatio + T3SIGMA
    _res = "%s %4.1f %3d %3d %4.3f %4.3f %4.1f %4.1f" %(_cover,_maxOff_dtw2,_maxOffnx2,_maxOffnx2-_maxOffset2,_maxSignal2,_maxBkgd2,T3SIGMA,metric)
    _strRes = "results_DTW_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.txt"%(float(arguments['--vetoThickR']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))
    # _strRes = 'res.txt'
    with open(_strRes,'a') as file:
        file.write(_res+'\n')



    print('\n\nWriting histograms to file',_str)
    f_root = TFile(_str,"recreate")
    h.Write()
    hn.Write()
    hne.Write()
    hee.Write()
    for offset in offsets_nx:
        _histograms["sOverB_%d"%(offset)].Write()
        _histograms["signal_%d"%(offset)].Write()
        _histograms["background_%d"%(offset)].Write()
    f_root.Close()

    print('\n\nall done.')



def findRate():


    #d,proc,coverage = loadSimulationParametersNew()
    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)

    print('\nLoading Shotcrete activity:')
    shotmass,act_238u,act_232th,act_40k = loadShotcreteActivity()
    print('done.')

    print('\nLoading ROCK activity:')
    rockmass,act_238u,act_232th,act_40k = loadRockActivity()
    print('done.')

    print('\nLoading CONCRETE activity:')
    concmass,concact_238u,concact_232th,concact_40k = loadConcreteActivity()
    print('done.')

    print('\nLoading TANK activity:')
    tankmass,tankact_60co,tankact_137cs = loadTankActivity()
    print('done.')

    print('\nLoading PSUP activity:')
    psupmass,psupact_60co,psupact_137cs = loadPSUPActivity()
    print('done.')

    print('\nLoading VETO activity:')
    mVETOs,mVETOsU238,mVETOsTh232,mVETOsK40 = loadVETOActivity()
    print('done.')

    print('\nLoading PMT activity:')
    mPMTs,mPMTsU238,mPMTsTh232,mPMTsK40 = loadPMTActivity()
    print('done.')

    tankRadius  = float(arguments["--tankRadius"])-float(arguments['--steelThick'])
    tankHeight  = float(arguments["--halfHeight"])-float(arguments['--steelThick'])
    nKiloTons = pi*pow(tankRadius/1000.,2)*(2.*tankHeight/1000.)
    rRn222 = float(arguments["--Rn222"])*nKiloTons
    print('\nLoaded Rn-222 activity of ',rRn222,'Bq per water volume, assumed a rate of %4.3e Bq/m^3'%(float(arguments["--Rn222"])))

    FreeProtons = 0.668559
    TNU         = FreeProtons* nKiloTons /1000./365./24./3600.
    boulbyIBDRate   = 721.173964*TNU
    print('\nLoaded an IBD rate of ',boulbyIBDRate,' events per water volume per second, assumed a rate of %4.3e TNU'%(boulbyIBDRate/TNU))

    innerRad = 12.5 #meters
    outerRad = 13.5 #meters
    rockMass = (pi*pow(outerRad,2)*(2.*outerRad)-pi*pow(innerRad,2)*(2.*innerRad))*power(100.,3)*2.39
    # volumeR         = power(22.,3)-power(20.,3)# Rock cavern e.g. (22m x 22m x 22m) - (20m x 20m x 20m)
    # density         = 2.39 #from McGrath
    # rockMass        = volumeR*power(100.,3)*density
    #Mass of rock evalyated
    avgMuon         = npa([180.,264.])
    avgMuonNC       = power(avgMuon,0.849)
    avgNFluxMag     = 1e-6
    muonRate        = npa([7.06e-7,4.09e-8]) # mu/cm2/s
    tenMeVRatio     = npa([7.51/34.1,1.11/4.86])
    fastNeutrons    = rockMass*avgMuonNC*avgNFluxMag*muonRate*tenMeVRatio
    FN_boulby = fastNeutrons[1]

    avgRNYieldRC    = power(avgMuon,0.73)
    skRNRate        = 0.5e-7 # 1/mu/g cm2
    avgMuonSK       = power(219.,0.73)
    skMuFlux        = 1.58e-7 #mu/cm2/sec
    radionuclideRate= (skRNRate*avgRNYieldRC/avgMuonSK)*muonRate*nKiloTons*1e9
    RN_boulby        = radionuclideRate[1]
    print('\nLoaded mass of rock %e g. Fast Neutron Yield %e per sec; radionuclide yield %e per sec'%(rockMass,FN_boulby,RN_boulby))
    _strRes = "rate_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM_%s.txt"%(float(arguments['--vetoThickR']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]),arguments["-C"])
    _res = "%e %e %d %e %e %e"%(boulbyIBDRate,rRn222,mPMTs[0],mPMTsU238[0],mPMTsTh232[0],mPMTsK40[0])
    with open(_strRes,'a') as file:
        file.write(_res+'\n')

    return boulbyIBDRate,rRn222,mPMTsU238,mPMTsTh232,mPMTsK40
