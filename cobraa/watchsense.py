#import watchmakers as PR

from .load import *
from .globals import *
from .singles import *

from decimal import *
setcontext(ExtendedContext)

from numpy import max



def readEfficiencyHistogramLassen():

    hist = {}

    additionalString,additionalCommands,additionalMacStr,additionalMacOpt = testEnabledCondition(arguments)
    _rateHart = 0.
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
                    hist[_tag] =  outfile.Get(_tag)
                    #print('Extraction efficiency histogram from file ',_tag,'and applying rate: ',,'evts/s')
                    _str = "%s_%s_%s"%(_element,_loc,_p)
                    _str = _str.replace(" ","")
                    _scale = rates[_str][0]
                    if 'FASTNEUTRONS' in _str:
                        _fnrate = _scale
                    if 'A_Z' in _str and 'li9' in _str:
                        _li9rate = _scale
                    if 'A_Z' in _str and 'n17' in _str:
                        _n17rate = _scale
                    if 'hartlepool' in _str:
                        _rateHart += _scale
                    if 'boulby_geo' in _str:
                        _rateGeo = _scale
                    if 'boulby_world' in _str:
                        _rateWorld = _scale
                        

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
                        elif '238U' in _tag or '232Th' in _tag or '40K' in _tag or '222Rn' in _tag or 'rock_neutrons' in _tag:
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
    offBin  = int(arguments["--offBin"])
    if  arguments["--positiveScan"]:
        offsets_nx = np.arange(0,scanMax+1,offBin,dtype=int)## Correction to add last value
    elif  arguments["--negativeScan"]:
        offsets_nx = np.arange(-scanMax,0+1,offBin,dtype=int)## Correction to add last value
    else:
        offsets_nx = np.arange(-scanMax,scanMax+1,offBin,dtype=int)## Correction to add last value
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
    #_strFN = "bonsai_root_files%s/fastneutron_results_%s.root"%(additionalMacStr,_cover)
    #fn_root = TFile(_strFN,"recreate")

    _fnrateDaily    = _fnrate*3600.*24.
    _li9rateDaily   = _li9rate*3600.*24.
    _n17rateDaily   = _n17rate*3600.*24.
    _rateHartDaily  = _rateHart*3600.*24.
    _rateWorldDaily = _rateWorld*3600.*24.
    _rateGeoDaily   = _rateGeo*3600.*24.

    med = arguments["--detectMedia"]
    if arguments["--cylinderSize"]:
        config = arguments["--cylinderSize"]
        conPct = arguments["--cylinderPct"]
    elif arguments["--letterboxSize"]:
        config = arguments["--letterboxSize"]
        conPct = arguments["--letterboxPct"]
    else:
        config = "unknown "
        conPct = "unknown "
    _scal  = 180.0
    print("""\\begin{sidewaystable}
\\centering""")
    fnsys   = arguments["--fnSys"]
    fnscale = arguments["--fnScale"]
    print("\\"+f"caption{{Events for configuration {config}m/{conPct}pct {med} in a six month period; Fast-neutron rate: {_fnrate} per sec or {'{0:.2f}'.format(_fnrateDaily)} per day using a systematic uncertainty of {fnsys} and a rate correction of {fnscale}. Li-9 rate of {_li9rate} per sec or {'{0:.2f}'.format(_li9rateDaily)} per day. N-17 rate of {_n17rate} per sec or {'{0:.2f}'.format(_n17rateDaily)} per day.}}")
    print("\\begin{tabular}{l c c c c c c c c c c c}")
    print("\\hline")
    print(f"offset & ibds-1  & ibds-2 & world-neutrino & geo-neutrino & accidental & fast-neutron & Li-9 & N-17 & S/$\sqrt{{S+B}}$ & S/$\sqrt{{S+B+\sigma_{{sys}}^{2}}}$ & S/$\sqrt{{B+\sigma_{{sys}}^{2}}}$ \\\\ \\hline")
    for offset in offsets_nx:
        for fv_offset in offsets_dtw:
            _offset = str(offset)
            _offset = _offset.replace("-","Minus")
            _histograms["sOverB_%s"%(_offset)]= h.Clone()
            _histograms["sOverB_%s"%(_offset)].SetZTitle('signal/sqrt(signal+background)')
            _histograms["sOverB_%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
            _histograms["sOverB_%s"%(_offset)].SetTitle('hSoB - offset %2d'%(offset))
            _histograms["sOverB_%s"%(_offset)].SetName('hSoB%s'%(_offset))
            _histograms["sOverB_%s"%(_offset)].Reset()

            _histograms["signal_%s"%(_offset)]= h.Clone()
            _histograms["signal_%s"%(_offset)].SetZTitle('evts/day')
            _histograms["signal_%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
            _histograms["signal_%s"%(_offset)].SetTitle('Signal - offset %2d'%(offset))
            _histograms["signal_%s"%(_offset)].SetName('hSignal%s'%(_offset))
            _histograms["signal_%s"%(_offset)].Reset()

            _histograms["background_%s"%(_offset)]= h.Clone()
            _histograms["background_%s"%(_offset)].SetZTitle('evts/day')
            _histograms["background_%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
            _histograms["background_%s"%(_offset)].SetTitle('accidentals - offset %2d'%(offset))
            _histograms["background_%s"%(_offset)].SetName('hBackground%s'%(_offset))
            _histograms["background_%s"%(_offset)].Reset() 
            
            _str_FNtest =  f"histWatchman_ROCK_2_fast_neutrons_FASTNEUTRONS_{_offset}"
            try:
                _histograms["fastneutron_%s"%(_offset)]= outfile.Get( _str_FNtest )
                _histograms["fastneutron_%s"%(_offset)].SetZTitle('evts/day')
                _histograms["fastneutron_%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
                _histograms["fastneutron_%s"%(_offset)].SetTitle('fast-neutron - offset %2d'%(offset))
                _histograms["fastneutron_%s"%(_offset)].SetName('hBackgroundFN%s'%(_offset))
            except:
                print('Could not  extract:',_str_FNtest)
                _histograms["fastneutron_%s"%(_offset)]=h.Clone()
                _histograms["fastneutron_%s"%(_offset)].Reset()

            _str_Li9test =  f"histWatchman_LIQUID_li9_A_Z_{_offset}"
            try:
                _histograms["Li9_%s"%(_offset)]= outfile.Get( _str_Li9test )
                _histograms["Li9_%s"%(_offset)].SetZTitle('evts/day')
                _histograms["Li9_%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
                _histograms["Li9_%s"%(_offset)].SetTitle('Li-9 - offset %2d'%(offset))
                _histograms["Li9_%s"%(_offset)].SetName('hBackgroundLi9%s'%(_offset))
            except:
                print('Could not  extract:',_str_Li9test)
                _histograms["Li9_%s"%(_offset)]=h.Clone()
                _histograms["Li9_%s"%(_offset)].Reset()

            _str_N17test =  f"histWatchman_LIQUID_n17_A_Z_{_offset}"

            try:
                _histograms["N17_%s"%(_offset)]= outfile.Get( _str_N17test )
                _histograms["N17_%s"%(_offset)].SetZTitle('evts/day')
                _histograms["N17_%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
                _histograms["N17_%s"%(_offset)].SetTitle('N17 - offset %2d'%(offset))
                _histograms["N17_%s"%(_offset)].SetName('hBackgroundN17%s'%(_offset))
            except:
                print('Could not  extract:',_str_N17test)
                _histograms["N17_%s"%(_offset)]=h.Clone()
                _histograms["N17_%s"%(_offset)].Reset()


            _str_big_hartlepool_test   =  f"histWatchman_LIQUID_big_hartlepool_pn_ibd_{_offset}"
            _str_boulby_geo_test       =  f"histWatchman_LIQUID_boulby_geo_pn_ibd_{_offset}"
            _str_boulby_world_test     =  f"histWatchman_LIQUID_boulby_world_pn_ibd_{_offset}"
            _str_heysham_background_test =  f"histWatchman_LIQUID_heysham_background_pn_ibd_{_offset}"
            _str_heysham_signal_test   =  f"histWatchman_LIQUID_heysham_signal_pn_ibd_{_offset}"
            _str_heysham_signal_test   =  f"histWatchman_LIQUID_small_hartlepool_pn_ibd_{_offset}"
            
            try:
                _histograms["BH_%s"%(_offset)]= outfile.Get( _str_big_hartlepool_test )
                _histograms["BH_%s"%(_offset)].SetZTitle('evts/day')
                _histograms["BH_%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
                _histograms["BH_%s"%(_offset)].SetTitle('Big Hartelepool - offset %2d'%(offset))
                _histograms["BH_%s"%(_offset)].SetName('hBackgroundBH%s'%(_offset))
            except:
                print('Could not  extract:',_str_big_hartlepool_test)
                _histograms["BH_%s"%(_offset)]=h.Clone()
                _histograms["BH_%s"%(_offset)].Reset()


            try:
                _histograms["GEO_%s"%(_offset)]= outfile.Get(_str_boulby_geo_test)
                _histograms["GEO_%s"%(_offset)].SetZTitle('evts/day')
                _histograms["GEO_%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
                _histograms["GEO_%s"%(_offset)].SetTitle('Geoneutrinos - offset %2d'%(offset))
                _histograms["GEO_%s"%(_offset)].SetName('hBackgroundGeo%s'%(_offset))
            except:
                print('Could not  extract:',_str_boulby_geo_test)
                _histograms["GEO_%s"%(_offset)]=h.Clone()
                _histograms["GEO_%s"%(_offset)].Reset()

            try:
                _histograms["WORLD_%s"%(_offset)]= outfile.Get(_str_boulby_world_test)
                _histograms["WORLD_%s"%(_offset)].SetZTitle('evts/day')
                _histograms["WORLD_%s"%(_offset)].SetYTitle('%s cut on prompt'%(_energyEstimator))
                _histograms["WORLD_%s"%(_offset)].SetTitle('World neutrinos - offset %2d'%(offset))
                _histograms["WORLD_%s"%(_offset)].SetName('hBackgroundWorld%s'%(_offset))
            except:
                print('Could not  extract:',_str_boulby_world_test)
                _histograms["WORLD_%s"%(_offset)]=h.Clone()
                _histograms["WORLD_%s"%(_offset)].Reset()





            _proc = '_%d_%d_%s'%(offset,fv_offset,_cov)

            #print(f"offset & ibds  & accidental & fast-neutron & S/$\sqrt{{S+B}}$ \\\\")

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

                    _fn         = _histograms["fastneutron_%s"%(_offset)].GetBinContent(_db+fv_offset,_nb+offset)*_fnrateDaily *  float(arguments["--fnScale"])
                    _fn_err     = _histograms["fastneutron_%s"%(_offset)].GetBinError(_db+fv_offset,_nb+offset)*_fnrateDaily  *  float(arguments["--fnScale"])                
                    _li9        = _histograms["Li9_%s"%(_offset)].GetBinContent(_db+fv_offset,_nb+offset)*_li9rateDaily
                    _li9_err    = _histograms["Li9_%s"%(_offset)].GetBinError(_db+fv_offset,_nb+offset)*_li9rateDaily
                    _n17        = _histograms["N17_%s"%(_offset)].GetBinContent(_db+fv_offset,_nb+offset)*_n17rateDaily
                    _n17_err    = _histograms["N17_%s"%(_offset)].GetBinError(_db+fv_offset,_nb+offset)*_n17rateDaily
                    _geo        = _histograms["GEO_%s"%(_offset)].GetBinContent(_db+fv_offset,_nb+offset)*_rateGeoDaily
                    _geo_err    = _histograms["GEO_%s"%(_offset)].GetBinError(_db+fv_offset,_nb+offset)*_rateGeoDaily
                    _world        = _histograms["WORLD_%s"%(_offset)].GetBinContent(_db+fv_offset,_nb+offset)*_rateWorldDaily
                    _world_err    = _histograms["WORLD_%s"%(_offset)].GetBinError(_db+fv_offset,_nb+offset)*_rateWorldDaily

                    _bh        = _histograms["BH_%s"%(_offset)].GetBinContent(_db+fv_offset,_nb+offset)*_rateHartDaily
                    _bh_err    = _histograms["BH_%s"%(_offset)].GetBinError(_db+fv_offset,_nb+offset)*_rateHartDaily                    
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

                    _background = _p_v_1*_n_v_1*timeAcc*accidentalContamination
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


                    #sob = _signal/sqrt(_signal+_background + _fn +_li9 + _n17 + _world + _geo)
                    sob = _bh / sqrt(_background + _fn + _fn*float(arguments["--fnSys"])*_fn*float(arguments["--fnSys"]) +   _li9 + _n17 + _world + _geo)
                    sob_err = sqrt( pow(_signal_err*(2.*_background+_signal)/2.0/pow(_signal+_background,3./2.),2) + pow(_background_err*_signal/2.0/pow(_signal+_background,3./2.),2)  ) * sob

                    _histograms["sOverB_%s"%(_offset)].SetBinContent(_db,_nb,sob)
                    _histograms["signal_%s"%(_offset)].SetBinContent(_db,_nb,_signal)
                    _histograms["background_%s"%(_offset)].SetBinContent(_db,_nb,_background)

                    if sob >_maxSoverB:
                        _maxSoverB         = sob
                        _maxSoverBErr      = sob_err
                        _maxSoverB6m       = sob*sqrt(_scal)
                        _maxSoverBErr6m    = sob_err*sqrt(_scal)

                        _maxSignal         = _signal*_scal
                        _maxSignalErr      = _signal_err*_scal
                        _maxBkgd           = _background*_scal
                        _maxBkgdErr        = _background_err*_scal
                        _maxOffnx          = _n_nx
                        _maxOff_dtw        = _n_d
                        _maxFN             = _fn*_scal
                        _maxFNErr          = _fn_err*_scal 
                        _maxFNSys          = _fn*float(arguments["--fnSys"])*_scal 
                        _maxLi9            = _li9*_scal
                        _maxLi9Err         = _li9_err*_scal
                        _maxN17            = _n17*_scal
                        _maxN17Err         = _n17_err*_scal
                        _maxGeo            = _geo*_scal
                        _maxGeoErr         = _geo_err*_scal
                        _maxWorld          = _world*_scal
                        _maxWorldErr       = _world_err*_scal

                        _maxBH             = _bh*_scal
                        _maxBHErr          = _bh_err*_scal 
                        
                        _maxSoverB6m       = _maxBH/sqrt(_maxBH + _maxBkgd + _maxFN + _maxLi9 + _maxN17 + _maxGeo + _maxWorld)
                        _maxSoverBErr6m    = _maxSoverB6m*sqrt(_scal)
                        _maxSoverB6mSys    = _maxBH/sqrt(_maxBH + _maxBkgd + _maxFN + _maxFNSys*_maxFNSys + _maxLi9Err + _maxN17Err + _maxGeoErr + _maxWorldErr  )
                        _maxSoverBErr6mSys = _maxSoverB6mSys*sqrt(_scal)
                        _maxOverB6mSys     = _maxBH/sqrt(              _maxBkgd + _maxFN + _maxFNSys*_maxFNSys + _maxLi9Err + _maxN17Err + _maxGeoErr + _maxWorldErr  )
                        _maxOverBErr6mSys  = _maxOverB6mSys*sqrt(_scal)



                        _line = ("Offset:%3d, beta/n: Wall-dist (%4.1f,%4.1f) m, %s cut (%d,%d), rel. (e+,n) eff (%4.2f,%4.2f), neutron rate :(%4.2f per day), combined rate : %4.2f per day;"\
                            %(offset,_p_d,_n_d\
                            ,_energyEstimator\
                            ,_p_nx,_n_nx\
                            ,_p_v,_n_v\
                            ,_rate_v*86400.\
                            ,_rate_v*_p_v*86400.),)
                        _line2 =    ("acc. rel. (e+,n) rate (%5.3f,%5.3f): acc. combined rate: %4.3f per day (all cuts)"\
                         %(_p_v_1,_n_v_1,_p_v_1*_n_v_1*timeAcc*accidentalContamination),)
                        _line3  = f"{offset}/{'{0:.1f}'.format(_p_d)}/{'{0:.0f}'.format(_p_nx)}/{'{0:.0f}'.format(_n_nx)} & {'{0:.1f}'.format(_maxSignal)} $\pm$ {'{0:.1f}'.format(_maxSignalErr)} & {'{0:.1f}'.format(_maxBH)} $\pm$ {'{0:.1f}'.format(_maxBHErr)} & {'{0:.1f}'.format(_maxWorld)} $\pm$ {'{0:.1f}'.format(_maxWorldErr)} &  {'{0:.1f}'.format(_maxGeo)} $\pm$ {'{0:.1f}'.format(_maxGeoErr)} & {'{0:.1f}'.format(_maxBkgd)} $\pm$ {'{0:.1f}'.format(_maxBkgdErr)} & {'{0:.1f}'.format(_maxFN)} $\pm$ {'{0:.1f}'.format(_maxFNErr)}  $\pm$ {'{0:.1f}'.format(_maxFNSys)} & {'{0:.1e}'.format(_maxLi9)} $\pm$ {'{0:.1e}'.format(_maxLi9Err)} & {'{0:.1e}'.format(_maxN17)} $\pm$ {'{0:.1e}'.format(_maxN17Err)} &{'{0:.2f}'.format(_maxSoverB6m)}  & {'{0:.2f}'.format(_maxSoverB6mSys)}  & {'{0:.2f}'.format(_maxOverB6mSys)}    \\\\"
                        #_lineMax  = _line3
                        if arguments['--Heysham']:
                          _line2 =    ("\nacc. rel. (e+,n) rate (%5.3f,%5.3f): acc. combined rate: %4.3f per day (all cuts)"\
                            %(_p_v_1,_n_v_1,_p_v_1*_n_v_1*timeAcc*accidentalContamination),) +\
                            ("\nIBD background (e+, n) eff (%4.2f,%4.2f), neutron rate: %4.2f per day, combined rate: %4.2f per day"\
                                      %(_p_v_bkg, _n_v, _rate_v_bkg*86400*hb_ratio, _rate_v_bkg*_p_v_bkg*86400.),)

                        if _maxOverB6mSys >_maxSoverB2:
                            _maxSoverB2    = _maxOverB6mSys
                            _maxSoverBErr2 = sob_err
                            _maxSignal2    = _signal
                            _maxSignalErr2 = _signal_err
                            _maxBkgd2      = _background
                            _maxBkgdErr2   = _background_err
                            _maxOffnx2     = _n_nx
                            _maxOff_dtw2   = _n_d
                            _maxOffset2    = offset
                            _maxFN2      = _fn
                            _maxFNErr2  = _fn_err
                            _maxLi9     = _li9
                            _maxLi9Err  = _li9_err
                            _lineMax    = f"{config}m_{conPct}pct_{med}_{fnsys}_{fnscale} & " + _line3
                            #f"{offset} & {'{0:.1f}'.format(_maxSignal)} $\pm$ {'{0:.1f}'.format(_maxSignalErr)} & {'{0:.1f}'.format(_maxBH)} $\pm$ {'{0:.1f}'.format(_maxBHErr)} & {'{0:.1f}'.format(_maxWorld)} $\pm$ {'{0:.1f}'.format(_maxWorldErr)} &  {'{0:.1f}'.format(_maxGeo)} $\pm$ {'{0:.1f}'.format(_maxGeoErr)} & {'{0:.1f}'.format(_maxBkgd)} $\pm$ {'{0:.1f}'.format(_maxBkgdErr)} & {'{0:.1f}'.format(_maxFN)} $\pm$ {'{0:.1f}'.format(_maxFNErr)}  $\pm$ {'{0:.1f}'.format(_maxFNSys)} & {'{0:.1e}'.format(_maxLi9)} $\pm$ {'{0:.1e}'.format(_maxLi9Err)} & {'{0:.1e}'.format(_maxN17)} $\pm$ {'{0:.1e}'.format(_maxN17Err)} &{'{0:.2f}'.format(_maxSoverB6m)}  & {'{0:.2f}'.format(_maxSoverB6mSys)}  & {'{0:.2f}'.format(_maxOverB6mSys)}    \\\\"
            #print('Offset:',str(offset).rjust(3,' '),',Found max S/sqrt(S+B)','{0:.4f}'.format(_maxSoverB),',(S,B,%s,dtw):('%(_energyEstimator),'{0:.4f}'.format(_maxSignal),'{0:.4f}'.format(_maxBkgd),_maxOffnx,'{0:.1f}'.format(_maxOff_dtw),')')
            line += (_line + _line2,)
            print(_line3)
    print("\\hline") 
    print("\\end{tabular}")
    print("\\end{sidewaystable}")
    # print line

    for _l in line:
        for i in range(len(_l)):
            print(_l[i], end=' ')
        print('')


    print('Max sob:\n',_lineMax)
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
    _res = "%s %4.1f %3d %3d %4.3f %4.3f %4.3f %4.3f %4.1f %4.1f %4.2f %4.2f" %(_cover,_maxOff_dtw2,_maxOffnx2,_maxOffnx2-_maxOffset2,_maxSignal2,_maxSignalErr2,_maxBkgd2,_maxBkgdErr2,_maxFN2,_maxFNErr2,_maxSoverB2*13.416,_maxSoverBErr2*13.416)
    _strRes = "results_DTW_%dmm_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.txt"%(float(arguments['--vetoThickR']),float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))
    _strRes2 = "%s/results_U238_%4.3fPPM_Th232_%4.3fPPM_K_%4.3fPPM.txt"%(home,float(arguments["--U238_PPM"]),float(arguments["--Th232_PPM"]),float(arguments["--K_PPM"]))

    # _strRes = 'res.txt'
    with open(_strRes,'a') as file:
        file.write(_res+'\n')

    with open(_strRes2,'a') as file:
        #file.write(_res+'\n')
        file.write(_lineMax+'\n')

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
        _offset = str(offset)
        _offset = _offset.replace("-","Minus")
        _histograms["sOverB_%s"%(_offset)].Write()
        _histograms["signal_%s"%(_offset)].Write()
        _histograms["background_%s"%(_offset)].Write()
        _histograms["fastneutron_%s"%(_offset)].Write()
    f_root.Close()


    print('\n\n')


