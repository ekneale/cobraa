import os.path

from stat import S_IRWXG,S_IRWXU
from shutil import rmtree
import warnings

import numpy as np
from numpy import power,absolute,logical_and,column_stack,zeros,empty,append,\
sqrt,absolute,recarray
from math import pow,exp,log10,pi


# This reads in flags and arguments from the command line and loads rates for the specified detector design.
# Author Marc Bergevin
# Adapted by Liz Kneale (2021)


docstring = """
    Usage: watchmakers.py [options]

    Arguments:

    Options:

    ## System options

    --force                Forcing the recreation of the root_file,bonsai_root_file and log folders
    -v                     Verbose. Allow print out of additional information.
    --cluster=<_clus>      Specify cluster to use (lassen/sheffield/edinburgh/glasgow,etc) [Default: local]

    ## Create macros and job scripts for a user defined detector configuration

    -m                     Generate rat-pac macro files
    -j                     Create rat-pac/bonsai submision scripts for option above. Can be run with -m.
    --jobTime=<_jt>        Length of job in minutes for LASSEN [Default: 200]
    --energyEst=<_EE>      Default energy estimator (n9,n100,n400,nX) [Default: n9]
    -N=<_N>                Number of runs to simulate [Default: 40]
    --reduced              Simulate a reduced number of processes/decays

    ## Define detector geometry and other features

    --muMetal=<_MM>        Implement muMetal [Default: 0]
    --lightCone=<_LC>      Implement lightConcentrators [Default: 0]
    --depth=<depthD>       Depth of detector (for fast neutron spectra) [Default: 2805.]
    --cylinderSize=<LS>    Height of cylinder (m)  [Default: 16]
    --cylinderPct=<LP>     Photocoverage in cylinder design (percent) [Default: 20]
    --detectMedia=<_dM>    Detector media (doped_water,...) [Default: doped_water]
    --collectionEff=<CE>   Collection efficiency (e.g.: 0.85,0.67,0.475)
    --pmtModel=<_PMTM>     PMT Model (r7081_lqe/r7081_hqe for 10inch or r11780_lqe/r11780_hqe for 12inch)
    --vetoModel=<_vPMTM>   Veto PMT Model (ETEL, r7081_lqe/r7081_hqe for 10inch or r11780_lqe/r11780_hqe for 12inch)
    --photocath =<_PC>     PMT photocathode (R7081HQE)
    --pmtCtrPoint          Point inner PMTs of Watchman geometry to detector center

    # Options to allow scaling of PMT numbers/activity
    --rPMT=<_rpmt>         Inner PMT radius in mm [Default: 5700]
    --rU238_IP=<_ruip>     Relative U238 Inner PMTs level [Default: 1.0]
    --rT232_IP=<_rtip>     Relative Th232 Inner PMTs level [Default: 1.0]
    --rK40_IP=<_rkip>      Relative K40 Inner PMTs level [Default: 1.0]

    --lightSim             Option to run a light simulation for testing purposes
    -e=<runBeamEntry>      Number of events to be simulated per macro [Default: 25000]
    --watchmakers          Option to run the simulation/analysis of individual event types

    ## Perform efficiency and sensitivity evaluation (after simulation and reconstruction).

    -M                     Merge result files from trial ntuples. Step one.
    --mergeRATFiles        Merge raw ratds files (off by default)
    --coincidences         Map the efficiencies of events which pass the cuts (analysis step 1)
    --evtype=<_ev>         Set process to evaluate for coincidences

    --sensitivity          Calculate the rates for final optimisation of signal significance (analysis step 2)
    --triggers             Get the number of triggers for singles processes
    --backgrounds          Plot backgrounds as a function of distance from rPMT

    ## Define the cuts/ranges over which to optimise

    --minNXprompt=<_minNXp>      Minimum threshold number of direct hits [Default: 9.]
    --maxNXprompt=<_maxNXp>      Maximum threshold number of direct hits [Default: 30.]
    --minNXdelayed=<_minNXd>     Minimum threshold number of direct hits [Default: 9.]
    --maxNXdelayed=<_maxNXd>     Maximum threshold number of direct hits [Default: 30.]
    --minFid =<_minFid>          Minimum fiducial cut (distance from PMT radius) [Default: 0.5]
    --maxFid =<_maxFid>          Maximum fiducial cut (distance from PMT radius) [Default: 2.0]
    --binwidthNX=<_binNX>        Bin width for scan over nX cut range [Default: 1]
    --binwidthFid=<_binFid>      Bin width for scan over closest PMT cut range [Default: 0.1]
    --binwidthdT=<binT>          Bin width for scan over dT cut [Default: 10]
    --binwidthdR=<binR>          Bin width for scan over dR cut [Default: 0.1]
    --dTmin=<_dTmin>             Minimum value for dt cut [Default: 100]
    --dTmax=<_dTmax>             Maximum value for dt cut [Default: 100]
    --dRmin=<_dRmin>             Minimum value for dR cut [Default: 2.0]
    --dRmax=<_dRmax>             Maximum value for dR cut [Default: 2.0]
    -g=<goodness>                Bonsai position goodness parameter [Default: 0.1]
    -G=<Goodness>                Bonsai direction goodness parameter [Default: 0.1]
    --se=<_se>                   Default signal efficiency [Default: 1.00]
    --Heysham                    Uses Heysham reactor spectrum in place of Hartlepool

    ## Original Watchmakers flags to be incorporated

    --histograms        Apply nx-position cuts and create efficiency histograms. Step two.
    --procType=<_pt>    Allows to single out a single process (ibd_p,ibd_n,CHAIN_238U_NA,CHAIN_232Th_NA,CHAIN_235U_NA,...)
    --procEle=<_pi>     Allows to single out a single element ()
    --procLoc=<_pl>     Allows to single out a single location for single process (LIQUID, PMT, ...)
    --procNew           Create a new sensitivity result file
    --PMTAccHists       Apply nx/good_pos/good_dir cuts and Look at Acceptance for events in liquid within PMT volume. Optional.
    --evalRate          Evaluate the rates and combine to efficiency histrograms Step three.
     --fnSys=<_fns>        Fast neutron systematics [Default: 0.27]
     --fnScale=<_fnsc>     Fast neutron scaling [Default: 1.0]
     --QFIT=<_QF>          Switch from bonsai to qfit [Default: 0]
     --positiveScan        Only vary the delayed cut as positive offset of the prompt cut
     --negativeScan        Only vary the delayed cut as negative offset of the prompt cut
     --scanMaxValue=<_sm>  Maximal offset value (etiher negative or positive) to scan over [Default: 20]
     --offBin=<_ob>        Binning in offset values [Default: 1]


    """

try:
    import docopt
    arguments = docopt.docopt(docstring)
    print('\nUsing docopt as the user control interface\n')
except ImportError:
    print('docopt is not a recognized module, it is required to run this module')

if (arguments['--Heysham']):
        print("Using Heysham spectrum (all 4 cores) and assuming Hartlepool is off")


def loadSimulationParameters():
    #Chain and subsequent isotopes
    d = {}
    d['CHAIN_238U_NA'] = {'LIQUID':['234Pa','214Pb','214Bi','210Bi','210Tl'],\
                         'PMT':['234Pa','214Pb','214Bi','210Bi','210Tl'],\
                         'TANK':['234Pa','214Pb','214Bi','210Bi','210Tl'],\
                         'ROCK_2':['234Pa','214Pb','214Bi','210Bi','210Tl'],\
                         'IBEAM':['234Pa','214Pb','214Bi','210Bi','210Tl'],\
                         'PSUP':['234Pa','214Pb','214Bi','210Bi','210Tl']}
    d['CHAIN_232Th_NA'] = {'LIQUID':['228Ac','212Pb','212Bi','208Tl'],\
            'PMT':['228Ac','212Pb','212Bi','208Tl'],\
            'TANK':['228Ac','212Pb','212Bi','208Tl'],\
            'ROCK_2':['228Ac','212Pb','212Bi','208Tl'],\
            'IBEAM':['228Ac','212Pb','212Bi','208Tl'],\
            'PSUP':['228Ac','212Pb','212Bi','208Tl']}
    d['CHAIN_235U_NA'] = {'LIQUID':['231Th','223Fr','211Pb','211Bi','207Tl'],\
            'TANK':['231Th','223Fr','211Pb','211Bi','207Tl'],\
            'IBEAM':['231Th','223Fr','211Pb','211Bi','207Tl'],\
            'PSUP':['231Th','223Fr','211Pb','211Bi','207Tl']}

    d['40K_NA'] = {'LIQUID':['40K'],\
            'TANK':['40K'],\
            'IBEAM':['40K'],\
            'PSUP':['40K'],\
            'PMT':['40K'],\
            'ROCK_2':['40K']}
    d['60Co_NA'] = {'PSUP':['60Co'],\
            'TANK':['60Co'],\
            'IBEAM':['60Co']}
    d['137Cs_NA'] = {'PSUP':['137Cs'],\
            'TANK':['137Cs'],\
            'IBEAM':['137Cs']}
    d['ibd_p'] = {'LIQUID':['IBDPositron']}
    d['ibd_p_hs'] = {'LIQUID':['IBDPositronHeyshamSig']}
    d['ibd_p_hb'] = {'LIQUID':['IBDPositronHeyshamBkg']}
    d['ibd_n'] = {'LIQUID':['IBDNeutron']}
    d['pn_ibd'] = {'LIQUID':['boulby_geo','big_hartlepool','small_hartlepool','boulby_world','heysham_signal','heysham_background']}

    d['singles'] = {'ALL':['singles']}
    d['A_Z'] = {'LIQUID':['li 9','n 17']}

    d['FASTNEUTRONS'] = {'ROCK_2':['fast_neutrons']}
 
    d['RADIOGENIC'] = {'ROCK_1':['rock_neutrons']}
    d['RADIOGENIC'] = {'ROCK_2':['rock_neutrons']}
    d['mono'] = {'LIQUID':['e+','e-','gamma']}

   
    if arguments['--reduced']:
        # Define which components are associated with each physical process
        # Removing negligible processes
        # (particularly useful for full coincidence evaluation)
        '''
        PMT             232Th: 208Tl, 212Bi,228Ac;      238U: 210Bi, 210Tl, 214Bi, 234Pa;           40K 
        PSUP            232Th: 208Tl, 212Bi,228Ac;      238U: 210Tl, 214Bi, 234Pa;                  40K;    60Co: TODO (U235/Cs: none)
        IBEAM           232Th: 208Tl;                   238U: 210Tl, 214Bi;                         40K; (U235/Co/Cs: none)
        TANK            232Th: 208Tl, 212Bi;            238U: 210Tl, 214Bi;                         40K;    60Co;
        GD-WATER        232Th: 208Tl, 212Bi,228Ac;      238U: 210Bi, 210Tl, 214Bi, 214Pb, 234Pa;    40K;       U235: TODO
        ROCK (inner)    232Th: 208Tl;                   238U: 210Tl,214Bi;                          Radiogenic neutrons
        ROCK (outer)                                                                                Radiogenic neutrons
        '''

        d.clear()
        d['CHAIN_238U_NA'] = {'LIQUID':['210Bi', '210Tl', '214Bi', '214Pb', '234Pa'],\
                'PMT':['210Bi', '210Tl', '214Bi', '234Pa'],\
                'TANK':['210Tl', '214Bi'],\
                'ROCK_2':['210Tl', '214Bi'],\
                'IBEAM':['210Tl', '214Bi'],\
                'PSUP':['210Tl', '214Bi', '234Pa']}

        d['CHAIN_232Th_NA'] = {'LIQUID':['208Tl', '212Bi','228Ac'],\
                'PSUP':['208Tl', '212Bi','228Ac'],\
                'PMT':['208Tl', '212Bi','228Ac'],\
                'TANK':['208Tl', '212Bi'],\
                'IBEAM':['208Tl'],\
                'ROCK_2':['208Tl']}

        d['CHAIN_235U_NA'] = {'LIQUID':['207Tl','211Pb','223Fr','231Th']}

        d['40K_NA'] = {'LIQUID':['40K'],\
                'PMT':['40K'],\
                'TANK':['40K'],\
                'IBEAM':['40K'],\
                'PSUP':['40K']}

        d['60Co_NA'] = {'PSUP':['60Co'],\
                'TANK':['60Co']}
        d['RADIOGENIC'] = {'ROCK_2':['rock_neutrons']}


        d['pn_ibd'] = {'LIQUID':['boulby_geo','big_hartlepool','small_hartlepool','boulby_world','heysham_signal','heysham_background']}

        d['singles'] = {'ALL':['singles']}
        d['A_Z'] = {'LIQUID':['li 9','n 17']}
        d['FASTNEUTRONS'] = {'ROCK_2':['fast_neutrons']}

    if arguments['--reduced'] and arguments['--lightSim']==0:
        # Define what components are associated with each physical process
        # Reduced number of processes for a full simulation
        # (useful for full coincidence evaluation) 
        process = {
        'CHAIN_238U_NA': ['PMT','PSUP','LIQUID','TANK','IBEAM','ROCK_2'],\
        'CHAIN_232Th_NA':['PMT','PSUP','LIQUID','TANK','IBEAM','ROCK_2'],\
        'CHAIN_235U_NA':['LIQUID'],\
        '40K_NA':        ['PMT','PSUP','LIQUID','TANK','IBEAM'],\
        '60Co_NA':       ['TANK','PSUP'],\
        'RADIOGENIC':    ['ROCK_2'],\
        'pn_ibd':        ['LIQUID'],\
        'A_Z':           ['LIQUID'],\
        'singles':       ['ALL'],\
        'FASTNEUTRONS':  ['ROCK_2']
        }

    if arguments['--reduced'] and arguments['--lightSim']==1:
        # Define what components are associated with each physical process
        # Reduced number of processes for a full simulation
        # (useful for full coincidence evaluation) 
        process = {
        'CHAIN_238U_NA':['PMT','LIQUID'],\
        'CHAIN_232Th_NA':['PMT','LIQUID'],\
        'CHAIN_235U_NA':['LIQUID'],\
        '40K_NA':['LIQUID','PMT'],\
        'RADIOGENIC':['ROCK_2'],\
        'pn_ibd':['LIQUID'],\
        'A_Z':['LIQUID'],\
        'singles':['ALL'],\
        'FASTNEUTRONS':['ROCK_2']}

    elif arguments['--lightSim']:
        # Define what components are associated with each physical process
        # Reduced number of processes for a test simulation 
        # (not for full simulation purposes) 
        process = {
        'CHAIN_238U_NA':['PMT','LIQUID'],\
        'CHAIN_232Th_NA':['PMT','LIQUID'],\
        'CHAIN_235U_NA':['LIQUID'],\
        '40K_NA':['LIQUID','PMT'],\
        'RADIOGENIC':['ROCK_2'],\
        'pn_ibd':['LIQUID'],\
        'ibd_p':['LIQUID'],\
        'ibd_p_hs':['LIQUID'],\
        'ibd_p_hb':['LIQUID'],\
        'ibd_n':['LIQUID'],\
        'A_Z':['LIQUID'],\
        'singles':['ALL'],\
        'FASTNEUTRONS':['ROCK_2']}

                
    else:
        # Define what components are associated with each physical process
        # (all processes included, some may not trigger a detector response)
        process = { 
        'CHAIN_238U_NA':['PMT','PSUP','IBEAM','TANK','ROCK_2','LIQUID'],\
        'CHAIN_232Th_NA':['PMT','PSUP','IBEAM','TANK','ROCK_2','LIQUID'],\
        'CHAIN_235U_NA':['TANK','PSUP','LIQUID','IBEAM'],\
        '40K_NA':['LIQUID','PMT','PSUP', 'IBEAM','TANK','ROCK_2'],\
        '60Co_NA':['TANK','PSUP','IBEAM'],\
        '137Cs_NA':['TANK','PSUP','IBEAM'],\
        'pn_ibd':['LIQUID'],\
        'ibd_p':['LIQUID'],\
        'ibd_p_hs':['LIQUID'],\
        'ibd_p_hb':['LIQUID'],\
        'ibd_n':['LIQUID'],\
        'A_Z':['LIQUID'],\
        'singles':['ALL'],\
        'mono':['LIQUID'],\
        'RADIOGENIC':['ROCK_2'],\
        'FASTNEUTRONS':['ROCK_2']}
        # removed ROCK_1 radiogenic due to strange key error
    print(d,"\n\n\n\n")
    print(process)
    ## First column is the production rate per second of the process, second column is the fractional changes to the event generation.

    uip   = float(arguments["--rU238_IP"])
    tip    = float(arguments["--rT232_IP"])
    kip    = float(arguments["--rK40_IP"])
    pmtVolCorr = 1.

    if int(arguments['--cylinderSize'])==12:
        if int(arguments['--cylinderPct'])==10:
            iPMTs = 1128./2238.
            pmtVolCorr = 1352385.722/1347679.758
        elif int(arguments['--cylinderPct'])==15:
            iPMTs = 1616./2238.
            pmtVolCorr = 1350316.794/1347679.758
    elif int(arguments['--cylinderSize'])==16 and int(arguments['--rPMT'])==5700:
        if int(arguments['--cylinderPct'])==10:
            iPMTs = 1614./3248.
            pmtVolCorr = 3210.15/3203.22
        elif int(arguments['--cylinderPct'])==15:
            iPMTs = 2552./3248.
            pmtVolCorr = 3206.17/3203.22
    elif int(arguments['--cylinderPct'])==16 and int(arguments['--rPMT'])==6700:
        if int(arguments['--cylinderPct'])==10:
            iPMTs = 1232./2230.
            pmtVolCorr = 3211.77/3207.11
            print('not currently correcting for PMT volume')
        elif int(arguments['--cylinderPct'])==15:
            iPMTs = 1808./2230.
            pmtVolCorr = 3209.33/3207.11
            print('not currenty correcting for PMT volume')


    if arguments['--cylinderSize'] and int(arguments['--cylinderSize'])==12:
        print('Using rates for 12m cylinder')
        jobRate = {\
'IBDPositron_LIQUID_ibd_p': [ 2.968e-05*pmtVolCorr , 1], \
'IBDPositronHeyshamSig_LIQUID_ibd_p_hs': [1.977e-06*pmtVolCorr , 1], \
'IBDPositronHeyshamBkg_LIQUID_ibd_p_hb': [3.880e-06 *pmtVolCorr, 1], \
'IBDNeutron_LIQUID_ibd_n': [ 2.968e-05*pmtVolCorr, 1], \
'boulby_geo_LIQUID_pn_ibd': [1.02e-06*pmtVolCorr , 1],\
'big_hartlepool_LIQUID_pn_ibd': [1.697e-05*pmtVolCorr , 1],\
'small_hartlepool_LIQUID_pn_ibd': [1.271e-05 *pmtVolCorr, 1],\
'boulby_world_LIQUID_pn_ibd': [5.857e-06 *pmtVolCorr, 1],\
'heysham_signal_LIQUID_pn_ibd': [1.977e-06 *pmtVolCorr, 1],\
'heysham_background_LIQUID_pn_ibd': [3.880e-06 *pmtVolCorr, 1],\
'40K_LIQUID_40K_NA': [5.39*pmtVolCorr , 1], \
'40K_PMT_40K_NA': [5.08E+03  * kip, 1], \
'40K_VETO_40K_NA': [1.33 , 1], \
'40K_IBEAM_40K_NA': [1.02E+04, 500], \
'40K_PSUP_40K_NA': [1.63E+03 , 1], \
'40K_TANK_40K_NA': [1.44E+04, 50], \
'40K_ROCK_2_40K_NA': [1.33E+04, 1000], \
'234Pa_PMT_CHAIN_238U_NA': [2.49E+03  * uip, 1], \
'214Pb_PMT_CHAIN_238U_NA': [2.49E+03  * uip, 1], \
'214Bi_PMT_CHAIN_238U_NA': [2.49E+03  * uip, 1], \
'210Bi_PMT_CHAIN_238U_NA': [2.49E+03  * uip, 1], \
'210Tl_PMT_CHAIN_238U_NA': [4.98E-01  * uip, 1], \
'234Pa_VETO_CHAIN_238U_NA': [3.33E-01, 1], \
'214Pb_VETO_CHAIN_238U_NA': [3.33E-01, 1], \
'214Bi_VETO_CHAIN_238U_NA': [3.33E-01, 1], \
'210Bi_VETO_CHAIN_238U_NA': [3.33E-01, 1], \
'210Tl_VETO_CHAIN_238U_NA': [3.33E-01*0.0002, 1], \
'234Pa_IBEAM_CHAIN_238U_NA': [2.28E+04, 1], \
'214Pb_IBEAM_CHAIN_238U_NA': [2.28E+04, 1], \
'214Bi_IBEAM_CHAIN_238U_NA': [2.28E+04, 1], \
'210Bi_IBEAM_CHAIN_238U_NA': [2.28E+04, 1], \
'210Tl_IBEAM_CHAIN_238U_NA': [2.28E+04*0.0002, 1], \
'234Pa_PSUP_CHAIN_238U_NA': [3.63E+03, 50], \
'214Pb_PSUP_CHAIN_238U_NA': [3.63E+03, 50], \
'214Bi_PSUP_CHAIN_238U_NA': [3.63E+03, 50], \
'210Bi_PSUP_CHAIN_238U_NA': [3.63E+03, 50], \
'210Tl_PSUP_CHAIN_238U_NA': [3.63E+03*0.0002, 50], \
'234Pa_TANK_CHAIN_238U_NA': [3.20E+04, 50], \
'214Pb_TANK_CHAIN_238U_NA': [3.20E+04, 50], \
'214Bi_TANK_CHAIN_238U_NA': [3.20E+04, 50], \
'210Bi_TANK_CHAIN_238U_NA': [3.20E+04, 50], \
'210Tl_TANK_CHAIN_238U_NA': [3.20E+04*0.0002, 50], \
'234Pa_ROCK_2_CHAIN_238U_NA': [2.76E+06, 1], \
'214Pb_ROCK_2_CHAIN_238U_NA': [2.76E+06, 1], \
'214Bi_ROCK_2_CHAIN_238U_NA': [2.76E+06, 1], \
'210Bi_ROCK_2_CHAIN_238U_NA': [2.76E+06, 1], \
'210Tl_ROCK_2_CHAIN_238U_NA': [2.76E+06*0.0002, 1], \
'234Pa_GD_CHAIN_238U_NA': [1.01E-01*pmtVolCorr , 1], \
'214Pb_GD_CHAIN_238U_NA': [1.01E-01*pmtVolCorr , 1], \
'214Bi_GD_CHAIN_238U_NA': [1.01E-01*pmtVolCorr , 1], \
'210Bi_GD_CHAIN_238U_NA': [1.01E-01*pmtVolCorr , 1], \
'210Tl_GD_CHAIN_238U_NA': [1.01E-01*0.0002*pmtVolCorr , 1], \
'228Ac_PMT_CHAIN_232Th_NA': [2.54E+03  * tip, 1], \
'212Pb_PMT_CHAIN_232Th_NA': [2.54E+03  * tip, 1], \
'212Bi_PMT_CHAIN_232Th_NA': [1.62E+03  * tip, 1], \
'208Tl_PMT_CHAIN_232Th_NA': [9.13E+02  * tip, 1], \
'228Ac_VETO_CHAIN_232Th_NA': [3.32E-02, 1], \
'212Pb_VETO_CHAIN_232Th_NA': [3.32E-02, 1], \
'212Bi_VETO_CHAIN_232Th_NA': [3.32E-02*0.64, 1], \
'208Tl_VETO_CHAIN_232Th_NA': [3.32E-02*0.36, 1], \
'228Ac_IBEAM_CHAIN_232Th_NA': [3.28E+03, 1], \
'212Pb_IBEAM_CHAIN_232Th_NA': [3.28E+03, 1], \
'212Bi_IBEAM_CHAIN_232Th_NA': [3.28E+03*0.64, 1], \
'208Tl_IBEAM_CHAIN_232Th_NA': [3.28E+03*0.36, 1], \
'228Ac_PSUP_CHAIN_232Th_NA': [5.22E+02, 50], \
'212Pb_PSUP_CHAIN_232Th_NA': [5.22E+02, 50], \
'212Bi_PSUP_CHAIN_232Th_NA': [5.22E+02*0.64, 50], \
'208Tl_PSUP_CHAIN_232Th_NA': [5.22E+02*0.36, 50], \
'228Ac_TANK_CHAIN_232Th_NA': [4.61E+03, 50], \
'212Pb_TANK_CHAIN_232Th_NA': [4.61E+03, 50], \
'212Bi_TANK_CHAIN_232Th_NA': [4.61E+03*0.64, 50], \
'208Tl_TANK_CHAIN_232Th_NA': [4.61E+03*0.36, 50], \
'228Ac_ROCK_2_CHAIN_232Th_NA': [4.27E+05, 1000], \
'212Pb_ROCK_2_CHAIN_232Th_NA': [4.27E+05, 1000], \
'212Bi_ROCK_2_CHAIN_232Th_NA': [4.27E+05*0.64, 1000], \
'208Tl_ROCK_2_CHAIN_232Th_NA': [4.27E+05*0.36, 1000], \
'228Ac_GD_CHAIN_232Th_NA': [5.04E-02*pmtVolCorr , 1], \
'212Pb_GD_CHAIN_232Th_NA': [5.04E-02*pmtVolCorr , 1], \
'212Bi_GD_CHAIN_232Th_NA': [5.04E-02*0.64*pmtVolCorr , 1], \
'208Tl_GD_CHAIN_232Th_NA': [5.04E-02*0.36*pmtVolCorr , 1], \
'228Ac_LIQUID_CHAIN_232Th_NA': [1.35E-01*pmtVolCorr , 1], \
'212Pb_LIQUID_CHAIN_232Th_NA': [1.35E-01*pmtVolCorr , 1], \
'212Bi_LIQUID_CHAIN_232Th_NA': [1.35E-01*0.64*pmtVolCorr , 1], \
'208Tl_LIQUID_CHAIN_232Th_NA': [1.35E-01*0.36*pmtVolCorr , 1], \
'231Th_IBEAM_CHAIN_235U_NA': [1.29e+03, 50], \
'223Fr_IBEAM_CHAIN_235U_NA': [1.29e+03*0.0138, 50], \
'211Pb_IBEAM_CHAIN_235U_NA': [1.29e+03, 50], \
'211Bi_IBEAM_CHAIN_235U_NA': [1.29e+03*0.00270, 50], \
'207Tl_IBEAM_CHAIN_235U_NA': [1.29e+03, 50], \
'231Th_PSUP_CHAIN_235U_NA': [2.06E+02, 50], \
'223Fr_PSUP_CHAIN_235U_NA': [2.06E+02*0.0138, 50], \
'211Pb_PSUP_CHAIN_235U_NA': [2.06E+02, 50], \
'211Bi_PSUP_CHAIN_235U_NA': [2.06E+02*0.00270, 50], \
'207Tl_PSUP_CHAIN_235U_NA': [2.06E+02, 50], \
'231Th_TANK_CHAIN_235U_NA': [1.82E+03, 50], \
'223Fr_TANK_CHAIN_235U_NA': [1.82E+03*0.0138, 50], \
'211Pb_TANK_CHAIN_235U_NA': [1.82E+03, 50], \
'211Bi_TANK_CHAIN_235U_NA': [1.82E+03*0.00270, 50], \
'207Tl_TANK_CHAIN_235U_NA': [1.82E+03, 50], \
'231Th_GD_CHAIN_235U_NA': [4.70E-03*pmtVolCorr , 1], \
'223Fr_GD_CHAIN_235U_NA': [4.70E-03*0.0138*pmtVolCorr, 1], \
'211Pb_GD_CHAIN_235U_NA': [4.70E-03*pmtVolCorr , 1], \
'211Bi_GD_CHAIN_235U_NA': [4.70E-03*0.00270*pmtVolCorr , 1], \
'207Tl_GD_CHAIN_235U_NA': [4.70E-03*pmtVolCorr , 1], \
'234Pa_LIQUID_CHAIN_238U_NA': [1.35E+00*pmtVolCorr , 1], \
'214Pb_LIQUID_CHAIN_238U_NA': [1.35E+00*pmtVolCorr , 1], \
'214Bi_LIQUID_CHAIN_238U_NA': [1.35E+00*pmtVolCorr , 1], \
'210Bi_LIQUID_CHAIN_238U_NA': [1.35E+00 *pmtVolCorr, 1], \
'210Tl_LIQUID_CHAIN_238U_NA': [1.35E+00*0.0002*pmtVolCorr , 1], \
'60Co_IBEAM_60Co_NA': [1.45e+04, 50], \
'137Cs_IBEAM_137Cs_NA': [1.52e+04, 50], \
'60Co_TANK_60Co_NA': [2.03E+04, 50], \
'137Cs_TANK_137Cs_NA': [2.14E+04, 50], \
'60Co_PSUP_60Co_NA': [2.30E+03, 50], \
'137Cs_PSUP_137Cs_NA': [2.43E+03, 50], \
'li9_LIQUID_A_Z': [1.705E-06 *pmtVolCorr, 1], \
'n17_LIQUID_A_Z': [1.713E-06 *pmtVolCorr, 1],\
'singles_ALL_singles': [1,1000],\
'rock_neutrons_ROCK_2_RADIOGENIC': [8.30E+00, 1],\
'rock_neutrons_ROCK_1_RADIOGENIC': [2.04E+02, 1],\
'fast_neutrons_ROCK_2_FASTNEUTRONS': [1.18E-02, 0.5],\
'mono_LIQUID_e-':[1,1],\
'mono_LIQUID_e+':[1,1],\
'mono_LIQUID_gamma':[1,1]}
# singles rate for lightSim option with rock neutrons
# NB veto rates are incorrect

    elif arguments['--cylinderSize']==16 and arguments['--rPMT']==6700:
        #print('Using rates for 16m cylinder with 6.7m inner PMT radius')
        jobRate = {\
'big_hartlepool_LIQUID_pn_ibd': [4.845e-05*pmtVolCorr , 1],\
'small_hartlepool_LIQUID_pn_ibd': [3.360e-05*pmtVolCorr , 1],\
'boulby_geo_LIQUID_pn_ibd': [3.565e-07*pmtVolCorr , 1],\
'boulby_world_LIQUID_pn_ibd': [2.227e-06*pmtVolCorr , 1],\
'heysham_signal_LIQUID_pn_ibd': [4.585e-06*pmtVolCorr , 1],\
'heysham_background_LIQUID_pn_ibd': [1.263e-05*pmtVolCorr , 1],\
'40K_LIQUID_40K_NA': [1.28*pmtVolCorr , 1], \
'40K_PMT_40K_NA': [3.58E+04  * kip, 1], \
'40K_VETO_40K_NA': [2.61e+02 * kip, 1], \
'40K_IBEAM_40K_NA': [1.70E+04, 500], \
'40K_PSUP_40K_NA': [2.39E+03 , 1], \
'40K_TANK_40K_NA': [2.55E+04, 50], \
'40K_ROCK_2_40K_NA': [2.14E+04, 1000], \
'234Pa_PMT_CHAIN_238U_NA': [3.62E+03  * uip, 1], \
'214Pb_PMT_CHAIN_238U_NA': [3.62E+03  * uip, 1], \
'214Bi_PMT_CHAIN_238U_NA': [3.62E+03  * uip, 1], \
'210Bi_PMT_CHAIN_238U_NA': [3.62E+03  * uip, 1], \
'210Tl_PMT_CHAIN_238U_NA': [3.62E+03*0.002  * uip, 1], \
'234Pa_VETO_CHAIN_238U_NA': [8.02E+01, 1], \
'214Pb_VETO_CHAIN_238U_NA': [8.02E+01, 1], \
'214Bi_VETO_CHAIN_238U_NA': [8.02E+01, 1], \
'210Bi_VETO_CHAIN_238U_NA': [8.02E+01, 1], \
'210Tl_VETO_CHAIN_238U_NA': [8.02E+01*0.0002, 1], \
'234Pa_IBEAM_CHAIN_238U_NA': [3.78e+04, 1], \
'214Pb_IBEAM_CHAIN_238U_NA': [3.78E+04, 1], \
'214Bi_IBEAM_CHAIN_238U_NA': [3.78E+04, 1], \
'210Bi_IBEAM_CHAIN_238U_NA': [3.78E+04, 1], \
'210Tl_IBEAM_CHAIN_238U_NA': [3.78E+04*0.0002, 1], \
'234Pa_PSUP_CHAIN_238U_NA': [5.32E+03, 50], \
'214Pb_PSUP_CHAIN_238U_NA': [5.32E+03, 50], \
'214Bi_PSUP_CHAIN_238U_NA': [5.32E+03, 50], \
'210Bi_PSUP_CHAIN_238U_NA': [5.32E+03, 50], \
'210Tl_PSUP_CHAIN_238U_NA': [5.32E+03*0.0002, 50], \
'234Pa_TANK_CHAIN_238U_NA': [5.68E+04, 50], \
'214Pb_TANK_CHAIN_238U_NA': [5.68E+04, 50], \
'214Bi_TANK_CHAIN_238U_NA': [5.68E+04, 50], \
'210Bi_TANK_CHAIN_238U_NA': [5.68E+04, 50], \
'210Tl_TANK_CHAIN_238U_NA': [5.68E+04*0.0002, 50], \
'234Pa_ROCK_2_CHAIN_238U_NA': [4.44E+06, 1], \
'214Pb_ROCK_2_CHAIN_238U_NA': [4.44E+06, 1], \
'214Bi_ROCK_2_CHAIN_238U_NA': [4.44E+06, 1], \
'210Bi_ROCK_2_CHAIN_238U_NA': [4.44E+06, 1], \
'210Tl_ROCK_2_CHAIN_238U_NA': [4.44E+06*0.0002, 1], \
'234Pa_GD_CHAIN_238U_NA': [3.17e-01 *pmtVolCorr, 1], \
'214Pb_GD_CHAIN_238U_NA': [3.17e-01 *pmtVolCorr, 1], \
'214Bi_GD_CHAIN_238U_NA': [3.17e-01 *pmtVolCorr, 1], \
'210Bi_GD_CHAIN_238U_NA': [3.17e-01 *pmtVolCorr, 1], \
'210Tl_GD_CHAIN_238U_NA': [3.17e-01*0.0002 *pmtVolCorr, 1], \
'234Pa_LIQUID_CHAIN_238U_NA': [3.20*pmtVolCorr , 1], \
'214Pb_LIQUID_CHAIN_238U_NA': [3.20*pmtVolCorr , 1], \
'214Bi_LIQUID_CHAIN_238U_NA': [3.20*pmtVolCorr , 1], \
'210Bi_LIQUID_CHAIN_238U_NA': [3.20*pmtVolCorr , 1], \
'210Tl_LIQUID_CHAIN_238U_NA': [3.20*0.0002*pmtVolCorr , 1], \
'228Ac_PMT_CHAIN_232Th_NA': [3.17E+03  * tip, 1], \
'212Pb_PMT_CHAIN_232Th_NA': [3.17E+03  * tip, 1], \
'212Bi_PMT_CHAIN_232Th_NA': [3.17E+03*0.64  * tip, 1], \
'208Tl_PMT_CHAIN_232Th_NA': [3.17E+03*0.36  * tip, 1], \
'228Ac_VETO_CHAIN_232Th_NA': [7.04E+01, 1], \
'212Pb_VETO_CHAIN_232Th_NA': [7.04E+01, 1], \
'212Bi_VETO_CHAIN_232Th_NA': [7.04E+01*0.64, 1], \
'208Tl_VETO_CHAIN_232Th_NA': [7.04E+01*0.36, 1], \
'228Ac_IBEAM_CHAIN_232Th_NA': [5.45E+03, 1], \
'212Pb_IBEAM_CHAIN_232Th_NA': [5.45E+03, 1], \
'212Bi_IBEAM_CHAIN_232Th_NA': [5.45E+03*0.64, 1], \
'208Tl_IBEAM_CHAIN_232Th_NA': [5.45E+03*0.36, 1], \
'228Ac_PSUP_CHAIN_232Th_NA': [7.66E+02, 50], \
'212Pb_PSUP_CHAIN_232Th_NA': [7.66E+02, 50], \
'212Bi_PSUP_CHAIN_232Th_NA': [7.66E+02*0.64, 50], \
'208Tl_PSUP_CHAIN_232Th_NA': [7.66E+02*0.36, 50], \
'228Ac_TANK_CHAIN_232Th_NA': [8.18E+03, 50], \
'212Pb_TANK_CHAIN_232Th_NA': [8.18E+03, 50], \
'212Bi_TANK_CHAIN_232Th_NA': [8.18E+03*0.64, 50], \
'208Tl_TANK_CHAIN_232Th_NA': [8.18E+03*0.36, 50], \
'228Ac_ROCK_2_CHAIN_232Th_NA': [6.87E+05, 1000], \
'212Pb_ROCK_2_CHAIN_232Th_NA': [6.87E+05, 1000], \
'212Bi_ROCK_2_CHAIN_232Th_NA': [6.87E+05*0.64, 1000], \
'208Tl_ROCK_2_CHAIN_232Th_NA': [6.87E+05*0.36, 1000], \
'228Ac_GD_CHAIN_232Th_NA': [1.59E-01*pmtVolCorr , 1], \
'212Pb_GD_CHAIN_232Th_NA': [1.59E-01*pmtVolCorr , 1], \
'212Bi_GD_CHAIN_232Th_NA': [1.59E-01*0.64*pmtVolCorr , 1], \
'208Tl_GD_CHAIN_232Th_NA': [1.59E-01*0.36*pmtVolCorr , 1], \
'228Ac_LIQUID_CHAIN_232Th_NA': [3.20E-01*pmtVolCorr , 1], \
'212Pb_LIQUID_CHAIN_232Th_NA': [3.20E-01*pmtVolCorr , 1], \
'212Bi_LIQUID_CHAIN_232Th_NA': [3.20E-01*0.64*pmtVolCorr , 1], \
'208Tl_LIQUID_CHAIN_232Th_NA': [3.20E-01*0.36*pmtVolCorr , 1], \
'231Th_IBEAM_CHAIN_235U_NA': [2.15e+03, 50], \
'223Fr_IBEAM_CHAIN_235U_NA': [2.15e+03*0.0138, 50], \
'211Pb_IBEAM_CHAIN_235U_NA': [2.15e+03, 50], \
'211Bi_IBEAM_CHAIN_235U_NA': [2.15e+03*0.00270, 50], \
'207Tl_IBEAM_CHAIN_235U_NA': [2.15e+03, 50], \
'231Th_PSUP_CHAIN_235U_NA': [3.02E+02, 50], \
'223Fr_PSUP_CHAIN_235U_NA': [3.02E+02*0.0138, 50], \
'211Pb_PSUP_CHAIN_235U_NA': [3.02E+02, 50], \
'211Bi_PSUP_CHAIN_235U_NA': [3.02E+02*0.00270, 50], \
'207Tl_PSUP_CHAIN_235U_NA': [3.02E+02, 50], \
'231Th_TANK_CHAIN_235U_NA': [3.22E+03, 50], \
'223Fr_TANK_CHAIN_235U_NA': [3.22E+03*0.0138, 50], \
'211Pb_TANK_CHAIN_235U_NA': [3.22E+03, 50], \
'211Bi_TANK_CHAIN_235U_NA': [3.22E+03*0.00270, 50], \
'207Tl_TANK_CHAIN_235U_NA': [3.22E+03, 50], \
'231Th_LIQUID_CHAIN_235U_NA': [1.49E-01*pmtVolCorr , 1], \
'223Fr_LIQUID_CHAIN_235U_NA': [1.49E-01*0.0138*pmtVolCorr, 1], \
'211Pb_LIQUID_CHAIN_235U_NA': [1.49E-01*pmtVolCorr , 1], \
'211Bi_LIQUID_CHAIN_235U_NA': [1.49E-01*0.00270*pmtVolCorr , 1], \
'207Tl_LIQUID_CHAIN_235U_NA': [1.49E-01*pmtVolCorr , 1], \
'60Co_IBEAM_60Co_NA': [2.40e+04, 50], \
'137Cs_IBEAM_137Cs_NA': [2.53e+04, 50], \
'60Co_TANK_60Co_NA': [3.61E+04, 50], \
'137Cs_TANK_137Cs_NA': [3.80E+04, 50], \
'60Co_PSUP_60Co_NA': [3.38E+03, 50], \
'137Cs_PSUP_137Cs_NA': [3.56E+03, 50], \
'li9_LIQUID_A_Z': [4.051E-06*pmtVolCorr , 1], \
'n17_LIQUID_A_Z': [4.072E-06*pmtVolCorr , 1],\
'singles_ALL_singles': [1,1000],\
'mono_LIQUID_e-':[1,1],\
'mono_LIQUID_e+':[1,1],\
'mono_LIQUID_gamma':[1,1],\
'rock_neutrons_ROCK_2_RADIOGENIC': [1.34e01,1],\
'rock_neutrons_ROCK_1_RADIOGENIC': [3.11e02, 1],\
'fast_neutrons_ROCK_2_FASTNEUTRONS': [1.85e-2, 0.5]}
# singles rate for lightSim option with rock neutrons
# NB veto rates are incorrect


    else:
        print('Using rates for 16m tank with 5.7m inner PMT radius')
        jobRate = {\
'big_hartlepool_LIQUID_pn_ibd': [4.845e-05*pmtVolCorr , 1],\
'small_hartlepool_LIQUID_pn_ibd': [3.360e-05*pmtVolCorr , 1],\
'boulby_geo_LIQUID_pn_ibd': [3.565e-07*pmtVolCorr , 1],\
'boulby_world_LIQUID_pn_ibd': [2.227e-06*pmtVolCorr , 1],\
'heysham_signal_LIQUID_pn_ibd': [4.585e-06*pmtVolCorr , 1],\
'heysham_background_LIQUID_pn_ibd': [1.263e-05*pmtVolCorr , 1],\
'40K_LIQUID_40K_NA': [1.28*pmtVolCorr , 1], \
'40K_PMT_40K_NA': [8.45E+03  * kip, 1], \
'40K_VETO_40K_NA': [2.61e+02 * kip, 1], \
'40K_IBEAM_40K_NA': [1.70E+04, 500], \
'40K_PSUP_40K_NA': [2.39E+03 , 1], \
'40K_TANK_40K_NA': [2.55E+04, 50], \
'40K_ROCK_2_40K_NA': [2.14E+04, 1000], \
'234Pa_PMT_CHAIN_238U_NA': [2.60E+03  * uip, 1], \
'214Pb_PMT_CHAIN_238U_NA': [2.60E+03  * uip, 1], \
'214Bi_PMT_CHAIN_238U_NA': [2.60E+03  * uip, 1], \
'210Bi_PMT_CHAIN_238U_NA': [2.60E+03  * uip, 1], \
'210Tl_PMT_CHAIN_238U_NA': [2.60E+03*0.002  * uip, 1], \
'234Pa_VETO_CHAIN_238U_NA': [8.02E+01, 1], \
'214Pb_VETO_CHAIN_238U_NA': [8.02E+01, 1], \
'214Bi_VETO_CHAIN_238U_NA': [8.02E+01, 1], \
'210Bi_VETO_CHAIN_238U_NA': [8.02E+01, 1], \
'210Tl_VETO_CHAIN_238U_NA': [8.02E+01*0.0002, 1], \
'234Pa_IBEAM_CHAIN_238U_NA': [3.78e+04, 1], \
'214Pb_IBEAM_CHAIN_238U_NA': [3.78E+04, 1], \
'214Bi_IBEAM_CHAIN_238U_NA': [3.78E+04, 1], \
'210Bi_IBEAM_CHAIN_238U_NA': [3.78E+04, 1], \
'210Tl_IBEAM_CHAIN_238U_NA': [3.78E+04*0.0002, 1], \
'234Pa_PSUP_CHAIN_238U_NA': [5.32E+03, 50], \
'214Pb_PSUP_CHAIN_238U_NA': [5.32E+03, 50], \
'214Bi_PSUP_CHAIN_238U_NA': [5.32E+03, 50], \
'210Bi_PSUP_CHAIN_238U_NA': [5.32E+03, 50], \
'210Tl_PSUP_CHAIN_238U_NA': [5.32E+03*0.0002, 50], \
'234Pa_TANK_CHAIN_238U_NA': [5.68E+04, 50], \
'214Pb_TANK_CHAIN_238U_NA': [5.68E+04, 50], \
'214Bi_TANK_CHAIN_238U_NA': [5.68E+04, 50], \
'210Bi_TANK_CHAIN_238U_NA': [5.68E+04, 50], \
'210Tl_TANK_CHAIN_238U_NA': [5.68E+04*0.0002, 50], \
'234Pa_ROCK_2_CHAIN_238U_NA': [4.44E+06, 1], \
'214Pb_ROCK_2_CHAIN_238U_NA': [4.44E+06, 1], \
'214Bi_ROCK_2_CHAIN_238U_NA': [4.44E+06, 1], \
'210Bi_ROCK_2_CHAIN_238U_NA': [4.44E+06, 1], \
'210Tl_ROCK_2_CHAIN_238U_NA': [4.44E+06*0.0002, 1], \
'234Pa_GD_CHAIN_238U_NA': [3.17e-01 *pmtVolCorr, 1], \
'214Pb_GD_CHAIN_238U_NA': [3.17e-01 *pmtVolCorr, 1], \
'214Bi_GD_CHAIN_238U_NA': [3.17e-01 *pmtVolCorr, 1], \
'210Bi_GD_CHAIN_238U_NA': [3.17e-01 *pmtVolCorr, 1], \
'210Tl_GD_CHAIN_238U_NA': [3.17e-01*0.0002 *pmtVolCorr, 1], \
'234Pa_LIQUID_CHAIN_238U_NA': [3.20*pmtVolCorr , 1], \
'214Pb_LIQUID_CHAIN_238U_NA': [3.20*pmtVolCorr , 1], \
'214Bi_LIQUID_CHAIN_238U_NA': [3.20*pmtVolCorr , 1], \
'210Bi_LIQUID_CHAIN_238U_NA': [3.20*pmtVolCorr , 1], \
'210Tl_LIQUID_CHAIN_238U_NA': [3.20*0.0002*pmtVolCorr , 1], \
'228Ac_PMT_CHAIN_232Th_NA': [2.28E+03  * tip, 1], \
'212Pb_PMT_CHAIN_232Th_NA': [2.28E+03  * tip, 1], \
'212Bi_PMT_CHAIN_232Th_NA': [2.28E+03*0.64  * tip, 1], \
'208Tl_PMT_CHAIN_232Th_NA': [2.28E+03*0.36  * tip, 1], \
'228Ac_VETO_CHAIN_232Th_NA': [7.04E+01, 1], \
'212Pb_VETO_CHAIN_232Th_NA': [7.04E+01, 1], \
'212Bi_VETO_CHAIN_232Th_NA': [7.04E+01*0.64, 1], \
'208Tl_VETO_CHAIN_232Th_NA': [7.04E+01*0.36, 1], \
'228Ac_IBEAM_CHAIN_232Th_NA': [5.45E+03, 1], \
'212Pb_IBEAM_CHAIN_232Th_NA': [5.45E+03, 1], \
'212Bi_IBEAM_CHAIN_232Th_NA': [5.45E+03*0.64, 1], \
'208Tl_IBEAM_CHAIN_232Th_NA': [5.45E+03*0.36, 1], \
'228Ac_PSUP_CHAIN_232Th_NA': [7.66E+02, 50], \
'212Pb_PSUP_CHAIN_232Th_NA': [7.66E+02, 50], \
'212Bi_PSUP_CHAIN_232Th_NA': [7.66E+02*0.64, 50], \
'208Tl_PSUP_CHAIN_232Th_NA': [7.66E+02*0.36, 50], \
'228Ac_TANK_CHAIN_232Th_NA': [8.18E+03, 50], \
'212Pb_TANK_CHAIN_232Th_NA': [8.18E+03, 50], \
'212Bi_TANK_CHAIN_232Th_NA': [8.18E+03*0.64, 50], \
'208Tl_TANK_CHAIN_232Th_NA': [8.18E+03*0.36, 50], \
'228Ac_ROCK_2_CHAIN_232Th_NA': [6.87E+05, 1000], \
'212Pb_ROCK_2_CHAIN_232Th_NA': [6.87E+05, 1000], \
'212Bi_ROCK_2_CHAIN_232Th_NA': [6.87E+05*0.64, 1000], \
'208Tl_ROCK_2_CHAIN_232Th_NA': [6.87E+05*0.36, 1000], \
'228Ac_GD_CHAIN_232Th_NA': [1.59E-01*pmtVolCorr , 1], \
'212Pb_GD_CHAIN_232Th_NA': [1.59E-01*pmtVolCorr , 1], \
'212Bi_GD_CHAIN_232Th_NA': [1.59E-01*0.64*pmtVolCorr , 1], \
'208Tl_GD_CHAIN_232Th_NA': [1.59E-01*0.36*pmtVolCorr , 1], \
'228Ac_LIQUID_CHAIN_232Th_NA': [3.20E-01*pmtVolCorr , 1], \
'212Pb_LIQUID_CHAIN_232Th_NA': [3.20E-01*pmtVolCorr , 1], \
'212Bi_LIQUID_CHAIN_232Th_NA': [3.20E-01*0.64*pmtVolCorr , 1], \
'208Tl_LIQUID_CHAIN_232Th_NA': [3.20E-01*0.36*pmtVolCorr , 1], \
'231Th_IBEAM_CHAIN_235U_NA': [2.15e+03, 50], \
'223Fr_IBEAM_CHAIN_235U_NA': [2.15e+03*0.0138, 50], \
'211Pb_IBEAM_CHAIN_235U_NA': [2.15e+03, 50], \
'211Bi_IBEAM_CHAIN_235U_NA': [2.15e+03*0.00270, 50], \
'207Tl_IBEAM_CHAIN_235U_NA': [2.15e+03, 50], \
'231Th_PSUP_CHAIN_235U_NA': [3.02E+02, 50], \
'223Fr_PSUP_CHAIN_235U_NA': [3.02E+02*0.0138, 50], \
'211Pb_PSUP_CHAIN_235U_NA': [3.02E+02, 50], \
'211Bi_PSUP_CHAIN_235U_NA': [3.02E+02*0.00270, 50], \
'207Tl_PSUP_CHAIN_235U_NA': [3.02E+02, 50], \
'231Th_TANK_CHAIN_235U_NA': [3.22E+03, 50], \
'223Fr_TANK_CHAIN_235U_NA': [3.22E+03*0.0138, 50], \
'211Pb_TANK_CHAIN_235U_NA': [3.22E+03, 50], \
'211Bi_TANK_CHAIN_235U_NA': [3.22E+03*0.00270, 50], \
'207Tl_TANK_CHAIN_235U_NA': [3.22E+03, 50], \
'231Th_LIQUID_CHAIN_235U_NA': [1.49E-01*pmtVolCorr , 1], \
'223Fr_LIQUID_CHAIN_235U_NA': [1.49E-01*0.0138*pmtVolCorr, 1], \
'211Pb_LIQUID_CHAIN_235U_NA': [1.49E-01*pmtVolCorr , 1], \
'211Bi_LIQUID_CHAIN_235U_NA': [1.49E-01*0.00270*pmtVolCorr , 1], \
'207Tl_LIQUID_CHAIN_235U_NA': [1.49E-01*pmtVolCorr , 1], \
'60Co_IBEAM_60Co_NA': [2.40e+04, 50], \
'137Cs_IBEAM_137Cs_NA': [2.53e+04, 50], \
'60Co_TANK_60Co_NA': [3.61E+04, 50], \
'137Cs_TANK_137Cs_NA': [3.80E+04, 50], \
'60Co_PSUP_60Co_NA': [3.38E+03, 50], \
'137Cs_PSUP_137Cs_NA': [3.56E+03, 50], \
'li9_LIQUID_A_Z': [4.051E-06*pmtVolCorr , 1], \
'n17_LIQUID_A_Z': [4.072E-06*pmtVolCorr , 1],\
'singles_ALL_singles': [1,1000],\
'mono_LIQUID_e-':[1,1],\
'mono_LIQUID_e+':[1,1],\
'mono_LIQUID_gamma':[1,1],\
'rock_neutrons_ROCK_2_RADIOGENIC': [1.34e01,1],\
'rock_neutrons_ROCK_1_RADIOGENIC': [3.11e02, 1],\
'fast_neutrons_ROCK_2_FASTNEUTRONS': [1.85e-2, 0.5]}

    return d,process,jobRate

