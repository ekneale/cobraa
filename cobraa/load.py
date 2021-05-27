import os.path

from stat import S_IRWXG,S_IRWXU
from shutil import rmtree
import warnings

import numpy as np
from numpy import sqrt
from numpy import array as npa
from numpy import power,absolute,logical_and,column_stack,zeros,empty,append,\
sqrt,absolute,recarray


from math import pow,exp,log10,pi

try:
    from rootpy.plotting import Canvas,Hist,Hist2D,Graph
    from rootpy.plotting.style import set_style
    from rootpy.io import root_open


    warnings.simplefilter("ignore")
except:
    print("Historical package rootpy is not loaded. This package is not needed for the core watchmakers tasks.")

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
    --cluster=<_clus>      Specify cluster to use [Default: local]

    ## Create macros and job scripts for a user defined detector configuration

    -m                     Generate rat-pac macro files
    -j                     Create rat-pac/bonsai submision scripts for option above. Can be run with -m.
    --jobTime=<_jt>        Length of job in minutes for LASSEN [Default: 200]
    --energyEst=<_EE>      Default energy estimator (n9,n100,n400,nX) [Default: n9]
    -N=<_N>                Number of runs to simulate [Default: 40]

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
    --rPMT=<_rpmt>         Inner PMT radius in mm [Default: 6700]
    --rU238_IP=<_ruip>     Relative U238 Inner PMTs level [Default: 1.0]
    --rT232_IP=<_rtip>     Relative Th232 Inner PMTs level [Default: 1.0]
    --rK40_IP=<_rkip>      Relative K40 Inner PMTs level [Default: 1.0]

    --lightSim             Option to run a light simulation for testing purposes
    -e=<runBeamEntry>      Number of events to be simulated per macro [Default: 25000]
    
    ## Perform efficiency and sensitivity evaluation (after simulation and reconstruction).

    -M                     Merge result files from trial ntuples. Step one.
    --mergeRATFiles        Merge raw ratds files (off by default)
    --coincidences         Map the efficiencies of events which pass the cuts (analysis step 1)
    --evtype=<_ev>         Set process to evaluate for coincidences
    --sensitivity          Calculate the rates for final optimisation of signal significance (analysis step 2)
    --triggers             Get the number of triggers for singles processes

    ## Define the cuts/ranges over which to optimise

    --minNXprompt=<_minNXp>      Minimum threshold number of direct hits [Default: 9.]
    --maxNXprompt=<_maxNXp>      Maximum threshold number of direct hits [Default: 30.]
    --minNXdelayed=<_minNXd>     Minimum threshold number of direct hits [Default: 9.]
    --maxNXdelayed=<_maxNXd>     Maximum threshold number of direct hits [Default: 30.]
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
    """

try:
    import docopt
    arguments = docopt.docopt(docstring)
    print('\nUsing docopt as the user control interface\n')
except ImportError:
    print('docopt is not a recognized module, it is required to run this module')

if (arguments['--Heysham']):
        print("Using Heysham spectrum (all 4 cores) and assuming Hartlepool is off")

def loadSimulationParametersCoincidence():
    
    location = {
        'LIQUID':['pn_ibd','40K_NA','CHAIN_238U_NA','CHAIN_232Th_NA','CHAIN_235U_NA','A_Z'],\
        'PMT':['40K_NA','CHAIN_238U_NA','CHAIN_232Th_NA'],\
        'ROCK_1':[],\
        'ROCK_2':['40K_NA','CHAIN_238U_NA','CHAIN_232Th_NA','RADIOGENIC','FASTNEUTRONS'],\
        'TANK':[],\
        'PSUP':[],\
        'IBEAM':[],\
        'ALL':[],\
    }

    

def loadSimulationParameters():
    #Chain and subsequent isotopes
    d = {}

    d['CHAIN_238U_NA'] =['234Pa','214Pb','214Bi','210Bi','210Tl']
    d['CHAIN_232Th_NA'] = ['228Ac','212Pb','212Bi','208Tl']#['228Ac','212Bi','208Tl']
    d['CHAIN_235U_NA'] = ['231Th','223Fr','211Pb','211Bi','207Tl']#['231Th','223Fr','207Tl']
    d['40K_NA']         = ['40K']
    d['STEEL_ACTIVITY'] = ['60Co','137Cs']

    d['ibd_p'] = ['IBDPositron']
    d['ibd_p_hs'] = ['IBDPositronHeyshamSig']
    d['ibd_p_hb'] = ['IBDPositronHeyshamBkg']
    d['ibd_n'] = ['IBDNeutron']
    d['pn_ibd']   = ['boulby_geo','core1_hartlepool','core2_hartlepool','boulby_world','heysham_signal','heysham_background']

    d['singles'] = ['singles']

    d['A_Z'] =  ['li 9','n 17']
    d['RADIOGENIC'] = ['rock_neutrons']
    d['FASTNEUTRONS'] = ['fast_neutrons']

    if arguments['--lightSim']:
        ## Define what components are associated with each physical process
        process = {
        'pn_ibd':['LIQUID'],\
        '40K_NA':['LIQUID','PMT'],\
        'CHAIN_238U_NA':['PMT','GD','LIQUID'],\
        'CHAIN_232Th_NA':['PMT','GD','LIQUID'],\
        'CHAIN_235U_NA':['GD'],\
        'A_Z':['LIQUID'],\
        'singles':['SINGLES'],\
        'RADIOGENIC':['ROCK_2'],\
        'FASTNEUTRONS':['ROCK_2']}

    else:
        ## Define what components are associated with each physical process
        process = { 
        'pn_ibd':['LIQUID'],
        '40K_NA':['LIQUID','PMT','PSUP', 'IBEAM','TANK','ROCK_2'],\
        'CHAIN_238U_NA':['PMT','PSUP','IBEAM','TANK','ROCK_2','LIQUID'],\
        'CHAIN_232Th_NA':['PMT','PSUP','IBEAM','TANK','ROCK_2','LIQUID'],\
        'CHAIN_235U_NA':['TANK','PSUP','LIQUID','IBEAM'],\
        'STEEL_ACTIVITY':['TANK','PSUP','IBEAM'],\
        'A_Z':['LIQUID'],\
        'singles':['SINGLES'],\
        'RADIOGENIC':['ROCK_2','ROCK_1'],\
        'FASTNEUTRONS':['ROCK_2']}

    ## First column is the production rate per second of the process, second column is the fractional changes to the event generation.

    uip   = float(arguments["--rU238_IP"])
    tip    = float(arguments["--rT232_IP"])
    kip    = float(arguments["--rK40_IP"])
    pmtVolCorr = 1.

    if arguments['--cylinderSize'] and int(arguments['--cylinderSize'])==12:
        if int(arguments['--cylinderPct'])==10:
            iPMTs = 1128./2238.
            pmtVolCorr = 1352385.722/1347679.758
        elif int(arguments['--cylinderPct'])==15:
            iPMTs = 1616./2238.
            pmtVolCorr = 1350316.794/1347679.758
    else:
        if int(arguments['--cylinderPct'])==10:
            iPMTs = 1614./3248.
#            pmtVolCorr = 782347.4824/779082.9845
            print('not currently correcting for PMT volume')
        elif int(arguments['--cylinderPct'])==15:
            iPMTs = 1128./3248.
            #pmtVolCorr = 780617.7225/779082.9845
            print('not currenty correcting for PMT volume')


    if arguments['--cylinderSize'] and int(arguments['--cylinderSize'])==12:
        print('Using rates for 12m cylinder')
        jobRate = {\
'IBDPositron_LIQUID_ibd_p': [ 2.968e-05*pmtVolCorr , 1], \
'IBDPositronHeyshamSig_LIQUID_ibd_p_hs': [1.977e-06*pmtVolCorr , 1], \
'IBDPositronHeyshamBkg_LIQUID_ibd_p_hb': [3.880e-06 *pmtVolCorr, 1], \
'IBDNeutron_LIQUID_ibd_n': [ 2.968e-05*pmtVolCorr, 1], \
'boulby_geo_LIQUID_pn_ibd': [1.02e-06*pmtVolCorr , 1],\
'core1_hartlepool_LIQUID_pn_ibd': [1.697e-05*pmtVolCorr , 1],\
'core2_hartlepool_LIQUID_pn_ibd': [1.271e-05 *pmtVolCorr, 1],\
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
'231Th_LIQUID_CHAIN_235U_NA': [4.70E-03*pmtVolCorr , 1], \
'223Fr_LIQUID_CHAIN_235U_NA': [4.70E-03*0.0138*pmtVolCorr, 1], \
'211Pb_LIQUID_CHAIN_235U_NA': [4.70E-03*pmtVolCorr , 1], \
'211Bi_LIQUID_CHAIN_235U_NA': [4.70E-03*0.00270*pmtVolCorr , 1], \
'207Tl_LIQUID_CHAIN_235U_NA': [4.70E-03*pmtVolCorr , 1], \
'234Pa_LIQUID_CHAIN_238U_NA': [1.35E+00*pmtVolCorr , 1], \
'214Pb_LIQUID_CHAIN_238U_NA': [1.35E+00*pmtVolCorr , 1], \
'214Bi_LIQUID_CHAIN_238U_NA': [1.35E+00*pmtVolCorr , 1], \
'210Bi_LIQUID_CHAIN_238U_NA': [1.35E+00 *pmtVolCorr, 1], \
'210Tl_LIQUID_CHAIN_238U_NA': [1.35E+00*0.0002*pmtVolCorr , 1], \
'60Co_IBEAM_STEEL_ACTIVITY': [1.45e+04, 50], \
'137Cs_IBEAM_STEEL_ACTIVITY': [1.52e+04, 50], \
'60Co_TANK_STEEL_ACTIVITY': [2.03E+04, 50], \
'137Cs_TANK_STEEL_ACTIVITY': [2.14E+04, 50], \
'60Co_PSUP_STEEL_ACTIVITY': [2.30E+03, 50], \
'137Cs_PSUP_STEEL_ACTIVITY': [2.43E+03, 50], \
'li9_LIQUID_A_Z': [1.705E-06 *pmtVolCorr, 1], \
'n17_LIQUID_A_Z': [1.713E-06 *pmtVolCorr, 1],\
'singles_SINGLES_singles': [1,40],\
'rock_neutrons_ROCK_2_RADIOGENIC': [8.30E+00, 1],\
'rock_neutrons_ROCK_1_RADIOGENIC': [2.04E+02, 1],\
'fast_neutrons_ROCK_2_FASTNEUTRONS': [1.18E-02, 0.5]}
# singles rate for lightSim option with rock neutrons
# NB veto rates are incorrect

    else:
        #print('Using rates for 16m cylinder with 6.7m inner PMT radius')
        jobRate = {\
'core1_hartlepool_LIQUID_pn_ibd': [4.845e-05*pmtVolCorr , 1],\
'core2_hartlepool_LIQUID_pn_ibd': [3.360e-05*pmtVolCorr , 1],\
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
'231Th_LIQUID_CHAIN_235U_NA': [1.48E-02*pmtVolCorr , 1], \
'223Fr_LIQUID_CHAIN_235U_NA': [1.48E-02*0.0138*pmtVolCorr, 1], \
'211Pb_LIQUID_CHAIN_235U_NA': [1.48E-02*pmtVolCorr , 1], \
'211Bi_LIQUID_CHAIN_235U_NA': [1.48E-02*0.00270*pmtVolCorr , 1], \
'207Tl_LIQUID_CHAIN_235U_NA': [1.48E-02*pmtVolCorr , 1], \
'60Co_IBEAM_STEEL_ACTIVITY': [2.40e+04, 50], \
'137Cs_IBEAM_STEEL_ACTIVITY': [2.53e+04, 50], \
'60Co_TANK_STEEL_ACTIVITY': [3.61E+04, 50], \
'137Cs_TANK_STEEL_ACTIVITY': [3.80E+04, 50], \
'60Co_PSUP_STEEL_ACTIVITY': [3.38E+03, 50], \
'137Cs_PSUP_STEEL_ACTIVITY': [3.56E+03, 50], \
'li9_LIQUID_A_Z': [4.051E-06*pmtVolCorr , 1], \
'n17_LIQUID_A_Z': [4.072E-06*pmtVolCorr , 1],\
'singles_SINGLES_singles': [1,40],\
'rock_neutrons_ROCK_2_RADIOGENIC': [1.34e01,1],\
'rock_neutrons_ROCK_1_RADIOGENIC': [3.11e02, 1],\
'fast_neutrons_ROCK_2_FASTNEUTRONS': [1.85e-2, 0.5]}
# singles rate for lightSim option with rock neutrons
# NB veto rates are incorrect


    return d,process,jobRate

