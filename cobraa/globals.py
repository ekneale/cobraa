from .load import *
from .style import setStyle

# Define global variables and settings

warnings.simplefilter("ignore")

#Define style settings

Tol_bright = ['#4477AA','#66CCEE','#228833','#CCBB44','#EE6677','#AA3377','#BBBBBB','k']
setStyle()

def drange(start, stop, step):
    rii= start
    while rii<stop:
        yield rii
        rii+=step

d,proc,rates = loadSimulationParameters()
singlespersec = 0
for _p in proc:
    for _loc in proc[_p]:
        for _element in d[_p][_loc]:
            if 'NA' in _p or 'RADIOGENIC' in _p:
                singlespersec+=rates['%s_%s_%s'%(_element,_loc,_p)][0]

nruns = int(arguments['-N'])
nsetSingles = 100
detectorStr = f"Watchman_rightcylinder_{arguments['--cylinderSize']}m_{arguments['--cylinderSize']}m_{arguments['--cylinderPct']}pct_{arguments['--rPMT']}mm"


# detector dimensions

if arguments["--cylinderSize"]=='12':
    detectorRadius = 6000.
    detectorHeight = 6000.
    pmtRadius      = 5465.
    pmtHeight      = 5465.
else:
    detectorRadius = 8000.
    detectorHeight = 8000.
    pmtRadius      = 6700.
    pmtHeight      = 6700.

### cut paramaters
posGood = float(arguments['-g'])
dirGood = float(arguments['-G'])
energyEstimator = arguments['--energyEst']

### search parameters
# energy (n9, n100, etc)
minNXprompt   = float(arguments["--minNXprompt"])
maxNXprompt   = float(arguments["--maxNXprompt"])
minNXdelayed  = float(arguments["--minNXdelayed"])
maxNXdelayed  = float(arguments["--maxNXdelayed"])
binwidthNX    = float(arguments["--binwidthNX"])
rangeNXpmin,rangeNXpmax = minNXprompt-binwidthNX/2.,maxNXprompt+binwidthNX/2.
rangeNXdmin,rangeNXdmax = minNXdelayed-binwidthNX/2.,maxNXdelayed+binwidthNX/2.
binNX = int((rangeNXpmax-rangeNXpmin)/binwidthNX)

# coincidence
binwidthdT = float(arguments["--binwidthdT"])
binwidthdR = float(arguments["--binwidthdR"])
dTmin,dTmax = float(arguments["--dTmin"]),float(arguments["--dTmax"])
dRmin,dRmax = float(arguments["--dRmin"]),float(arguments["--dRmax"])
rangedTmin,rangedTmax = dTmin-binwidthdT/2.,dTmax+binwidthdT/2.
rangedRmin,rangedRmax = dRmin-binwidthdR/2.,dRmax+binwidthdR/2.


# fiducial volume
binwidthFid = float(arguments["--binwidthFid"])
minFid = float(arguments["--minFid"])
maxFid = float(arguments["--maxFid"])
rangeFidmin,rangeFidmax = minFid-binwidthFid/2.,maxFid+binwidthFid/2.
binFid = int(round((rangeFidmax-rangeFidmin)/binwidthFid))


def testEnabledCondition(arguments):

    # defines directory/file naming protocol and macro commands

    additionalString      = ""
    additionalMacOpt      = ""

    # file naming
    if arguments['--cylinderSize']:
       additionalString    = "_%sm_%spct"%(arguments['--cylinderSize'],arguments['--cylinderPct'])

    if (arguments['--lightSim']):
        additionalString += "_lightSim"

    if (arguments['--rPMT']):
        additionalString += "_rPMT_%smm"%(arguments['--rPMT'])

    # additional macro commands
    if (arguments['--detectMedia']):
        additionalMacOpt +="/rat/db/set GEO[detector_veto1] material \"%s\"\n"%(arguments['--detectMedia'])
        additionalMacOpt += "/rat/db/set GEO[detector_veto2] material \"%s\"\n"%(arguments['--detectMedia'])
        additionalMacOpt += "/rat/db/set GEO[detector_target_gb] material \"%s\"\n"%(arguments['--detectMedia'])
        additionalMacOpt += "/rat/db/set GEO[detector_target_fv] material \"%s\"\n"%(arguments['--detectMedia'])
        additionalString += "_detectorMedia_%s" %(arguments['--detectMedia'])

    if (arguments['--collectionEff']):
        additionalMacOpt += "/rat/db/set GEO[inner_pmts] efficiency_correction %f\n" %(float(arguments['--collectionEff']))
        additionalString += "_collectionEfficiency_%f" %(float(arguments['--collectionEff']))

    if (arguments['--pmtModel']):
        additionalMacOpt += "/rat/db/set GEO[inner_pmts] pmt_model \"%s\"\n" %((arguments['--pmtModel']))
        additionalString += "_pmtModel_%s" %((arguments['--pmtModel']))

    if (arguments['--vetoModel']):
        additionalMacOpt += "/rat/db/set GEO[veto_pmts] pmt_model \"%s\"\n" %((arguments['--vetoModel']))
        additionalString += "_vetoModel_%s" %((arguments['--vetoModel']))

    if (int(arguments['--muMetal'])==1):
        additionalMacOpt += "/rat/db/set GEO[inner_pmts] mu_metal %s\n" %((arguments['--muMetal']))
        additionalString += "_muMetal_%s" %((arguments['--muMetal']))

    if (int(arguments['--lightCone'])==1):
        additionalMacOpt += "/rat/db/set GEO[inner_pmts] light_cone %s\n" %((arguments['--lightCone']))
        additionalString += "_lightCone_%s" %((arguments['--lightCone']))

    if arguments['--pmtCtrPoint']:
        additionalMacOpt += '/rat/db/set GEO[inner_pmts] orientation "point"\n'
        additionalMacOpt += '/rat/db/set GEO[shield] orientation_inner "point"\n'
        additionalString += "_pmtCtrPoint_"

    return  additionalString,additionalMacOpt

additionalString,additionalMacOpt = testEnabledCondition(arguments)


