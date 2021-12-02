**Cobraa Coincident-background reactor antineutrino analysis** 

Coincidence analysis for reactor discovery

Adapted from the WATCHMAN Watchmakers package (Author: Marc Bergevin, LLNL).

Toolkit which handles full-detector simulation through to sensitivity analysis
for the WATCHMAN project.

Constitutes a significant overhaul of the original Watchmakers code. Retains the convenience of the directory organisation and macro/job production but incorporates streamlining and a full analysis update.

To be used alongside the FRED (BONSAI) reconstruction or CoRe (BONSAI) pair reconstruction.


Analysis performs full evaluation of coincidences for both signal and background and 
optimises sensitivity to a reactor as a function of cuts in up to 7 dimensions: 

prompt energy threshold

delayed energy threshold

fiducial volume

dT time between triggers

dR distance between triggers (FRED only)

Fit quality (BONSAI timing goodness)

Maximum prompt energy


To install:

    $git clone https://github.com/ekneale/cobraa

    $cd cobraa

    $./configure

    $source env_wm.sh


**Basic use:**

1. create macros and simulation jobs

     ```cobraa -m -j ```

2. run jobs

3. run CoRe reconstruction

4. Map coincidences after cuts

     ```cobraa --coincidences```

5. Calculate rates and optimise signal significance

     ```cobraa --sensitivity```


**Macro generation**

RAT-PAC macros are generated for the following elements of the rat simulation:

1. Detector: detector*.mac sets the detector geometry. Relevant flags are:

```
    --muMetal 1               Turn on the mu metal shields (default is off)

    --lightCone 1             Turn on the light cones (default is off)

    --cylinderSize 22         Set the tank diameter/height to 22m (default is 16m)

    --cylinderPct 10          Set the photocoverage to 10% (default is 20%)

    --cylinderPct 15          Set the photocoverage to 15% (default is 20%)
    
    --rPMT 9000               Set the inner PMT radius to 6.7m (default is 5.7m)

    --detectMedia [index]     Set the detector medium (default is doped_water)
```
   *specify index as per rat-pac/data/Watchman/OPTICS_Watchman.ratdb*
```
    --collectionEff           Set the PMT collection efficiency to a different value (default is 0.9)

    --pmtModel [index]        Set the PMT model (default is r7081_lqe)
```
   *specify index as per rat-pac/data/PMT.ratdb*
```
    --vetoModel [index]       Set the veto PMT model (default is r7081_lqe)

    --photocath               Set the photocathode surface (default is R7081HQE)
```
   *specify index as per rat-pac/data/Watchman/OPTICS_Watchman.ratdb*
```
    --pmtCtrPoint             Point inner PMTs towards centre of detector
```
2. Processors: process.mac sets the daq and features for the log output

3. Physics: phys*.mac sets the rat-pac generator for a given physics process

4. Geo: geo*.mac sets the location in the detector where a given physics process takes place

5. Rates: rates*.mac sets the rate of the given physics process for the number of days being simulated. Relevant flags:
```
   -e                       Sets the number of events to be simulated per macro (default is 25000 for full simulation)
  ```
   For testing purposes, it is recommended to run with a smaller number of events e.g. -e 2500
   
   The number of events simulated per macro run for each event type is then:

   e * rate factor

   where rate factor is a reduction/increase applied where we can justify simulating fewer/more than the number of events in that time period (e.g. where there are many events or very high trigger efficiency).

Macros are found in the cobraa/mac directory.

**Job generation**

Cobraa produces jobs/scripts for the following types of events by default:

Hartlepool cores (called 'hartlepool_1' and 'hartlepool_2')

Heysham cores - ('heysham_full' 4-core signal and 'heysham_2' 2-core signal)

Torness cores ('torness_full')

More distant reactor cores ('gravelines_full','sizewell_B','hinkley_C') for risk mitigation study

Reactor background from more distant reactors at Boulby (boulby_worldbg).

Geoneutrino background at Boulby (boulby_geo)

Fast neutron background

Radionuclide background

Accidental background (singles)

Relevant flags are:
```
    --lightSimWater           Runs simulations for only the decays which contribute significantly to accidentals in water
    
    --lightSimWbLS            Runs simulations for only the decays which contribute significantly to accidentals in WbLS (TODO)
    
    --cluster                 Used to specify the job submission script header options (default is to run locally)
                              Options: lassen, sheffield, edinburgh, warwick so far.

    -N                        Sets the number of times to run each script (default is 40)
```
Outputs are:

1. Scripts: script*.sh scripts to run the rat simulation for each type of event.

   Sources environment variables; runs rat with detector + process + geo + phys + rates macros.


2. Jobs: job*.sh jobs to either run the scripts locally (default) or on Lassen/Sheff clusters
   
   Defines cluster-dependent settings if required
  
   Runs the script 40 times
   
   For singles events, a number of jobs equivalent to ```nsetSingles``` (currently fixed in globals.py) are created. This is to make it possible to run enough statistics at the same time as keeping run times within cluster limits. Should really be increased. 



***Reconstruction***

RAT-PAC output must be passed to the FRED (https://github.com/AIT-WATCHMAN/FRED) or CoRe (currently at https://github.com/ekneale/CoRe -contact eskneale1@sheffield.ac.uk for access) reconstruction prior to the next stages in the process.
   


***Statistics evaluation***

This is an additional function which is useful to use after the reconstruction stage for ensuring that you have sufficient statistics.

To run: ```cobraa --triggers [detector options]```

Reads in the fred output for all event types available and finds:

 1. Number of events simulated

 2. Length of time simulated

 3. Trigger efficiency (trigger = interaction with PMT hits > 6 in 800 microseconds)

 4. Reconstruction efficiency (event passes minimal post-reconstruction cuts: inner PMT hits > 10, n9 (unscattered light) hits > 9, reconstructed vertex > 0.5m from inner PMT radius)

The information is saved to a triggerdata.txt in latex formatting (if using in LaTex, please use ```\usepackage{siunitx}``` in the preamble). Where a root file does not exist, the details are written in simsmissing.txt. Where an event is found to contribute significantly to singles rates, the details are writting in simsrequired.txt.

***Backgrounds plot***

This is another additional function which is useful for comparing the relative contributions of the components to the radioactive backgrounds. Plots the singles rates as a function of fiducial cut. NB this option only works where the backgrounds have been simulated individually.

To run: ```cobraa --backgrounds [detector options] [fiducial cut options]```

**Analysis**

This is the post-reconstruction stage. 

The analysis ultimately optimises reactor sensitivity in 4 dimensions:

```nx_p>minNX``` threshold number of hits n in x ns detected from the prompt event in the pair
```nx_d>minNX``` threshold number of hits n in x ns detected from the delayed event in the pair
```nx_p<Epmax``` maximum number of hits n in x ns detected from the prompt event in the pair
```closestPMT``` perpendicular distance to the nearest PMT, as measured from the PMT radius ```rPMT``` (negative-valued for events occurring in the buffer/veto region) - cut applied to both prompt and delayed event
```dT``` time between the two interactions in a pair
```dR``` distance between the two interactions in a pair (FRED only)
```g``` BONSAI goodness

There are two steps to the analysis:
 1. Evaluate the number of coincidences and scale to per day rate ```--coincidences```
 2. Evaluate the sensitivity based on the coincidence rates ```--sensitivity```

Relevant flags for the analysis steps define the ranges over which to optimise (_optimisation options_):

```
--core                                 run using output from CoRe pair reconstruction

--minNXprompt/--maxNXprompt            define nx_p range

--minNXdelayed/--maxNXdelayed          define nx_d range

--dTmin/--dTmax                        define dT range (us) - by default (130,160)

--dRmin/--dRmax                        define dR range (m) - by default (1.8,2.2)

--minFid/--maxFid                      define fiducial cut range (defined as distance from rPMT) - by default (0.5,2.5)

--gmin/--gmax                          define range for bonsai timing goodness cut - by default (0.1,0.3

--minEpmax/--maxEpmax                  define range for maximum prompt energy cut-off

--binwidthNX                           
--binwidthFid                          
--binwidthdT                           set optimisation granularity
--binwidthdR
--binwidthg
--binwidthEpmax

[--positiveScan                         only evaluate coincidences/sensitivity for delayed_nxcut>=prompt_nxcut (for Gd-water)]

[--negativeScan                         only evaluate coincidences/sensitivity for prompt_nxcut>=delayed_nxcut (for WbLS)]

```


***Step 1 - coincidences***

To run: ```cobraa --coincidences [detector options] [optimisation options]```

Additional option to evaluate individual part of the simulation:

```--evtype [element]		Runs individual signal or background (options are hartlepool_1, hartlepool_2, boulby_geo, singles, fasteneutrons etc).

Iterates over the other ranges to produce 2D coincidence maps in ```nx_p``` and ```closestPMT```.

Additional cuts:

```inner_hit``` > 4 applied to both prompt and delayed

```veto_hit``` < 4 applied to both prompt and delayed

All coincidences passing above cuts are then passed through a fast-neutron multiplicity cut which checks for further coincidences before and after pair.

Coincidence maps are then scaled to day rate and saved to file in ```fred_root_files*/coincidence_results.root``` (or ```core_root_files*/coincidence_results.root```).

***Step 2 - sensitivity***

To run: ```cobraa --sensitivity [detector options] [optimisation options] [sensitivity metric options]```

Additional _sensitivity metric options_:

```--poisson```           calculate significance and anomaly dwell time using poisson (gaussian-distributed bg)

```--poissonpoisson```    calculate significance and anomaly dwell time using poisson (poisson-distributed bg)

```--knoll```             calculate dwell time to 3 sigma detection at 95% confidence (gaussian sig and bg)

```--2sigma```            calculate dwell time to 2 sigma detection at 90% confidence (gaussian sig and bg)

```--optimiseSoB```       optimise signal over background rather than dwell time (useful in stats-limited regime)

Reads in the coincidence maps from Step 1, also adding together radionuclides and IBD backgrounds.

By default, performs optimisation of the anomaly dwell time from Gaussian significance over all optimisation variables in parallel.

Saves histograms of rates per day and significance for each combination of values in the nx_d and dT ranges to e.g. ```core_root_files*/sensitivity_results.root```.

Outputs details of rates, statistical and systematic errors plus dwell time to results_*.txt.
