**Cobraa Coincident-background reactor antineutrino analysis** 

Coincidence analysis for reactor discovery

Adapted from the WATCHMAN Watchmakers package (Author: Marc Bergevin, LLNL).

Toolkit which handles full-detector simulation through to sensitivity analysis
for the WATCHMAN project.

Constitutes a significant overhaul of the original Watchmakers code. Retains the convenience of the directory organisation and macro/job production but incorporates streamlining and a full analysis update.

To be used alongside the FRED reconstruction.


Analysis performs full evaluation of coincidences for both signal and background and 
optimises sensitivity to a reactor as a function of cuts in 4 dimensions: 

prompt energy threshold

delayed energy threshold

fiducial volume

dT time between triggers

dR distance between triggers

Additional cuts applied:

Fit quality (BONSAI timing goodness)


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

    --cylinderSize 12         Set the tank diameter/height to 12m (default is 16m)

    --cylinderPct 10          Set the photocoverage to 10% (default is 20%)

    --cylinderPct 15          Set the photocoverage to 15% (default is 20%)
    
    --rPMT 6700               Set the inner PMT radius to 6.7m (default is 5.7m)

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

Macros are found in the watch-core/mac directory.

**Job generation**

Cobraa produces jobs/scripts for the following types of events by default:

Hartlepool cores (called 'big' and 'small', which relates to the currently reported power output)

Reactor background from more distant reactors at Boulby (boulby_world)

Geoneutrino background at Boulby (boulby_geo)

Fast neutron background

Radionuclide background

Accidental background (singles)

Relevant flags are:
```

    --reduced                 Runs only the processes and decays which tend to trigger/reconstruct in the fiducial. (Recommended for the coincidence                                     analysis and can be combined with --lightSim)
    
    --watchmakers             Generates jobs/scripts for individual singles simulations (all radioactive decays, IBD positron, IBD neutron)
    
    --lightSim                Runs the singles only in the components which are the main contributors to accidentals (can be combined with --reduced)

    --cluster                 Used to specify the job submission script header options (default is to run locally)
                              Options: lassen, sheffield, edinburgh, glasgow, etc (only lassen and sheffield so far).

    -N                        Sets the number of times to run each script (default is 40)
    
    --Heysham                 Uses Heysham IBD signal and background instead of Hartlepool and Boulby geo/world
```
Outputs are:

1. Scripts: script*.sh scripts to run the rat simulation for each type of event.

   Sources environment variables; runs rat with detector + process + geo + phys + rates macros.


2. Jobs: job*.sh jobs to either run the scripts locally (default) or on Lassen/Sheff clusters
   
   Defines cluster-dependent settings if required
  
   Runs the script 40 times
   
   For singles events, a number of jobs equivalent to ```nsetSingles``` (currently fixed in globals.py) are created. This is to make it possible to run enough statistics at the same time as keeping run times within cluster limits. 



***Reconstruction***

RAT-PAC output must be passed to the FRED reconstruction prior to the next stages in the process.
   


***Statistics evaluation***

This is an additional function which is useful to use after the reconstruction stage for ensuring that you have sufficient statistics.

To run: ```cobraa --triggers [detector options]```

Reads in the fred output for all event types available and finds:

 1. Number of events simulated

 2. Length of time simulated

 3. Trigger efficiency (trigger = interaction with PMT hits > 6 in 800 microseconds)

 4. Reconstruction efficiency (event passes minimal post-reconstruction cuts: inner PMT hits > 10, n9 (unscattered light) hits > 9, reconstructed vertex > 0.5m from inner PMT radius)

The information is saved to a triggerdata.txt in latex formatting (if using in LaTex, please use ```\usepackage{siunitx}``` in the preamble). Where a root file does not exist, the details are written in simsrequired.txt.

***Backgrounds plot***

This is another additional function which is useful for comparing the relative contributions of the components to the radioactive backgrounds. Plots the singles rates as a function of fiducial cut. NB this option only works where the backgrounds have been simulated individually.

To run: ```cobraa --backgrounds [detector options] [fiducial cut options]```

**Analysis**

This is the post-reconstruction stage. To run the reconstruction, use fred (https://github.com/AIT-WATCHMAN/FRED).

The analysis ultimately optimises reactor sensitivity in 4 dimensions:

```nx_p``` number of hits n in x ns detected from the prompt event in the pair
```nx_d``` number of hits n in x ns detected from the delayed event in the pair
```closestPMT``` perpendicular distance to the nearest PMT, as measured from the PMT radius ```rPMT``` (negative-valued for events occurring in the buffer/veto region)
```dT``` time between the two interactions in a pair
```dR``` distance between the two interactions in a pair

 There are two steps to the analysis:
 1. Evaluate the number of coincidences and scale to per day rate ```--coincidences```
 2. Evaluate the sensitivity based on the coincidence rates ```--sensitivity```

Relevant flags for the analysis steps define the ranges over which to optimise:

```
--minNXprompt/--maxNXprompt            define nx_p range

--minNXdelayed/--maxNXdelayed          define nx_d range

--dTmin/--dTmax                        define dT range (by default, uses a fixed cut in dT of 100us)

--dRmin/--dRmax                        define dR range (by default, uses a fixed cut in dR of 2m)

--minFid/--maxFid                      define fiducial cut range (defined as distance from rPMT)

--binwidthNX                           
--binwidthFid                          set optimisation granularity
--binwidthdT

-g                                     set bonsai timing goodness cut (between 0 and 1)

```


***Step 1 - coincidences***

To run: ```cobraa --coincidences [detector options]```

Iterates over the ```nx_d``` and ```dT``` ranges to produce 2D coincidence maps in ```nx_p``` and ```closestPMT```.

```closestPMT``` cut applied to both prompt and delayed event

```nx_p``` cut applied to prompt event in pair

```nx_d``` cut applied to delayed event in pair

Additional cuts:

```inner_hit``` > 4 applied to both prompt and delayed

```veto_hit``` < 4 applied to both prompt and delayed

```good_pos``` cut applied to both prompt and delayed

All coincidences passing above cuts are then passed through a fast-neutron multiplicity cut which checks for further coincidences before and after pair.

Coincidence maps are then scaled to day rate and saved to file in ```core_root_files*/coincidence_results.root```.

***Step 2 - sensitivity***

To run: ```cobraa --sensitivity [detector options]```

Reads in the coincidence maps from Step 1, also adding together radionuclides and IBD backgrounds.

Performs optimisation of the signal significance in terms of s/sqrt(b) in all of the 4 dimensions in parallel.

Saves histograms of rates per day and significance for each combination of values in the nx_d and dT ranges to ```core_root_files*/sensitivity_results.root```.

