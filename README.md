**Cobraa Coincident-background reactor antineutrino analysis** is an analysis for reactor discovery via the simulation and evaluation of coincident events. It can also optionally use reconstruction of coincident events.

It's a toolchain which handles full-detector simulation through to sensitivity analysis
for the current NEO detector design options: 16m/22m with Gd-water or Gd-WbLS.

Adapted from the WATCHMAN Watchmakers package (Author: Marc Bergevin, LLNL), it constitutes a significant overhaul of the original Watchmakers code. Retains the convenience of the directory organisation and macro/job production but incorporates streamlining and a full analysis update.

To be used alongside the FRED (BONSAI) reconstruction or CoRe (BONSAI) pair reconstruction.


Analysis performs full evaluation of coincidences for both signal and background and 
optimises sensitivity to a reactor as a function of cuts in up to 7 dimensions: 

1. Prompt energy threshold
2. Delayed energy threshold
3. Fiducial volume
4. Time between triggers dT
5. Distance between triggers dR (FRED only)
6. Fit quality (BONSAI timing goodness)
7. Maximum prompt energy


To install:

    $git clone https://github.com/ait-watchman/cobraa

    $cd cobraa

    $./configure

    $source env_wm.sh


**Basic use:**

1. Create macros and simulation jobs for a small-ish number of events (16m/Gd-water/20% photocoverage by default)

     ```cobraa -m -j -e 2500 --lightSimWater```

2. Run jobs (locally by default)

   ```source job/job*.sh```

3. Run FRED or CoRe reconstruction (this stage is not incorporated into Cobraa)

4. Merge ntuple root files from reconstruction

   ```cobraa -M```   

4. Map coincidences after cuts

     ```cobraa --coincidences [--core]```

5. Calculate rates and optimise signal significance

     ```cobraa --sensitivity [--core]```


See Wiki for more details
