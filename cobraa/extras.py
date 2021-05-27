'''
This is a place for additional tools which can be useful in the 
simulation-reconstruction-analysis process for verification 
and other purposes.

Author: Liz Kneale (May 2021)
'''

from ROOT import TFile
import pandas as pd
from .globals import *


def triggers():

    # Outputs details of the number of events/time simulated, plus
    # trigger and minimal-cut (reconstruction) efficiencies for 
    # all event types.
    # This is useful for checking that you have produced sufficient
    # statistics.
    pd.options.display.float_format = '{:.10f}'.format
    triggerdata = open("triggerdata.txt","w+")
    simsrequired = open("simsrequired.txt","w+")
    loclist=[]
    decaylist=[]
    isolist=[]
    timelist=[]
    eventlist=[]
    triggerlist=[]
    singleslist =[]
    for _p in proc:
        for _loc in proc[_p]:
            for _element in d[_p]:
                _tag = "%s_%s_%s"%(_element,_loc,_p)
                _tag = _tag.replace(" ","")
                _file = "fred_root_files%s/new_merged_Watchman_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
                _file = _file.replace(" ","")
                print(_tag," from ",_file)
                fredfile = TFile(_file)
                # check the root file is valid
                try:
                    runSummary = fredfile.Get('runSummary')
                    Entries = runSummary.GetEntries()
                    totalEvents = 0
                    for i in range(Entries):
                        runSummary.GetEntry(i)
                        totalEvents += runSummary.nEvents
                    data = fredfile.Get('data')
                    triggers = data.GetEntries()
                    singles = data.Draw("","n9>9 && closestPMT>0.5")
                    triggerrate = triggers/totalEvents*rates[_tag][0]
                    singlesrate = singles/totalEvents*rates[_tag][0]
                    rate = rates[_tag][0]
                    simtime = totalEvents/rates[_tag][0]/60./60./24.
                    loclist.append(_loc)
                    decaylist.append(_p)
                    isolist.append(_element)
                    eventlist.append(totalEvents)
                    timelist.append(simtime)
                    triggerlist.append(triggerrate)
                    singleslist.append(singlesrate)
                except:
                    simsrequired.writelines(f"{_tag}\n")

    # create a pandas dataframe with all the information
    df = pd.DataFrame({"Component":loclist, "Origin":decaylist,"Isotope":isolist,"Events":eventlist,"Time (days)":timelist,"Trigger rate (Hz)":triggerlist,"Singles rate (Hz)":singleslist})
    # format the names to make them more presentation-friendly
    df = df.replace("CHAIN_","",regex=True)
    df = df.replace("_NA","",regex=True)
    df = df.replace("LIQUID","GD-WATER")
    df = df.replace("STEEL_ACTIVITY","Co/Cs")
    df = df.replace("ROCK_1","ROCK (outer)")
    df = df.replace("ROCK_2","ROCK (inner)")
    df = df.replace("rock_neutrons","Neutrons")
    df = df.replace("FASTNEUTRONS","COSMOGENIC")
    df = df.replace("fast_neutrons","Fast neutrons")
    df = df.replace("pn_ibd","IBD")
    df = df.replace("A_Z","COSMOGENIC")
    df = df.replace("singles","Radioactivity")
    df = df.replace("SINGLES","All")
    df = df.replace("_hartlepool","",regex=True)
    df = df.sort_values(by=["Component","Origin","Isotope"])
    # convert to LaTex and do some formatting to work with siunitx
    triggerdata.writelines(df.to_latex(index=False,escape=False).replace('\\toprule', '\\hline').replace('\\midrule', '\\hline').replace('\\bottomrule','\\hline').replace('lllrrrr','|l|l|l|S|S|S|S|').replace("Component","{Component}").replace("Origin","{Origin}").replace("Isotope","{Isotope}").replace("Events","{Events}").replace("Time (days)","{Time (days)}").replace("Trigger rate (Hz)","{Trigger rate (Hz)}").replace("Singles rate (Hz)","{Singles rate (Hz)}"))
    return 0

