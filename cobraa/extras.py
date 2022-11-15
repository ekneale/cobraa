
from ROOT import TFile,TH1D,TCanvas,THStack,TLegend,gPad,TColor
import pandas as pd
from .globals import *

# This is a place for additional tools which can be useful in the 
# simulation-reconstruction-analysis process for verification 
# and other purposes.
# Author: Liz Kneale (May 2021)


def triggers():

    # Outputs to 'triggerdata.txt' details of the number of events/time
    # simulated, plus daq-trigger and minimal-cut (energy and fiducial)
    # rates for all event types. Gives 90% upper-confidence limit on 
    # singles rate.
    # This is useful for checking that you have produced sufficient
    # statistics.
    # Outputs to 'simsmissing.txt' any root files not present.
    # Outputs to 'simsrequired.txt' any decay with rate > 1e-4Hz
    # (should agree with relevant lightSim option where complete)
    pd.options.display.float_format = '{:.10f}'.format
    triggerdata = open("triggerdata.txt","w+")
    simsmissing = open("simsmissing.txt","w+")
    simsrequired = open("simsrequired.txt","w+")
    loclist=[]
    decaylist=[]
    isolist=[]
    timelist=[]
    eventlist=[]
    triggerlist=[]
    singleslist =[]
    uclist = []
    for _p in proc:
        for _loc in proc[_p]:
            for _element in d[_p][_loc]:
                _tag = "%s_%s_%s"%(_element,_loc,_p)
                _tag = _tag.replace(" ","")
                _file = "fred_root_files%s/merged_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
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
                    singles = data.Draw("","n9>0")
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
                    # calculate 90% upper-confidence limit on singles count (normal)
                    uc = singles+1.645*sqrt(singles/totalEvents)
                    # calculate 90% upper-confidence limit on singles count (binomial)
                    #uc = singles+1.645/totalEvents*(singles*(1-singles/totalEvents))
                    # convert upper-confidence limit to singles rate
                    uc *= 1/totalEvents*rates[_tag][0]
                    uclist.append(uc)
                except:
                    simsmissing.writelines(f"{_tag}\n")
                if singlesrate>1e-3 or '210Tl' in _tag:
                    simsrequired.writelines(f"{_tag}\n")

    # create a pandas dataframe with all the information
    df = pd.DataFrame({"{Component}":loclist, "{Origin}":decaylist,"{Isotope}":isolist,"{Events}":eventlist,"{Time (days)}":timelist,"{Trigger rate (Hz)}":triggerlist,"{Singles rate (Hz)}":singleslist,"{90% UCL}":uclist})
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
    df = df.sort_values(by=["{Component}","{Origin}","{Isotope}"])
    # convert to LaTex and do some formatting to work with siunitx
    triggerdata.writelines(df.to_latex(index=False,escape=False).replace('\\toprule', '\\hline').replace('\\midrule', '\\hline').replace('\\bottomrule','\\hline').replace('lllrrrrr','|l|l|l|S|S|S|S|S|'))
    return 0


def backgrounds():

    hPMT = TH1D('hPMT','hPMT',binFid,rangeFidmin,rangeFidmax)
    hPSUP = TH1D('hPSUP','hPSUP',binFid,rangeFidmin,rangeFidmax)
    hTANK = TH1D('hTANK','hTANK',binFid,rangeFidmin,rangeFidmax)
    hLIQUID = TH1D('hLIQUID','hLIQUID',binFid,rangeFidmin,rangeFidmax)
    hFN = TH1D('hFN','hFN',binFid,rangeFidmin,rangeFidmax)
    hIBEAM = TH1D('hIBEAM','hIBEAM',binFid,rangeFidmin,rangeFidmax)
    hROCK = TH1D('hROCK','hROCK',binFid,rangeFidmin,rangeFidmax)
    hRN = TH1D('hRN','hRN',binFid,rangeFidmin,rangeFidmax)
    hIBD = TH1D('hIBD','hIBD',binFid,rangeFidmin,rangeFidmax)

    for _p in proc:
        for _loc in proc[_p]:
            for _element in d[_p]:
                _tag = "%s_%s_%s"%(_element,_loc,_p)
                _tag = _tag.replace(" ","")
                _file = "fred_root_files%s/merged_%s_%s_%s.root"%(additionalString,_element,_loc,_p)
                _file = _file.replace(" ","")
                if 'hartlepool' in _tag or 'heysham' in _tag or 'mono' in _tag or 'boulby' in _tag:
                    continue
                fredfile = TFile(_file)
                # check the root file is valid
                totalEvents = 0
                try:
                    runSummary = fredfile.Get('runSummary')
                    Entries = runSummary.GetEntries()
                    for i in range(Entries):
                        runSummary.GetEntry(i)
                        totalEvents += runSummary.nEvents
                    data = fredfile.Get('data')

                except:
                    print("cannot open ",_tag,"from ",_file)
                    continue
                print("opening ",_tag,"from ",_file)

                for fidcut in drange(minFid,rangeFidmax,binwidthFid):
                    nevts = data.Draw("","n9>0 && closestPMT/1000.>%f"%(fidcut),"goff")
                    rate = nevts/totalEvents*rates[_tag][0]
                    if 'PMT' in _tag:
                        hPMT.Fill(fidcut,rate)
                    elif 'PSUP' in _tag:
                        hPSUP.Fill(fidcut,rate)
                    elif 'TANK' in _tag:
                        hTANK.Fill(fidcut,rate)
                    elif 'LIQUID' in _tag and 'NA' in _tag:
                        hLIQUID.Fill(fidcut,rate)
                    elif 'FASTNEUTRONS' in _tag:
                        hFN.Fill(fidcut,rate)
                    elif 'IBEAM' in _tag:
                        hIBEAM.Fill(fidcut,rate)
                    elif 'ROCK'in _tag and 'NA' in _tag:
                        hROCK.Fill(fidcut,rate)
                    elif 'RADIOGENICS' in _tag:
                        hROCK.Fill(fidcut,rate)
                    elif 'A_Z' in _tag:
                        hRN.Fill(fidcut,rate)
#                    elif 'pn_ibd' in _tag:
#                        hIBD.Fill(fidcut,rate)
    

    hPMT.SetLineColor(TColor.GetColor(Tol_bright[0]))
    hPSUP.SetLineColor(TColor.GetColor(Tol_bright[1]))
    hTANK.SetLineColor(TColor.GetColor(Tol_bright[2]))
    hLIQUID.SetLineColor(TColor.GetColor(Tol_bright[3]))
#    hFN.SetLineColor(TColor.GetColor(Tol_bright[4]))
    hIBEAM.SetLineColor(TColor.GetColor(Tol_bright[5]))
    hROCK.SetLineColor(TColor.GetColor(Tol_bright[6]))
#    hRN.SetLineColor(TColor.GetColor(Tol_bright[7]))
#    hIBD.SetLineColor(1)

    c1 = TCanvas()
    hs = THStack("hs","")
    hs.Add(hPMT,"hist")
    hs.Add(hPSUP,"hist")
    hs.Add(hTANK,"hist")
    hs.Add(hLIQUID,"hist")
    hs.Add(hIBEAM,"hist")
    hs.Add(hROCK,"hist")
    hs.Draw("nostack")
    hs.GetXaxis().SetTitle("Fiducial cut - distance from inner PMT radius (m)")
    hs.GetYaxis().SetTitle("Singles rate (Hz)")
    gPad.SetLogy()

    leg = TLegend(0.68,0.68,0.98,0.98)
    leg.AddEntry(hPMT,"Inner PMT","l")
    leg.AddEntry(hPSUP,"PSUP","l")
    leg.AddEntry(hTANK,"Tank","l")
    leg.AddEntry(hLIQUID,"Detector medium","l")
    leg.AddEntry(hIBEAM,"I beam","l")
    leg.AddEntry(hROCK,"Rock","l")
    leg.Draw()
    c1.SaveAs("backgrounds.png")
    c1.SaveAs("backgrounds.C")

    return hs
