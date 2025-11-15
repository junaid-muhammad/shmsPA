#! /usr/bin/python
#
# Description:
# ================================================================
# Time-stamp: "2025-11-14 01:29:19 junaid"
# ================================================================
#
# Author:  Muhammad Junaid <mjo147@uregina.ca>
#
# Copyright (c) junaid
#
###################################################################################################################################################

# Calculation of missed triggers correction
# Script takes two input argument <ROOT Filename> <central_momentum GeV>

# Import relevant packages
import sys, math, os, subprocess
import array
import ROOT
from ROOT import TCanvas, TList, TPaveLabel, TColor, TGaxis, TH1F, TH2F, TPad, TStyle, gStyle, gPad, TLegend, TGaxis, TLine, TMath, TLatex, TPaveText, TArc, TGraphPolar, TText, TString
from ROOT import kBlack, kCyan, kRed, kGreen, kMagenta, kBlue

##################################################################################################################################################

# Input filename
Inputfile = sys.argv[1]
central_mom = sys.argv[2]

# Extract physics setting from Inputfile and rootfile name
physics_setting = "_".join(Inputfile.split("_")[:3])
center_mom = Inputfile.split("_")[2]
Input_ROOTfile = Inputfile + ".root"

# import file, tree
f1 = ROOT.TFile(Input_ROOTfile)
t1 = f1.Get("PA")

# define variables and assign to tree
s1x_time = array.array('d', [0.])
s1x_energy = array.array('d', [0.])
s1y_time = array.array('d', [0.])
s1y_energy = array.array('d', [0.])
s2x_time = array.array('d', [0.])
s2x_energy = array.array('d', [0.])
s2y_time = array.array('d', [0.])
s2y_energy = array.array('d', [0.])
NPE = array.array('i', [0])
copyNo = array.array('i', [0])
AGC_NPE = array.array('i', [0])
NGC_NPE = array.array('i', [0])
HGC_NPE = array.array('i', [0])
cal_energy = array.array('d', [0.])

t1.SetBranchAddress("S1XEnergy", s1x_energy)
t1.SetBranchAddress("S1YEnergy", s1y_energy)
t1.SetBranchAddress("S2XEnergy", s2x_energy)
t1.SetBranchAddress("S2YEnergy", s2y_energy)
t1.SetBranchAddress("S1XTime", s1x_time)
t1.SetBranchAddress("S1YTime", s1y_time)
t1.SetBranchAddress("S2XTime", s2x_time)
t1.SetBranchAddress("S2YTime", s2y_time)
t1.SetBranchAddress("S2YNPE", NPE)
t1.SetBranchAddress("CopyNo", copyNo)
t1.SetBranchAddress("AGCNPE", AGC_NPE)
t1.SetBranchAddress("HGCNPE", HGC_NPE)
t1.SetBranchAddress("NGCNPE", NGC_NPE)
t1.SetBranchAddress("CalEnergy", cal_energy)

# other things needed:
energy_threshold = 0.5  # MeV
time_window = 20        # ns
npe_threshold = 100     # unitless
nEvents = t1.GetEntries()

# initialize counters to zero
nMissed = 0
nProton = 0
nPion = 0
nPositron = 0
nKaon = 0
nContam = 0
nStopped = 0
P = float(central_mom) * 1000  # convert to MeV

# Define PID cuts as functions
def is_proton(AGC_NPE, HGC_NPE):
    # proton: no signal in HGC, AGC, hadron in Cal
    return AGC_NPE < 20 and HGC_NPE < 50

def is_kaon(AGC_NPE, HGC_NPE):
    # kaon: signal in AGC, hadron in Cal
    return AGC_NPE > 20 and HGC_NPE < 50

def is_positron(AGC_NPE, HGC_NPE, NGC_NPE, cal_energy, P):
    # electron: signals in all Cherenkovs, lepton in Cal
    return AGC_NPE > 20 and HGC_NPE > 50 and NGC_NPE > 5 and cal_energy/P > 0.7

def is_pion(AGC_NPE, HGC_NPE):
    # pion: signals in HGC, AGC, hadron in Cal
    return AGC_NPE > 20 and HGC_NPE > 50

# loop through tree to simulate 3/4 trigger
for i in range(nEvents):
    t1.GetEntry(i)
    nTrig = 0
    if s1x_energy[0] > energy_threshold:
        nTrig += 1
    dt1 = s1y_time[0] - s1x_time[0]
    dt2 = s2x_time[0] - s1x_time[0]
    dt3 = s2y_time[0] - s1x_time[0]
    if s1y_energy[0] > energy_threshold and abs(dt1) < time_window:
        nTrig += 1
    if s2x_energy[0] > energy_threshold and abs(dt2) < time_window:
        nTrig += 1
    if NPE[0] > npe_threshold and abs(dt3) < time_window:
        nTrig += 1
    if nTrig < 3:
        nMissed += 1
    if copyNo[0] != 0 and copyNo[0] != 4:
        nStopped += 1
    if nTrig > 2:
        if is_proton(AGC_NPE[0], HGC_NPE[0]):
            nProton += 1
        elif is_kaon(AGC_NPE[0], HGC_NPE[0]):
            nKaon += 1
        elif is_positron(AGC_NPE[0], HGC_NPE[0], NGC_NPE[0], cal_energy[0], P):
            nPositron += 1
        elif is_pion(AGC_NPE[0], HGC_NPE[0]):
            nPion += 1
        else:
            nContam += 1

# calculate absorption corrections
Proton_Absorption_Corr = (nEvents - nProton) * 100 / nEvents
Kaon_Absorption_Corr = (nEvents - nKaon) * 100 / nEvents
Pion_Absorption_Corr = (nEvents - nPion) * 100 / nEvents
Positron_Absorption_Corr = (nEvents - nPositron) * 100 / nEvents

# print results
print("=" * 50)
print(f"Total events: {nEvents}")
print(f"Missed 3/4 triggers: {nMissed} or {nMissed*100/nEvents:.2f}%")
print(f"Total stopped tracks: {nStopped} or {nStopped*100/nEvents:.2f}%")
print(f"3/4 was triggered by: {nProton} p+, {nPion} pi+, {nKaon} K+, {nPositron} positrons, and {nContam} contaminated tracks")
print("-" * 50)
print("Absorption correction in percentage (by incident particle type):")
print("Physics Setting:", physics_setting, " at Central Momentum:", central_mom, "GeV")
print(f"p+:  {Proton_Absorption_Corr:.2f}%")
print(f"K+:  {Kaon_Absorption_Corr:.2f}%")
print(f"pi+: {Pion_Absorption_Corr:.2f}%")
print(f"e+:  {Positron_Absorption_Corr:.2f}%")
print("=" * 50)

# Create output CSV file with new header and format
output_filename = f"{physics_setting}_pion_absorption_corrections_{central_mom}GeV.csv"
with open(output_filename, 'w') as outfile:
    outfile.write("physics_setting,central_mom,pion_absorption_corr(%)\n")
    outfile.write(f"{physics_setting},{central_mom},{Pion_Absorption_Corr:.2f}\n")
print(f"Absorption corrections saved to {output_filename}")

##################################################################################################################################################
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Set ROOT to batch mode explicitly, does not splash anything to screen
###############################################################################################################################################

# Histograms for plotting
HGCNPE_pion_uncut = ROOT.TH1D("HGCNPE_pion_uncut", "HGC NPE for Pions (no cut); HGC_NPE; Events", 100, 0, 800)
AeroNPE_pion_uncut = ROOT.TH1D("AeroNPE_pion_uncut", "AGC NPE for Pions (no cut); AGC_NPE; Events", 100, 0, 300)
HGCNPE_pion_cut = ROOT.TH1D("HGCNPE_pion_cut", "HGC NPE for Pions (AGC and HGC cuts); HGC_NPE; Events", 100, 0, 800)
AeroNPE_pion_cut = ROOT.TH1D("AeroNPE_pion_cut", "AGC NPE for Pions (AGC and HGC cuts); AGC_NPE; Events", 100, 0, 300)

# Filling histograms

for i in range(nEvents):
    t1.GetEntry(i)
    nTrig = 0
    if s1x_energy[0] > energy_threshold:
        nTrig += 1
    dt1 = s1y_time[0] - s1x_time[0]
    dt2 = s2x_time[0] - s1x_time[0]
    dt3 = s2y_time[0] - s1x_time[0]
    if s1y_energy[0] > energy_threshold and abs(dt1) < time_window:
        nTrig += 1
    if s2x_energy[0] > energy_threshold and abs(dt2) < time_window:
        nTrig += 1
    if NPE[0] > npe_threshold and abs(dt3) < time_window:
        nTrig += 1
    if nTrig > 2:
        HGCNPE_pion_uncut.Fill(HGC_NPE[0])
        AeroNPE_pion_uncut.Fill(AGC_NPE[0])
        if is_pion(AGC_NPE[0], HGC_NPE[0]):
            HGCNPE_pion_cut.Fill(HGC_NPE[0])
            AeroNPE_pion_cut.Fill(AGC_NPE[0])

# Plotting section
c1 = TCanvas("c1", "Absorption Corrections", 1000, 1000)
c1.Divide(2,2)
c1.cd(1)
HGCNPE_pion_uncut.SetLineColor(kRed)
HGCNPE_pion_uncut.SetTitle("HGC NPE for Pions (no cut)")
HGCNPE_pion_uncut.Draw()
c1.cd(2)
HGCNPE_pion_cut.SetLineColor(kBlue)
HGCNPE_pion_cut.SetTitle("HGC NPE for Pions (AGC and HGC cuts)")
HGCNPE_pion_cut.Draw()
c1.cd(3)
AeroNPE_pion_uncut.SetLineColor(kRed)
AeroNPE_pion_uncut.SetTitle("AGC NPE for Pions (no cut)")
AeroNPE_pion_uncut.Draw()
c1.cd(4)
AeroNPE_pion_cut.SetLineColor(kBlue)
AeroNPE_pion_cut.SetTitle("AGC NPE for Pions (AGC and HGC cuts)")
AeroNPE_pion_cut.Draw()
c1.SaveAs(f"Pion_Absorption_Corr_NPE_Distributions_{physics_setting}_{central_mom}GeV.pdf")

print("Process completed successfully.")