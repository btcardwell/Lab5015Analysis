#! /usr/bin/env python
import os
import shutil
import glob
import math
import array
import sys
import time
import argparse
import json


parser = argparse.ArgumentParser(description='Module characterization summary plots')
parser.add_argument("-i",  "--inputLabels",   required=True, type=str, help="comma-separated list of input labels")
parser.add_argument("-m",  "--resMode",       required=True, type=int, help="resolution mode: 2 - tDiff, 1 - tAve")
parser.add_argument("-o",  "--outFolder",     required=True, type=str, help="out folder")
parser.add_argument("-l",  "--plotLabel",     required=False, type=str, help="single label for all plots")
args = parser.parse_args()


import ROOT
import CMS_lumi, tdrstyle

#set the tdr style
tdrstyle.setTDRStyle()
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetLabelSize(0.055,'X')
ROOT.gStyle.SetLabelSize(0.055,'Y')
ROOT.gStyle.SetTitleSize(0.07,'X')
ROOT.gStyle.SetTitleSize(0.07,'Y')
ROOT.gStyle.SetTitleOffset(1.05,'X')
ROOT.gStyle.SetTitleOffset(1.1,'Y')
ROOT.gStyle.SetLegendFont(42)
ROOT.gStyle.SetLegendTextSize(0.040)
ROOT.gStyle.SetPadTopMargin(0.07)
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kWarning


def getTimeResolution(h1_deltaT):
   
   tRes = [-1,-1]

   h1_deltaT.GetXaxis().SetRangeUser(h1_deltaT.GetMean() - 5*h1_deltaT.GetRMS(), h1_deltaT.GetMean() + 5*h1_deltaT.GetRMS())
                    
   fitFunc = ROOT.TF1('fitFunc','gaus',-10000, 10000)
   fitFunc.SetLineColor(ROOT.kGreen+3)
   fitFunc.SetLineWidth(2)
   fitFunc.SetParameters(h1_deltaT.GetMaximum(),h1_deltaT.GetMean(), h1_deltaT.GetRMS())
   
   fitXMin = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) - 200
   fitXMax = h1_deltaT.GetBinCenter(h1_deltaT.GetMaximumBin()) + 200.
   #fitXMin = h1_deltaT.GetMean() - 3*h1_deltaT.GetRMS()
   #fitXMax = h1_deltaT.GetMean() + 3*h1_deltaT.GetRMS()
   fitFunc.SetRange(fitXMin, fitXMax)
   h1_deltaT.Fit('fitFunc','QNRL','', fitXMin, fitXMax)
   #fitFunc.SetRange(fitFunc.GetParameter(1) - 3.0*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 3.0*fitFunc.GetParameter(2))
   fitFunc.SetRange(fitFunc.GetParameter(1) - 1.0*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 1.0*fitFunc.GetParameter(2))
   h1_deltaT.Fit('fitFunc','QNRL')
   fitFunc.SetRange(fitFunc.GetParameter(1) - 2.5*fitFunc.GetParameter(2), fitFunc.GetParameter(1) + 2.5*fitFunc.GetParameter(2))
   h1_deltaT.Fit('fitFunc','QRSL+')

   #if (fitFunc==None): continue                    
   #if (fitFunc.GetParameter(2) > 1000): continue
   #if (fitFunc.GetParameter(2) < 20): continue
   #if (fitFunc.GetParError(2) > 200): continue
   tRes = [ fitFunc.GetParameter(2),fitFunc.GetParError(2)]
   return tRes


# ====================================


# INPUT
inputdir = '/afs/cern.ch/user/r/rmenonra/bTB/Lab5015Analysis/plots/TOFHIR2C/'
#source = 'Laser'
source = 'TB'


# OUTPUT
outdir  = '/eos/user/r/rmenonra/www/MTD/TOPHIR2C/plotssummary'
outdir=outdir+args.outFolder
outFileName = inputdir+'/summaryPlots_'+args.outFolder+'.root'
print 'Saving root file ', outFileName
print 'Saving plots in ', outdir
outfile = ROOT.TFile(outFileName, 'RECREATE' )



# ranges for plots
tResMin = 0
tResMax = 140
tResMaxTh = 200
vovMax = 4.0


# create files list
label_list = (args.inputLabels.split(','))
print label_list

# resolution mode : 0 : /1;  1: /sqrt(2) if CTR, 2: /2 if TDiff
kscale = 2.
if (args.resMode == 0): kscale = 1
if (args.resMode == 1): kscale = math.sqrt(2)


peaks = []
enBins = []
refPeak = 0
if (source == 'Na22'):
    #peaks = [511, 1275]
    #enBins = [5, 7] ##??? non me li ricordo
    #refPeak = 511
    peaks = [1275]
    enBins = [1] ##??? non me li ricordo
    refPeak = 1275
if (source == 'Laser'):
    peaks = [0]
    enBins = [1]
if (source == 'TB'):
    peaks = [0]
    enBins = [1]


# --- colors
cols = { 0.50 : 51,
         0.60 : 51+8,
         0.80 : 51+16,
         1.00 : 51+24,
         1.25 : 51+32,
         1.50 : 51+40,
         2.00 : 51+44,
         2.50 : 51+46,
         3.00 : 51+47,
         3.50 : 51+48
}
#module colors
mols={ 813: 06,
       820: 51+5,
       815: 46,
       200: 50}
marks={ 813: 8, 820: 33 }

       
# --- prepare output dir
if (os.path.isdir(outdir) == False): 
    os.system('mkdir %s'%outdir)
os.system('mkdir %s/summaryPlots/'%outdir)
os.system('mkdir %s/summaryPlots/tot/'%outdir)
os.system('mkdir %s/summaryPlots/energy/'%outdir)
os.system('mkdir %s/summaryPlots/timeResolution/'%outdir)
os.system('mkdir %s/summaryPlots/timeResolution/fits/'%outdir)
    

# -- ref threhsold
thRef = 11

# -- get list of bars, Vovs, thresholds to be analyzed
bars = []
thresholds = []
Vovs = []
Module = [] 
for label in label_list:
    module =int(label.split('_')[2][4:7])
    Module.append(module)
    inputFile = None
    inputFile = ROOT.TFile.Open(inputdir+'/moduleCharacterization_step2_%s.root'%label)
    listOfKeys = [key.GetName().replace('h1_deltaT_totRatioCorr_','') for key in ROOT.gDirectory.GetListOfKeys() if key.GetName().startswith('h1_deltaT_totRatioCorr_bar')]
    for k in listOfKeys:
       # print k
        barNum = int (k.split('_')[0][3:5])
        bars.append(barNum)
        vov = float (k.split('_')[1][3:7])
        Vovs.append(vov)
        thr = int (k.split('_')[2][2:4])
        thresholds.append(thr)
# remove duplicates
bars = [i for n, i in enumerate(bars) if i not in bars[:n]]
Vovs = [i for n, i in enumerate(Vovs) if i not in Vovs[:n]]
thresholds = [i for n, i in enumerate(thresholds) if i not in thresholds[:n]]
Module = [i for n, i in enumerate(Module) if i not in Module[:n]]
bars.sort()
Vovs.sort()
thresholds.sort()
Module.sort()
goodBars = {}
goodThrs= {}
VovsEff = {}
plots_label = ''

print('Vovs:',Vovs)
print('Bars:',bars)

if ('LYSO528' in args.outFolder):
    plots_label = 'HPK (15 #mum, type2) + LYSO528 (prod5, type2)'
    for vov in Vovs:
        VovsEff[vov] = vov 
    goodBars[3.50] = [6,7,8,9,10,11,12,13]
    goodBars[2.00] = [6,7,8,9,10,11,12,13]
    goodBars[1.50] = [6,7,8,9,10,11,12,13]
    goodBars[1.25] = [6,7,8,9,10,11,12,13]
    goodBars[1.00] = [6,7,8,9,10,11,12,13]
    goodBars[0.80] = [6,7,8,9,10,11,12,13]
    goodBars[0.60] = [6,7,8,9,10,11,12,13]

elif ('LYSO813' in args.outFolder):
    plots_label = 'HPK (25 #mum, type2) + LYSO813 (prod1, type2)'
    for vov in Vovs:
        VovsEff[vov] = vov 
    goodBars[3.50] = [6,7,8,9,10,11,12,13]
    goodBars[2.00] = [6,7,8,9,10,11,12,13]
    goodBars[1.50] = [6,7,8,9,10,11,12,13]
    goodBars[1.25] = [6,7,8,9,10,11,12,13]
    goodBars[1.00] = [6,7,8,9,10,11,12,13]
    goodBars[0.80] = [6,7,8,9,10,11,12,13]
    goodBars[0.60] = [6,7,8,9,10,11,12,13]
    goodBars[0.50] = [6,7,8,9,10,11,12,13]

elif ('LYSO818' in args.outFolder):
    plots_label = 'HPK (25 #mum, type1) + LYSO818 (prod1, type1)'
    for vov in Vovs:
        VovsEff[vov] = vov 
    goodBars[3.50] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13]
    goodBars[2.00] = [0,1,2,3,4,5,7,8,9,10,11,12,13]
    goodBars[1.50] = [0,1,2,3,4,5,7,8,9,10,11,12,13]
    goodBars[1.00] = [0,1,2,3,4,5,7,8,9,10,11,12,13]
    goodBars[0.80] = [0,1,2,3,4,5,7,8,9,10,11,12,13]
    goodBars[0.50] = [0,2,3,4,5,7,8,9,10,11,12,13]

if ('LYSO820' in args.outFolder):
   # plots_label = 'HPK (30 #mum, type2) + LYSO820 (prod1, type2)'
    for vov in Vovs:
        VovsEff[vov] = vov 
    goodBars[3.50] = [6,7,8,9,10,11,12,13]
   # goodBars[2.00] = [6,7,8,9,10,11,12,13]
    goodBars[1.50] = [6,7,8,9,10,11,12,13]
    goodBars[1.25] = [6,7,8,9,10,11,12,13]
    goodBars[1.00] = [6,7,8,9,10,11,12,13]
    goodBars[0.80] = [6,7,8,9,10,11,12,13]
    goodBars[0.60] = [6,7,8,9,10,11,12,13]

   # goodThrs[2.00]=[5,7,11,15,20,25]
    goodThrs[1.50]=[5,7,11,15,20,25]
    goodThrs[1.25]=[5,7,11,15,20,25]
    goodThrs[1.00]=[5,7,11,15,20,25]
    goodThrs[0.80]=[5,7,11,15,20,25]
    goodThrs[0.60]=[5,7,11,15,20,25]
    goodThrs[3.50]=[15,20,25]



   
   
elif ('LYSO528' in args.outFolder):
   plots_label = 'HPK (15#mum) + LYSO528 (prod5, type2)'
   for vov in Vovs:
      VovsEff[vov] = vov 
      if (args.versionTOFHIR=='TOFHIR2C'): 
         goodBars[vov] = bars
      else:
         goodBars[3.50] = [2,3,4,5,7,8,9,10,11,12,13] 
         goodBars[2.00] = [2,3,4,5,7,8,9,10,11,12,13] 
         goodBars[1.50] = [2,3,4,5,7,8,9,10,11,12,13] 
         goodBars[1.00] = [2,3,4,5,7,8,9,10,11,12,13] 

elif ('HPK_nonIrr_LYSO813_T-30C' in args.outFolder):
   plots_label = 'HPK (25 #mum) + LYSO813 (prod1, type2) T=-30#circC '
   for vov in Vovs:
      VovsEff[vov] = vov 
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[3.50] = [0,3,4,5,7,8,9,10,11,12,13,15] 
         goodBars[1.50] = [0,3,4,5,7,8,9,10,11,12,13,15] 
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,8,12,13,15]

elif ('HPK_nonIrr_LYSO813_T-15C' in args.outFolder):
   plots_label = 'HPK (25 #mum) + LYSO813 (prod1, type2) T=-15#circC '
   for vov in Vovs:
      VovsEff[vov] = vov 
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[3.50] = [0,2,3,4,5,7,8,9,10,11,12,13,15] 
         goodBars[1.50] = [0,3,4,5,7,8,9,10,11,12,13,15] 
         goodBars[0.80] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,8,11,12,13,15]

elif ('HPK_nonIrr_LYSO813_T0C' in args.outFolder):
   plots_label = 'HPK (25 #mum) + LYSO813 (prod1, type2) T=0#circC '
   for vov in Vovs:
      VovsEff[vov] = vov 
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[3.50] = [0,1,2,3,4,5,7,8,9,10,11,12,13,15] 
         goodBars[1.50] = [0,2,3,4,5,7,8,9,10,11,12,13,15] 
         goodBars[0.80] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,8,9,12,13,15]

elif ('HPK_nonIrr_LYSO813_T5C' in args.outFolder):
   plots_label = 'HPK (25 #mum) + LYSO813 (prod1, type2) T=5#circC '
   for vov in Vovs:
      VovsEff[vov] = vov 
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[3.50] = [0,1,2,3,4,5,7,8,9,10,11,12,13,15] 
         goodBars[1.50] = [0,2,3,4,5,7,8,9,10,11,12,13,15] 
         goodBars[0.80] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[0.50] = [0,3,4,5,7,8,9,12,13,15]

elif ('HPK_nonIrr_LYSO813_T15C' in args.outFolder):
   plots_label = 'HPK (25 #mum) + LYSO813 (prod1, type2) T=15#circC '
   for vov in Vovs:
      VovsEff[vov] = vov 
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[3.50] = [0,1,2,3,4,5,7,8,9,10,11,12,13,14,15] 
         goodBars[1.50] = [0,1,2,3,4,5,7,8,9,10,11,12,13,15] 
         goodBars[0.80] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,8,9,11,12,13,15]

elif ('LYSO814' in args.outFolder):
   plots_label = 'HPK (20 #mum) + LYSO814 (prod1, type2)'
   for vov in Vovs:
      VovsEff[vov] = vov 
   if (args.versionTOFHIR=='TOFHIR2C'):
      goodBars[vov] = bars
   else:
      goodBars[3.50] = [0,2,3,4,5,7,8,9,10,11,12,13] 
      goodBars[1.50] = [0,2,3,4,5,7,8,9,10,11,12,13]
      goodBars[1.00] = [0,2,3,4,5,7,8,9,10,11,12,13]
      goodBars[0.80] = [0,3,4,5,7,8,9,10,11,12,13]
      goodBars[0.50] = [0,3,4,5,7,8,9,10,11,12,13]

elif ('LYSO824' in args.outFolder):
   plots_label = 'HPK (25 #mum, low Cgrid) + LYSO824 (prod5, type2)'
   for vov in Vovs:
      VovsEff[vov] = vov 
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[3.50] = [0,1,3,4,5,7,8,9,10,11,12,13] 
         goodBars[1.50] = [0,1,3,4,5,7,8,9,10,11,12,13] 
         goodBars[0.80] = [0,2,3,4,5,7,8,9,10,11,12,13]
         goodBars[0.50] = [0,3,4,5,7,8,9,11,12,13]

elif ('HPK_2E14_LYSO815_T-40C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 2E14) + LYSO815 (prod1,type2) T=-40#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data,'HPK_2E14_LYSO815_T-40C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:      
         goodBars[2.00] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.50] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.25] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,12,13]

elif ('HPK_2E14_LYSO815_T-35C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 2E14) + LYSO815 (prod1,type2) T=-35#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data,'HPK_2E14_LYSO815_T-35C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = [ b for b in bars if b!=3] 
      else:
         goodBars[2.00] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.50] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.25] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,12,13]

elif ('HPK_2E14_LYSO815_T-30C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 2E14) + LYSO815 (prod1,type2) T=-30#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data, 'HPK_2E14_LYSO815_T-30C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = [ b for b in bars if b!=3] 
      else:
         goodBars[2.00] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.50] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.25] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,12,13]

elif ('HPK_2E14_LYSO825_T-40C' in args.outFolder):
   plots_label = 'HPK (20 #mum, 2E14) + LYSO825 (prod5,type2) T=-40#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data,'HPK_2E14_LYSO825_T-40C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[2.00] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.50] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[1.25] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,12,13]
         goodBars[0.60] = [4,5,7]

elif ('HPK_2E14_LYSO825_T-35C' in args.outFolder):
   plots_label = 'HPK (20 #mum, 2E14) + LYSO825 (prod5,type2) T=-35#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data,'HPK_2E14_LYSO825_T-35C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov]  = [0,1,2,4,5,7,8,9,10,11,12,13,14]
         goodBars[0.60] = [1,2,4,5,7,8,9,10,11]
      else:
         goodBars[2.00] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.50] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[1.25] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,12,13]
         goodBars[0.60] = [4,5,7]

elif ('HPK_2E14_LYSO825_T-30C' in args.outFolder):
   plots_label = 'HPK (20 #mum, 2E14) + LYSO825 (prod5,type2) T=-30#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data, 'HPK_2E14_LYSO825_T-30C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[2.00] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.50] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[1.25] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,12,13]
         goodBars[0.60] = [4,5,7]

elif ('HPK_1E14_LYSO819_T-37C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 1E14) + LYSO819 (prod1,type1) T=-37#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data,'HPK_1E14_LYSO819_T-37C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[2.50] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[2.00] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[1.50] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.25] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,12,13]
         
elif ('HPK_1E14_LYSO819_T-32C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 1E14) + LYSO819 (prod1,type1) T=-32#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data, 'HPK_1E14_LYSO819_T-32C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[2.50] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[2.00] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[1.50] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.25] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,12,13]

elif ('HPK_1E14_LYSO819_T-27C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 1E14) + LYSO819 (prod1,type1) T=-27#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data,'HPK_1E14_LYSO819_T-27C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[2.50] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,15]
         goodBars[2.00] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.50] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.25] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,12,13]

elif ('HPK_1E14_LYSO819_T-22C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 1E14) + LYSO819 (prod1,type1) T=-22#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data, 'HPK_1E14_LYSO819_T-22C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[2.50] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,15]
         goodBars[2.00] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.50] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.25] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[1.00] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,12,13]

elif ('HPK_1E13_LYSO829_T-19C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 1E13) + LYSO829 (prod5,type1) T=-19#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data, 'HPK_1E13_LYSO829_T-19C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[3.00] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[2.00] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[1.50] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         
elif ('HPK_1E13_LYSO829_T-32C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 1E13) + LYSO829 (prod5,type1) T=-32#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data, 'HPK_1E13_LYSO829_T-32C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:      
         goodBars[3.00] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[2.00] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[1.50] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13]
      
elif ('HPK_1E13_LYSO829_T0C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 1E13) + LYSO829 (prod5,type1) T=0#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data, 'HPK_1E13_LYSO829_T0C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[2.50] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[2.00] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[1.50] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[1.25] = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
         goodBars[0.80] = [0,3,4,5,7,8,9,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,8,9,11,12,13,15]

elif ('HPK_1E13_LYSO829_T12C' in args.outFolder):
   plots_label = 'HPK (25 #mum, 1E13) + LYSO829 (prod5,type1) T=12#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data, 'HPK_1E13_LYSO829_T12C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
      else:
         goodBars[1.25] = [0,2,3,4,5,7,8,9,10,11,12,13,14,15]
         goodBars[1.00] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[0.80] = [0,2,3,4,5,7,8,9,10,11,12,13,15]
         goodBars[0.60] = [0,3,4,5,7,8,9,11,12,13,15]

elif ('HPK_1E14_LYSO844_T-30C' in args.outFolder):
   plots_label = 'HPK (15 #mum, 1E14) + LYSO844 (prod10,type2) T=-30#circC'
   for vov in Vovs:
      VovsEff[vov] = getVovEffDCR(data, 'HPK_1E14_LYSO844_T-30C', ('%.02f'%vov))[0]
      if (args.versionTOFHIR=='TOFHIR2C'):
         goodBars[vov] = bars
if('LYSO200' in args.outFolder):
   for vov in Vovs:
        VovsEff[vov] = vov

   goodBars[2.00] = [6,7,8,9,10,11,12,13]
   goodBars[1.50] = [5,6,7,8,9,10,11,12,13]
   goodBars[1.25] = [5,6,7,8,9,10,12,13]
   goodBars[1.00] = [3,4,5,6,7,8,9,10,11,12,13,14,15]
   goodBars[0.80] = [3,4,5,6,7,8,9,10,11,12,13,14,15]
   goodBars[0.60] = [3,4,5,6,7,8,9,10,11,12,13,14,15]

   goodThrs[2.00]=[15,20,25]
   goodThrs[1.50]=[11,15,20,25]
   goodThrs[1.25]=[11,15,20,25]
   goodThrs[1.00]=[5,7,11,15,20,25]
   goodThrs[0.80]=[5,7,11,15,20,25]
   goodThrs[0.60]=[5,7,11,15,20,25]
   
 
else:
   for vov in Vovs:
      VovsEff[vov] = vov
      goodBars[vov] = bars 


print 'bars:', bars
print 'good bars:', goodBars
print 'Vovs:',Vovs
print 'thresholds:', thresholds
print 'good thresholds: ', goodThrs


# --- Summary graphs
g_tot_vs_th  = {} # g [bar, l, vov] 
g_tot_vs_vov  = {} # g [bar, l, th] 
g_tot_vs_bar = {} # g [l, vov] 

g_energy_vs_th  = {} # g [bar, l, vov, peak] 
g_energy_vs_vov  = {} # g [bar, l, th, peak] 
g_energy_vs_bar = {} # g [l, vov, peak,module] 

g_deltaT_energyRatioCorr_vs_th  = {} # g [bar, vov, energyBin] 
g_deltaT_energyRatioCorr_vs_vov  = {} # g [bar, th, energyBin] 
g_deltaT_energyRatioCorr_vs_bar = {} # g [viv, thr, energyBin,module] 

g_deltaT_totRatioCorr_vs_th  = {} # g [bar, vov, energyBin] 
g_deltaT_totRatioCorr_vs_vov  = {} # g [bar, th, energyBin] 
g_deltaT_totRatioCorr_vs_bar = {} # g [viv, thr, energyBin,module] 

g_deltaT_energyRatioCorr_totRatioCorr_vs_th  = {} # g [bar, vov, energyBin] 
g_deltaT_energyRatioCorr_totRatioCorr_vs_vov  = {} # g [bar, th, energyBin] 
g_deltaT_energyRatioCorr_totRatioCorr_vs_bar = {} # g [viv, thr, energyBin] 

g_deltaT_energyRatioCorr_bestTh_vs_vov = {} # g [bar, energyBin] 
g_deltaT_energyRatioCorr_bestTh_vs_bar = {} # g [vov, energyBin,module] 

g_deltaT_totRatioCorr_bestTh_vs_vov = {} # g [bar, energyBin] 
g_deltaT_totRatioCorr_bestTh_vs_bar = {} # g [vov, energyBin,module] 

g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov = {} # g [bar, energyBin] 
g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar = {} # g [vov, energyBin,module] 

g_deltaT_totRatioCorr_bestTh_vs_vov_average = {} # g [energyBin]
g_deltaT_energyRatioCorr_bestTh_vs_vov_average = {} # g [energyBin]
g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average = {} # g [energyBin]


for bar in bars:
   for l in ['L','R','L-R']:
      for vov in Vovs:
          for module in Module:
            # bars= goodBars[vov]
            # print "Added=",vov, bars
             g_tot_vs_th[bar, l, vov] = ROOT.TGraphErrors()
             for peak in peaks: g_energy_vs_th[bar, l, vov, peak,module] = ROOT.TGraphErrors()
             for enBin in enBins: 
                g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module] = ROOT.TGraphErrors()
                g_deltaT_totRatioCorr_vs_th[bar, vov, enBin] = ROOT.TGraphErrors()
                g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin] = ROOT.TGraphErrors()
             for thr in thresholds:
                g_tot_vs_vov[bar, l, thr] = ROOT.TGraphErrors()
                for peak in peaks: g_energy_vs_vov[bar, l, thr, peak,module] = ROOT.TGraphErrors()
                for enBin in enBins: 
                   g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin] = ROOT.TGraphErrors()
                   g_deltaT_totRatioCorr_vs_vov[bar, thr, enBin] = ROOT.TGraphErrors()
                   g_deltaT_energyRatioCorr_totRatioCorr_vs_vov[bar, thr, enBin] = ROOT.TGraphErrors()

for l in ['L','R','L-R']:
    for vov in Vovs:
          for thr in thresholds:
             for module in Module:
                g_tot_vs_bar[l, vov, thr,module] = ROOT.TGraphErrors()                          
                for peak in peaks: g_energy_vs_bar[l, vov, thr, peak,module] = ROOT.TGraphErrors()

              
                if (l=='L-R'):
                   for enBin in enBins: 
                      g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin,module] = ROOT.TGraphErrors()
                      g_deltaT_totRatioCorr_vs_bar[vov, thr, enBin,module] = ROOT.TGraphErrors()
                      g_deltaT_energyRatioCorr_totRatioCorr_vs_bar[vov, thr, enBin,module] = ROOT.TGraphErrors()
                      
for bar in bars:
   for enBin in enBins: 
      g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin] = ROOT.TGraphErrors()
      g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin] = ROOT.TGraphErrors()
      g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin] = ROOT.TGraphErrors()
      
for module in Module:
 for vov in Vovs:
    for enBin in enBins: 
        g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module] = ROOT.TGraphErrors()
        g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module] = ROOT.TGraphErrors()
        g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module] = ROOT.TGraphErrors()
        
for enBin in enBins: 
   g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin] = ROOT.TGraphErrors()
   g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin] = ROOT.TGraphErrors()
   g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin] = ROOT.TGraphErrors()


# --- Read the histograms from moduleCharacterization_step2 file:
for label in label_list:
       module =int(label.split('_')[2][4:7])
       print label
       inputFile == None
       inputFile = ROOT.TFile.Open(inputdir+'/moduleCharacterization_step2_%s.root'%label)

       for vov in Vovs:
          
         for bar in bars:
           if (bar not in goodBars[vov]): continue
           
           for thr in thresholds:
              if (thr not in goodThrs[vov]): continue
              for l in ['L','R','L-R']:
                                                         
                   if ( inputFile.GetListOfKeys().Contains('h1_tot_bar%02d%s_Vov%.02f_th%02d'%(bar, l, vov, thr)) ): 
                      h1_tot = inputFile.Get('h1_tot_bar%02d%s_Vov%.02f_th%02d'%(bar, l, vov, thr))
                      if (h1_tot==None or h1_tot.GetEntries()==0): continue
                                                 
                      max1 = h1_tot.GetBinCenter(h1_tot.GetMaximumBin())
                      h1_tot.GetXaxis().SetRangeUser(0.25*max1,2.*max1)
                      fitFunc = ROOT.TF1('fitFunc','gaus',max1-0.05*max1,max1+0.05*max1)
                      h1_tot.Fit('fitFunc','QNRS+')
                      h1_tot.Fit('fitFunc','QNS+','', fitFunc.GetParameter(1)-fitFunc.GetParameter(2), fitFunc.GetParameter(1)+fitFunc.GetParameter(2))
                      
                      g_tot_vs_th[bar, l, vov].SetPoint(g_tot_vs_th[bar, l, vov].GetN(), thr, fitFunc.GetParameter(1) )    
                      g_tot_vs_th[bar, l, vov].SetPointError(g_tot_vs_th[bar, l, vov].GetN()-1, 0, fitFunc.GetParError(1) )    
                      
                      g_tot_vs_vov[bar, l, thr].SetPoint(g_tot_vs_vov[bar, l, thr].GetN(), VovsEff[vov], fitFunc.GetParameter(1) )    
                      g_tot_vs_vov[bar, l, thr].SetPointError(g_tot_vs_vov[bar, l, thr].GetN()-1, 0, fitFunc.GetParError(1) )    
                      
                      if (thr == thRef):
                         g_tot_vs_bar[l, vov, thr,module].SetPoint(g_tot_vs_bar[l, vov, thr,module].GetN(), bar, fitFunc.GetParameter(1) )    
                         g_tot_vs_bar[l, vov, thr,module].SetPointError(g_tot_vs_bar[l, vov, thr,module].GetN()-1, 0, fitFunc.GetParError(1) )    
                         # -- energy vs thr, Vov, bar
                      if ( inputFile.GetListOfKeys().Contains('h1_energy_bar%02d%s_Vov%.02f_th%02d'%(bar, l, vov, thr)) ): 
                         h1_energy = inputFile.Get('h1_energy_bar%02d%s_Vov%.02f_th%02d'%(bar, l, vov, thr))
                         if ( h1_energy == None or h1_energy.GetEntries()==0): continue
                         energyPeak = {}
                         if (source == 'Na22'):
                            for peak in peaks:
                               energyPeak[peak] = [h1_energy.GetFunction('fitFunc_%d'%peak).GetParameter(1), h1_energy.GetFunction('fitFunc_%d'%peak).GetParError(1)]
                         elif ( source == 'Laser'):
                            fitFunc = ROOT.TF1('fitFunc','gaus',0, 1000)
                            h1_energy.Fit('fitFunc','QNRS+')
                            h1_energy.Fit('fitFunc','QNS+','', fitFunc.GetParameter(1)-fitFunc.GetParameter(2), fitFunc.GetParameter(1)+fitFunc.GetParameter(2))
                            for peak in peaks:
                               energyPeak[peak] = [ fitFunc.GetParameter(1), fitFunc.GetParError(1)]
                         elif ( source == 'TB'):
                            #fitFunc = h1_energy.GetFunction('f_landau_bar%02d%s_Vov%.02f_vth1_%02d'%(bar, l, vov, thr))
                            #if (fitFunc==None): continue
                            fitFunc = ROOT.TF1('f_landau_bar%02d%s_Vov%.02f_vth1_%02d'%(bar, l, vov, thr), '[0]*TMath::Landau(x,[1],[2])', 0,1000.)
                            if (vov==0.6 or vov==0.8):h1_energy.GetXaxis().SetRangeUser(50,800)
                            else:h1_energy.GetXaxis().SetRangeUser(100,800)
                            emax = h1_energy.GetBinCenter(h1_energy.GetMaximumBin())
                            fitFunc.SetParameters(10, emax, 30)
                            fitFunc.SetRange(0.8*emax, 1.5*emax)
                            h1_energy.Fit(fitFunc,'QR')
                            #print bar, vov, thr, emax, fitFunc.GetParameter(1), fitFunc.GetParError(1)
                            for peak in peaks:
                               energyPeak[peak] = [fitFunc.GetParameter(1), fitFunc.GetParError(1)] 
                                    

                         for peak in peaks:
                            g_energy_vs_th[bar, l, vov, peak,module].SetPoint(g_energy_vs_th[bar, l, vov, peak,module].GetN(), thr, energyPeak[peak][0] )
                            g_energy_vs_th[bar, l, vov, peak,module].SetPointError(g_energy_vs_th[bar, l, vov, peak,module].GetN()-1, 0, energyPeak[peak][1])
                           # print "My print:", module, l, vov, thr, bar 
                            g_energy_vs_vov[bar, l, thr, peak,module].SetPoint(g_energy_vs_vov[bar, l, thr, peak,module].GetN(), VovsEff[vov], energyPeak[peak][0] )
                            g_energy_vs_vov[bar, l, thr, peak,module].SetPointError(g_energy_vs_vov[bar, l, thr, peak,module].GetN()-1, 0, energyPeak[peak][1] )
                            
                            if (thr == thRef):
                              # print g_energy_vs_bar
                              # print "New test happening:", module, vov, l, thr, bar 
                               g_energy_vs_bar[l, vov, thr, peak,module].SetPoint(g_energy_vs_bar[l, vov, thr, peak,module].GetN(), bar, energyPeak[peak][0] )
                               g_energy_vs_bar[l, vov, thr, peak,module].SetPointError(g_energy_vs_bar[l, vov, thr, peak,module].GetN()-1, 0, energyPeak[peak][1] )
   
   # -- tRes vs Vov, thr, bar
for label in label_list:
       module =int(label.split('_')[2][4:7])
       
       inputFile == None
       inputFile = ROOT.TFile.Open(inputdir+'/moduleCharacterization_step2_%s.root'%label)

    
       
       bestRes_totCorr = {}
       bestRes_energyCorr = {}
       bestRes_energyCorr_totCorr = {}
           
       for bar in bars:
          for vov in Vovs:
            for enBin in enBins:    
                    bestRes_totCorr[bar, vov, enBin] = [9999, 9999]
                    bestRes_energyCorr[bar, vov, enBin] = [9999, 9999]
                    bestRes_energyCorr_totCorr[bar, vov, enBin] = [9999, 9999]
                    
            
       for bar in bars:
          for vov in Vovs:
                 if (bar not in goodBars[vov]): continue
                 for thr in thresholds:
                    if (thr not in goodThrs[vov]): continue
                    tRes_energyCorr = {}
                    tRes_totCorr = {}
                    tRes_energyCorr_totCorr = {}
                    for enBin in enBins:
                       h1_deltaT_totCorr    = inputFile.Get('h1_deltaT_totRatioPhaseCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
                       h1_deltaT_energyCorr = inputFile.Get('h1_deltaT_energyRatioPhaseCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
                       h1_deltaT_energyCorr_totCorr = inputFile.Get('h1_deltaT_energyRatioCorr_totRatioCorr_phaseCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
    
                       #h1_deltaT_totCorr    = inputFile.Get('h1_deltaT_totRatioCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
                       #h1_deltaT_energyCorr = inputFile.Get('h1_deltaT_energyRatioCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
                       #h1_deltaT_energyCorr_totCorr = inputFile.Get('h1_deltaT_energyRatioCorr_totRatioCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d'%(bar, vov, thr, enBin))
    
                       # totRatio + phase corr
                       if (h1_deltaT_totCorr == None): continue
                       if (h1_deltaT_totCorr.GetEntries() < 200 ): continue
                       tRes_totCorr[enBin] = getTimeResolution(h1_deltaT_totCorr)
                       if ( tRes_totCorr[enBin][0] < bestRes_totCorr[bar, vov, enBin][0]):
                          bestRes_totCorr[bar, vov, enBin] = tRes_totCorr[enBin]
                       ctemp = ROOT.TCanvas()
                       h1_deltaT_totCorr.GetYaxis().SetRangeUser(0, h1_deltaT_totCorr.GetBinContent(h1_deltaT_totCorr.GetMaximumBin())*1.2)                
                       h1_deltaT_totCorr.Draw()                
                       #ctemp.SaveAs(outdir+'/summaryPlots/timeResolution/fits/'+'/c_h1_deltaT_totRatioCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d.png'%(bar, vov, thr, enBin))
    
    
                       # energyRatio + phase corr
                       if (h1_deltaT_energyCorr == None): continue
                       if (h1_deltaT_energyCorr.GetEntries() < 200 ): continue
                       tRes_energyCorr[enBin] = getTimeResolution(h1_deltaT_energyCorr)
                       if ( tRes_energyCorr[enBin][0] < bestRes_energyCorr[bar, vov, enBin][0] ):
                          #if ( ('1E14' in args.outFolder or  '2E14' in args.outFolder ) and VovsEff[vov] <= 1.50 and thr >= 13): continue
                          bestRes_energyCorr[bar, vov, enBin] = tRes_energyCorr[enBin]
                       ctemp = ROOT.TCanvas()
                       h1_deltaT_energyCorr.GetYaxis().SetRangeUser(0, h1_deltaT_energyCorr.GetBinContent(h1_deltaT_energyCorr.GetMaximumBin())*1.2)                
                       h1_deltaT_energyCorr.Draw()                
                       #ctemp.SaveAs(outdir+'/summaryPlots/timeResolution/fits/'+'/c_h1_deltaT_energyRatioCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d.png'%(bar, vov, thr, enBin))
    
                       # energyRatio + totRatio + phase corr
                       if (h1_deltaT_energyCorr_totCorr == None): continue
                       if (h1_deltaT_energyCorr_totCorr.GetEntries() < 200 ): continue
                       tRes_energyCorr_totCorr[enBin] = getTimeResolution(h1_deltaT_energyCorr_totCorr)
                       if ( tRes_energyCorr_totCorr[enBin][0] < bestRes_energyCorr_totCorr[bar, vov, enBin][0]):
                          #if ( ('1E14' in args.outFolder or  '2E14' in args.outFolder ) and VovsEff[vov] <= 1.50 and thr >= 13): continue
                          bestRes_energyCorr_totCorr[bar, vov, enBin] = tRes_energyCorr_totCorr[enBin]
                       ctemp = ROOT.TCanvas()
                       h1_deltaT_energyCorr_totCorr.GetYaxis().SetRangeUser(0, h1_deltaT_energyCorr_totCorr.GetBinContent(h1_deltaT_energyCorr_totCorr.GetMaximumBin())*1.2)
                       h1_deltaT_energyCorr_totCorr.Draw()                
                       #ctemp.SaveAs(outdir+'/summaryPlots/timeResolution/fits/'+'/c_h1_deltaT_energyRatioCorr_totRatioCorr_bar%02dL-R_Vov%.02f_th%02d_energyBin%02d.png'%(bar, vov, thr, enBin))
    
    
                       # graphs vs threshold
                       g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module].SetPoint(g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module].GetN(), thr, tRes_energyCorr[enBin][0]/kscale )
                       g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module].SetPointError(g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module].GetN()-1, 0, tRes_energyCorr[enBin][1]/kscale)
    
                       g_deltaT_totRatioCorr_vs_th[bar, vov, enBin].SetPoint(g_deltaT_totRatioCorr_vs_th[bar, vov, enBin].GetN(), thr, tRes_totCorr[enBin][0]/kscale )
                       g_deltaT_totRatioCorr_vs_th[bar, vov, enBin].SetPointError(g_deltaT_totRatioCorr_vs_th[bar, vov, enBin].GetN()-1, 0, tRes_totCorr[enBin][1]/kscale)
    
                       g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin].SetPoint(g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin].GetN(), thr, tRes_energyCorr_totCorr[enBin][0]/kscale )
                       g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin].SetPointError(g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin].GetN()-1, 0, tRes_energyCorr_totCorr[enBin][1]/kscale)               
                       
                       # graphs vs Vov
                       g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetPoint(g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].GetN(), VovsEff[vov], tRes_energyCorr[enBin][0]/kscale )
                       g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetPointError(g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].GetN()-1, 0, tRes_energyCorr[enBin][1]/kscale)
    
                       g_deltaT_totRatioCorr_vs_vov[bar, thr, enBin].SetPoint(g_deltaT_totRatioCorr_vs_vov[bar, thr, enBin].GetN(), VovsEff[vov], tRes_totCorr[enBin][0]/kscale )
                       g_deltaT_totRatioCorr_vs_vov[bar, thr, enBin].SetPointError(g_deltaT_totRatioCorr_vs_vov[bar, thr, enBin].GetN()-1, 0, tRes_totCorr[enBin][1]/kscale)
    
                       g_deltaT_energyRatioCorr_totRatioCorr_vs_vov[bar, thr, enBin].SetPoint(g_deltaT_energyRatioCorr_totRatioCorr_vs_vov[bar, thr, enBin].GetN(), VovsEff[vov], tRes_energyCorr_totCorr[enBin][0]/kscale )
                       g_deltaT_energyRatioCorr_totRatioCorr_vs_vov[bar, thr, enBin].SetPointError(g_deltaT_energyRatioCorr_totRatioCorr_vs_vov[bar, thr, enBin].GetN()-1, 0, tRes_energyCorr_totCorr[enBin][1]/kscale)
                       
                       if (thr == thRef):
                         
                          g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin,module].SetPoint(g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin,module].GetN(), bar, tRes_energyCorr[enBin][0]/kscale )
                          g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin,module].SetPointError(g_deltaT_energyRatioCorr_vs_bar[vov, thr, enBin,module].GetN()-1, 0, tRes_energyCorr[enBin][1]/kscale)
    
                          g_deltaT_totRatioCorr_vs_bar[vov, thr, enBin,module].SetPoint(g_deltaT_totRatioCorr_vs_bar[vov, thr, enBin,module].GetN(), bar, tRes_totCorr[enBin][0]/kscale )
                          g_deltaT_totRatioCorr_vs_bar[vov, thr, enBin,module].SetPointError(g_deltaT_totRatioCorr_vs_bar[vov, thr, enBin,module].GetN()-1, 0, tRes_totCorr[enBin][1]/kscale)
                        
                 # graphs at best Th
                 for enBin in enBins:
                       
                    # -- tRes  vs Vov at the best threshold    
                    if (bestRes_energyCorr[bar, vov, enBin][0]!= 9999):
                       g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetPoint(g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].GetN(), VovsEff[vov], (bestRes_energyCorr[bar, vov, enBin][0])/kscale )
                       g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetPointError(g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].GetN()-1, 0, (bestRes_energyCorr[bar, vov, enBin][1])/kscale)
    
                    if (bestRes_totCorr[bar, vov, enBin][0]!= 9999): 
                       g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin].SetPoint(g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin].GetN(), VovsEff[vov], (bestRes_totCorr[bar, vov, enBin][0])/kscale )
                       g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin].SetPointError(g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin].GetN()-1, 0, (bestRes_totCorr[bar, vov, enBin][1])/kscale)
                       
                    if (bestRes_energyCorr_totCorr[bar, vov, enBin][0]!= 9999):              
                       g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin].SetPoint(g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin].GetN(), VovsEff[vov], (bestRes_energyCorr_totCorr[bar, vov, enBin][0])/kscale )
                       g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin].SetPointError(g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin].GetN()-1, 0, (bestRes_energyCorr_totCorr[bar, vov, enBin][1])/kscale)
    
                    # -- tRes  vs bar at the best threshold
                    if (bestRes_energyCorr[bar, vov, enBin][0]!= 9999): 
                       g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module].SetPoint(g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module].GetN(), bar, (bestRes_energyCorr[bar, vov, enBin][0])/kscale )
                       g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module].SetPointError(g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module].GetN()-1, 0, (bestRes_energyCorr[bar, vov, enBin][1])/kscale)            
    
                    if (bestRes_totCorr[bar, vov, enBin][0]!= 9999): 
                       g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetPoint(g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].GetN(), bar, (bestRes_totCorr[bar, vov, enBin][0])/kscale )
                       g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetPointError(g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].GetN()-1, 0, (bestRes_totCorr[bar, vov, enBin][1])/kscale)            
                        
                    if (bestRes_energyCorr_totCorr[bar, vov, enBin][0]!= 9999):
                       g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetPoint(g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].GetN(), bar, (bestRes_energyCorr_totCorr[bar, vov, enBin][0])/kscale )
                       g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetPointError(g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].GetN()-1, 0, (bestRes_energyCorr_totCorr[bar, vov, enBin][1])/kscale)            
                       

# fill graphs with average tRes vs OV
for module in Module:
 for enBin in enBins:
   for vov in Vovs:
      fitpol0 = ROOT.TF1('fitpol0','pol0',-1,16) 
      g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].Fit(fitpol0,'QSN')
      ave, err = [fitpol0.GetParameter(0),fitpol0.GetParError(0)]
      g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin].SetPoint(g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin].GetN(), VovsEff[vov], ave)
      g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin].SetPointError(g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin].GetN()-1, 0, err)

      g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module].Fit(fitpol0,'QSN')
      ave, err = [fitpol0.GetParameter(0),fitpol0.GetParError(0)]
      g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin].SetPoint(g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin].GetN(), VovsEff[vov], ave)
      g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin].SetPointError(g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin].GetN()-1, 0, err)

      g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].Fit(fitpol0,'QSN')
      ave, err = [fitpol0.GetParameter(0),fitpol0.GetParError(0)]
      g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin].SetPoint(g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin].GetN(), VovsEff[vov], ave)
      g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin].SetPointError(g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin].GetN()-1, 0, err)
      

# -- Draw
latex = ROOT.TLatex(0.18,0.94,'%s'%(plots_label))
latex.SetNDC()
latex.SetTextSize(0.05)
latex.SetTextFont(42)

for bar in bars:
    ctot1 = ROOT.TCanvas('c_tot_vs_th_bar%.02d'%bar)
    hPad1 = ROOT.TH2F('hPad1','', 100, -1., 64.,40, 0.,40.)
    hPad1.SetTitle(";threshold [DAC];ToT [ns]")
    hPad1.Draw()
    ctot1.SetGridy()
    leg = ROOT.TLegend(0.70, 0.70, 0.89, 0.89)  
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
   
    for i, vov in enumerate(Vovs):

                       
            for l in ['L','R']:
                g_tot_vs_th[bar, l, vov].Sort()
                g_tot_vs_th[bar, l, vov].SetMarkerStyle(20)
                if (l == 'R'): g_tot_vs_th[bar, l, vov].SetMarkerStyle(24)
                g_tot_vs_th[bar, l, vov].SetMarkerColor(i+1)
                g_tot_vs_th[bar, l, vov].SetLineColor(i+1)
                g_tot_vs_th[bar, l, vov].Draw('plsame')
            leg.AddEntry(g_tot_vs_th[bar,'L', vov], 'V_{OV}^{eff} = %.02f V'%VovsEff[vov], 'PL')
    leg.Draw()    
    latex.Draw('same')
    ctot1.SaveAs(outdir+'/summaryPlots/tot/'+ctot1.GetName()+'.png')
    ctot1.SaveAs(outdir+'/summaryPlots/tot/'+ctot1.GetName()+'.pdf')
    hPad1.Delete()
    
    # -- tot vs Vov
    ctot2 = ROOT.TCanvas('c_tot_vs_Vov_bar%.02d'%bar)
    hPad2 = ROOT.TH2F('hPad2','', 10, 0., vovMax, 40, 0.,40.)
    hPad2.SetTitle(";V_{OV}^{eff} [V];ToT [ns]")
    hPad2.Draw()
    ctot2.SetGridy()
    leg = ROOT.TLegend(0.70, 0.70, 0.89, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, thr in enumerate(thresholds):
        for l in ['L','R']:
            g_tot_vs_vov[bar,l, thr].Sort()
            g_tot_vs_vov[bar,l, thr].SetMarkerStyle(20)
            if (l == 'R'): g_tot_vs_vov[bar, l, thr].SetMarkerStyle(24)
            g_tot_vs_vov[bar,l, thr].SetMarkerColor(i+1)
            g_tot_vs_vov[bar,l, thr].SetLineColor(i+1)
            g_tot_vs_vov[bar,l, thr].Draw('plsame')
        leg.AddEntry(g_tot_vs_vov[bar,'L', thr], 'th. = %d'%thr, 'PL')
    leg.Draw()    
    ctot2.SaveAs(outdir+'/summaryPlots/tot/'+ctot2.GetName()+'.png')
    ctot2.SaveAs(outdir+'/summaryPlots/tot/'+ctot2.GetName()+'.pdf')
    hPad2.Delete()
    
    # -- energy vs threshold
for bar in bars:
    cen1 = ROOT.TCanvas('c_energy_vs_th_bar%.02d'%bar)
    hPadEn1 = ROOT.TH2F('hPadEn1','Bar', 100, -1., 30.,500, 0.,800.)
    hPadEn1.GetYaxis().SetTitle('energy')
    hPadEn1.GetXaxis().SetTitle('threshold [DAC]')
    hPadEn1.Draw()
    cen1.SetGridy()
 
    leg = ROOT.TLegend(0.50, 0.50, 0.99, 0.99)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
   # leg.SetNColumns(4) 
    for i, vov in enumerate(Vovs):
        for module in Module:
            for l in ['L','R']:
                g_energy_vs_th[bar, l, vov, refPeak, module].Sort()
                g_energy_vs_th[bar, l, vov, refPeak, module].SetMarkerColor(mols[module]) 
                if (l == 'R'): g_energy_vs_th[bar, l, vov, refPeak,module].SetMarkerColor(mols[module]+1) 
                g_energy_vs_th[bar, l, vov, refPeak,module].SetLineColor(1)
                
                leg.AddEntry(g_energy_vs_th[bar,l, vov, refPeak,module], '%s,%s, V_{OV}^{eff} = %.02f V'%(module,l,VovsEff[vov]), 'PL')
    leg.SetHeader('Bar %s'%bar)
    leg.Draw()
    cen1.SaveAs(outdir+'/summaryPlots/energy/'+cen1.GetName()+'.png')
    cen1.SaveAs(outdir+'/summaryPlots/energy/'+cen1.GetName()+'.pdf')
hPadEn1.Delete()

# -- energy vs Vov    (needs to be fixed)
for bar in bars:
    cen2 = ROOT.TCanvas('c_energy_vs_Vov_bar%.02d'%bar)
    hPadEn2 = ROOT.TH2F('hPadEn2','', 10, 0., vovMax,50, 0.,1000.)
    hPadEn2.SetTitle(";V_{OV}^{eff} [V]; energy")
    hPadEn2.Draw()
    cen2.SetGridy()
    leg = ROOT.TLegend(0.60, 0.50, 0.99, 0.99)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for i, thr in enumerate(thresholds):
        for module in Module:
            for l in ['L','R']:
                g_energy_vs_vov[bar,l, thr, refPeak,module].Sort()
                g_energy_vs_vov[bar,l, thr, refPeak,module].SetMarkerStyle(20)
                if (l == 'R'): g_energy_vs_vov[bar, l, thr, refPeak,module].SetMarkerStyle(24)
                g_energy_vs_vov[bar,l, thr, refPeak,module].SetMarkerColor(i+1)
                g_energy_vs_vov[bar,l, thr, refPeak,module].SetLineColor(i+1)
                g_energy_vs_vov[bar,l, thr, refPeak,module].Draw('plsame')
                leg.AddEntry(g_energy_vs_vov[bar,l, thr, refPeak,module], '%s, %s, thr= %.02f DAC'%(module,l,thr), 'PL')
    leg.SetHeader('Bar %s'%bar)
           
    leg.Draw()
    cen2.SaveAs(outdir+'/summaryPlots/energy/'+cen2.GetName()+'.png')
    cen2.SaveAs(outdir+'/summaryPlots/energy/'+cen2.GetName()+'.pdf') 
hPadEn2.Delete()   


    # -- time resolution vs threshold(next!)
for bar in bars:
    for enBin in enBins:
        ctres1 = ROOT.TCanvas('c_tRes_energyRatioCorr_vs_th_bar%.02d_enBin%02d'%(bar,enBin))
        hPadT1 = ROOT.TH2F('hPadT1','', 100, -1., 64.,100, tResMin,tResMaxTh)
        hPadT1.SetTitle(";threshold [DAC]; #sigma_{t} [ps]")
        hPadT1.Draw()
        ctres1.SetGridy()
        leg = ROOT.TLegend(0.43, 0.55, 0.89, 0.99)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        for module in Module:
            for i, vov in enumerate(Vovs):
                
               # g_deltaT_totRatioCorr_vs_th[bar, vov, enBin].SetMarkerStyle(mols[module])
               # g_deltaT_totRatioCorr_vs_th[bar, vov, enBin].SetMarkerColor(cols[vov])
               # g_deltaT_totRatioCorr_vs_th[bar, vov, enBin].SetLineColor(1)
               # g_deltaT_totRatioCorr_vs_th[bar, vov, enBin].Draw('plsame')
                g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module].SetMarkerStyle(cols[vov])
                g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module].SetMarkerColor(mols[module])
               # g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin].SetLineColor(2)
                g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module].Draw('plsame')
               # g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin].SetMarkerStyle(mols[module])
               # g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin].SetMarkerColor(cols[vov])
               # g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin].SetLineColor(3)
               # g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin].Draw('plsame')
               # leg.AddEntry(g_deltaT_totRatioCorr_vs_th[bar, vov, enBin], '%s,%s, V_{OV}^{eff} = %.02f V'%(module,"totRatioCorr",VovsEff[vov]), 'PL')
                leg.SetHeader('Bar %s'%bar)
                leg.AddEntry(g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module], '%s,%s, V_{OV}^{eff} = %.02f V'%(module,"energyRatioCorr",VovsEff[vov]), 'PL')
               # leg.AddEntry(g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin], '%s, V_{OV}^{eff} = %.02f V'%("energytotRatioCorr",VovsEff[vov]), 'PL')
                
                
    
                outfile.cd()
                g_deltaT_energyRatioCorr_vs_th[bar, vov, enBin,module].Write('g_deltaT_energyRatioCorr_vs_th_bar%02d_Vov%.02f_enBin%02d'%(bar, vov, enBin))
               # g_deltaT_totRatioCorr_vs_th[bar, vov, enBin].Write('g_deltaT_totRatioCorr_vs_th_bar%02d_Vov%.02f_enBin%02d'%(bar, vov, enBin))
               # g_deltaT_energyRatioCorr_totRatioCorr_vs_th[bar, vov, enBin].Write('g_deltaT_energyRatioCorr_totRatioCorr_vs_th_bar%02d_Vov%.02f_enBin%02d'%(bar, vov, enBin))
        leg.Draw()
        latex.Draw('same')
        ctres1.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres1.GetName()+'.png')
        ctres1.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres1.GetName()+'.pdf')
        hPadT1.Delete()
    
    
            # -- time resolution vs Vov
for bar in bars:
        ctres2 = ROOT.TCanvas('c_tRes_energyRatioCorr_vs_Vov_bar%.02d_enBin%02d'%(bar,enBin))
        hPadT2 = ROOT.TH2F('hPadT2','', 6, 0., vovMax, 10, tResMin,tResMax)
        hPadT2.SetTitle(";V_{OV}^{eff} [V];#sigma_{t} [ps]")
        hPadT2.Draw()
        ctres2.SetGridy()
        leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        for i, thr in enumerate(thresholds):
            g_deltaT_totRatioCorr_vs_vov[bar, thr, enBin].SetMarkerStyle(20)
            g_deltaT_totRatioCorr_vs_vov[bar, thr, enBin].SetMarkerColor(i+1)
            g_deltaT_totRatioCorr_vs_vov[bar, thr, enBin].SetLineColor(i+1)
            g_deltaT_totRatioCorr_vs_vov[bar, thr, enBin].Draw('plsame')
            g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetMarkerStyle(24)
            g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetMarkerColor(i+1)
            g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].SetLineColor(i+1)
            g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin].Draw('plsame')
            g_deltaT_energyRatioCorr_totRatioCorr_vs_vov[bar, thr, enBin].SetMarkerStyle(34)
            g_deltaT_energyRatioCorr_totRatioCorr_vs_vov[bar, thr, enBin].SetMarkerColor(i+1)
            g_deltaT_energyRatioCorr_totRatioCorr_vs_vov[bar, thr, enBin].SetLineColor(i+1)
            g_deltaT_energyRatioCorr_totRatioCorr_vs_vov[bar, thr, enBin].Draw('plsame')
            leg.AddEntry(g_deltaT_energyRatioCorr_vs_vov[bar, thr, enBin], 'th. = %d'%thr, 'PL')
            outfile.cd()
            g_deltaT_totRatioCorr_vs_vov[bar,thr, enBin].Write('g_deltaT_totRatioCorr_vs_vov_bar%02d_th%02d_enBin%02d'%(bar, thr, enBin)) 
            g_deltaT_energyRatioCorr_vs_vov[bar,thr, enBin].Write('g_deltaT_energyRatioCorr_vs_vov_bar%02d_th%02d_enBin%02d'%(bar, thr, enBin)) 
            g_deltaT_energyRatioCorr_vs_vov[bar,thr, enBin].Write('g_deltaT_energyRatioCorr_totRatioCorr_vs_vov_bar%02d_th%02d_enBin%02d'%(bar, thr, enBin)) 
        leg.SetHeader('Bar %s'%bar)

        leg.Draw()
        latex.Draw('same')
        ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.png')
        ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.pdf')
        hPadT2.Delete()   

        # -- time resolution vs Vov at the best Th
        ctres2 = ROOT.TCanvas('c_tRes_energyRatioCorr_bestTh_vs_Vov_bar%.02d_enBin%02d'%(bar,enBin))
        hPadT2 = ROOT.TH2F('hPadT2','', 6, 0., vovMax,10, tResMin,tResMax)
        hPadT2.SetTitle(";V_{OV}^{eff} [V];#sigma_{t} [ps]")
        hPadT2.Draw()
        leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        ctres2.SetGridy()
       # g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin].SetMarkerStyle(20)
       # g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin].SetMarkerColor(1)
       # g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin].SetLineColor(1)
       # g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin].Draw('plsame')
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetMarkerStyle(24)
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetMarkerColor(2)
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].SetLineColor(2)
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].Draw('plsame')
       # g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin].SetMarkerStyle(34)
       # g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin].SetMarkerColor(4)
       # g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin].SetLineColor(4)
       # g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin].Draw('plsame')
       # leg.AddEntry(g_deltaT_totRatioCorr_bestTh_vs_vov[bar,enBin], '%s'%"totRatioCorr")
        leg.AddEntry(g_deltaT_energyRatioCorr_bestTh_vs_vov[bar,enBin], '%s'%"energyRatioCorr")
       # leg.AddEntry(g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar,enBin], '%s'%"energytotRatioCorr")
        leg.SetHeader('Bar %s'%bar)
        leg.Draw()
        latex.Draw('same')
        outfile.cd() 
        g_deltaT_energyRatioCorr_bestTh_vs_vov[bar, enBin].Write('g_deltaT_energyRatioCorr_bestTh_vs_vov_bar%02d_enBin%02d'%(bar, enBin))
       # g_deltaT_totRatioCorr_bestTh_vs_vov[bar, enBin].Write('g_deltaT_totRatioCorr_bestTh_vs_vov_bar%02d_enBin%02d'%(bar, enBin))
       # g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov[bar, enBin].Write('g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_bar%02d_enBin%02d'%(bar, enBin))
        ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.png')
        ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.pdf')
        hPadT2.Delete()   

# average time resolution vs OV
for enBin in enBins:
   ctres2 = ROOT.TCanvas('c_tRes_totRatioCorr_bestTh_vs_Vov_enBin%02d_average'%(enBin))
   hPadT2 = ROOT.TH2F('hPadT2','', 6, 0.0, vovMax,10, tResMin,tResMax)
   hPadT2.SetTitle(";V_{OV}^{eff} [V];#sigma_{t} [ps]")
   hPadT2.Draw()
   ctres2.SetGridy()
   g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin].SetMarkerStyle(20)
   g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin].SetMarkerColor(1)
   g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin].SetLineColor(1)
   g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin].Draw('plsame')
   g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin].SetMarkerStyle(24)
   g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin].SetMarkerColor(1)
   g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin].SetLineColor(1)
   g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin].Draw('plsame')
   g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin].SetMarkerStyle(34)
   g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin].SetMarkerColor(1)
   g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin].SetLineColor(1)
   g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin].Draw('plsame')
   latex.Draw('same')
   outfile.cd() 
   g_deltaT_energyRatioCorr_bestTh_vs_vov_average[enBin].Write('g_deltaT_energyRatioCorr_bestTh_vs_vov_enBin%02d_average'%(enBin))
   g_deltaT_totRatioCorr_bestTh_vs_vov_average[enBin].Write('g_deltaT_totRatioCorr_bestTh_vs_vov_enBin%02d_average'%(enBin))
   g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_average[enBin].Write('g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_vov_enBin%02d_average'%(enBin))
   ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.png')
   ctres2.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres2.GetName()+'.pdf')
   hPadT2.Delete()   
   
# -- plots vs bar(done!)

# -- tot vs bar(not focussing)
for module in Module: 
 for i, vov in enumerate(Vovs):
   ctot3 = ROOT.TCanvas('c_tot_vs_bar_Vov%.02f'%vov)
   hPad3 = ROOT.TH2F('hPad3','', 100, -0.5, 15.5,40, 0.,40.)
   hPad3.SetTitle("; bar; ToT [ns]")
   hPad3.Draw()
   ctot3.SetGridy()
   leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
   leg.SetBorderSize(0)
   leg.SetFillStyle(0)
   for l in ['L','R']:
      g_tot_vs_bar[l, vov, thRef,module].Sort()
      g_tot_vs_bar[l, vov, thRef,module].SetMarkerStyle(20)
      if (l == 'R'): g_tot_vs_bar[l, vov, thRef,module].SetMarkerStyle(24)
      g_tot_vs_bar[l, vov, thRef,module].SetMarkerColor(cols[vov])
      g_tot_vs_bar[l, vov, thRef,module].SetLineColor(cols[vov])
      g_tot_vs_bar[l, vov, thRef,module].Draw('plsame')
      leg.AddEntry(g_tot_vs_bar['L', vov, thRef,module], '%s, V_{OV}^{eff} = %.02f V'%(module,VovsEff[vov]), 'PL')
   leg.Draw()
   ctot3.SaveAs(outdir+'/summaryPlots/tot/'+ctot3.GetName()+'.png')
   ctot3.SaveAs(outdir+'/summaryPlots/tot/'+ctot3.GetName()+'.pdf')    
   hPad3.Delete()

# -- energy vs bar

for i, vov in enumerate(Vovs):
   cen3 = ROOT.TCanvas('c_energy_vs_bar_Vov%.02f'%vov)
   hPadEn3 = ROOT.TH2F('hPadEn3','', 100, -0.5, 15.5, 40, 0.,1000.)
   hPadEn3.SetTitle("; bar; energy")
   hPadEn3.Draw()
   cen3.SetGridy()
   leg = ROOT.TLegend(0.65, 0.50, 0.89, 0.89)
   leg.SetBorderSize(0)
   leg.SetFillStyle(0)
   for module in Module:
    
     for l in ['L','R']:
       g_energy_vs_bar[l, vov, thRef, refPeak,module].Sort() 
       g_energy_vs_bar[l, vov, thRef, refPeak,module].SetMarkerColor(mols[module])
       if (l == 'R'): g_energy_vs_bar[l, vov, thRef, refPeak,module].SetMarkerColor(mols[module]+20)
       g_energy_vs_bar[l, vov, thRef, refPeak,module].Draw('plsame')
       outfile.cd()  
       g_energy_vs_bar[l, vov, thRef, refPeak,module].Write('g_energy%s_vs_bar_Vov%.02f_th%02d'%(l,vov,thRef))
       leg.AddEntry(g_energy_vs_bar[l, vov, thRef, refPeak,module],'%s, %s, V_{OV}^{eff} = %.02f V'%(module,l,VovsEff[vov]),'PL')
   
   leg.Draw()
   cen3.SaveAs(outdir+'/summaryPlots/energy/'+cen3.GetName()+'.png')
   cen3.SaveAs(outdir+'/summaryPlots/energy/'+cen3.GetName()+'.pdf')    
   hPadEn3.Delete()
    
# -- time resolution vs bar at the ref threshold (done!)
for vov in Vovs:
 for enBin in enBins:
   ctres31 = ROOT.TCanvas('c_tRes_energyRatioCorr_refTh_vs_bar_Vov%.02f_enBin%02d'%(vov,enBin))
   hPadt31 = ROOT.TH2F('hPadT3','', 100, -0.5, 15.5,100, tResMin,tResMax)
   hPadt31.SetTitle("; bar; #sigma_{t}[ps]")
   hPadt31.Draw()
   ctres31.SetGridy()
   leg = ROOT.TLegend(0.70, 0.96, 0.90, 0.63)
   leg.SetBorderSize(0)
   leg.SetFillStyle(0)
   for module in Module:
      g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin,module].SetMarkerColor(mols[module])
      g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin,module].Draw('psame')
      outfile.cd() 
      g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin,module].Write('g_deltaT_energyRatioCorr_vs_bar__Vov%.02f_th%02d'%(vov,thRef))
      leg.AddEntry(g_deltaT_energyRatioCorr_vs_bar[vov, thRef, enBin,module],'%s , V_{OV}^{eff} = %.02f V'%(module,VovsEff[vov]), 'PL')
      leg.SetHeader('refVth=%s DAC'%thRef)
      leg.Draw()
      latex.Draw('same')
      ctres31.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres31.GetName()+'.png')
      ctres31.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres31.GetName()+'.pdf')    
 hPadt31.Delete()  

# -- time resolution vs bar at the best threshold
for enBin in enBins:
    
       ctres3 = ROOT.TCanvas('c_tRes_energyRatioCorr_bestTh_vs_bar_enBin%02d'%(enBin))
       hPadT3 = ROOT.TH2F('hPadT3','', 100, -0.5, 15.5,100, tResMin,tResMax)
       hPadT3.SetTitle("; bar; #sigma_{t} [ps]")
       hPadT3.Draw()
       ctres3.SetGridy()
       #leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
       leg = ROOT.TLegend(0.70, 0.90, 0.60, 0.70)
       leg.SetBorderSize(0)
       leg.SetFillStyle(0)
       if (len(Vovs)>4): 
          leg.SetNColumns(2);
          leg.SetColumnSeparation(0.2);
       for module in Module:
                                
           for i, vov in enumerate(Vovs):
             
             # g_deltaT_energyRatioCorr_bestTh_vs_bar[vov,enBin,module].SetMarkerColor(8)
              g_deltaT_energyRatioCorr_bestTh_vs_bar[vov,enBin,module].SetMarkerStyle(cols[vov])
              g_deltaT_energyRatioCorr_bestTh_vs_bar[vov,enBin,module].SetMarkerColor(mols[module])
             
            
              g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module].Draw('psame')
              fitRes = ROOT.TF1('fitRes','pol0',0,16)
              g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module].Fit(fitRes,'QRN')
              
              leg.AddEntry(g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module],'%s ,V_{OV}^{eff} = %.02f V'%(module,VovsEff[vov]), 'PL')
              outfile.cd()
              g_deltaT_energyRatioCorr_bestTh_vs_bar[vov, enBin,module].Write('g_deltaT_energyRatioCorr_bestTh_vs_bar_Vov%.02f_enBin%02d'%(vov, enBin))
          
           leg.Draw()
           latex.Draw('same')
           ctres3.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.png')
           ctres3.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.pdf')    
       hPadT3.Delete()
   #continue... 
        
for enBin in enBins: 
       ctres30 = ROOT.TCanvas('c_tRes_totRatioCorr_bestTh_vs_bar_enBin%02d'%(enBin))
       hPadT4 = ROOT.TH2F('hPadT4','', 100, -0.5, 15.5,100, tResMin,tResMax)
       hPadT4.SetTitle("; bar; #sigma_{t} [ps]")
       hPadT4.Draw()
       ctres30.SetGridy()
       #leg = ROOT.TLegend(0.70, 0.50, 0.89, 0.89)
       leg = ROOT.TLegend(0.20, 0.90, 0.60, 0.70)
       leg.SetBorderSize(0)
       leg.SetFillStyle(0)
       if (len(Vovs)>4): 
          leg.SetNColumns(2);
          leg.SetColumnSeparation(0.2);
       for module in Module:
                                
           for i, vov in enumerate(Vovs):
              g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetMarkerStyle(cols[vov])
              g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetMarkerColor(mols[module])
             # g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetLineColor(cols[vov])
              g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].Draw('psame')
              fitRes = ROOT.TF1('fitRes','pol0',0,16)
              g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].Fit(fitRes,'QRN')
             
              leg.AddEntry(g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module], '%s, V_{OV}^{eff} = %.02f V'%(module,VovsEff[vov]), 'PL')
              outfile.cd()
              g_deltaT_totRatioCorr_bestTh_vs_bar[vov, enBin,module].Write('g_deltaT_totRatioCorr_bestTh_vs_bar_Vov%.02f_enBin%02d'%(vov, enBin))
           leg.Draw()
           latex.Draw('same')
           ctres30.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.png')
           ctres30.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.pdf')    
       hPadT4.Delete()  
        
for enBin in enBins: 
       ctres3 = ROOT.TCanvas('c_tRes_energyRatioCorr_totRatioCorr_bestTh_vs_bar_enBin%02d'%(enBin))
       hPadT5 = ROOT.TH2F('hPadT5','', 100, -0.5, 15.5,100, tResMin,tResMax)
       hPadT5.SetTitle("; bar; #sigma_{t} [ps]")
       hPadT5.Draw()
       ctres3.SetGridy()
       
       leg = ROOT.TLegend(0.20, 0.90, 0.60, 0.70)
       leg.SetBorderSize(0)
       leg.SetFillStyle(0)
       if (len(Vovs)>4): 
          leg.SetNColumns(2);
          leg.SetColumnSeparation(0.2);
       for module in Module:
                                
           for i, vov in enumerate(Vovs):
              g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetMarkerStyle(cols[vov])
              g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetMarkerColor(mols[module])
             # g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].SetLineColor(cols[vov])
              g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].Draw('psame')
              fitRes = ROOT.TF1('fitRes','pol0',0,16)
              g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].Fit(fitRes,'QRN')
              
              leg.AddEntry(g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module], '%s, V_{OV}^{eff} = %.02f V'%(module,VovsEff[vov]), 'PL')
              outfile.cd()
              g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar[vov, enBin,module].Write('g_deltaT_energyRatioCorr_totRatioCorr_bestTh_vs_bar_Vov%.02f_enBin%02d'%(vov, enBin))
           leg.Draw()
           latex.Draw('same')
           ctres3.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.png')
           ctres3.SaveAs(outdir+'/summaryPlots/timeResolution/'+ctres3.GetName()+'.pdf')    
       hPadT5.Delete()  
        
outfile.Close()
    
    #raw_input('OK?')
