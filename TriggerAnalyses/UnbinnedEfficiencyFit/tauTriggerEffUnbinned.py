#!/usr/bin/env python

import imp
import sys
import os
import subprocess

import math
from array import array
import ROOT

from tdrstyle import *
from turnOnHelperFunctions import *

# Settings:
startFit = 20
endFit = 1000
#filename = "/user/kargoll/results/TurnOns/forMoriond2013/DoubleTauJet/turnOnUnbinned_Run2012A_rerunH2TauHLTPaths_MisomvaMelecrej_L3tionProng4_ptThr-44-64-25-25_dZ0.2_MT40_smB.root"
periods = ["2012A","2012B","2012Cv3","2012Cv4","2012Dpart1","2012Dpart2"]
levels = ["L1","L2","L2p5","L3","L1L2","L1L2L2p5","L1L2L2p5L3"]
etaBins = ["Total","Central","Forward"]

# RUN THE MACRO
drawGrid = False
setTDRStyle(drawGrid)

canvas=ROOT.TCanvas("canvas","myCanvas")

resultString = ""

for period in periods:
    filename = "/user/kargoll/results/TurnOns/forMoriond2013/DoubleTauJet/Unbinned/turnOnUnbinned_Run"+period+".root"
    file = ROOT.TFile(filename,"READ")
    
    resultString = resultString + "\nHere are the results for "+period+":\n"
    
    for etaBin in etaBins:
        for level in levels:
            graph = ROOT.TGraph(file.Get("eff"+level+"_"+etaBin))
            graph.SetMarkerStyle(2)
            graph.SetTitle("")
            graph.GetYaxis().SetRangeUser(-0.1,1.1)
            graph.GetYaxis().SetTitle(level + " Efficiency")
            graph.GetXaxis().SetTitle("PFTau E_{T} / GeV")
            graph.Draw("AP")
            
            text = ROOT.TPaveText(0.5,0.3,0.9,0.5,"NDC")
            text.SetFillColor(0)
            text.SetShadowColor(0)
            text.SetTextAlign(11)
            text.SetTextFont(42)
            text.AddText(etaBin+" "+level+" eff. in "+period+":")
            
            # fit using chi2 method
#            fit = ROOT.TF1("fit",errorfit,startFit,endFit,3)
#            fit.SetParameters(90.,35.,0.8)
#            fit.SetParLimits(0,0.,1.)
#            graph.Fit("fit","RNWSQ")
#            fit.SetLineColor(ROOT.kRed)
#            fit.Draw("same")
            
            # fit using unbinned maximum likelihood
            myUnbinnedFit = unbinnedFit(graph,startFit)
            fitResults = myUnbinnedFit.fitUnbinned()
            fit = ROOT.TF1("fit",errorfit,startFit,endFit,3)
            fit.SetParameters(fitResults[0],fitResults[1],fitResults[2])
            fit.SetParError(0,fitResults[3])
            fit.SetParError(1,fitResults[4])
            fit.SetParError(2,fitResults[5])
            fit.SetLineColor(ROOT.kRed)
            fit.Draw("same")
            
            #text.AddText('{0:10} = {1:.4f}'.format("chi2",fit.GetChisquare()))
            #text.AddText('{0:10} = {1}'.format("ndf",fit.GetNDF()))
            text.AddText('{0:10} = {1:.4f} +/- {2:.4f}'.format("plateau",fit.GetParameter(0),fit.GetParError(0)))
            text.AddText('{0:10} = {1:.4f} +/- {2:.4f}'.format("m0",fit.GetParameter(1),fit.GetParError(1)))
            text.AddText('{0:10} = {1:.4f} +/- {2:.4f}'.format("sigma",fit.GetParameter(2),fit.GetParError(2)))
            graph.GetListOfFunctions().Add(text)
            
            graph.Draw("P")
            
            canvas.Update()
            
            #canvas.Print("test.png")
            
            output = "/user/kargoll/CMSSW/CMSSW_5_3_3_patch3/src/H2TauTauJet/DoubleTauJetResults/Moriond2013/Unbinned/turnOnTau_UnbinnedlogL_"\
                        +period+"_"+level+"_"+etaBin+"_fitFrom"+str(startFit)
            epsout = output+".eps"
            pdfout = output+".pdf"
            canvas.Print(epsout)
            subprocess.check_call("epstopdf --outfile="+pdfout+" "+epsout,shell=True)
            
    
            resultString = resultString + \
            '{0:7} : {1:8} : {2}*0.5*(TMath::Erf((pt-{3})/2./{4}/sqrt(pt))+1.))\n'\
            .format(etaBin,level,fit.GetParameter(0),fit.GetParameter(1),fit.GetParameter(2))

    
print resultString
