#include "TH1F.h"
#include "TFile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include <iostream>

Double_t fcn1(Double_t* x, Double_t* par)
  {
    const double sqrtPiOver2 = 1.2533141373;
    const double sqrt2 = 1.4142135624;
    double sig = fabs((double) par[1]);
    double t = (x[0] - par[0])/sig;
    if(par[2] < 0)
      t = -t;
    double absAlpha = fabs(par[2]/sig);
    double a = TMath::Power(par[3]/absAlpha,par[3])*exp(-0.5*absAlpha*absAlpha);
    double b = absAlpha - par[3]/absAlpha;
    double ApproxErf;
    double arg = absAlpha / sqrt2;
    if (arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    double leftArea = (1 + ApproxErf) * sqrtPiOver2;
    double rightArea = ( a * 1/TMath::Power(absAlpha - b,par[3]-1)) / (par[3] - 1);
    double area = leftArea + rightArea;
    if( t <= absAlpha ){
      arg = t / sqrt2;
      if(arg > 5.) ApproxErf = 1;
      else if (arg < -5.) ApproxErf = -1;
      else ApproxErf = TMath::Erf(arg);
      return par[4] * (1 + ApproxErf) * sqrtPiOver2 / area;
    }
    else{
      return par[4] * (leftArea + a * (1/TMath::Power(t-b,par[3]-1) -
				       1/TMath::Power(absAlpha - b,par[3]-1)) / (1 - par[3])) / area;
    }
  }
  



void jetTriggerEff(TString inputFile, TString saveName){
  gROOT->LoadMacro("./tdrstyle.C");
  setTDRStyle();
  tdrGrid(1);
  gStyle->SetOptFit(0);
  //load fitter 

  gROOT->LoadMacro("./logFitTest.C");
  TFile* theFile = new TFile(inputFile);
  TH1F* all = (TH1F*)theFile->Get("denominatorJet");
  TH1F* passedL1 = (TH1F*)theFile->Get("numeratorJet");
  
  TCanvas* theCanvas = new TCanvas();
  //cout << " *************************** Fit L1 *************************** " << endl;
  TGraphAsymmErrors* theL1Eff = new TGraphAsymmErrors(passedL1, all, "cl=0.683 b(1,1) mode");
  theL1Eff->SetTitle(0);
  theL1Eff->GetYaxis()->SetTitle("L2 Efficiency");
  theL1Eff->GetXaxis()->SetTitle("PFJet p_{T}");
  theL1Eff->SetMinimum(0);
  theL1Eff->SetMaximum(1.125);
  theL1Eff->Draw("AP");
  
  TF1* l1Fit = logFitTest(2,passedL1,all);
  l1Fit->SetLineWidth(2);
  l1Fit->Draw("same");
  
//  TF1* myfitfunc = new TF1("myfitfunc",fcn1,0,250,5);
//  //myfitfunc->SetParameters(30,1.5,0.5,10,1.0);
//  myfitfunc->SetParameters(30.,6.5,3.5,100,1.0);
//  theL1Eff->Fit("myfitfunc");
//  myfitfunc->SetLineColor(3);
//  myfitfunc->SetLineWidth(2);
//  myfitfunc->Draw("same");

  theL1Eff->SetMarkerStyle(20);
  theL1Eff->Draw("P");

  theCanvas->SaveAs(saveName);
   
  //cout << " *************************** Parameters *************************** " << endl;
    cout << "\n Fit:\t [0]:" << l1Fit->GetParameter(0) 
           << "\t [1]:" << l1Fit->GetParameter(1) 
           << "\t [2]:" << l1Fit->GetParameter(2)
     <<endl;
}
