
#include "plotCommon.hpp"
#include "littleHelpers.hpp"

#include <map>
#include <vector>
#include <iostream>

#include <TH1F.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include "TMultiGraph.h"

#include <sstream> 
#include <THStack.h>
#include <TCanvas.h>
#include "TRandom.h"
#include "TDatime.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TUnuran.h"

//======================================================================
void
Wahrscheinlichkeiten(
		      )
{
    InitgStyle();

    TCanvas * c1=new TCanvas("c1","c1");

    TH1F * hist1=new TH1F("hist1","hist1",100,6000,7000);
    TH1F * hist2=new TH1F("hist2","hist2",100,6000,7000);

    TUnuran poiss, binom;
    // Initialize unuran to generate normal random numbers from the Poisson
    // distribution with parameter mu
    poiss.InitPoisson(6483.);
    binom.InitBinomial(24733, 6483./24733);//24733
    
    int N = 1000000;
    // Sample distributions N times (generate N random numbers)
    for (int i = 0; i<N; ++i){
      hist1->Fill(poiss.SampleDiscr());
      hist2->Fill(binom.SampleDiscr());
    }
    
    hist1->Scale(1./N);
    hist2->Scale(1./N);

    hist1->SetLineColor(1);
    hist2->SetLineColor(2);

    hist2->Draw();
    hist1->Draw("sames");

//     //Layout
//     TString xtitle = "p_{T}(Jet1) [GeV]";
//     TString ytitle = "Anzahl / %.1f GeV";
//     hist1->GetXaxis()->SetTitleOffset(xtitleOffset);
//     hist1->GetYaxis()->SetTitleOffset(ytitleOffset);
//     hist1->GetXaxis()->SetTitle(xtitle);
//     hist1->GetYaxis()->SetTitle(Form(ytitle.Data(), hist1->GetBinWidth(1)));

//     TH1F * hist2=new TH1F("hist2","hist2",50,0,200);
 
//       //Layout
//     hist2->GetXaxis()->SetTitleOffset(xtitleOffset);
//     hist2->GetYaxis()->SetTitleOffset(ytitleOffset);
//     hist2->GetXaxis()->SetTitle(xtitle);
//     hist2->GetYaxis()->SetTitle(Form(ytitle.Data(), hist2->GetBinWidth(1)));

 
}


