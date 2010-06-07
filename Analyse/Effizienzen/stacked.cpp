
//======================================================================
/**
 * @file
 *
 * @author      RWTH Aachen IIIb Top
 *
 * @version     0.0
 *
 * @brief
 *              Plot macro for N-1-plots to investigate selection cuts.
 *              Execute by typing .x plotMacro_N_minus_One.cpp+("../testHistos.root")
 *
 */

//======================================================================


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

//======================================================================
void
stack(
		      )
{
    InitgStyle();

    TTree *tree             = fileToScaledTree(inputfile , 385.           );
    TTree *tree2            = fileToScaledTree(inputfile2, 1944.          );

//     double scaleTtbar             = scale10pb(tree            ,385.               );
//     tree->SetWeight(scaleTtbar);

    TCanvas * c1=new TCanvas("c1","c1");
    c1->Divide(1,3);
    c1->cd(1);

    //Complete Selection: "p_Objects & p_Trigger & p_TwoLeptons & p_TwoJets & p_OppositeCharge & p_SameFlavour"
    //for Ttbar: add "p_Dilepton" or "!p_Dilepton"
    //alternative jet-collection: "p_TwoJetsSC5"
    TH1F * hist1=new TH1F("hist1","hist1",50,0,200);
    hist1->Sumw2();

    tree->Draw("jet1_pt >> hist1","p_Dilepton & p_Objects & p_Trigger & p_TwoLeptons & p_TwoJets & p_OppositeCharge & p_SameFlavour"/*add all filters except the one you like to inspect*/,"HIST"); // to do: Paths am besten auch als Uebergabewerte 

    //Layout
    TString xtitle = "p_{T}(Jet1) [GeV]";
    TString ytitle = "Anzahl / %.1f GeV";
    hist1->GetXaxis()->SetTitleOffset(xtitleOffset);
    hist1->GetYaxis()->SetTitleOffset(ytitleOffset);
    hist1->GetXaxis()->SetTitle(xtitle);
    hist1->GetYaxis()->SetTitle(Form(ytitle.Data(), hist1->GetBinWidth(1)));

    c1->cd(2);
    TH1F * hist2=new TH1F("hist2","hist2",50,0,200);
    hist2->Sumw2();

    tree2->Draw("jet1_pt >> hist2","p_Objects & p_Trigger & p_TwoLeptons & p_TwoJets & p_OppositeCharge & p_SameFlavour"/*add all filters except the one you like to inspect*/,"HIST"); // to do: Paths am besten auch als Uebergabewerte 

    //Layout
    hist2->GetXaxis()->SetTitleOffset(xtitleOffset);
    hist2->GetYaxis()->SetTitleOffset(ytitleOffset);
    hist2->GetXaxis()->SetTitle(xtitle);
    hist2->GetYaxis()->SetTitle(Form(ytitle.Data(), hist2->GetBinWidth(1)));

    // Histogramme stacken
    c1->cd(3);
    THStack *stack = new THStack("stack", "stack");
    stack->Add(hist1);
    stack->Add(hist2);
    stack->Draw("HIST");

    c1->Print("N-1Plots.pdf");

}


