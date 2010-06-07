
//======================================================================
/**
 * @file
 *
 * @author      RWTH Aachen IIIb Top
 *
 * @version     0.0
 *
 * @brief
 *              compare SC5 and anti-kt
 *              
 *
 */

//======================================================================

#include "/home/home2/institut_3b/kargoll/Analyse/plotCommon.hpp"
#include "/home/home2/institut_3b/kargoll/Analyse/includes.hpp"
#include "/home/home2/institut_3b/kargoll/Analyse/littleHelpers.hpp"

//======================================================================

void
compareAlgorithms( TString variable , TString varName = "replaceMe", TString unit = "", int Nbins = 50, double lower=0., double upper=50., bool ylog = false)
{
    using namespace std;
    InitgStyle();
    //gStyle->SetPadRightMargin(0.2); // line commented: boxes inside, line not commented: boxes outside
    //    gStyle->SetTitleYOffset(1.0);

    TCanvas * c1=new TCanvas("c1","c1");
//     c1->Divide(1,3);
//     c1->cd(1);

    cout << "Tree einlesen:" << endl;
    TTree * Tree = fileToScaledTree("/user/kargoll/results/allObjectsAfterTrigger/10TeV/Ttbar/allObjectsAfterTrigger-10TeV-Ttbar.root",385.);

    cout << "Histogramme anlegen:" << endl;
    Double_t maximumBin,minimumBin;
    TH1D* Antikt = histoFromTree("anti-k_{T}", Tree, "jet_"+variable      , ""/*Selektion*/, Nbins, lower, upper, "goff"/*options*/,  varName);
    Antikt->SetLineColor(1);
    Antikt->SetLineStyle(1);
    //   Antikt->Scale(1/Antikt->Integral());
    if(unit = ""){
      Antikt->GetXaxis()->SetTitle(varName);
      Antikt->GetYaxis()->SetTitle("Anzahl");
    }
    else{
      Antikt->GetXaxis()->SetTitle(varName+" / "+unit);
      TString ytitle = "Anzahl / %.1f "+unit;
      Antikt->GetYaxis()->SetTitle(Form(ytitle.Data(), Antikt->GetBinWidth(1)));
    }
    maximumBin = Antikt->GetMaximum();
    minimumBin = Antikt->GetMinimum();
    Antikt->Draw("HIST");
    c1->Update();

    TH1D* SC5    = histoFromTree("SIScone5"  , Tree, "jetSC5_"+variable   , ""/*Selektion*/, Nbins, lower, upper, "goff"/*options*/,  varName);
    SC5->SetLineColor(2);
    SC5->SetLineStyle(2);
    SC5->GetXaxis()->SetTitle(varName);
    //    SC5->Scale(1/SC5->Integral());
    if (SC5->GetMaximum() > maximumBin ) maximumBin = SC5->GetMaximum();
    if (SC5->GetMinimum() < minimumBin ) minimumBin = SC5->GetMinimum();
    gStyle->SetStatY(0.75);
    SC5->Draw("HIST,sames");
    c1->Update();
    TPaveStats* statBox = dynamic_cast<TPaveStats*>( SC5->GetListOfFunctions()->FindObject("stats") );
    statBox->SetLineColor(2);
    statBox->SetLineStyle(2);
    statBox->Draw();
    //    c1->Update();

    c1->SetLogy(ylog);
//     if(ylog) Antikt->GetYaxis()->SetRangeUser((Antikt->GetMinimum())*0.9,maximumBin*1.1);
//     else     Antikt->GetYaxis()->SetRangeUser(0.0,maximumBin*1.1);
    Antikt->GetYaxis()->SetRangeUser(minimumBin*0.9,maximumBin*1.1);
    c1->Update();

    c1->Print(varName+".pdf");
    c1->Print(varName+".ps");
    c1->Print(varName+".eps");
    c1->Print(varName+".root");

}

