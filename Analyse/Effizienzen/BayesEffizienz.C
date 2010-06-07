#include "TROOT.h"
#include "TCanvas.h"
#include "fehlerrechnungtools.h"
#include "TTree.h"
#include "TFile.h"
#include "littleHelpers.hpp"
#include "plotCommon.hpp"

void BayesEffizienz(){


  InitgStyle();
//   gSystem->CompileMacro("fehlerrechnung.h");
//   gSystem->CompileMacro("fehlerrechnungtools.h");

  using namespace std;
  TCanvas * c1=new TCanvas("c1","c1");
  
  cout << "File einlesen:" << endl;
  TTree *Ttbar             = fileToTree( "/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-Ttbar.root"           );

  cout << "Versuche Methode ohne Histogramme:" << endl;
  double eff =0., tief =0., hoch =0.;
  TGraphAsymmErrors * dummy = new TGraphAsymmErrors();
  dummy->Efficiency(Ttbar->GetEntries("p_Dilepton & p_Objects & p_Trigger & p_TwoLeptons & p_TwoJets & p_OppositeCharge & p_SameFlavour"),Ttbar->GetEntries("p_Dilepton") ,0.6827 ,eff , tief, hoch);
  cout << eff << " mit Intervall [" << tief << "," << hoch << "] entsprechend +" << hoch-eff << " -" << eff-tief << endl;

  

//   TH1F * alle=new TH1F("alle","alle",20,30.,230.);
//   Ttbar->Draw("jet1_pt >> alle","p_Dilepton","HIST,goff"); 

//   TH1F * sel=new TH1F("sel","sel",20,30.,230.);
//   Ttbar->Draw("jet1_pt >> sel","p_Dilepton & p_Objects & p_Trigger & p_TwoLeptons & p_TwoJets & p_OppositeCharge & p_SameFlavour","HIST,goff"); 

//   TGraphAsymmErrors * effizienz = new TGraphAsymmErrors(sel,alle);

//   //Layout
// //   hist2->GetXaxis()->SetTitleOffset(xtitleOffset);
// //   hist2->GetYaxis()->SetTitleOffset(ytitleOffset);
//   effizienz->GetXaxis()->SetTitle("p_{T} (Jet1) [GeV]");
//   effizienz->GetYaxis()->SetTitle("Effizienz");

  TH1F * alle=new TH1F("alle","alle",12,-2.4,2.4);
  Ttbar->Draw("jet1_eta >> alle","p_Dilepton","HIST,goff"); 

  TH1F * sel=new TH1F("sel","sel",12,-2.4,2.4);
  Ttbar->Draw("jet1_eta >> sel","p_Dilepton & p_Objects & p_Trigger","HIST,goff"); //& p_TwoLeptons & p_TwoJets & p_OppositeCharge & p_SameFlavour","HIST,goff"); 

  TGraphAsymmErrors * effizienz = new TGraphAsymmErrors(sel,alle);

  //Layout
//   hist2->GetXaxis()->SetTitleOffset(xtitleOffset);
//   hist2->GetYaxis()->SetTitleOffset(ytitleOffset);
  effizienz->GetXaxis()->SetTitle("#eta_{T} (Jet1)");
  effizienz->GetYaxis()->SetTitle("Triggereffizienz");
  
  effizienz->Draw("ALP");

  //cout << "Die Effizienz ist " << (double)sel->GetEntries() / alle->GetEntries() << " mit unterem Fehler " << effizienz->GetErrorYlow(0) << " und oberem Fehler " << effizienz->GetErrorYhigh(0) << endl;


}
