#include "plotCommon.hpp"
#include "littleHelpers.hpp"
#include "includes.hpp"

//-----------------------------------------------------------
struct sample{
  TString name;
  TTree* tree;
  double scale;
  int color;
  int style;
}

//-----------------------------------------------------------

void DYBackground(){
  InitgStyle();
  gStyle->SetPadRightMargin(0.2);
  gStyle->SetTitleYOffset(1.0);
  
  TCanvas * c1=new TCanvas("c1","c1");
  //     c1->Divide(1,3);
  //     c1->cd(1);
  
  cout << "Files einlesen:" << endl;

  TFile *fileZmumu = new TFile("CheckSelectionSaveFilters_Zmumu.root");
  fileZmumu->cd("Analyze");
  gDirectory->pwd();
  TTree *TreeZmumu=0; gDirectory->GetObject("Event",TreeZmumu);	

  TFile *fileZee = new TFile("CheckSelectionSaveFilters_Zee.root");
  fileZmumu->cd("Analyze");
  gDirectory->pwd();
  TTree *TreeZee=0; gDirectory->GetObject("Event",TreeZee);	

  TFile *fileTtbar = new TFile("CheckSelectionSaveFilters_Ttbar.root");
  fileTtbar->cd("Analyze");
  gDirectory->pwd();
  TTree *TreeTtbar=0; gDirectory->GetObject("Event",TreeTtbar);	

  cout << "Histogramme anlegen:" << endl;
  int binNumber = 30;
  TH1D * Histo_Zmumu = new TH1D("Histo_Zmumu","Histo_Zmumu",binNumber,50,150);
  TH1D * Histo_Zee = new TH1D("Histo_Zee","Histo_Zee",binNumber,50,150);
  TH1D * Histo_Signal = new TH1D("Histo_Signal","Histo_Signal",binNumber,50,150);
  TH1D * Histo_TtbarBG = new TH1D("Histo_TtbarBG","Histo_TtbarBG",binNumber,50,150);

  cout << "Histogramme fuellen:" << endl;
  cout << "Zusaetzlicher MET > 35 - Cut auf Different Flavour Events um Kurvenverlaeufe moeglichst gleich zu haben." << endl;
  std::string And = " && ";
  std::string SelectionString = "p_Objects && p_OppositeCharge && p_Trigger && p_TwoJets && p_TwoLeptons";
  std::string SignalString = "p_Dilepton";
  std::string TtbarBGString = "!p_Dilepton";
  std::string SameFlavString = "abs(sameFlavour) == 2 && Met > 35";
  std::string DiffFlavString = "sameFlavour == 0 && Met > 35";
  std::string OutString = "(InvM > 106. || InvM < 76.)";
  std::string InString = "(InvM < 106. && InvM > 76.)";
  TreeZmumu->Draw("InvM>>Histo_Zmumu", (SelectionString + And + SameFlavString).c_str(),"goff");
  TreeZee->Draw("InvM>>Histo_Zee",  (SelectionString + And + SameFlavString).c_str(),"goff");
  TreeTtbar->Draw("InvM>>Histo_Signal",  (SignalString + And + SelectionString + And + SameFlavString).c_str(),"goff");
  TreeTtbar->Draw("InvM>>Histo_TtbarBG",  (TtbarBGString + And + SelectionString + And + SameFlavString).c_str(),"goff");

  cout << "Normierung auf 10pb-1:" << endl;
  double ScaleZmumu = 1.056e-2;
  double ScaleZee = 7.247e-3;
  double ScaleTtbar = 7.268e-3;
  Histo_Zmumu->Scale(ScaleZmumu);
  Histo_Zee->Scale(ScaleZee);
  Histo_Signal->Scale(ScaleTtbar);
  Histo_TtbarBG->Scale(ScaleTtbar);

  cout << "Different Flavour Histogramme anlegen und fuellen:" << endl;
  TH1D * Histo_Zmumu_DF = new TH1D("Histo_Zmumu_DF","Histo_Zmumu_DF",binNumber,50,150);
  TH1D * Histo_Zee_DF = new TH1D("Histo_Zee_DF","Histo_Zee_DF",binNumber,50,150);
  TH1D * Histo_Signal_DF = new TH1D("Histo_Signal_DF","Histo_Signal_DF",binNumber,50,150);
  TH1D * Histo_TtbarBG_DF = new TH1D("Histo_TtbarBG_DF","Histo_TtbarBG_DF",binNumber,50,150);
  TreeZmumu->Draw("InvM>>Histo_Zmumu_DF",  (SelectionString + And + DiffFlavString).c_str(),"goff");
  TreeZee->Draw("InvM>>Histo_Zee_DF",  (SelectionString + And + DiffFlavString).c_str(),"goff");
  TreeTtbar->Draw("InvM>>Histo_Signal_DF",  (SignalString + And + SelectionString + And + DiffFlavString).c_str(),"goff");
  TreeTtbar->Draw("InvM>>Histo_TtbarBG_DF",  (TtbarBGString + And + SelectionString + And + DiffFlavString).c_str(),"goff");

  cout << "Histogramme malen:" << endl;
  c1->SetLogy(true);
  Histo_Zmumu->Add(Histo_Zee);

  int LineWidth = 3;
  Histo_Zmumu->SetXTitle("Invariante Masse [GeV]"); Histo_Zmumu->SetYTitle("Ereignisse pro 10pb^{-1}");
  int upperLimit = 25.;
  Histo_Zmumu->GetYaxis()->SetRangeUser(0.01,upperLimit);
  Histo_Zmumu->GetXaxis()->SetTitleOffset(1.1);
  Histo_Zmumu->GetYaxis()->SetTitleOffset(1.1);
  Histo_Zmumu->SetLineWidth(LineWidth);
  Histo_Zmumu->SetLineColor(kRed);
  Histo_Zmumu->Draw("");
  Histo_Signal->SetLineWidth(LineWidth);
  Histo_Signal->SetLineColor(kGreen+2);
  Histo_Signal->Draw("same");
  Histo_TtbarBG->SetLineWidth(LineWidth);
  Histo_TtbarBG->SetLineColor(kBlue);
  Histo_TtbarBG->Draw("same");

  cout << "Different Flavour Histogramme malen:" << endl;
  Histo_Zmumu_DF->Scale(ScaleZmumu);
  Histo_Zee_DF->Scale(ScaleZee);
  Histo_Signal_DF->Scale(ScaleTtbar);
  Histo_TtbarBG_DF->Scale(ScaleTtbar);
  Histo_Zmumu_DF->Add(Histo_Zee_DF);
  Histo_Zmumu_DF->SetLineWidth(LineWidth);
  Histo_Zmumu_DF->SetLineColor(kRed);
  Histo_Zmumu_DF->SetLineStyle(2);
  Histo_Zmumu_DF->Draw("same");
  Histo_Signal_DF->SetLineWidth(LineWidth);
  Histo_Signal_DF->SetLineColor(kGreen+2);
  Histo_Signal_DF->SetLineStyle(2);
  Histo_Signal_DF->Draw("same");
  Histo_TtbarBG_DF->SetLineWidth(LineWidth);
  Histo_TtbarBG_DF->SetLineColor(kBlue);
  Histo_TtbarBG_DF->SetLineStyle(2);
  Histo_TtbarBG_DF->Draw("same");
	
  left = new TLine(76.,-0.5,76.,upperLimit);
  right = new TLine(106.,-0.5,106.,upperLimit);
  left->SetLineWidth(LineWidth);
  left->SetLineColor(34);
  left->SetLineStyle(2);
  right->SetLineWidth(LineWidth);
  right->SetLineColor(34);
  right->SetLineStyle(2);
  left->Draw();
  right->Draw();

  //Legende
  TH1D * Histo_Dummy1 = new TH1D("Histo_Dummy1","Histo_Dummy1",binNumber,50,150);
  Histo_Dummy1->SetLineWidth(LineWidth);
  Histo_Dummy1->SetLineColor(kBlack);
  Histo_Dummy1->SetLineStyle(1);
  TH1D * Histo_Dummy2 = new TH1D("Histo_Dummy2","Histo_Dummy2",binNumber,50,150);
  Histo_Dummy2->SetLineWidth(LineWidth);
  Histo_Dummy2->SetLineColor(kBlack);
  Histo_Dummy2->SetLineStyle(2);

  TLegend *legende1 = new TLegend(0.53,0.76,0.88,0.89);
  legende1->AddEntry(Histo_Signal,"Signal","l");
  legende1->AddEntry(Histo_Zmumu,"Drell-Yan","l");
  legende1->AddEntry(Histo_TtbarBG,"andere t #bar{t}","l");
  legende1->SetEntrySeparation(0.1);
  legende1->Draw();
  TLegend *legende2 = new TLegend(0.53,0.65,0.88,0.75);
  legende2->AddEntry(Histo_Dummy1,"ee/#mu#mu-Events","l");
  legende2->AddEntry(Histo_Dummy2,"e#mu-Events","l");
  legende2->SetEntrySeparation(0.1);
  legende2->Draw();
	
  cout << "Berechne die Zahlenwerte:" << endl;
  // 	double Zmumu_out = TreeZmumu->GetEntries((SelectionString + And + SameFlavString + And + OutString).c_str()) * ScaleZmumu;
  // 	double Zee_out = TreeZee->GetEntries((SelectionString + And + SameFlavString + And + OutString).c_str()) * ScaleZee;
  // 	double Zmumu_in = TreeZmumu->GetEntries((SelectionString + And + SameFlavString + And + InString).c_str()) * ScaleZmumu;
  // 	double Zee_in = TreeZee->GetEntries((SelectionString + And + SameFlavString + And + InString).c_str()) * ScaleZee;
  // 	double Routin = (Zmumu_out + Zee_out) / (Zmumu_in + Zee_in);

  // 	cout << "Zmumu: There are " << Zmumu_out << " Events outside and " << Zmumu_in << " Events inside the Z-Peak." <<  endl;
  // 	cout << "Zee: There are " << Zee_out << " Events outside and " << Zee_in << " Events inside the Z-Peak." <<  endl;
  // 	cout << "Sum: There are " << Zee_out + Zmumu_out << " Events outside and " << Zee_in + Zmumu_in << " Events inside the Z-Peak." <<  endl;
  // 	cout << "Routin = " << Routin << endl;

  cout << "MonteCarlo:" << endl;
  fehlerrechnung Zmumu_out, Zmumu_in, Zee_out, Zee_in, ROutIn;
  Zmumu_out = fehlerrechnung(TreeZmumu->GetEntries((SelectionString + And + SameFlavString + And + OutString).c_str()),"poisson") * ScaleZmumu;
  Zee_out = fehlerrechnung(TreeZee->GetEntries((SelectionString + And + SameFlavString + And + OutString).c_str()),"poisson") * ScaleZee;
  Zmumu_in = fehlerrechnung( TreeZmumu->GetEntries((SelectionString + And + SameFlavString + And + InString).c_str()),"poisson") * ScaleZmumu;
  Zee_in = fehlerrechnung( TreeZee->GetEntries((SelectionString + And + SameFlavString + And + InString).c_str()),"poisson") * ScaleZee;
  ROutIn = (Zmumu_out + Zee_out) / (Zmumu_in + Zee_in);
  //cprint(" DY-Events im Datenbereich: ", (Zmumu_out + Zee_out));
  cprint("R_outin = ",ROutIn, 3);
  cprint("DY innerhalb unskaliert: ", (Zmumu_in/ScaleZmumu)+(Zee_in/ScaleZee));
  cprint("DY ausserhalb unskaliert: ", (Zmumu_out/ScaleZmumu)+(Zee_out/ScaleZee));

  fehlerrechnung Ttbar_in, Total_in, Ttbar_in_DF, Zmumu_in_DF, Zee_in_DF, Total_in_DF,DYinPeak, DYinData;
  Ttbar_in = fehlerrechnung(TreeTtbar->GetEntries((SelectionString + And + SameFlavString + And + InString).c_str()),"poisson") * ScaleTtbar;
  cprint("Ttbar-SameFlavor-Events im Z-Peak unskaliert: ", Ttbar_in /ScaleTtbar);
  cprint("Ttbar-SameFlavor-Events im Z-Peak: ", Ttbar_in);
  Total_in = Ttbar_in + Zmumu_in + Zee_in;
  Ttbar_in_DF = fehlerrechnung(TreeTtbar->GetEntries((SelectionString + And + DiffFlavString + And + InString).c_str()),"poisson") * ScaleTtbar;
  cprint("Ttbar-DiffFlavor-Events im Z-Peak unskaliert: ", Ttbar_in_DF /ScaleTtbar);
  cprint("Ttbar-DiffFlavor-Events im Z-Peak: ", Ttbar_in_DF);
  cprint("Quotient (sollte 1 sein): ", Ttbar_in/Ttbar_in_DF);
  Zmumu_in_DF = fehlerrechnung(TreeZmumu->GetEntries((SelectionString + And + DiffFlavString + And + InString).c_str()),"poisson") * ScaleZmumu;
  Zee_in_DF   = fehlerrechnung(TreeZee->GetEntries((SelectionString + And + DiffFlavString + And + InString).c_str()),"poisson") * ScaleZee;
  Total_in_DF = Ttbar_in_DF + Zmumu_in_DF + Zee_in_DF;
  DYinPeak = Total_in - Total_in_DF;
  cprint("In MC-Daten wirklich vorhandene Events innerhalb des Peaks: ", (Zmumu_in + Zee_in));
  cprint("Aus \"Daten\": DY-Ereignisse innerhalb des Peaks:             ", DYinPeak);
  DYinData = ROutIn * DYinPeak;
  cprint("Vorhersage der Methode fuer DY-Events im Datenbereich:      ", DYinData);
  cprint("In MC-Daten wirklich vorhandene Events im Datenbereich:     ", (Zmumu_out + Zee_out));
  cprint("Differenz Rekonstruktion - Wahrheit:                        ", (DYinData - Zmumu_out - Zee_out));	

  fehlerrechnung Ttbar_out, Ttbar_out_DF;
  Ttbar_out = fehlerrechnung(TreeTtbar->GetEntries((SelectionString + And + SameFlavString + And + OutString).c_str()),"poisson") * ScaleTtbar;
  cprint("Ttbar-SameFlavor-Events !im Z-Peak unskaliert: ", Ttbar_out /ScaleTtbar);
  cprint("Ttbar-SameFlavor-Events !im Z-Peak: ", Ttbar_out);
  Ttbar_out_DF = fehlerrechnung(TreeTtbar->GetEntries((SelectionString + And + DiffFlavString + And + OutString).c_str()),"poisson") * ScaleTtbar;
  cprint("Ttbar-DiffFlavor-Events !im Z-Peak unskaliert: ", Ttbar_out_DF /ScaleTtbar);
  cprint("Ttbar-DiffFlavor-Events !im Z-Peak: ", Ttbar_out_DF);
  cprint("Quotient (sollte 1 sein): ", Ttbar_out/Ttbar_out_DF);



  // 	TCanvas *c2 = new TCanvas("c2","c2",100,100,900,900);
  // 	//	TH1D * Histo_Diff = new TH1D("Histo_Diff","Histo_Diff",binNumber,50,150);
  // 	Histo_Signal->Add(Histo_TtbarBG);
  // 	Histo_Signal_DF->Add(Histo_TtbarBG_DF);
  // 	//	Histo_Signal->Add(Histo_Signal_DF,-1.);
  // 	Histo_Signal->Divide(Histo_Signal_DF);
  // 	Histo_Signal->Draw("");
  // // 	Histo_Diff->(Histo_Signal,Histo_Signal_DF,1,-1);
  // // 	Histo_Diff->Draw("");
  // 	left->Draw();
  // 	right->Draw();

}
