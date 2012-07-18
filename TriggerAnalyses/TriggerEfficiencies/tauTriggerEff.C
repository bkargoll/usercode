{ 
  gROOT->LoadMacro("./tdrstyle.C");
  setTDRStyle();
  tdrGrid(1);
  gStyle->SetOptFit(0);
  //load fitter 
  
  gROOT->LoadMacro("./logFitTest.C");
  TFile* theFile = new TFile("/user/kargoll/results/turnOn_Run2012A_hltDoubleTauJetofflineSelectedMuTau_offlineSelectedTausXXMisodbMelecrej_tauPt30_TrackPt5MediumIso_dZ0.2_nProngs5_noWJetsMT40_allOfflineProngs_maxInvM9999999.root");
  TH1F* all = (TH1F*)theFile->Get("denominator");
  TH1F* passedL1 = (TH1F*)theFile->Get("numeratorL1");
  TH1F* passedL2 = (TH1F*)theFile->Get("numeratorL2p5");
  //TH1F* passedL2 = (TH1F*)theFile->Get("numeratorL2");
  TH1F* passedHLT = (TH1F*)theFile->Get("numeratorL3");
  TH1F* passedL1_HLT = (TH1F*)theFile->Get("numeratorL1_HLT");
  
  TCanvas* theCanvas = new TCanvas("theCanvas","theCanvas",800,600);
  theCanvas->Divide(2,2);
  //cout << " *************************** Fit L1 *************************** " << endl;
  theCanvas->cd(1);
  TGraphAsymmErrors* theL1Eff = new TGraphAsymmErrors(passedL1, all, "cl=0.683 b(1,1) mode");
  theL1Eff->SetTitle("L1 Efficiency");
  theL1Eff->GetYaxis()->SetTitle("L1 Efficiency");
  theL1Eff->GetXaxis()->SetTitle("PFTau E_{T}");
  theL1Eff->SetMinimum(0);
  theL1Eff->SetMaximum(1.125);
  theL1Eff->Draw("AP");
  TF1* l1Fit = logFitTest(2,passedL1,all);
  l1Fit->SetLineWidth(2);
  l1Fit->Draw("same");
  theL1Eff->Draw("P");
  
  //cout << " *************************** Fit L2 *************************** " << endl;
  theCanvas->cd(2);
  //TGraphAsymmErrors* theL2Eff = new TGraphAsymmErrors(passedL2, passedL1, "cl=0.683 b(1,1) mode");
  TGraphAsymmErrors* theL2Eff = new TGraphAsymmErrors(passedL2, all, "cl=0.683 b(1,1) mode");
  theL2Eff->SetMinimum(0);
  theL2Eff->SetMaximum(1.125);
  theL2Eff->SetTitle("L2.5 Efficiency");
  theL2Eff->GetYaxis()->SetTitle("L2.5 Efficiency");
  theL2Eff->GetXaxis()->SetTitle("PFTau E_{T}");
  theL2Eff->Draw("AP");
  TF1* l2Fit = logFitTest(3,passedL2,all);
  l2Fit->SetLineWidth(2);
  l2Fit->Draw("same");
  theL2Eff->Draw("P");
  
  //cout << " *************************** Fit HLT *************************** " << endl;
  theCanvas->cd(3);
  TGraphAsymmErrors* theHLTEff = new TGraphAsymmErrors(passedHLT, all, "cl=0.683 b(1,1) mode");
  theHLTEff->SetTitle("hltPFTaus Efficiency");
  theHLTEff->GetYaxis()->SetTitle("L3 Efficiency");
  theHLTEff->GetXaxis()->SetTitle("PFTau E_{T}");
  theHLTEff->SetMaximum(1.125);
  theHLTEff->SetMinimum(0);
  theHLTEff->Draw("AP");
  TF1* l3Fit = logFitTest(3,passedHLT,all);
  l3Fit->SetLineWidth(2);
  l3Fit->Draw("same");  
  theHLTEff->Draw("P");
  /*
  TF1* hltFit = new TF1("hltFit", "[0] - exp(([1]-x)*[2])", 15, 200);
  //TF1* hltFit = new TF1("hltFit", "[0]", 15, 160);
  hltFit->SetParameter(0,0.85);
  hltFit->SetLineWidth(2);
  theHLTEff->Fit("hltFit", "R+");
  */
  TGraphAsymmErrors* theTauEff = new TGraphAsymmErrors(passedL1_HLT, all, "cl=0.683 b(1,1) mode");
  theTauEff->SetTitle("L1+HLT Efficiency");
  theTauEff->GetYaxis()->SetTitle("L1+HLT Efficiency");
  theTauEff->GetXaxis()->SetTitle("PFTau E_{T}");
  theTauEff->SetMaximum(1.125);
  theTauEff->SetMinimum(0);
  //get number of bins from all histo
  int nBins = all->GetNbinsX();
  
  //cout << " *************************** Full Trigger *************************** " << endl;
  theCanvas->cd(4);
  //create combined function
  
  string l1Func = "([0]*0.5*(TMath::Erf((x-[1])/2./[2]/sqrt(x))+1.))";
  string l2Func = "([3]*0.5*(TMath::Erf((x-[4])/2./[5]/sqrt(x))+1.))";
  string l3Func = "([6]*0.5*(TMath::Erf((x-[7])/2./[8]/sqrt(x))+1.))";
  //string hltFunc = "[6]";
  string mult = "*";
  string theFcn = l1Func + mult + l2Func + mult + l3Func;
  TF1* tauEffFunc = new TF1("tauEffFunc",theFcn.c_str(), 0, 200);
  tauEffFunc->SetParameter(0,l1Fit->GetParameter(0));
  tauEffFunc->SetParameter(1,l1Fit->GetParameter(1));
  tauEffFunc->SetParameter(2,l1Fit->GetParameter(2));
  tauEffFunc->SetParameter(3,l2Fit->GetParameter(0));
  tauEffFunc->SetParameter(4,l2Fit->GetParameter(1));
  tauEffFunc->SetParameter(5,l2Fit->GetParameter(2));
  tauEffFunc->SetParameter(6,l3Fit->GetParameter(0));
  tauEffFunc->SetParameter(7,l3Fit->GetParameter(1));
  tauEffFunc->SetParameter(8,l3Fit->GetParameter(2));
  tauEffFunc->SetLineColor(2);
  tauEffFunc->SetLineWidth(2);
  theTauEff->Draw("ap");
  tauEffFunc->Draw("same");
  TF1* fullFit = logFitTest(3, passedL1_HLT, all);
  fullFit->SetLineColor(4);
  fullFit->SetLineWidth(2);
  fullFit->Draw("same");
  theTauEff->Draw("P");
  
  //cout << " *************************** Parameters *************************** " << endl;
  cout << "\n L1Fit:\t [0]:" << l1Fit->GetParameter(0) 
       << "\t [1]:" << l1Fit->GetParameter(1) 
       << "\t [2]:" << l1Fit->GetParameter(2)
       << "\n L2Fit:\t [0]:" << l2Fit->GetParameter(0) 
       << "\t [1]:" << l2Fit->GetParameter(1) 
       << "\t [2]:" << l2Fit->GetParameter(2)
       << "\n L3Fit:\t [0]:" << l3Fit->GetParameter(0) 
       << "\t [1]:" << l3Fit->GetParameter(1) 
       << "\t [2]:" << l3Fit->GetParameter(2)
       << "\n Full function: [0]:"<<fullFit->GetParameter(0)
       << "\t [1]:" << fullFit->GetParameter(1) 
       << "\t [2]:" << fullFit->GetParameter(2)
       << endl;

}
