{ 
  gROOT->LoadMacro("./tdrstyle.C");
  setTDRStyle();
  tdrGrid(1);
  gStyle->SetOptFit(0);
  //load fitter 
  
  gROOT->LoadMacro("./logFitTest.C");
  TFile* theFile = new TFile("/user/kargoll/results/TurnOns/DoubleTauForMichal/turnOn_Whole2012_hltDiTauForMichalofflineSelectedMuTau_offlineSelectedTausXXMisodbMelecrej_l2Pt35_l3Pt35_TrackPt1MediumIso_dZ0.2_Prong4_noWJetsMT40_allOfflineProngs.root");
  TString dataset = "Whole2012_Prong4_allOfflineProngs";
  //TString isolation = "combinedIso";
  
  
  TH1F* all = (TH1F*)theFile->Get("denominator");
  TH1F* passedL1 = (TH1F*)theFile->Get("numeratorL1");
  TH1F* passedL1L2 = (TH1F*)theFile->Get("numeratorL1L2");
  TH1F* passedL1L2L2p5 = (TH1F*)theFile->Get("numeratorL1L2L2p5");
  TH1F* passedL1_HLTFull = (TH1F*)theFile->Get("numeratorL1_HLTFull");
  TH1F* passedL1_HLTReg = (TH1F*)theFile->Get("numeratorL1_HLTReg");
  
  
  cout << " *************************** Fit L1 *************************** " << endl;  
  TCanvas* canvasL1 = new TCanvas();
  canvasL1->cd();
  TGraphAsymmErrors* theL1Eff = new TGraphAsymmErrors(passedL1, all, "cl=0.683 b(1,1) mode");
  theL1Eff->SetTitle(0);
  theL1Eff->GetYaxis()->SetTitle("L1 Efficiency");
  theL1Eff->GetXaxis()->SetTitle("PFTau E_{T}");
  theL1Eff->SetMinimum(0);
  theL1Eff->SetMaximum(1.125);
  theL1Eff->Draw("AP");
  TF1* l1Fit = logFitTest(2,passedL1,all);
  l1Fit->SetLineWidth(2);
  l1Fit->Draw("same");
  theL1Eff->Draw("P");
  canvasL1->SaveAs("turnOnTau_L1_"+dataset+".pdf");
  
  cout << " *************************** Fit L1L2 *************************** " << endl;
  TCanvas* canvasL1L2 = new TCanvas();
  canvasL1L2->cd();
  TGraphAsymmErrors* theL2Eff = new TGraphAsymmErrors(passedL1L2, all, "cl=0.683 b(1,1) mode");
  theL2Eff->SetMinimum(0);
  theL2Eff->SetMaximum(1.125);
  theL2Eff->SetTitle(0);
  theL2Eff->GetYaxis()->SetTitle("L1+L2 Efficiency");
  theL2Eff->GetXaxis()->SetTitle("PFTau E_{T}");
  theL2Eff->Draw("AP");
  TF1* l2Fit = logFitTest(3,passedL1L2,all);
  l2Fit->SetLineWidth(2);
  l2Fit->Draw("same");
  theL2Eff->Draw("P");
  canvasL1L2->SaveAs("turnOnTau_L1L2_"+dataset+".pdf");
  
  cout << " *************************** Fit L1L2L2p5 *************************** " << endl;
  TCanvas* canvasL1L2L2p5 = new TCanvas();
  canvasL1L2L2p5->cd();
  TGraphAsymmErrors* theL2p5Eff = new TGraphAsymmErrors(passedL1L2L2p5, all, "cl=0.683 b(1,1) mode");
  theL2p5Eff->SetTitle(0);
  theL2p5Eff->GetYaxis()->SetTitle("L1+L2+L2.5 Efficiency");
  theL2p5Eff->GetXaxis()->SetTitle("PFTau E_{T}");
  theL2p5Eff->SetMaximum(1.125);
  theL2p5Eff->SetMinimum(0);
  theL2p5Eff->Draw("AP");
  TF1* l2p5Fit = logFitTest(3,passedL1L2L2p5,all);
  l2p5Fit->SetLineWidth(2);
  l2p5Fit->Draw("same");  
  theL2p5Eff->Draw("P");
  canvasL1L2L2p5->SaveAs("turnOnTau_L1L2L2p5_"+dataset+".pdf"); 
  

  cout << " *************************** Full Trigger *************************** " << endl;
  TCanvas* canvasL1_HLT = new TCanvas();
  canvasL1_HLT->cd();
  TGraphAsymmErrors* theTauEffFull = new TGraphAsymmErrors(passedL1_HLTFull, all, "cl=0.683 b(1,1) mode");
  theTauEffFull->SetTitle(0);
  theTauEffFull->GetYaxis()->SetTitle("L1+L2+L2.5+L3 Efficiency");
  theTauEffFull->GetXaxis()->SetTitle("PFTau E_{T}");
  theTauEffFull->SetMaximum(1.125);
  theTauEffFull->SetMinimum(0);
  theTauEffFull->SetMarkerColor(kBlack);
  theTauEffFull->Draw("ap");
  TF1* fullFitFull = logFitTest(3, passedL1_HLTFull, all);
  fullFitFull->SetLineWidth(2);
  fullFitFull->SetLineColor(kGray+1);
  fullFitFull->Draw("same");
  
  TGraphAsymmErrors* theTauEffReg = new TGraphAsymmErrors(passedL1_HLTReg, all, "cl=0.683 b(1,1) mode");
  theTauEffReg->SetMarkerStyle(22);
  theTauEffReg->SetMarkerColor(kOrange+10);
  theTauEffReg->SetLineColor(kOrange+10);
  theTauEffReg->Draw("P");
  TF1* fullFitReg = logFitTest(3, passedL1_HLTReg, all);
  fullFitReg->SetLineWidth(2);
  fullFitReg->SetLineColor(kOrange);
  fullFitReg->Draw("same");
  
  theTauEffFull->Draw("P");
  theTauEffReg->Draw("P");
  
  leg = new TLegend(0.5,0.2,0.9,0.4);
  leg->SetNColumns(2);
  leg->SetFillColor(0);
  leg->AddEntry(theTauEffFull,"","P");
  leg->AddEntry(fullFitFull,"Full Tracking","L");
  leg->AddEntry(theTauEffReg,"","P");
  leg->AddEntry(fullFitReg,"Regional Tracking","L");
  leg->Draw();
  
  canvasL1_HLT->SaveAs("turnOnTau_L1L2L2p5L3_DoubleTauForMichal_"+dataset+".pdf");  
  
  //cout << " *************************** Parameters *************************** " << endl;
//  cout << "\n L1Fit:\t [0]:" << l1Fit->GetParameter(0) 
//      << "\t [1]:" << l1Fit->GetParameter(1) 
//       << "\t [2]:" << l1Fit->GetParameter(2)
//       << "\n L2Fit:\t [0]:" << l2Fit->GetParameter(0) 
//       << "\t [1]:" << l2Fit->GetParameter(1) 
//       << "\t [2]:" << l2Fit->GetParameter(2)
//       << "\n L2.5Fit:\t [0]:" << l2p5Fit->GetParameter(0) 
//       << "\t [1]:" << l2p5Fit->GetParameter(1) 
//       << "\t [2]:" << l2p5Fit->GetParameter(2)
   cout<< "\n Full function full tracking: [0]:"<<fullFitFull->GetParameter(0)
       << "\t [1]:" << fullFitFull->GetParameter(1) 
       << "\t [2]:" << fullFitFull->GetParameter(2)
       << "\n Full function reg. tracking: [0]:"<<fullFitReg->GetParameter(0)
       << "\t [1]:" << fullFitReg->GetParameter(1) 
       << "\t [2]:" << fullFitReg->GetParameter(2)
       << endl;

}

