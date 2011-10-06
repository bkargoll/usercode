#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TText.h"
#include "TFile.h"
//#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

void plotDeltaRMatching(TString name = "", bool log = true, int nEvents = 0) {
//	gStyle->SetPaintTextFormat("2.1f");
//	gStyle->SetOptStat(0);
//	gStyle->SetPalette(8);

//	TCanvas * c1 = new TCanvas("c1", "c1", 1000, 800);
//	TCanvas * c2 = new TCanvas("c2", "c2", 1000, 800);
//	TCanvas * c3 = new TCanvas("c3", "c3", 1000, 800);
//	TCanvas * c4 = new TCanvas("c4", "c4", 1000, 800);

	// read in file
	TFile *f = 0;
	f = new TFile("/user/kargoll/results/MatchingEfficiency_DYToTauTau_chargedMatchingdR002dPt02Charge.root");
	f->cd("matchEff");
	gDirectory->pwd();
	TTree *Tree = 0;
	gDirectory->GetObject("Tree", Tree);

	if (nEvents == 0)
		nEvents = Tree->GetEntries();

	//Declaration of leaves types
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_p4;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_p4Matching;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_genp4;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_genp4Matching;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_genp4Vis;
	std::vector<bool> * recoTau_matched;
//	std::vector<unsigned int> * recoTau_decayMode;
	std::vector<int> * recoTau_HPSdecayMode;
//	std::vector<float> * recoTau_discr_DecayModeFinding;

	// set branch adresses
	Tree->SetBranchAddress("recoTau_p4", &recoTau_p4);
	Tree->SetBranchAddress("recoTau_p4Matching", &recoTau_p4Matching);
	Tree->SetBranchAddress("recoTau_genp4", &recoTau_genp4);
	Tree->SetBranchAddress("recoTau_genp4Matching", &recoTau_genp4Matching);
	Tree->SetBranchAddress("recoTau_genp4Vis", &recoTau_genp4Vis);
	Tree->SetBranchAddress("recoTau_matched", &recoTau_matched);
//	Tree->SetBranchAddress("recoTau_decayMode", &recoTau_decayMode);
	Tree->SetBranchAddress("recoTau_HPSdecayMode", &recoTau_HPSdecayMode);
//	Tree->SetBranchAddress("recoTau_discr_DecayModeFinding", &recoTau_discr_DecayModeFinding);

	TH1D * dR = new TH1D("deltaR", "deltaR"+name, 1000, 0., 0.1);
	TH1D * dPhi = new TH1D("deltaPhi", "deltaPhi"+name, 1000, 0., 0.1);
	TH1D * dEta = new TH1D("deltaEta", "deltaEta"+name, 1000, 0., 0.1);
	TH1D * ratioPt = new TH1D("ratioPt", "ratioPt"+name, 1000, -1.0, 1.0);
	// event loop
	for (int event = 0; event < nEvents; event++) {
		Tree->GetEntry(event);
		int nTaus = recoTau_p4->size();
		// tau loop
		for (int tau = 0; tau < nTaus; tau++) {
			if (recoTau_matched->at(tau) && recoTau_HPSdecayMode->at(tau) == 1 ) {
				ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > first = recoTau_p4->at(tau);
				ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > second = recoTau_genp4Vis->at(tau);
					dR->Fill(ROOT::Math::VectorUtil::DeltaR(first,second));
					dPhi->Fill(ROOT::Math::VectorUtil::DeltaPhi(first,second));
					dEta->Fill(fabs(first.Eta() - second.Eta()));
					ratioPt->Fill((first.Pt() - second.Pt()) / first.Pt());
			}
		}
	}

//	c1->SetLogy(log);
//	c2->SetLogy(log);
//	c3->SetLogy(log);
//
//	c1->cd();
//	dR->Draw();
//	if (name != "") {
////		c1->Print(name+"_deltaR_chargedMatchingdR002dPt02Charge.gif");
//		c1->Print(name+"_deltaR_chargedMatchingdR002dPt02Charge.root");
//	}
//	c2->cd();
//	dPhi->Draw();
//	if (name != "") {
////		c2->Print(name+"_deltaPhi__chargedMatchingdR002dPt02Charge.gif");
//		c2->Print(name+"_deltaPhi__chargedMatchingdR002dPt02Charge.root");
//	}
//	c3->cd();
//	dEta->Draw();
//	if (name != "") {
////		c3->Print(name+"_deltaEta__chargedMatchingdR002dPt02Charge.gif");
//		c3->Print(name+"_deltaEta__chargedMatchingdR002dPt02Charge.root");
//	}
//	c4->cd();
//	ratioPt->Draw();
//	if (name != "") {
////		c4->Print(name+"_ratioPt__chargedMatchingdR002dPt02Charge.gif");
//		c4->Print(name+"_ratioPt__chargedMatchingdR002dPt02Charge.root");
//	}

	if (name == "") name = "test";
	TString name2 = "_Histos_chargedMatchingdR002dPt02Charge.root";
	TString n = name+name2;
	TFile * file = new TFile(n, "recreate");
	file->cd();
	gDirectory->pwd();
	dR->Write();
	dPhi->Write();
	dEta->Write();
	ratioPt->Write();
	file->Write();
//	dR->Write();

}
