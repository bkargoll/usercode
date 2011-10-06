#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"

void efficiency(TString variable, TString xTitle, int nbins, double lower, double upper,
		bool safe = false, bool logY = false, TString style = "histo,text00",
		double ratioLower = 0.0, double ratioUpper = 1.1) {

	gStyle->SetPaintTextFormat("5.2f");
	gStyle->SetOptStat(0);

	TCanvas * c1 = new TCanvas("c1", "c1", 1000, 800);
	c1->Divide(1, 2);

	TFile *f = 0;
	f = new TFile("/user/kargoll/results/MatchingEfficiency_DYToTauTau_chargedMatchingdR002dPt02Charge.root");
	f->cd("matchEff");
	gDirectory->pwd();
	TTree *Tree = 0;
	gDirectory->GetObject("Tree", Tree);

	c1->cd(1);
	c1->cd(1)->SetLogy(logY);
	TH1D * all = new TH1D("all", "all", nbins, lower, upper);
	all->GetXaxis()->SetTitle(xTitle);
//	all->GetXaxis()->SetNdivisions(nbins);
	all->GetYaxis()->SetTitle("a.u.");
	Tree->Draw(variable + " >> all", "abs(genTau_p4Matching.Eta()) < 2.5", "");
	TH1D * matched = new TH1D("matched", "matched", nbins, lower, upper);
	matched->SetLineStyle(2);
	Tree->Draw(variable + " >> matched",
			"abs(genTau_p4Matching.Eta()) < 2.5 && genTau_matched", "sames");

	c1->cd(2);
	TH1D * ratio = new TH1D("ratio", "ratio", nbins, lower, upper);
	for (int i = 0; i < nbins+1; ++i) {
		if (all->GetBinContent(i) == 0)
			ratio->SetBinContent(i, 0);
		else if (matched->GetBinContent(i) == 0) {
			std::cout << "WARNING: Zero matching eff" << std::endl;
			ratio->SetBinContent(i, 0);
		} else
			ratio->SetBinContent(i, matched->GetBinContent(i)
					/ all->GetBinContent(i));
	}
	ratio->GetYaxis()->SetRangeUser(ratioLower, ratioUpper);
	ratio->SetMarkerSize(2.0);
	ratio->GetXaxis()->SetTitle(xTitle);
//	ratio->GetXaxis()->SetNdivisions(nbins);
	ratio->GetYaxis()->SetTitle("matched/all");
	ratio->Draw(style);

	if (safe)
		c1->Print("matchingEff_chargedMatchingdR002dPt02Charge" + variable + ".gif");

}
