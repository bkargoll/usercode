#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TLine.h"

void analyze(TString matching, TString variable, int nbins, double lower, double upper,
		bool safe, bool logY, TString style,
		double ratioLower, double ratioUpper, TFile* f, TCanvas* c1, TTree* Tree){
	f->cd(matching);
	gDirectory->pwd();
	gDirectory->GetObject("Tree", Tree);

	c1->cd(1);
	c1->cd(1)->SetLogy(logY);
	TH1D * all = new TH1D("all", "all", nbins, lower, upper);
	all->GetXaxis()->SetTitle(variable);
//	all->GetXaxis()->SetNdivisions(nbins);
	all->GetYaxis()->SetTitle("a.u.");
	Tree->Draw(variable + " >> all", "abs(genTau_p4Vis.Eta()) < 2.5", "");
	TH1D * matched = new TH1D("matched", "matched", nbins, lower, upper);
	matched->SetLineStyle(2);
	Tree->Draw(variable + " >> matched",
			"abs(genTau_p4Vis.Eta()) < 2.5 && genTau_matched", "sames");

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
	ratio->GetXaxis()->SetTitle(variable);
//	ratio->GetXaxis()->SetNdivisions(nbins);
	ratio->GetYaxis()->SetTitle("matched/all");
	ratio->Draw(style);

	double mean = double(Tree->GetEntries("genTau_matched")) / Tree->GetEntries("");
	std::cout << mean << std::endl;
	TLine * meanLine = new TLine(lower,mean,upper,mean);
	meanLine->SetLineStyle(3);
	meanLine->Draw();

	if (safe)
		c1->Print(matching + "_" + variable + ".gif");
}


void efficiency(TString variable, int nbins, double lower, double upper,
		bool safe = false, bool logY = false, TString style = "histo,text00",
		double ratioLower = 0.0, double ratioUpper = 1.1) {

	gStyle->SetPaintTextFormat("5.2f");
	gStyle->SetOptStat(0);

	TCanvas * c1 = new TCanvas("c1", "c1", 1000, 800);
	c1->Divide(1, 2);

	TFile *f = 0;
	f = new TFile("/user/kargoll/results/EffMatrixMatchEff_DYToTauTau_multipleMatching.root");
	TTree *Tree = 0;
	analyze("matchEffCfdR03dPt1000", variable, nbins, lower, upper, safe, logY, style, ratioLower, ratioUpper, f, c1, Tree);
	analyze("matchEffCtdR03dPt1000", variable, nbins, lower, upper, safe, logY, style, ratioLower, ratioUpper, f, c1, Tree);
	analyze("matchEffCfdR003dPt1000", variable, nbins, lower, upper, safe, logY, style, ratioLower, ratioUpper, f, c1, Tree);
	analyze("matchEffCfdR03dPt02", variable, nbins, lower, upper, safe, logY, style, ratioLower, ratioUpper, f, c1, Tree);
}
