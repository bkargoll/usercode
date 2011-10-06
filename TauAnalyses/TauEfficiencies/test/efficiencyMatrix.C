#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TText.h"
#include "Math/Vector4D.h"

void efficiencyMatrix(int nEvents = 0, int genMax = 35, int recMax = 35,
		int verbose = 0, TString drawStyle = "text", const char* discrNames =
				"none", bool safe = false) {
	//		TString variable, int nbins, double lower, double upper,
	//		bool safe = false, bool logY = false, TString style = "histo,text00",
	//		double ratioLower = 0.0, double ratioUpper = 1.1) {

	gStyle->SetPaintTextFormat("2.1f");
	gStyle->SetOptStat(0);
	gStyle->SetPalette(8);

	TCanvas * c1 = new TCanvas("c1", "c1", 1000, 800);

	// read in file
	TFile *f = 0;
	f = new TFile("/user/kargoll/results/EfficiencyMatrix_DYToTauTau_chargedMatchingdR002dPt02Charge.root");
//		f = new TFile("EfficiencyMatrix_DYToTauTau_test.root");
	f->cd("effMatrix");
	gDirectory->pwd();
	TTree *Tree = 0;
	gDirectory->GetObject("Tree", Tree);

	if (nEvents == 0)
		nEvents = Tree->GetEntries();

	//Declaration of leaves types
	std::vector<bool>* genTau_goodMother;
	std::vector<int>* recoTau_mode;
	std::vector<int>* recoTau_matchedGenMode;
	std::vector<bool>* recoTau_matchedGenGoodMother;
	std::vector<bool>* recoTau_discr_DecayModeFinding;
	std::vector<bool>* recoTau_discr_LooseElectronRejection;
	std::vector<bool>* recoTau_discr_LooseIsolation;
	std::vector<bool>* recoTau_discr_LooseMuonRejection;
	std::vector<bool>* recoTau_discr_MediumElectronRejection;
	std::vector<bool>* recoTau_discr_MediumIsolation;
	std::vector<bool>* recoTau_discr_TightElectronRejection;
	std::vector<bool>* recoTau_discr_TightIsolation;
	std::vector<bool>* recoTau_discr_TightMuonRejection;
	std::vector<bool>* recoTau_discr_VLooseIsolation;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_p4;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* genTau_p4Vis;

	// set branch adresses
	Tree->SetBranchAddress("genTau_goodMother", &genTau_goodMother);
	Tree->SetBranchAddress("recoTau_mode", &recoTau_mode);
	Tree->SetBranchAddress("recoTau_matchedGenMode", &recoTau_matchedGenMode);
	Tree->SetBranchAddress("recoTau_p4", &recoTau_p4);
	Tree->SetBranchAddress("genTau_p4Vis", &genTau_p4Vis);
	Tree->SetBranchAddress("recoTau_matchedGenGoodMother",
			&recoTau_matchedGenGoodMother);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByDecayModeFinding",
			&recoTau_discr_DecayModeFinding);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByLooseElectronRejection",
			&recoTau_discr_LooseElectronRejection);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByLooseIsolation",
			&recoTau_discr_LooseIsolation);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByLooseMuonRejection",
			&recoTau_discr_LooseMuonRejection);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByMediumElectronRejection",
			&recoTau_discr_MediumElectronRejection);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByMediumIsolation",
			&recoTau_discr_MediumIsolation);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByTightElectronRejection",
			&recoTau_discr_TightElectronRejection);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByTightIsolation",
			&recoTau_discr_TightIsolation);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByTightMuonRejection",
			&recoTau_discr_TightMuonRejection);
	Tree->SetBranchAddress("hpsPFTauDiscriminationByVLooseIsolation",
			&recoTau_discr_VLooseIsolation);

	// access generator decay modes
	double ptCut = 10.0;
	double ptCutHigh = 0.0;
	double etaCut = 2.5;
	TString cuts = "";
	if (ptCut != 0.0) {
		char cut[50];
		sprintf(cut, "genTau_p4Vis.Pt() >=  %f", ptCut);
		cuts += cut;
		if (ptCutHigh != 0.0) {
			char cut2[50];
			sprintf(cut2, " && genTau_p4Vis.Pt() <  %f", ptCutHigh);
			cuts += cut2;
		}
		if (etaCut != 0.0)
			cuts += " && ";
	}
	if (etaCut != 0.0) {
		char cut[50];
		sprintf(cut, "fabs(genTau_p4Vis.Eta()) <=  %f", etaCut);
		cuts += cut;
	}

	TH1D * genModes = new TH1D("genModes", "genModes", 37, -1.5, 35.5);
	Tree->Draw("genTau_mode >> genModes", cuts, "goff", nEvents);
	int cntGenSelected = genModes->GetEntries();
	TH1D * temp = new TH1D("temp", "temp", 37, -1.5, 35.5);
	Tree->Draw("genTau_mode >> temp", "", "goff", nEvents);
	int cntGenAll = temp->GetEntries();
	if (verbose > 1) {
		TCanvas * c2 = new TCanvas("c2", "c2", 1000, 800);
		c2->cd();
		genModes->Draw();
		c1->cd();
	}

	// create efficiency matrix histogram
	// x axis: reconstructed mode, y axis: generated mode
	TH2D * effM = new TH2D("effM", "reconstruction efficiency (in %)", 37,
			-1.5, 35.5, 38, -2.5, 35.5);
	effM->GetXaxis()->SetTitle("reconstructed decay mode");
	effM->GetYaxis()->SetTitle("generated decay mode");

	// event loop
	unsigned cntAcceptedTaus = 0, cntDeclinedTaus = 0, cntAllTaus = 0;
	for (int event = 0; event < nEvents; event++) {
		Tree->GetEntry(event);
		if (verbose > 0)
			std::cout << "Event No. " << event << " contains "
					<< genTau_goodMother->size() << " genTaus and "
					<< recoTau_matchedGenGoodMother->size() << " = "
					<< recoTau_matchedGenMode->size() << " = "
					<< recoTau_mode->size() << " recoTaus." << std::endl;

		int nRecoTaus = recoTau_mode->size();
		cntAllTaus += nRecoTaus;
		if (nRecoTaus != recoTau_matchedGenMode->size() || nRecoTaus
				!= recoTau_matchedGenGoodMother->size()) {
			std::cout << "ERROR: WRONG VECTOR SIZE" << std::endl;
			return;
		}
		for (int tau = 0; tau < nRecoTaus; tau++) {
			// discard taus with bad kinematics
			if (etaCut == 0.0)
				etaCut = 999999.9;
			if (ptCutHigh == 0.0)
				ptCutHigh = 999999.9;
			if (fabs(recoTau_p4->at(tau).Eta()) <= etaCut
					&& recoTau_p4->at(tau).Pt() >= ptCut
					&& recoTau_p4->at(tau).Pt() < ptCutHigh) {
				// check discriminators
				if ( true
											&& recoTau_discr_DecayModeFinding->at(tau)
											&& recoTau_discr_LooseElectronRejection->at(tau)
											&& recoTau_discr_LooseIsolation->at(tau)
											&& recoTau_discr_LooseMuonRejection->at(tau)
//											&& recoTau_discr_MediumElectronRejection->at(tau)
//											&& recoTau_discr_MediumIsolation->at(tau)
//											&& recoTau_discr_TightElectronRejection->at(tau)
//											&& recoTau_discr_TightIsolation->at(tau)
//											&& recoTau_discr_TightMuonRejection->at(tau)
//											&& recoTau_discr_VLooseIsolation->at(tau)
						) {
					++cntAcceptedTaus;
					effM->Fill(recoTau_mode->at(tau),
							recoTau_matchedGenMode->at(tau));
					if (verbose > 1)
						std::cout << "Fill bin (" << recoTau_mode->at(tau)
								<< "," << recoTau_matchedGenMode->at(tau)
								<< ")" << std::endl;
				} else {
					++cntDeclinedTaus;
				}
			}

		}
	}

	// calculate ratios
	std::vector<int> recModes(38);
	int genBin;
	int recBin;
	for (genBin = 3; genBin < 39; genBin++) {
		int nGen = genModes->GetBinContent(genBin - 1);
		for (recBin = 1; recBin < 38; recBin++) {
			int oldContent = effM->GetBinContent(recBin, genBin);
			recModes.at(recBin) += oldContent;
			double newContent = (nGen != 0) ? double(oldContent) / nGen * 100.
					: 0.0;
			effM->SetBinContent(recBin, genBin, newContent);
		}
	}

	genBin = 1; // not matched PFTaus
	for (recBin = 1; recBin < 38; recBin++) {
		int oldContent = effM->GetBinContent(recBin, genBin);
		double newContent = (recModes.at(recBin) != 0) ? double(oldContent)
				/ recModes.at(recBin) : 0.0;
		effM->SetBinContent(recBin, genBin, newContent);
	}
	//	genBin = 2; // decay mode not determined. this should never happen


	effM->GetXaxis()->SetRangeUser(-1.5, recMax + 0.5);
	effM->GetYaxis()->SetRangeUser(-2.5, genMax + 0.5);
	effM->GetXaxis()->SetNdivisions(200 + recMax + 4);
	effM->GetYaxis()->SetNdivisions(200 + genMax + 3);
	effM->Draw(drawStyle);
	TLine *hLine = new TLine(-1.5, -0.5, recMax + 1.5, -0.5);
	hLine->SetLineStyle(2);
	hLine->Draw();
	TLine *vLine = new TLine(-0.5, -2.5, -0.5, genMax + 1.5);
	vLine->SetLineStyle(2);
	vLine->Draw();
	TPaveText *pave = new TPaveText(4., genMax - 2.0, 8.5, genMax - 10.0);
	pave->SetFillColor(0);
	pave->SetLineColor(0);
	pave->SetShadowColor(0);
	pave->SetTextSize(0.03);
	pave->SetTextAlign(12);
	//	float textSize = 0.02;
	/*TText *t1=*/
	pave->AddText("Kinematic cuts:\n");
	pave->AddText("p_{T} > 10GeV, |#eta| < 2.5\n");
	char genkineff[50];
	sprintf(genkineff, "gen. cut eff.: %2.3f%%\n", float(cntGenSelected) / cntGenAll);
	pave->AddText(genkineff);
	char recokineff[50];
	sprintf(recokineff, "reco. cut eff.: %2.3f%%\n", float(cntAcceptedTaus + cntDeclinedTaus) / cntAllTaus);
	pave->AddText(recokineff);
	pave->AddText("\n");
	pave->AddText("Tau discriminators:\n");
	//	t1->SetTextSize(textSize);
	TString discrString = discrNames;
	pave->AddText(discrString);
	//	t2->SetTextSize(textSize);
	char disceff[50];
	sprintf(disceff, "discr. eff.: %2.1f%%", float(cntAcceptedTaus)
			/ (cntAcceptedTaus + cntDeclinedTaus) * 100.);
	pave->AddText(disceff);
	pave->Draw();
	c1->Update();

	//	discrString += "Discriminators";
	if (safe)
		c1->Print("efficiencyMatrix_DYToTauTau_finalMatching_" + discrString + "_pt10eta25.gif");

	printf("%i genTaus have been processed in total \n", cntGenAll);
	printf("%i (%2.1f%) passed kinematic cuts 10GeV < p_T < 20GeV, abs(eta)<2.5 \n",
				cntGenSelected, float(cntGenSelected) / cntGenAll);
	printf("%i PFTaus have been processed in total \n", cntAllTaus);
	printf("%i (%2.1f%) passed kinematic cuts 10GeV < p_T < 20GeV, abs(eta)<2.5 \n",
			cntAcceptedTaus + cntDeclinedTaus, float(cntAcceptedTaus + cntDeclinedTaus) / cntAllTaus);
	printf("%i (%2.1f%) passed discriminators, %i (%2.1f%) did not\n",
			cntAcceptedTaus, float(cntAcceptedTaus) / (cntAcceptedTaus + cntDeclinedTaus) * 100.,
			cntDeclinedTaus, float(cntDeclinedTaus) / (cntAcceptedTaus + cntDeclinedTaus) * 100.);

}
