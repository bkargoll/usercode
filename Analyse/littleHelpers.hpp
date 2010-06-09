//======================================================================
/**
 * @file        littleHelpers.hpp
 *
 * @author      Bastian Kargoll
 *
 * @version     2.3
 *
 * @brief       multiple small helper functions which are
 *              useful for various macros, generally.
 *
 */

//======================================================================

#include "includes.hpp"

// Skalierungsfaktor für verschiedene Sample berechnen, Skalierung auf 10pb-1
// crosssection in pb, efficiency = Effizienz der Generator-Selektion
inline double scale10pb(TTree* sample, double crosssection, double efficiency=1.){
	return efficiency * crosssection * 10. / sample->GetEntries() ;
}

// File einlesen und daraus Tree einlesen
TTree* fileToTree(TString file){
	TFile *f = new TFile(file);
	f->cd("Analyze");
	gDirectory->pwd();
	TTree *Tree=0;
	gDirectory->GetObject("Event",Tree);
	return Tree;
}

// File einlesen, daraus Tree einlesen und diesen skalieren
TTree* fileToScaledTree(TString file, double crosssection, double efficiency=1.){
	TFile *f = new TFile(file);
	f->cd("Analyze");
	gDirectory->pwd();
	TTree *Tree=0; gDirectory->GetObject("Event",Tree);
	double scale = scale10pb(Tree, crosssection , efficiency);
	Tree->SetWeight(scale);
	return Tree;
}

// Tree einlesen und daraus Histogramm erstellen
TH1D* histoFromTree(TString name, TTree* sample, TString variable, TString selection, unsigned int nBins, double low, double high, TString options="", TString xtitle = "test", double xtitleOffset=1.1, double ytitleOffset=1.3){
	TH1D * hist=new TH1D(name,name,nBins,low,high);
	hist->Sumw2();
	std::cout << name << ": " << sample->Draw(variable+" >> "+name,selection,options) << " Events gefuellt" << std::endl;
	//Layout
	if(xtitle = "test") xtitle = variable;
	TString ytitle = "Anzahl / %.1f GeV";
	hist->GetXaxis()->SetTitleOffset(xtitleOffset);
	hist->GetYaxis()->SetTitleOffset(ytitleOffset);
	hist->GetXaxis()->SetTitle(xtitle);
	hist->GetYaxis()->SetTitle(Form(ytitle.Data(), hist->GetBinWidth(1)));
	return hist;
}

// Statistik-Boxen an den Rand der Canvas malen, automatische Größenanpassung
void drawStatBox(TH1D* histo, int& step, int nSteps, double FrameSize, int color = -1, int style = 0){
	TPaveStats* statBox = dynamic_cast<TPaveStats*>( histo->GetListOfFunctions()->FindObject("stats") );

	double statboxSpacing = FrameSize / nSteps;    // gleichmaessiges Aufteilen der Statboxes ueber den Rand
	double statboxHeight  = 0.1 * statboxSpacing * (9.+1./nSteps);;  // 1/10 Abstand, (9/10+Abstand/nSteps) Statbox 
	if(color == -1) color = step+1;

	statBox->SetX1NDC(0.80);
	statBox->SetX2NDC(0.99);
	statBox->SetY2NDC(0.95-step*statboxSpacing);
	statBox->SetY1NDC(0.95-step*statboxSpacing-statboxHeight);
	statBox->SetFillColor(color);
	statBox->SetFillStyle(style);
	statBox->Draw();
	step++;
}
// überladen: Statbox-Größen manuell eingeben
void drawStatBox(int& step, TH1D* histo, int color = -1, double statboxHeight = 0.1,  double statboxSpacing = 0.15){
	TPaveStats* statBox = dynamic_cast<TPaveStats*>( histo->GetListOfFunctions()->FindObject("stats") );

	if(color == -1) color = step+1;
	statBox->SetX1NDC(0.80);
	statBox->SetX2NDC(0.99);
	statBox->SetY2NDC(0.95-step*statboxSpacing);
	statBox->SetY1NDC(0.95-step*statboxSpacing-statboxHeight);
	statBox->SetTextColor(color);
	statBox->Draw();
	step++;
}

// struct als Container für alle wichtigen Größen eines Datensatzes
struct sample{
	TTree* tree;
	TString file;
	TString subselection;
	double crosssection;
	double mcEff;
	double scale;
	int color;
	int style;
	TH1D* histo;
};

sample createSample( TString File, TString Subselection, double Crosssection, double McEff){
	sample Sample;
	Sample.file         = File;
	Sample.subselection = Subselection;
	Sample.crosssection = Crosssection;
	Sample.mcEff        = McEff;
	return Sample;
}
