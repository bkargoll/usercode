
//======================================================================
/**
 * @file
 *
 * @author      RWTH Aachen IIIb Top
 *
 * @version     0.0
 *
 * @brief
 *              Plot macro for N-1-plots to investigate selection cuts.
 *              Execute by typing .x plotMacro_N_minus_One.cpp+("../testHistos.root")
 *
 */

//======================================================================


#include "plotCommon.hpp"
#include "includes.hpp"
#include "littleHelpers.hpp"

//======================================================================

void
stackThem( TString variable , TString xTitle = "replaceMe", int Nbins = 50, double lower=0., double upper=50., bool ylog = false)
{
    using namespace std;
    InitgStyle();
    gStyle->SetPadRightMargin(0.2);
    gStyle->SetTitleYOffset(1.0);

    TCanvas * c1=new TCanvas("c1","c1");
//     c1->Divide(1,3);
//     c1->cd(1);

    cout << "!!!!!WARNUNG!!!!! " << endl;
    cout << "Da von den BCtoE- und EMenriched-Events ohnehin keines die Selektion passiert, wurden diese Files auskommentiert!" << endl;
    cout << "!!!!!WARNUNG!!!!! " << endl;

    cout << "Datensamples anmelden:" << endl;
    map<TString,sample> samples;
    samples["Signal"]            = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-Ttbar.root"            , "p_Dilepton &" , 385.     , 1.);
    samples["TtbarBG"]           = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-Ttbar.root"            , "!p_Dilepton &", 385.     , 1.);
    samples["WW"]                = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-WW.root"               , ""             , 74.      , 1.);
    samples["ZZ"]                = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-ZZ.root"               , ""             , 10.5     , 1.);
    samples["WZ"]                = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-WZ.root"               , ""             , 32.      , 1.);
    samples["Zee"]               = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-Zee.root"              , ""             , 1944.    , 1.);
    samples["Zmumu"]             = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-Zmumu.root"            , ""             , 1944.    , 1.);
    samples["Ztautau"]           = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-Ztautau.root"          , ""             , 1944.    , 1.);
    samples["InclusiveMu15"]     = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-InclusiveMu15.root"    , ""             , 0.5091e9 , 2.881e-4);
//     samples["BCtoE20to30"]       = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-BCtoE20to30.root"      , ""             , 0.4e9    , 4.8e-4);
//     samples["BCtoE30to80"]       = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-BCtoE30to80.root"      , ""             , 0.1e9    , 2.4e-3);
//     samples["BCtoE80to170"]      = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-BCtoE80to170.root"     , ""             , 1.9e6    , 0.012);
//     samples["EMenriched20to30"]  = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-EMenriched20to30.root" , ""             , 0.4e9    , 0.008);
//     samples["EMenriched30to80"]  = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-EMenriched30to80.root" , ""             , 0.1e9    , 0.047);
//     samples["EMenriched80to170"] = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-EMenriched80to170.root", ""             , 1.9e6    , 0.15);
    samples["Wenu"]              = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-Wenu.root"             , ""             , 42800./3 , 0.738);
    samples["Wmunu"]             = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-Wmunu.root"            , ""             , 42800./3 , 0.691);
    samples["SingleTops"]        = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-SingleTops.root"       , ""             , 5.       , 0.32442);
    samples["SingleTopt"]        = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-SingleTopt.root"       , ""             , 130.     , 0.32442);
    samples["SingleToptW"]       = createSample("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-SingleToptW.root"      , ""             , 29.      , 1.);

    cout << "Files einlesen, Selektion ausfuehren und Histogramme fuellen:" << endl;
    map<TString,sample>::iterator iter;
    TString selection = "11 == 11";//"p_Objects & p_Trigger & p_TwoLeptons & p_TwoJets & p_OppositeCharge & p_SameFlavour";
    for ( iter=samples.begin(); iter!=samples.end(); ++iter){
      (*iter).second.tree  = fileToScaledTree( (*iter).second.file, (*iter).second.crosssection, (*iter).second.mcEff );
      (*iter).second.histo = histoFromTree( "Hist"+(*iter).first, (*iter).second.tree, variable, (*iter).second.subselection+selection, Nbins, lower, upper, "goff");
    }

    cout << "Histogramme zusammenfassen:" << endl;
    TH1D *HistSignal = new TH1D(*samples["Signal"].histo);
    HistSignal->SetName("Signal");
    TH1D* HistZlike = new TH1D(*samples["ZZ"].histo);
    HistZlike->Add(samples["WZ"].histo);
    HistZlike->Add(samples["Zee"].histo);
    HistZlike->Add(samples["Zmumu"].histo);
    HistZlike->Add(samples["Ztautau"].histo);
    HistZlike->SetName("Z-like");
    TH1D* HistToplike = new TH1D(*samples["TtbarBG"].histo);
    HistToplike->Add(samples["SingleTops"].histo);
    HistToplike->Add(samples["SingleTopt"].histo);
    HistToplike->Add(samples["SingleToptW"].histo);
    HistToplike->SetName("Top-like");
    TH1D* HistOtherBG = new TH1D(*samples["WW"].histo);
    HistOtherBG->Add(samples["Wenu"].histo);
    HistOtherBG->Add(samples["Wmunu"].histo);
    HistOtherBG->Add(samples["InclusiveMu15"].histo);
//     HistOtherBG->Add(samples["BCtoE20to30"].histo);
//     HistOtherBG->Add(samples["BCtoE30to80"].histo);
//     HistOtherBG->Add(samples["BCtoE80to170"].histo);
//     HistOtherBG->Add(samples["EMenriched20to30"].histo);
//     HistOtherBG->Add(samples["EMenriched30to80"].histo);
//     HistOtherBG->Add(samples["EMenriched80to170"].histo);
    HistOtherBG->SetName("otherBackgrounds");

    cout << "Layout festlegen:" << endl;
    gStyle->SetStatX(2.0); gStyle->SetStatY(2.0); // hack: move the individual statbox out of the visible canvas

    HistSignal->SetLineColor(30);
    HistSignal->SetFillColor(30);
    HistSignal->SetFillStyle(3001);
    HistSignal->Draw();
    c1->Update();

    HistZlike->SetLineColor(46);
    HistZlike->SetFillColor(46);
    HistZlike->SetFillStyle(3002);
    HistZlike->Draw();
    c1->Update();

    HistToplike->SetLineColor(40);
    HistToplike->SetFillColor(40);
    HistToplike->SetFillStyle(3003);
    HistToplike->Draw();
    c1->Update();

    HistOtherBG->SetLineColor(28);
    HistOtherBG->SetFillColor(28);
    HistOtherBG->SetFillStyle(3004);
    HistOtherBG->Draw();
    c1->Update();

    cout << "Histogramme stacken:" << endl;
    THStack *stack = new THStack("stack", "stack");
//     stack->Add(HistOtherBG);
//     stack->Add(HistToplike);
//     stack->Add(HistZlike);
//     stack->Add(HistSignal);
 stack->Add(HistSignal);
 stack->Add(HistZlike);
 stack->Add(HistToplike);
 // stack->Add(HistOtherBG);

    stack->Draw("HIST");
    stack->GetYaxis()->SetTitle("Anzahl");
    stack->GetXaxis()->SetTitle( (xTitle == "replaceMe") ? variable : xTitle );
    if (ylog) c1->SetLogy();
    c1->Update();
    
    cout << "Stat-Boxen anlegen:" << endl;
    double frameSize = 1. - c1->GetTopMargin() - c1->GetBottomMargin();
    int step = 0;
    int steps = 4;
    drawStatBox(HistSignal   , step, steps, frameSize, 30,3001);
    drawStatBox(HistZlike    , step, steps, frameSize, 46,3002);
    drawStatBox(HistToplike  , step, steps, frameSize, 40,3003);
    drawStatBox(HistOtherBG  , step, steps, frameSize, 28,3004);

    c1->Print(variable+"_stacked.pdf");
    c1->Print(variable+"_stacked.ps");
    c1->Print(variable+"_stacked.eps");
    c1->Print(variable+"_stacked.root");


}

