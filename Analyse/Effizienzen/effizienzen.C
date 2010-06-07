#include "TROOT.h"
#include "TCanvas.h"
#include "fehlerrechnungtools.h"
#include "TTree.h"
#include "TFile.h"

// File einlesen und daraus Tree einlesen
TTree* fileToTree(TString sample, int SPE){
  TFile *f = 0;
  if(SPE == 10)      f = new TFile("/user/kargoll/results/MC-10TeV/nTuple/MC-10TeV-"+sample+".root");
  else if (SPE == 7) f = new TFile("/user/kargoll/results/MC-7TeV/nTuple/MC-7TeV-"+sample+".root");
  f->cd("Analyze");
  gDirectory->pwd();
  TTree *Tree=0; gDirectory->GetObject("Event",Tree);
  return Tree;
}

// Skalierungsfaktor für verschiedene Sample berechnen, Skalierung auf 10pb-1
// crosssection in pb, efficiency = Effizienz der Generator-Selektion
inline double scale10pb(TTree* sample, double crosssection, double efficiency=1.){
  return efficiency * crosssection * 10. / sample->GetEntries() ;
}

void cutFlow(TTree* sample, TString name, double scale, fehlerrechnung* N){
  TString step[6];
  step[0] = "ERROR";
  if (name == "Signal") step[0] = "p_Dilepton & p_Objects";
  else if (name == "TtbarBG") step[0] = "!p_Dilepton & p_Objects";
  else step[0] = "p_Objects";
  step[1] = step[0]+" & p_Trigger";
  step[2] = step[1]+" & p_TwoLeptons";
  step[3] = step[2]+" & p_TwoJets";
  step[4] = step[3]+" & p_OppositeCharge";
  step[5] = step[4]+" & p_SameFlavour";

  for (int i=0;i<6;i++){
    N[i] = scale * fehlerrechnung(sample->GetEntries(step[i]),"poisson");
    cprint(name,N[i],3);
  }
 }
  


void effizienzen(int SPE){

//   gSystem->CompileMacro("fehlerrechnung.h");
//   gSystem->CompileMacro("fehlerrechnungtools.h");

  using namespace std;
  if ( !(SPE == 10 || SPE == 7) ){
    cout << "Schwerpunktsenergie nicht gueltig." << endl;
    return;
  }
  cout << "Es werden die Effizienzen fuer eine Schwerpunktsenergie von " << SPE << " TeV berechnet." << endl;
  cout << "Soll der gesamte Cut-Flow berechnet werden? Sonst werden nur die Endergebnisse ausgegeben." << endl;
  cout << "\"Ja\" oder \"y\" fuer Cut-Flow" << endl;
  TString cutflow;
  cin >> cutflow;
  
  cout << "Files einlesen:" << endl;
  TTree *Ttbar             = fileToTree("Ttbar"            , SPE);
  TTree *WW                = fileToTree("WW"               , SPE);
  TTree *ZZ                = fileToTree("ZZ"               , SPE);
  TTree *WZ                = fileToTree("WZ"               , SPE);
  TTree *Zee               = fileToTree("Zee"              , SPE);
  TTree *Zmumu             = fileToTree("Zmumu"            , SPE);
  TTree *Ztautau           = fileToTree("Ztautau"          , SPE);
  TTree *InclusiveMu15     = fileToTree("InclusiveMu15"    , SPE);
  TTree *BCtoE20to30       = fileToTree("BCtoE20to30"      , SPE);
  TTree *BCtoE30to80       = fileToTree("BCtoE30to80"      , SPE);
  TTree *BCtoE80to170      = fileToTree("BCtoE80to170"     , SPE);
  TTree *EMenriched20to30  = fileToTree("EMenriched20to30" , SPE);
  TTree *EMenriched30to80  = fileToTree("EMenriched30to80" , SPE);
  TTree *EMenriched80to170 = fileToTree("EMenriched80to170", SPE);
  TTree *Wenu              = fileToTree("Wenu"             , SPE);
  TTree *Wmunu             = fileToTree("Wmunu"            , SPE);
  TTree *SingleTops        = fileToTree("SingleTops"       , SPE);
  TTree *SingleTopt        = fileToTree("SingleTopt"       , SPE);
  TTree *SingleToptW       = fileToTree("SingleToptW"      , SPE);
  

  double scaleTtbar=0.,scaleWW=0.,scaleZZ=0.,scaleWZ=0.,scaleZee=0.,scaleZmumu=0.,scaleZtautau=0.,scaleInclusiveMu15=0.,scaleBCtoE20to30=0.,scaleBCtoE30to80=0.,scaleBCtoE80to170=0.,scaleEMenriched20to30=0.,scaleEMenriched30to80=0.,scaleEMenriched80to170=0.,scaleWenu=0.,scaleWmunu=0.,scaleSingleTops=0.,scaleSingleTopt=0.,scaleSingleToptW=0.;
  if (SPE == 10){
    cout << "Skalierungsfaktoren berechnen fuer 10 TeV:" << endl;
    // Kommentare: Quellenangabe, mit
    // TopXSec = https://twiki.cern.ch/twiki/bin/view/CMS/CrossSections_3XSeries,   NLO-Angaben
    // Prod    = https://twiki.cern.ch/twiki/bin/viewauth/CMS/ProductionSummer2009, LO-Angaben
    // Effizienzen sind immer von der Production Page (ausser für Single Top)
    scaleTtbar             = scale10pb(Ttbar            ,385.               ); // TopXSec
    scaleWW                = scale10pb(WW               ,74.                ); // TopXSec
    scaleZZ                = scale10pb(ZZ               ,10.5               ); // TopXSec
    scaleWZ                = scale10pb(WZ               ,32.                ); // TopXSec
    scaleZee               = scale10pb(Zee              ,1944.              ); // Prod
    scaleZmumu             = scale10pb(Zmumu            ,1944.              ); // Prod
    scaleZtautau           = scale10pb(Ztautau          ,1944.              ); // Prod
    scaleInclusiveMu15     = scale10pb(InclusiveMu15    ,0.5091e9 ,2.881e-4 ); // Prod
    scaleBCtoE20to30       = scale10pb(BCtoE20to30      ,0.4e9    ,4.8e-4   ); // Prod
    scaleBCtoE30to80       = scale10pb(BCtoE30to80      ,0.1e9    ,2.4e-3   ); // Prod
    scaleBCtoE80to170      = scale10pb(BCtoE80to170     ,1.9e6    ,0.012    ); // Prod
    scaleEMenriched20to30  = scale10pb(EMenriched20to30 ,0.4e9    ,0.008    ); // Prod
    scaleEMenriched30to80  = scale10pb(EMenriched30to80 ,0.1e9    ,0.047    ); // Prod
    scaleEMenriched80to170 = scale10pb(EMenriched80to170,1.9e6    ,0.15     ); // Prod
    scaleWenu              = scale10pb(Wenu             ,42800./3 ,0.738    ); // TopXSec
    scaleWmunu             = scale10pb(Wmunu            ,42800./3 ,0.691    ); // TopXSec
    scaleSingleTops        = scale10pb(SingleTops       ,5.       ,0.32442  ); // TopXSec
    scaleSingleTopt        = scale10pb(SingleTopt       ,130.     ,0.32442  ); // TopXSec
    scaleSingleToptW       = scale10pb(SingleToptW      ,29.                ); // TopXSec
  }
  if (SPE == 7){
    cout << "Skalierungsfaktoren berechnen fuer 7 TeV:" << endl;
    // Kommentare: Quellenangabe, mit
    // TopXSec = https://twiki.cern.ch/twiki/bin/view/CMS/CrossSections_3XSeries,   NLO-Angaben               // TODO
    // Prod    = https://twiki.cern.ch/twiki/bin/viewauth/CMS/ProductionSummer2009, LO-Angaben
    // Effizienzen sind immer von der Production Page (ausser für Single Top)
    scaleTtbar             = scale10pb(Ttbar            , 165.               ); // TopXSec, approx. NNLO
    scaleWW                = scale10pb(WW               , 43.                ); // TopXSec
    scaleZZ                = scale10pb(ZZ               , 5.9                ); // TopXSec
    scaleWZ                = scale10pb(WZ               , 18.                ); // TopXSec
    scaleZee               = scale10pb(Zee              , 1300.              ); // Prod
    scaleZmumu             = scale10pb(Zmumu            , 1300.              ); // Prod
    scaleZtautau           = scale10pb(Ztautau          , 1300.              ); // Prod
    scaleInclusiveMu15     = scale10pb(InclusiveMu15    , 0.2969e9 , 0.00037 ); // Prod
    scaleBCtoE20to30       = scale10pb(BCtoE20to30      , 0.2355e9 , 0.00046 ); // Prod
    scaleBCtoE30to80       = scale10pb(BCtoE30to80      , 0.0593e9 , 0.00234 ); // Prod
    scaleBCtoE80to170      = scale10pb(BCtoE80to170     , 0.906e6  , 0.0104  ); // Prod
    scaleEMenriched20to30  = scale10pb(EMenriched20to30 , 0.2355e9  , 0.0073  ); // Prod
    scaleEMenriched30to80  = scale10pb(EMenriched30to80 , 0.0593e9 , 0.059   ); // Prod
    scaleEMenriched80to170 = scale10pb(EMenriched80to170, 0.906e6  , 0.148   ); // Prod
    scaleWenu              = scale10pb(Wenu             , 28000./3 , 0.779   ); // TopXSec
    scaleWmunu             = scale10pb(Wmunu            , 28000./3 , 0.742   ); // TopXSec
    scaleSingleTops        = scale10pb(SingleTops       , 4.6      , 0.32442 ); // TopXSec, NNNLO
    scaleSingleTopt        = scale10pb(SingleTopt       , 63.      , 0.32442 ); // TopXSec
    scaleSingleToptW       = scale10pb(SingleToptW      , 10.6               ); // TopXSec
  }

  
  if(cutflow == "Ja" || cutflow == "y"){
    cout << "Groessen des Cut-Flows werden bestimmt." << endl;
    cout << "Der Cut-Flow ist: Alle -- nach Trigger -- nach 2 Leptonen -- nach 2 Jets -- nach versch. Ladung -- nach SameFlavour-Schnitten" << endl;
    fehlerrechnung CFSignal[6], CFTtbarBG[6], CFWW[6], CFZZ[6], CFWZ[6], CFZee[6], CFZmumu[6], CFZtautau[6], CFInclusiveMu15[6], CFBCtoE20to30[6], CFBCtoE30to80[6], CFBCtoE80to170[6], CFEMenriched20to30[6], CFEMenriched30to80[6], CFEMenriched80to170[6], CFWenu[6], CFWmunu[6], CFSingleTops[6], CFSingleTopt[6], CFSingleToptW[6];
    cutFlow(Ttbar,"Signal",scaleTtbar,CFSignal);
    cutFlow(Ttbar,"TtbarBG",scaleTtbar,CFTtbarBG);
    cutFlow(WW,"WW",scaleWW,CFWW);
    cutFlow(ZZ,"ZZ",scaleZZ,CFZZ);
    cutFlow(WZ,"WZ",scaleWZ,CFWZ);
    cutFlow(Zee,"Zee",scaleZee,CFZee);
    cutFlow(Zmumu,"Zmumu",scaleZmumu,CFZmumu);
    cutFlow(Ztautau,"Ztautau",scaleZtautau,CFZtautau);
    cutFlow(InclusiveMu15,"InclusiveMu15",scaleInclusiveMu15,CFInclusiveMu15);
    cutFlow(BCtoE20to30,"BCtoE20to30",scaleBCtoE20to30,CFBCtoE20to30);
    cutFlow(BCtoE30to80,"BCtoE30to80",scaleBCtoE30to80,CFBCtoE30to80);
    cutFlow(BCtoE80to170,"BCtoE80to170",scaleBCtoE80to170,CFBCtoE80to170);
    cutFlow(EMenriched20to30,"EMenriched20to30",scaleEMenriched20to30,CFEMenriched20to30);
    cutFlow(EMenriched30to80,"EMenriched30to80",scaleEMenriched30to80,CFEMenriched30to80);
    cutFlow(EMenriched80to170,"EMenriched80to170",scaleEMenriched80to170,CFEMenriched80to170);
    cutFlow(Wenu,"Wenu",scaleWenu,CFWenu);
    cutFlow(Wmunu,"Wmunu",scaleWmunu,CFWmunu);
    cutFlow(SingleTops,"SingleTops",scaleSingleTops,CFSingleTops);
    cutFlow(SingleTopt,"SingleTopt",scaleSingleTopt,CFSingleTopt);
    cutFlow(SingleToptW,"SingleToptW",scaleSingleToptW,CFSingleToptW);
  }


  cout << "Selektion ausfuehren:" << endl;
  TString Selektion = "p_Objects & p_OppositeCharge & p_SameFlavour & p_Trigger & p_TwoJets & p_TwoLeptons";
  fehlerrechnung NTtbar,NSignal,NTtbarBG,NWW,NZZ,NWZ,NZee,NZmumu,NZtautau,NInclusiveMu15,NBCtoE20to30,NBCtoE30to80,NBCtoE80to170,NEMenriched20to30,NEMenriched30to80,NEMenriched80to170,NWenu,NWmunu,NSingleTops,NSingleTopt,NSingleToptW;
  NTtbar             = scaleTtbar             * fehlerrechnung(Ttbar->GetEntries(Selektion));//,"poisson");
  NSignal            = scaleTtbar             * fehlerrechnung(Ttbar->GetEntries("p_Dilepton & "+Selektion),"poisson");
  NTtbarBG           = scaleTtbar             * fehlerrechnung(Ttbar->GetEntries("!p_Dilepton & "+Selektion),"poisson");
  NWW                = scaleWW                * fehlerrechnung(WW->GetEntries(Selektion));//"poisson");
  NZZ                = scaleZZ                * fehlerrechnung(ZZ->GetEntries(Selektion));//"poisson");
  NWZ                = scaleWZ                * fehlerrechnung(WZ->GetEntries(Selektion));//"poisson");
  NZee               = scaleZee               * fehlerrechnung(Zee->GetEntries(Selektion));//"poisson");
  NZmumu             = scaleZmumu             * fehlerrechnung(Zmumu->GetEntries(Selektion));//"poisson");
  NZtautau           = scaleZtautau           * fehlerrechnung(Ztautau->GetEntries(Selektion));//"poisson");
  NInclusiveMu15     = scaleInclusiveMu15     * fehlerrechnung(InclusiveMu15->GetEntries(Selektion));//"poisson");
  NBCtoE20to30       = scaleBCtoE20to30       * fehlerrechnung(BCtoE20to30->GetEntries(Selektion));//"poisson");
  NBCtoE30to80       = scaleBCtoE30to80       * fehlerrechnung(BCtoE30to80->GetEntries(Selektion));//"poisson");
  NBCtoE80to170      = scaleBCtoE80to170      * fehlerrechnung(BCtoE80to170->GetEntries(Selektion));//"poisson");
  NEMenriched20to30  = scaleEMenriched20to30  * fehlerrechnung(EMenriched20to30->GetEntries(Selektion));//"poisson");
  NEMenriched30to80  = scaleEMenriched30to80  * fehlerrechnung(EMenriched30to80->GetEntries(Selektion));//"poisson");
  NEMenriched80to170 = scaleEMenriched80to170 * fehlerrechnung(EMenriched80to170->GetEntries(Selektion));//"poisson");
  NWenu              = scaleWenu              * fehlerrechnung(Wenu->GetEntries(Selektion));//"poisson");
  NWmunu             = scaleWmunu             * fehlerrechnung(Wmunu->GetEntries(Selektion));//"poisson");
  NSingleTops        = scaleSingleTops        * fehlerrechnung(SingleTops->GetEntries(Selektion));//"poisson");
  NSingleTopt        = scaleSingleTopt        * fehlerrechnung(SingleTopt->GetEntries(Selektion));//"poisson");
  NSingleToptW       = scaleSingleToptW       * fehlerrechnung(SingleToptW->GetEntries(Selektion));//"poisson");

  fehlerrechnung NDiboson,NDY,NQCD,NW,NSingleTop;
  NDiboson   = NWW + NZZ + NWZ;
  NDY        = NZee + NZmumu + NZtautau;
  NQCD       = NInclusiveMu15 + NBCtoE20to30 + NBCtoE30to80 + NBCtoE80to170 + /*NEMenriched20to30 +*/ NEMenriched30to80 + NEMenriched80to170;
  NW         = NWenu + NWmunu;
  NSingleTop = NSingleTops + NSingleTopt + NSingleToptW;

  fehlerrechnung NUntergrund;
  NUntergrund = NTtbarBG + NDiboson + NDY + NQCD + NW + NSingleTop;
  fehlerrechnung NAlle = NTtbar + NDiboson + NDY + NQCD + NW + NSingleTop;//NSignal + NUntergrund;
  fehlerrechnung NMessung = fehlerrechnung(NAlle.mittel(),"poisson");

  fehlerrechnung Effizienz = NSignal / (scaleTtbar * fehlerrechnung(Ttbar->GetEntries("p_Dilepton"))); //,"poisson"));
  fehlerrechnung Reinheit  = NSignal / NAlle;

  cout << "********Resultate:*********" << endl;
  cout << "Ereignisse in 10pb-1 nach Sample:" << endl;
  cprint("Signal:           ", NSignal, 4);
  cprint("ttbar BG:         ", NTtbarBG, 4);
  cprint("WW:               ", NWW, 4);
  cprint("ZZ:               ", NZZ, 4);
  cprint("WZ:               ", NWZ, 4);
  cprint("Zee:              ", NZee, 4);
  cprint("Zmumu:            ", NZmumu, 4);
  cprint("Ztautau:          ", NZtautau, 4);
  cprint("InclusiveMu15:    ", NInclusiveMu15, 4);
  cprint("BCtoE20to30:      ", NBCtoE20to30, 4);
  cprint("BCtoE30to80:      ", NBCtoE30to80, 4);
  cprint("BCtoE80to170:     ", NBCtoE80to170, 4);
  cprint("EMenriched20to30: ", NEMenriched20to30, 4);
  cprint("EMenriched30to80: ", NEMenriched30to80, 4);
  cprint("EMenriched80to170:", NEMenriched80to170, 4);
  cprint("Wenu:             ", NWenu, 4);
  cprint("Wmunu:            ", NWmunu, 4);
  cprint("SingleTops:       ", NSingleTops, 4);
  cprint("SingleTopt:       ", NSingleTopt, 4);
  cprint("SingleToptW:      ", NSingleToptW, 4);
  cout << "" << endl;
  cout << "Ereignisse in 10pb-1 nach Prozess:" << endl;
  cprint("Signal:           ", NSignal, 4);
  cprint("ttbar BG:         ", NTtbarBG, 4);
  cprint("Diboson:          ", NDiboson, 4);
  cprint("Drell-Yan:        ", NDY, 4);
  cprint("QCD:              ", NQCD, 4);
  cprint("W+Jets:           ", NW, 4);
  cprint("Single Top:       ", NSingleTop, 4);
  cout << "" << endl;
  cout << "Insgesamt:" << endl;
  cprint("Signal:           ", NSignal, 4);
  cprint("Untergrund:       ", NUntergrund, 4);
  cprint("Alle:             ", NAlle, 4);
  cprint("Messung:          ", NMessung, 4);
  cprint("Signaleffizienz:  ", Effizienz, 4);
  cprint("S/B:              ", NSignal / NUntergrund, 4);
  cprint("Reinheit S/(S+B): ", Reinheit , 4);
  cout << "Mit diesen Werten ergibt sich ein Wirkungsquerschnitt ttbar->dilep: sigma = (NMessung * Reinheit) / (Effizienz * Lumi) = " << endl;
  fehlerrechnung sigma = (NMessung * Reinheit) / (Effizienz * 10.);
  cprint(sigma,4);
  cout << "Damit ist der totale ttbar-WQ gegeben durch:" << endl;
  fehlerrechnung ttbarWQ = sigma / 0.0455;
  cprint("PYTHIA-Zerfallsbreite: 0,0468 , sigma(ttbar) = ", sigma / 0.0468, 4);
  cprint("PDG-Zerfallsbreite:    0,0455 , sigma(ttbar) = ", ttbarWQ, 4);
  cout << "Dies entspricht einer Signifikanz von sigma/delta(sigma) = " << ttbarWQ.signifikanz() << endl;
  
}
