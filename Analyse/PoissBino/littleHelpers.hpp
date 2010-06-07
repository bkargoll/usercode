//======================================================================
/**
 * @file        littleHelpers.hpp
 *
 * @author      Bastian Kargoll
 *
 * @version     0.0
 *
 * @brief       multiple small helper functions which are
 *              useful for various macros
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
  TTree *Tree=0; gDirectory->GetObject("Event",Tree);
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
