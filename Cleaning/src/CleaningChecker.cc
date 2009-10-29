// -*- C++ -*-
//
// 
/**\   Tools/Cleaning/src/CleaningChecker.cc

 Description: basic analyzer to test if cleaning works

 Note: 
*/
//
// Original Author:  Bastian Kargoll
//         Created:  Thu Oct  29 10:30 CEST 2009
//
//

// system include files
#include <memory>
#include <iostream>
#include <string>

// root include files
#include "TH1.h"

// CMSSW include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/Electron.h>
#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/Tau.h>
//
// class decleration
//
class CleaningChecker : public edm::EDAnalyzer {
   public:
      explicit CleaningChecker(const edm::ParameterSet&);
      ~CleaningChecker();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // map for histogramms
      std::map<std::string,TH1F*> histos_;

  // global counters
  int nElecClean;
  int nMuonClean;
  int nTauClean;
  int nJetClean;
  int nPhotonClean;

  // input tags  
  edm::InputTag photonSrc_;
  edm::InputTag elecSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag tauSrc_;
  edm::InputTag jetSrc_;
  edm::InputTag cleanphotonSrc_;
  edm::InputTag cleanelecSrc_;
  edm::InputTag cleanmuonSrc_;
  edm::InputTag cleantauSrc_;
  edm::InputTag cleanjetSrc_;
};

//
// constructors and destructor
//
CleaningChecker::CleaningChecker(const edm::ParameterSet& iConfig):
  //histos_(),
  	photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
  	elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
	muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
        tauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc" )),
        jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" )),
  	cleanphotonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("cleanphotonSrc")),
  	cleanelecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("cleanelectronSrc")),
	cleanmuonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("cleanmuonSrc")),
        cleantauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("cleantauSrc" )),
        cleanjetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("cleanjetSrc" ))
{

}


CleaningChecker::~CleaningChecker()
{
 
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
CleaningChecker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // get electron collections
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(elecSrc_,electrons);
  edm::Handle<edm::View<pat::Electron> > cleanelectrons;
  iEvent.getByLabel(cleanelecSrc_,cleanelectrons);

  // get muon collections
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);
  edm::Handle<edm::View<pat::Muon> > cleanmuons;
  iEvent.getByLabel(cleanmuonSrc_,cleanmuons);

  // get tau collections
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByLabel(tauSrc_,taus);
  edm::Handle<edm::View<pat::Tau> > cleantaus;
  iEvent.getByLabel(cleantauSrc_,cleantaus);

  // get jet collection
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);
  edm::Handle<edm::View<pat::Jet> > cleanjets;
  iEvent.getByLabel(cleanjetSrc_,cleanjets);
  
  // get photon collection  
  edm::Handle<edm::View<pat::Photon> > photons;
  iEvent.getByLabel(photonSrc_,photons);
  edm::Handle<edm::View<pat::Photon> > cleanphotons;
  iEvent.getByLabel(cleanphotonSrc_,cleanphotons);

  // fill histograms with size of collections
  histos_["nElectrons"]->Fill(electrons->size());
  histos_["nCleanElectrons"]->Fill(cleanelectrons->size());
  histos_["nMuons"]->Fill(muons->size());
  histos_["nCleanMuons"]->Fill(cleanmuons->size());
  histos_["nTaus"]->Fill(taus->size());
  histos_["nCleanTaus"]->Fill(cleantaus->size());
  histos_["nJets"]->Fill(jets->size());
  histos_["nCleanJets"]->Fill(cleanjets->size());
  histos_["nPhotons"]->Fill(photons->size());
  histos_["nCleanPhotons"]->Fill(cleanphotons->size());

  // print out if something went terribly wrong
  if ( electrons->size() < cleanelectrons->size() ) std::cout << "MISTAKE: more electrons after cleaning than before" << std::endl;
  if ( muons->size() < cleanmuons->size() ) std::cout << "MISTAKE: more muons after cleaning than before" << std::endl;
  if ( taus->size() < cleantaus->size() ) std::cout << "MISTAKE: more taus after cleaning than before" << std::endl;
  if ( jets->size() < cleanjets->size() ) std::cout << "MISTAKE: more jets after cleaning than before" << std::endl;
  if ( photons->size() < cleanphotons->size() ) std::cout << "MISTAKE: more photons after cleaning than before" << std::endl;

  // count in how many events particles where dropped during cleaning
  if ( electrons->size() > cleanelectrons->size() ) ++nElecClean;
  if ( muons->size() > cleanmuons->size() ) ++nMuonClean;
  if ( taus->size() > cleantaus->size() ) ++nTauClean;
  if ( jets->size() > cleanjets->size() ) ++nJetClean;
  if ( photons->size() > cleanphotons->size() ) ++nPhotonClean;

  // produce histograms showing distance to next particles
  for(edm::View<pat::Electron>::const_iterator elec =electrons->begin(); elec!=electrons->end(); ++elec){
    for(edm::View<pat::Muon>::const_iterator muon =muons->begin(); muon!=muons->end(); ++muon){
      histos_["dRElecMuon"]->Fill(  DeltaR<pat::Electron, pat::Muon>()(*elec, *muon));
    }
  }
  for(edm::View<pat::Electron>::const_iterator elec =cleanelectrons->begin(); elec!=cleanelectrons->end(); ++elec){
    for(edm::View<pat::Muon>::const_iterator muon =cleanmuons->begin(); muon!=cleanmuons->end(); ++muon){
      histos_["dRElecMuonClean"]->Fill(  DeltaR<pat::Electron, pat::Muon>()(*elec, *muon));
    }
  }
  for(edm::View<pat::Tau>::const_iterator tau =taus->begin(); tau!=taus->end(); ++tau){
    for(edm::View<pat::Muon>::const_iterator muon =muons->begin(); muon!=muons->end(); ++muon){
      histos_["dRTauMuon"]->Fill(  DeltaR<pat::Tau, pat::Muon>()(*tau, *muon));
    }
    for(edm::View<pat::Electron>::const_iterator elec =electrons->begin(); elec!=electrons->end(); ++elec){
      histos_["dRTauElec"]->Fill(  DeltaR<pat::Tau, pat::Electron>()(*tau, *elec));
    }
  }
  for(edm::View<pat::Tau>::const_iterator tau =cleantaus->begin(); tau!=cleantaus->end(); ++tau){
    for(edm::View<pat::Muon>::const_iterator muon =cleanmuons->begin(); muon!=cleanmuons->end(); ++muon){
      histos_["dRTauMuonClean"]->Fill(  DeltaR<pat::Tau, pat::Muon>()(*tau, *muon));
    }
    for(edm::View<pat::Electron>::const_iterator elec =cleanelectrons->begin(); elec!=cleanelectrons->end(); ++elec){
      histos_["dRTauElecClean"]->Fill(  DeltaR<pat::Tau, pat::Electron>()(*tau, *elec));
    }
  }
  for(edm::View<pat::Photon>::const_iterator phot =photons->begin(); phot!=photons->end(); ++phot){
    for(edm::View<pat::Electron>::const_iterator elec =electrons->begin(); elec!=electrons->end(); ++elec){
      histos_["dRPhotElec"]->Fill(  DeltaR<pat::Photon, pat::Electron>()(*phot, *elec));
    }
  }
  for(edm::View<pat::Photon>::const_iterator phot =cleanphotons->begin(); phot!=cleanphotons->end(); ++phot){
    for(edm::View<pat::Electron>::const_iterator elec =cleanelectrons->begin(); elec!=cleanelectrons->end(); ++elec){
      histos_["dRPhotElecClean"]->Fill(  DeltaR<pat::Photon, pat::Electron>()(*phot, *elec));
    }
  }
  for(edm::View<pat::Jet>::const_iterator jet =jets->begin(); jet!=jets->end(); ++jet){
    for(edm::View<pat::Muon>::const_iterator muon =muons->begin(); muon!=muons->end(); ++muon){
      histos_["dRJetMuon"]->Fill(  DeltaR<pat::Jet, pat::Muon>()(*jet, *muon));
    }
    for(edm::View<pat::Electron>::const_iterator elec =electrons->begin(); elec!=electrons->end(); ++elec){
      histos_["dRJetElec"]->Fill(  DeltaR<pat::Jet, pat::Electron>()(*jet, *elec));
    }
    for(edm::View<pat::Tau>::const_iterator tau =taus->begin(); tau!=taus->end(); ++tau){
      histos_["dRJetTau"]->Fill(  DeltaR<pat::Jet, pat::Tau>()(*jet, *tau));
    }
    for(edm::View<pat::Photon>::const_iterator phot =photons->begin(); phot!=photons->end(); ++phot){
      histos_["dRJetPhot"]->Fill(  DeltaR<pat::Jet, pat::Photon>()(*jet, *phot));
    }
  }
  for(edm::View<pat::Jet>::const_iterator jet =cleanjets->begin(); jet!=cleanjets->end(); ++jet){
    for(edm::View<pat::Muon>::const_iterator muon =cleanmuons->begin(); muon!=cleanmuons->end(); ++muon){
      histos_["dRJetMuonClean"]->Fill(  DeltaR<pat::Jet, pat::Muon>()(*jet, *muon));
    }
    for(edm::View<pat::Electron>::const_iterator elec =cleanelectrons->begin(); elec!=cleanelectrons->end(); ++elec){
      histos_["dRJetElecClean"]->Fill(  DeltaR<pat::Jet, pat::Electron>()(*jet, *elec));
    }
    for(edm::View<pat::Tau>::const_iterator tau =cleantaus->begin(); tau!=cleantaus->end(); ++tau){
      histos_["dRJetTauClean"]->Fill(  DeltaR<pat::Jet, pat::Tau>()(*jet, *tau));
    }
    for(edm::View<pat::Photon>::const_iterator phot =cleanphotons->begin(); phot!=cleanphotons->end(); ++phot){
      histos_["dRJetPhotClean"]->Fill(  DeltaR<pat::Jet, pat::Photon>()(*jet, *phot));
    }
  }

}



// ------------ method called once each job just before starting event loop  ------------
void 
CleaningChecker::beginJob(const edm::EventSetup&)
{
  // register to the TFileService
  edm::Service<TFileService> fs;

  // book histograms:
  histos_[  "nElectrons"  ]=fs->make<TH1F>( "nElectrons"   ,   "Number of electrons in event"       , 15, 0, 15);
  histos_["nCleanElectrons"]=fs->make<TH1F>( "nCleanElectrons"   ,   "Number of electrons in cleaned event"       , 15, 0, 15);
  histos_[  "nMuons"  ]=fs->make<TH1F>( "nMuons"   ,   "Number of muons in event"       , 15, 0, 15);
  histos_["nCleanMuons"]=fs->make<TH1F>( "nCleanMuons"   ,   "Number of muons in cleaned event"       , 15, 0, 15);
  histos_[  "nTaus"  ]=fs->make<TH1F>( "nTaus"   ,   "Number of taus in event"       , 15, 0, 15);
  histos_["nCleanTaus"]=fs->make<TH1F>( "nCleanTaus"   ,   "Number of taus in cleaned event"       , 15, 0, 15);
  histos_[  "nJets"  ]=fs->make<TH1F>( "nJets"   ,   "Number of jets in event"       , 15, 0, 15);
  histos_["nCleanJets"]=fs->make<TH1F>( "nCleanJets"   ,   "Number of jets in cleaned event"       , 15, 0, 15);
  histos_[  "nPhotons"  ]=fs->make<TH1F>( "nPhotons"   ,   "Number of photons in event"       , 15, 0, 15);
  histos_["nCleanPhotons"]=fs->make<TH1F>( "nCleanPhotons"   ,   "Number of photons in cleaned event"       , 15, 0, 15);
  histos_[  "dRElecMuon"  ]=fs->make<TH1F>( "dRElecMuon"   ,   "DeltaR between electrons and muons"       ,20 , 0.,5. );
  histos_["dRElecMuonClean"]=fs->make<TH1F>( "dRElecMuonClean"   , "DeltaR between clean electrons and muons"       , 20, 0., 5.);
  histos_[  "dRTauMuon"  ]=fs->make<TH1F>( "dRTauMuon"   ,   "DeltaR between taus and muons"       ,20 , 0.,5. );
  histos_["dRTauMuonClean"]=fs->make<TH1F>( "dRTauMuonClean"   , "DeltaR between clean taus and muons"       , 20, 0., 5.);
  histos_[  "dRTauElec"  ]=fs->make<TH1F>( "dRTauElec"   ,   "DeltaR between taus and electrons"       ,20 , 0.,5. );
  histos_["dRTauElecClean"]=fs->make<TH1F>( "dRTauElecClean"   , "DeltaR between clean taus and electrons"       , 20, 0., 5.);
  histos_[  "dRPhotElec"  ]=fs->make<TH1F>( "dRPhotElec"   ,   "DeltaR between photons and electrons"       ,20 , 0.,5. );
  histos_["dRPhotElecClean"]=fs->make<TH1F>( "dRPhotElecClean"   , "DeltaR between clean photons and electrons"       , 20, 0., 5.);
  histos_[  "dRJetElec"  ]=fs->make<TH1F>( "dRJetElec"   ,   "DeltaR between jets and electrons"       ,20 , 0.,5. );
  histos_["dRJetElecClean"]=fs->make<TH1F>( "dRJetElecClean"   , "DeltaR between clean jets and electrons"       , 20, 0., 5.);
  histos_[  "dRJetMuon"  ]=fs->make<TH1F>( "dRJetMuon"   ,   "DeltaR between jets and muons"       ,20 , 0.,5. );
  histos_["dRJetMuonClean"]=fs->make<TH1F>( "dRJetMuonClean"   , "DeltaR between clean jets and muons"       , 20, 0., 5.);
  histos_[  "dRJetTau"  ]=fs->make<TH1F>( "dRJetTau"   ,   "DeltaR between jets and taus"       ,20 , 0.,5. );
  histos_["dRJetTauClean"]=fs->make<TH1F>( "dRJetTauClean"   , "DeltaR between clean jets and taus"       , 20, 0., 5.);
  histos_[  "dRJetPhot"  ]=fs->make<TH1F>( "dRJetPhot"   ,   "DeltaR between jets and photons"       ,20 , 0.,5. );
  histos_["dRJetPhotClean"]=fs->make<TH1F>( "dRJetPhotClean"   , "DeltaR between clean jets and photons"       , 20, 0., 5.);

  // set global counters
  nElecClean = 0;
  nMuonClean = 0;
  nTauClean = 0;
  nJetClean = 0;
  nPhotonClean = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CleaningChecker::endJob() {
  using namespace std;
  // print out results
  cout << "Electrons were dropped during cleaning in " << nElecClean << " events. (Should be >0)" << endl;
  cout << "Muons were dropped during cleaning in " << nMuonClean << " events. (Should be 0)" << endl;
  cout << "Taus were dropped during cleaning in " << nTauClean << " events. (Should be >0)" << endl;
  cout << "Jets were dropped during cleaning in " << nJetClean << " events. (Should be >0)" << endl;
  cout << "Photons were dropped during cleaning in " << nPhotonClean << " events. (Should be 0, as photons have been dropped during PAT-production)" << endl;
  cout << "Some histograms have been produced and stored in testCleaning.root to compare collections before and after cleaning." << endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(CleaningChecker);
