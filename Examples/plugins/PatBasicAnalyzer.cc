#include <map>
#include <string>
#include <iostream>

#include "TH1.h"
#include "math.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include <DataFormats/Candidate/interface/Candidate.h>

class PatBasicAnalyzer : public edm::EDAnalyzer {

public:
  explicit PatBasicAnalyzer(const edm::ParameterSet&);
  ~PatBasicAnalyzer();

  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  double ET(const reco::Candidate*, const reco::Candidate*, const reco::Candidate*);
  
  // simple map to contain all histograms; 
  // histograms are booked in the beginJob() 
  // method
  std::map<std::string,TH1F*> histContainer_; 

  // input tags  
  edm::InputTag photonSrc_;
  edm::InputTag elecSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag tauSrc_;
  edm::InputTag jetSrc_;
  edm::InputTag metSrc_;

  int vollhad;
  int semilep;
  int dilep;
  int events;
  int selektion[100];
  int leptonen;
  int id;
  int winner;
  int loser;
};

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include <DataFormats/Candidate/interface/const_iterator.h>


PatBasicAnalyzer::PatBasicAnalyzer(const edm::ParameterSet& iConfig):
  histContainer_(),
  photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
  elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc")),
  muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
  tauSrc_(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc" )),
  jetSrc_(iConfig.getUntrackedParameter<edm::InputTag>("jetSrc" )),
  metSrc_(iConfig.getUntrackedParameter<edm::InputTag>("metSrc" ))
{
}

PatBasicAnalyzer::~PatBasicAnalyzer()
{
}

double
PatBasicAnalyzer::ET(const reco::Candidate *a, const reco::Candidate *b, const reco::Candidate *c)
{
double ergebnis = sqrt( (a->mass() + b->mass() + c->mass())*(a->mass() + b->mass() + c->mass()) + (a->px()+b->px()+c->px())*(a->px()+b->px()+c->px()) + (a->py()+b->py()+c->py())*(a->py()+b->py()+c->py()) + (a->pz()+b->pz()+c->pz())*(a->pz()+b->pz()+c->pz()) );

 return ergebnis;
}


void
PatBasicAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{ using namespace std; 

  // get generator particles
  edm::Handle<edm::View<reco::GenParticle> > genParticles;
  iEvent.getByLabel("genParticles",genParticles);

  // get electron collection
  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(elecSrc_,electrons);

  // get muon collection
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_,muons);

  // get tau collection  
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByLabel(tauSrc_,taus);

  // get jet collection
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_,jets);

  // get met collection  
  edm::Handle<edm::View<pat::MET> > mets;
  iEvent.getByLabel(metSrc_,mets);
  
  // get photon collection  
  edm::Handle<edm::View<pat::Photon> > photons;
  iEvent.getByLabel(photonSrc_,photons);
    

  //loop over genParticles
  int id=0;
  int id2=0;
  int id3=0;
  bool WpHadron = false;
  bool WmHadron = false;
  const reco::Candidate *b=0;
  const reco::Candidate *bbar=0;
  const reco::Candidate *Wp=0;
  const reco::Candidate *Wm=0;
  const reco::Candidate *qWp=0;
  const reco::Candidate *qWm=0;
  const reco::Candidate *qbarWp=0;
  const reco::Candidate *qbarWm=0;

  for(edm::View<reco::GenParticle>::const_iterator gen=genParticles->begin(); gen!=genParticles->end(); ++gen){
    id = gen->pdgId();
    if(id==6){
      cout << "Toechter des Tops" << endl;
      for(reco::candidate::const_iterator gen2=gen->begin(); gen2!=gen->end(); ++gen2){
	id2 = gen2->pdgId();
	cout << id2 << endl;
	if(id2==5)  b= &(*gen2);
	if(id2==24){
	  Wp = &(*gen2);
	  for(reco::candidate::const_iterator gen3=gen2->begin(); gen3!=gen2->end(); ++gen3){
	    id3 = gen3->pdgId();
	    if(id3==2 || id3==4) {qWp = &(*gen3); WpHadron = true;}
	    if(id3==-1||id3==-3) qbarWp= &(*gen3);	    
	  }
 	}
      }
    }
    if(id==-6){
      cout << "Toechter des Antitops" << endl;
      for(reco::candidate::const_iterator gen2=gen->begin(); gen2!=gen->end(); ++gen2){
	id2 = gen2->pdgId();
	cout << id2 << endl;
	if(id2==-5)  bbar= &(*gen2);
	if(id2==-24){
	  Wm = &(*gen2);
	  for(reco::candidate::const_iterator gen3=gen2->begin(); gen3!=gen2->end(); ++gen3){
	    id3 = gen3->pdgId();
	    if(id3==-2 || id3==-4) {qbarWm = &(*gen3); WmHadron = true;}
	    if( id3==1 || id3==3)  qWm= &(*gen3);	    
	  }
 	}
      }
    }
  }
    histContainer_["Benergy"]->Fill(b->et());
    histContainer_["Benergy"]->Fill(bbar->et());
    histContainer_["Wenergy"]->Fill(Wp->et());
    histContainer_["Wenergy"]->Fill(Wm->et());
  if( !WpHadron && !WmHadron ){
    cout << "Das Ereignis ist dileptonisch." << endl;
  }
  if( WpHadron && WmHadron ){
      cout << "Das Ereignis ist vollhadronisch." << endl;
      histContainer_["qenergy"]->Fill(qWp->et());
      histContainer_["qenergy"]->Fill(qbarWp->et());
      histContainer_["qenergy"]->Fill(qWm->et());
      histContainer_["qenergy"]->Fill(qbarWm->et());
      histContainer_["diff"]->Fill(qWp->et() - b->et());
      histContainer_["diff"]->Fill(qbarWp->et() - b->et());
      histContainer_["diff"]->Fill(qWm->et() - bbar->et());
      histContainer_["diff"]->Fill(qbarWm->et() - bbar->et());
      histContainer_["wrongdiff"]->Fill(qWp->et() - bbar->et());
      histContainer_["wrongdiff"]->Fill(qbarWp->et() - bbar->et());
      histContainer_["wrongdiff"]->Fill(qWm->et() - b->et());
      histContainer_["wrongdiff"]->Fill(qbarWm->et() - b->et());
      histContainer_["ratio"]->Fill(qWp->et() / b->et());
      histContainer_["ratio"]->Fill(qbarWp->et() / b->et());
      histContainer_["ratio"]->Fill(qWm->et() / bbar->et());
      histContainer_["ratio"]->Fill(qbarWm->et() / bbar->et());
      histContainer_["wrongratio"]->Fill(qWp->et() / bbar->et());
      histContainer_["wrongratio"]->Fill(qbarWp->et() / bbar->et());
      histContainer_["wrongratio"]->Fill(qWm->et() / b->et());
      histContainer_["wrongratio"]->Fill(qbarWm->et() / b->et());
  }
  if( WpHadron && !WmHadron){
    cout << "Das Ereignis ist semileptonisch (W+ -> qq)." << endl;
    histContainer_["qenergy"]->Fill(qWp->et());
    histContainer_["qenergy"]->Fill(qbarWp->et());
    histContainer_["ratio"]->Fill(qWp->et() / b->et());
    histContainer_["ratio"]->Fill(qbarWp->et() / b->et());
    histContainer_["wrongratio"]->Fill(qWp->et() / bbar->et());
    histContainer_["wrongratio"]->Fill(qbarWp->et() / bbar->et());
    histContainer_["diff"]->Fill(qWp->et() - b->et());
    histContainer_["diff"]->Fill(qbarWp->et() - b->et());
    histContainer_["wrongdiff"]->Fill(qWp->et() - bbar->et());
    histContainer_["wrongdiff"]->Fill(qbarWp->et() - bbar->et());
    double ETnoqbar = ET(b,bbar,qWp);
    double ETnoq = ET(b,bbar,qbarWp);
    double ETnobbar = ET(b,qWp,qbarWp);
    double ETnob = ET(bbar,qWp,qbarWp);
    if(ETnobbar > ETnoqbar && ETnobbar > ETnoq && ETnobbar > ETnob){
      cout << "Hier hat's geklappt!                                                                                 +++++++ :-) +++++++"<< endl;
      ++winner;
    }
        else{
      cout << "Hier ist es schiefgegangen.                                                                          ------- :-( -------"<< endl;
      ++loser;
    }
  }
  if( !WpHadron && WmHadron){
    cout << "Das Ereignis ist semileptonisch (W- -> qq)." << endl;
    histContainer_["qenergy"]->Fill(qWm->et());
    histContainer_["qenergy"]->Fill(qbarWm->et());
    histContainer_["ratio"]->Fill(qWm->et() / bbar->et());
    histContainer_["ratio"]->Fill(qbarWm->et() / bbar->et());
    histContainer_["wrongratio"]->Fill(qWm->et() / b->et());
    histContainer_["wrongratio"]->Fill(qbarWm->et() / b->et());
    histContainer_["diff"]->Fill(qWm->et() - bbar->et());
    histContainer_["diff"]->Fill(qbarWm->et() - bbar->et());
    histContainer_["wrongdiff"]->Fill(qWm->et() - b->et());
    histContainer_["wrongdiff"]->Fill(qbarWm->et() - b->et());
    double ETnoqbar = ET(qWm,b,bbar);
    double ETnoq = ET(qbarWm,b,bbar);
    double ETnobbar = ET(qWm,qbarWm,b);
    double ETnob = ET(qWm,qbarWm,bbar);
    if(ETnob > ETnoqbar && ETnob > ETnoq && ETnob > ETnobbar){
      cout << "Hier hat's geklappt!                                                                                 +++++++ :-) +++++++"<< endl;
      ++winner;
    }
        else{
      cout << "Hier ist es schiefgegangen.                                                                          ------- :-( -------"<< endl;
      ++loser;
      }
  }

//   // loop over genParticles                                    genParticle rumspielerei
//   int nTop = 0;
//   int nTopb= 0;
//   int nLepton =0;
//   int nLeptonpt =0;
//   int nElec = 0;
//   int nElecpt = 0;
//   int nMuon =0;
//   int nMuonpt =0;
//   int nTau =0;
//   int nTaupt =0;
//   for(edm::View<reco::GenParticle>::const_iterator gen=genParticles->begin(); gen!=genParticles->end(); ++gen){
//     id = gen->pdgId();
//     if(id == 6 ) ++nTop;
//     if(id == -6) ++nTopb;
//     if(id == 11 || id == -11){
//       ++nLepton;
//       ++nElec;
//       if(gen->pt() > 30){ ++nLeptonpt; ++nElecpt;}
//     }
//     if(id == 13 || id == -13){
//       ++nLepton;
//       ++nMuon;
//       if(gen->pt() > 30){ ++nLeptonpt; ++nMuonpt;}
//     }
//     if(id == 15 || id == -15){
//       ++nLepton;
//       ++nTau;
//       if(gen->pt() > 30){ ++nLeptonpt; ++nTaupt;}
//     }
//   }
//   histContainer_["tops"]->Fill(nTop);
//   histContainer_["antitops"]->Fill(nTopb);
//   histContainer_["tops+antitops"]->Fill(nTop + nTopb);
//   histContainer_["leptons"]->Fill(nLepton);
//   histContainer_["leptonspt"]->Fill(nLeptonpt);

//   // loop over leptons
//   int nLeptonreco = 0;
//   int nLeptonrecopt = 0;
//   int nElecreco =0;
//   int nElecrecopt =0;
//   int nMuonreco =0;
//   int nMuonrecopt =0;
//   int nTaureco =0;
//   int nTaurecopt =0;
//   for(edm::View<pat::Electron>::const_iterator elec=electrons->begin(); elec!=electrons->end(); ++elec){
//     ++nLeptonreco;
//     ++nElecreco;
//     if(elec->pt() > 30){ ++nLeptonrecopt; ++nElecrecopt;}
//   }
//   for(edm::View<pat::Muon>::const_iterator muon =muons->begin(); muon!=muons->end(); ++muon){
//     ++nLeptonreco;
//     ++nMuonreco;
//     if(muon->pt() > 30){ ++nLeptonrecopt; ++nMuonrecopt;}
//   }
//   for(edm::View<pat::Tau>::const_iterator tau=taus->begin(); tau!=taus->end(); ++tau){
//     ++nLeptonreco;
//     ++nTaureco;
//     if(tau->pt() > 30){ ++nLeptonrecopt; ++nTaurecopt;}
//   }
//   histContainer_["deltalept"]->Fill(nLeptonreco - nLepton);
//   histContainer_["deltaleptpt"]->Fill(nLeptonrecopt - nLeptonpt);
//   histContainer_["deltaelec"]->Fill(nElecreco - nElec);
//   histContainer_["deltaelecpt"]->Fill(nElecrecopt - nElecpt);
//   histContainer_["deltamuon"]->Fill(nMuonreco - nMuon);
//   histContainer_["deltamuonpt"]->Fill(nMuonrecopt - nMuonpt);
//   histContainer_["deltatau"]->Fill(nTaureco - nTau);
//   histContainer_["deltataupt"]->Fill(nTaurecopt - nTaupt);

  /*
  // loop over jets
  size_t nJets0=0;
  size_t nJets50=0;
  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
    ++nJets0;
    if(jet->pt()>50)
      ++nJets50;
  }
  histContainer_["jets0"]->Fill(nJets0);
  histContainer_["jets50"]->Fill(nJets50);

  // loops over leptons
  size_t lpt=0;
  size_t l=0;
  double drtemp=101.;
  double deltar=100.;
  double deltarpt=100.;
  double dphi = 0.;

  cout << "Elektronenzahl" << electrons->size() << endl;
  cout << "Myonenzahl" << muons->size() << endl;
  cout << "Tauzahl" << taus->size() << endl;
  // loop over electrons
  for(edm::View<pat::Electron>::const_iterator elec=electrons->begin(); elec!=electrons->end(); ++elec){
  deltar = 100.;
  deltarpt = 100.;
  drtemp = 101.;    
  ++l;
    for(edm::View<pat::Electron>::const_iterator elec2=electrons->begin(); elec2!=elec; ++elec2){
      dphi = elec->phi() - elec2->phi();
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      drtemp = sqrt( ((elec->eta() - elec2->eta())*(elec->eta() - elec2->eta())) + (dphi * dphi) );
      cout << "e-e drtemp = " << drtemp << endl;
      if(drtemp < deltar){
	deltar = drtemp;
	if(elec2->pt() >30) deltarpt = drtemp;
      }
    }
    for(edm::View<pat::Muon>::const_iterator muon2=muons->begin(); muon2!=muons->end(); ++muon2){
      dphi = elec->phi() - muon2->phi();
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      drtemp = sqrt( ((elec->eta() - muon2->eta())*(elec->eta() - muon2->eta())) + (dphi * dphi) );
      cout << "e-m drtemp = " << drtemp << endl;
      if(drtemp < deltar){
	deltar = drtemp;
	if(muon2->pt() >30) deltarpt = drtemp;
      }
    }
    for(edm::View<pat::Tau>::const_iterator tau2=taus->begin(); tau2!=taus->end(); ++tau2){
      dphi = elec->phi() - tau2->phi();
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      drtemp = sqrt( ((elec->eta() - tau2->eta())*(elec->eta() - tau2->eta())) + (dphi * dphi) );
      cout << "e-t drtemp = " << drtemp << endl;
      if(drtemp < deltar){
	deltar = drtemp;
	if(tau2->pt() >30) deltarpt = drtemp;
      }
    }
//     for(edm::View<pat::Photon>::const_iterator phot2=photons->begin(); phot2!=photons->end(); ++phot2){
//       dphi = elec->phi() - phot2->phi();
//       if(dphi > M_PI) dphi = 2*M_PI - dphi;
//       drtemp = sqrt( ((elec->eta() - phot2->eta())*(elec->eta() - phot2->eta())) + (dphi * dphi) );
//       if(drtemp < deltar){
// 	deltar = drtemp;
// 	if(phot2->pt() >30) deltarpt = drtemp;
//       }
//     }
//     for(edm::View<pat::Jet>::const_iterator jet2=jets->begin(); jet2!=jets->end(); ++jet2){
//       dphi = elec->phi() - jet2->phi();
//       if(dphi > M_PI) dphi = 2*M_PI - dphi;
//       drtemp = sqrt( ((elec->eta() - jet2->eta())*(elec->eta() - jet2->eta())) + (dphi * dphi) );
//       if(drtemp < deltar){
// 	deltar = drtemp;
// 	if(jet2->pt() >30) deltarpt = drtemp;
//       }
//     }
    if(elec->pt()>30){
      ++lpt;
      histContainer_["drelecpt"]->Fill(deltarpt);
      histContainer_["drallept"]->Fill(deltarpt);
      histContainer_["trackIsopt"]->Fill(elec->trackIso());
      if(deltarpt >80.) cout << "ALARMALARMALARM!!!!!! beim electron deltarpt = " << deltarpt << endl;
    }
    histContainer_["drelec"]->Fill(deltar);
    histContainer_["dralle"]->Fill(deltar);
    histContainer_["eta"]->Fill(elec->eta());
    histContainer_["phi"]->Fill(elec->phi());
    histContainer_["trackIso"]->Fill(elec->trackIso());
    if(deltar > 80.) cout << "ALARMALARMALARM!!!!!! beim electron deltar = " << deltar << endl;
  }
  // loop over muons
  for(edm::View<pat::Muon>::const_iterator muon =muons->begin(); muon!=muons->end(); ++muon){
  deltar = 100.;
  deltarpt = 100.;
  drtemp = 101.;    
  ++l;
    for(edm::View<pat::Electron>::const_iterator elec2=electrons->begin(); elec2!=electrons->end(); ++elec2){
      dphi = muon->phi() - elec2->phi();
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      drtemp = sqrt( ((muon->eta() - elec2->eta())*(muon->eta() - elec2->eta())) + (dphi * dphi) );
      cout << "m-e drtemp = " << drtemp << endl;
      if(drtemp < deltar){
	deltar = drtemp;
	if(elec2->pt() >30) deltarpt = drtemp;
      }
    }
    for(edm::View<pat::Muon>::const_iterator muon2=muons->begin(); muon2!=muon; ++muon2){
      dphi = muon->phi() - muon2->phi();
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      drtemp = sqrt( ((muon->eta() - muon2->eta())*(muon->eta() - muon2->eta())) + (dphi * dphi) );
      cout << "m-m drtemp = " << drtemp << endl;
      if(drtemp < deltar){
	deltar = drtemp;
	if(muon2->pt() >30) deltarpt = drtemp;
      }
    }
    for(edm::View<pat::Tau>::const_iterator tau2=taus->begin(); tau2!=taus->end(); ++tau2){
      dphi = muon->phi() - tau2->phi();
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      drtemp = sqrt( ((muon->eta() - tau2->eta())*(muon->eta() - tau2->eta())) + (dphi * dphi) );
      cout << "m-t drtemp = " << drtemp << endl;
      if(drtemp < deltar){
	deltar = drtemp;
	if(tau2->pt() >30) deltarpt = drtemp;
      }
    }
//     for(edm::View<pat::Photon>::const_iterator phot2=photons->begin(); phot2!=photons->end(); ++phot2){
//       dphi = muon->phi() - phot2->phi();
//       if(dphi > M_PI) dphi = 2*M_PI - dphi;
//       drtemp = sqrt( ((muon->eta() - phot2->eta())*(muon->eta() - phot2->eta())) + (dphi * dphi) );
//       if(drtemp < deltar){
// 	deltar = drtemp;
// 	if(phot2->pt() >30) deltarpt = drtemp;
//       }
//     }
//     for(edm::View<pat::Jet>::const_iterator jet2=jets->begin(); jet2!=jets->end(); ++jet2){
//       dphi = muon->phi() - jet2->phi();
//       if(dphi > M_PI) dphi = 2*M_PI - dphi;
//       drtemp = sqrt( ((muon->eta() - jet2->eta())*(muon->eta() - jet2->eta())) + (dphi * dphi) );
//       if(drtemp < deltar){
// 	deltar = drtemp;
// 	if(jet2->pt() >30) deltarpt = drtemp;
//       }
//     }
    if(muon->pt()>30){
      ++lpt;
      histContainer_["drmuonpt"]->Fill(deltarpt);
      histContainer_["drallept"]->Fill(deltarpt);
      histContainer_["trackIsopt"]->Fill(muon->trackIso());
      if(deltarpt >80.) cout << "ALARMALARMALARM!!!!!! beim muon deltarpt = " << deltarpt << endl;
    }
    histContainer_["drmuon"]->Fill(deltar);
    histContainer_["dralle"]->Fill(deltar);
    histContainer_["eta"]->Fill(muon->eta());
    histContainer_["phi"]->Fill(muon->phi());
    histContainer_["trackIso"]->Fill(muon->trackIso());
    if(deltar >80.) cout << "ALARMALARMALARM!!!!!! beim muon deltar = " << deltar << endl;
  }
  // loop over taus
  for(edm::View<pat::Tau>::const_iterator tau=taus->begin(); tau!=taus->end(); ++tau){
  drtemp = 101.;
  deltar = 100.;
  deltarpt = 100.;
  ++l;
    for(edm::View<pat::Electron>::const_iterator elec2=electrons->begin(); elec2!=electrons->end(); ++elec2){
      dphi = tau->phi() - elec2->phi();
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      drtemp = sqrt( ((tau->eta() - elec2->eta())*(tau->eta() - elec2->eta())) + (dphi * dphi) );
      cout << "t-e drtemp = " << drtemp << endl;
      if(drtemp < deltar){
	deltar = drtemp;
 	if(elec2->pt() >30) deltarpt = drtemp;
     }
    }    histContainer_["tops"]->Fill(nTop);
    histContainer_["antitops"]->Fill(nTopb);
    histContainer_["leptons"]->Fill(nLepton);
    histContainer_["antileptons"]->Fill(nLeptonb);
    for(edm::View<pat::Muon>::const_iterator muon2=muons->begin(); muon2!=muons->end(); ++muon2){
      dphi = tau->phi() - muon2->phi();
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      drtemp = sqrt( ((tau->eta() - muon2->eta())*(tau->eta() - muon2->eta())) + (dphi * dphi) );
      cout << "t-m drtemp = " << drtemp << endl;
      if(drtemp < deltar){
	deltar = drtemp;
	if(muon2->pt() >30) deltarpt = drtemp;
      }
    }
    for(edm::View<pat::Tau>::const_iterator tau2=taus->begin(); tau2!=tau; ++tau2){
      dphi = tau->phi() - tau2->phi();
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      drtemp = sqrt( ((tau->eta() - tau2->eta())*(tau->eta() - tau2->eta())) + (dphi * dphi) );
      cout << "t-t drtemp = " << drtemp << endl;
      if(drtemp < deltar){
	deltar = drtemp;
	if(tau2->pt() >30) deltarpt = drtemp;
      }
    }
//     for(edm::View<pat::Photon>::const_iterator phot2=photons->begin(); phot2!=photons->end(); ++phot2){
//       dphi = tau->phi() - phot2->phi();
//       if(dphi > M_PI) dphi = 2*M_PI - dphi;
//       drtemp = sqrt( ((tau->eta() - phot2->eta())*(tau->eta() - phot2->eta())) + (dphi * dphi) );
//       if(drtemp < deltar){
// 	deltar = drtemp;
// 	if(phot2->pt() >30) deltarpt = drtemp;
//       }
//     }
//     for(edm::View<pat::Jet>::const_iterator jet2=jets->begin(); jet2!=jets->end(); ++jet2){
//       dphi = tau->phi() - jet2->phi();
//       if(dphi > M_PI) dphi = 2*M_PI - dphi;
//       drtemp = sqrt( ((tau->eta() - jet2->eta())*(tau->eta() - jet2->eta())) + (dphi * dphi) );
//       if(drtemp < deltar){
// 	deltar = drtemp;
// 	if(jet2->pt() >30) deltarpt = drtemp;
//       }
//     }
    if(tau->pt()>30){
      ++lpt;
      histContainer_["drtaupt"]->Fill(deltarpt);
      histContainer_["drallept"]->Fill(deltarpt);
      histContainer_["trackIsopt"]->Fill(tau->trackIso());
      if(deltarpt >80.) cout << "ALARMALARMALARM!!!!!! beim tau deltarpt = " << deltarpt << endl;
    }
    histContainer_["drtau"]->Fill(deltar);
    histContainer_["dralle"]->Fill(deltar);
    histContainer_["eta"]->Fill(tau->eta());
    histContainer_["phi"]->Fill(tau->phi());
    histContainer_["trackIso"]->Fill(tau->trackIso());
    if(deltar >80.) cout << "ALARMALARMALARM!!!!!! beim tau deltar = " << deltarpt << endl;
  }
  leptonen=leptonen+l;
  histContainer_["leptons"]->Fill(l);
  histContainer_["leptons_pt"]->Fill(lpt);
*/
  
  // Events hochzählen
  ++events;
}

void 
PatBasicAnalyzer::beginJob(const edm::EventSetup&)
{
  // register to the TFileService
  edm::Service<TFileService> fs;
  
  // book histograms:
  histContainer_["Benergy"]=fs->make<TH1F>("Benergy",   "transverse energy of b quarks"  , 25, 0., 250.);
  histContainer_["Wenergy"]=fs->make<TH1F>("Wenergy",   "transverse energy of W bosons"  , 25, 0., 250.);
  histContainer_["qenergy"]=fs->make<TH1F>("qenergy",   "transverse energy of light quarks"  , 25, 0., 250.);
  histContainer_["ratio"]=fs->make<TH1F>("ratio",   "et(q)/et(b) ratio"  , 25, 0., 5.);
  histContainer_["wrongratio"]=fs->make<TH1F>("wrongratio",   "et(q)/et(b) ratio (different branches)"  , 25, 0., 5.);
  histContainer_["diff"]=fs->make<TH1F>("diff",   "et(q)-et(b)"  , 50, -200., 200.);
  histContainer_["wrongdiff"]=fs->make<TH1F>("wrongdiff",   "et(q)-et(b) (different branches)"  , 50, -200., 200.);
//    histContainer_["tops" ]=fs->make<TH1F>("tops",  "top  multiplicity",     10, 0,  10);
//    histContainer_["antitops" ]=fs->make<TH1F>("antitops",  "antitop  multiplicity",     10, 0,  10);
//    histContainer_["leptons" ]=fs->make<TH1F>("leptons",  "lepton  multiplicity",     10, 0,  10);
//    histContainer_["tops+antitops"]=fs->make<TH1F>("tops+antitops", "t+tbar multiplicity", 10,0,10);
//    histContainer_["deltalept"]=fs->make<TH1F>("deltalept", "#recoLeptons - #genLeptons", 20,-10,10);
//    histContainer_["deltaleptpt"]=fs->make<TH1F>("deltaleptpt", "#recoLeptons - #genLeptons (pt>30)", 20,-10,10);
//    histContainer_["leptonspt"]=fs->make<TH1F>("leptonspt", "l+lbar (pt>30) multiplicity", 20,0,20);
//    histContainer_["deltaelec"]=fs->make<TH1F>("deltaelec", "#recoElectrons - #genElectrons", 20,-10,10);
//    histContainer_["deltaelecpt"]=fs->make<TH1F>("deltaelecpt", "#recoElectrons - #genElectrons (pt>30)", 20,-10,10);
//    histContainer_["deltamuon"]=fs->make<TH1F>("deltamuon", "#recoMuons - #genMuons", 20,-10,10);
//    histContainer_["deltamuonpt"]=fs->make<TH1F>("deltamuonpt", "#recoMuons - #genMuons (pt>30)", 20,-10,10);
//    histContainer_["deltatau"]=fs->make<TH1F>("deltatau", "#recoTaus - #genTaus", 20,-10,10);
//    histContainer_["deltataupt"]=fs->make<TH1F>("deltataupt", "#recoTaus - #genTaus (pt>30)", 20,-10,10);
//   histContainer_["jets50" ]=fs->make<TH1F>("jets50",  "jet pt>50  multiplicity",     10, 0,  10);
//   histContainer_["jets0"  ]=fs->make<TH1F>("jets0",   "jet pt>0  multiplicity",      40, 0,  40);
//   histContainer_["leptons_pt" ]=fs->make<TH1F>("leptons_pt",  "leptons pt>30  multiplicity",     10, 0,  10);
//   histContainer_["leptons"  ]=fs->make<TH1F>("leptons",   "leptons pt>0  multiplicity",      40, 0,  40);
//   histContainer_["eta"    ]=fs->make<TH1F>("eta",     "leptons eta"           ,      30,  -3., 3. );
//   histContainer_["phi"    ]=fs->make<TH1F>("phi",     "leptons phi"           ,     100,-3.5 , 3.5);  
//   histContainer_["drelec" ]=fs->make<TH1F>("drelec",  "electrons delta r min" ,      50, 0., 5.0);
//   histContainer_["drmuon" ]=fs->make<TH1F>("drmuon",  "muons delta r min" ,      50, 0., 5.0);
//   histContainer_["drtau" ]=fs->make<TH1F>("drtau",  "taus delta r min" ,      50, 0., 5.0);
//   histContainer_["dralle" ]=fs->make<TH1F>("dralle",  "all leptons delta r min" ,      50,0., 5.0);
//   histContainer_["drelecpt" ]=fs->make<TH1F>("drelecpt",  "electrons (pt>30) delta r min" ,      50, 0., 5.0);
//   histContainer_["drmuonpt" ]=fs->make<TH1F>("drmuonpt",  "muons(pt>30) delta r min" ,      50, 0., 5.0);
//   histContainer_["drtaupt" ]=fs->make<TH1F>("drtaupt",  "taus(pt>30) delta r min" ,      50, 0., 5.0);
//   histContainer_["drallept" ]=fs->make<TH1F>("drallept",  "all leptons(pt>30) delta r min" ,      50,0., 5.0);
//   histContainer_["trackIso" ]=fs->make<TH1F>("trackIso",  "all leptons trackIso" ,      50,-2., 10.);
//   histContainer_["trackIsopt"]=fs->make<TH1F>("trackIsopt","all leptons (pt>30) trackIso", 50,-2.,10.);

  vollhad=0;
  semilep=0;
  dilep=0;
  events=0;
  leptonen=0;
  winner=0;
  loser=0;

}

void 
PatBasicAnalyzer::endJob() 
{
  using namespace std;
  cout << "Es wurden " << winner << " Mal die richtigen Quarks ausgewählt, aber auch " << loser << " Mal die falschen." << endl;
  cout << "Das macht eine Trefferquote von " << double(winner)/double(winner+loser)*100. << "% im Vergleich zu 25% durch raten." << endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatBasicAnalyzer);
