// -*- C++ -*-
//
// Package:    GenTauChecker
// Class:      GenTauChecker
// 
/**\class GenTauChecker GenTauChecker.cc GenLevelTools/GenTauChecker/src/GenTauChecker.cc

 Description: investigate taus on generator level

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Bastian Kargoll
//         Created:  Thu Jun 16 15:18:38 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/View.h"

#include "Math/Vector4D.h"
#include "TTree.h"
#include "TMath.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

//
// class declaration
//

class GenTauChecker: public edm::EDAnalyzer {
public:
	explicit GenTauChecker(const edm::ParameterSet&);
	~GenTauChecker();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	void checkDaughters(reco::Candidate const*, bool&, bool&, unsigned int&,
			unsigned int&, unsigned int&, unsigned int&, bool&);

private:
	virtual void beginJob();
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void endRun(edm::Run const&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock const&,
			edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock const&,
			edm::EventSetup const&);

	// ----------member data ---------------------------
	edm::InputTag genParticles_;
	TTree* tree;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* tau_p4;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* tau_p4Vis;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* tau_p4Invis;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* z_p4;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* z_p4_byDecay;
	std::vector<bool>* tau_hadronic;
	std::vector<bool>* tau_muon;
	std::vector<bool>* tau_electron;
	std::vector<unsigned int>* tau_nDaughters;
	std::vector<unsigned int>* tau_prongs;
	std::vector<unsigned int>* tau_piZeros;
	std::vector<unsigned int>* tau_KProngs;
	std::vector<unsigned int>* tau_KZeros;
	std::vector<int>* tau_mother;
	std::vector<int>* tau_grandma;
	std::vector<double>* deltaR_TrueVis;
	std::vector<double>* deltaR_VisInvis;
	std::vector<double>* deltaR_InvisTrue;
	double deltaR_taus;
	double deltaPhi_taus;
	double deltaEta_taus;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenTauChecker::GenTauChecker(const edm::ParameterSet& iConfig) :
	genParticles_(iConfig.getParameter<edm::InputTag> ("genParticles")) {

	// use the TFileService
	edm::Service<TFileService> fs;
	// create the Tree
	tree = fs->make<TTree> ("Tree", "Check Gen Taus");

	// create all the needed vectors
	tau_p4 = new std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<
			double> > >;
	tau_p4Vis = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	tau_p4Invis = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	z_p4 = new std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<
			double> > >;
	z_p4_byDecay = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	tau_hadronic = new std::vector<bool>;
	tau_muon = new std::vector<bool>;
	tau_electron = new std::vector<bool>;
	tau_nDaughters = new std::vector<unsigned int>;
	tau_prongs = new std::vector<unsigned int>;
	tau_piZeros = new std::vector<unsigned int>;
	tau_KProngs = new std::vector<unsigned int>;
	tau_KZeros = new std::vector<unsigned int>;
	tau_mother = new std::vector<int>;
	tau_grandma = new std::vector<int>;
	deltaR_TrueVis = new std::vector<double>;
	deltaR_VisInvis = new std::vector<double>;
	deltaR_InvisTrue = new std::vector<double>;

	// register the branches
	tree->Branch(
			"tau_p4",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&tau_p4);
	tree->Branch(
			"tau_p4Vis",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&tau_p4Vis);
	tree->Branch(
			"tau_p4Invis",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&tau_p4Invis);
	tree->Branch(
			"z_p4",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&z_p4);
	tree->Branch(
			"z_p4_byDecay",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&z_p4_byDecay);
	tree->Branch("tau_hadronic", "std::vector<bool>", &tau_hadronic);
	tree->Branch("tau_muon", "std::vector<bool>", &tau_muon);
	tree->Branch("tau_electron", "std::vector<bool>", &tau_electron);
	tree->Branch("tau_nDaughters", "std::vector<unsigned int>", &tau_nDaughters);
	tree->Branch("tau_prongs", "std::vector<unsigned int>", &tau_prongs);
	tree->Branch("tau_piZeros", "std::vector<unsigned int>", &tau_piZeros);
	tree->Branch("tau_KProngs", "std::vector<unsigned int>", &tau_KProngs);
	tree->Branch("tau_KZeros", "std::vector<unsigned int>", &tau_KZeros);
	tree->Branch("tau_mother", "std::vector<int>", &tau_mother);
	tree->Branch("tau_grandma", "std::vector<int>", &tau_grandma);
	tree->Branch("deltaR_TrueVis", "std::vector<double>", &deltaR_TrueVis);
	tree->Branch("deltaR_VisInvis", "std::vector<double>", &deltaR_VisInvis);
	tree->Branch("deltaR_InvisTrue", "std::vector<double>", &deltaR_InvisTrue);
	tree->Branch("deltaR_taus", &deltaR_taus, "deltaR_taus/D");
	tree->Branch("deltaPhi_taus", &deltaPhi_taus, "deltaPhi_taus/D");
	tree->Branch("deltaEta_taus", &deltaEta_taus, "deltaEta_taus/D");
}

GenTauChecker::~GenTauChecker() {

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void GenTauChecker::analyze(const edm::Event& iEvent,
		const edm::EventSetup& iSetup) {

	// clear vectors
	tau_p4->clear();
	tau_p4Vis->clear();
	tau_p4Invis->clear();
	z_p4->clear();
	z_p4_byDecay->clear();
	tau_hadronic->clear();
	tau_muon->clear();
	tau_electron->clear();
	tau_nDaughters->clear();
	tau_prongs->clear();
	tau_piZeros->clear();
	tau_KProngs->clear();
	tau_KZeros->clear();
	tau_mother->clear();
	tau_grandma->clear();
	deltaR_TrueVis->clear();
	deltaR_VisInvis->clear();
	deltaR_InvisTrue->clear();

	// read generator particles
	edm::Handle<edm::View<reco::GenParticle> > particles;
	iEvent.getByLabel(genParticles_, particles);

	// loop over particles
	for (edm::View<reco::GenParticle>::const_iterator p = particles->begin(); p
			!= particles->end(); ++p) {
		if (p->pdgId() == 23 && p->status() == 3) {
			z_p4->push_back(p->p4());
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > zFour;
			for (size_t d = 0; d < p->numberOfDaughters(); d++) {
				const reco::Candidate* dau = p->daughter(d);
				if (std::abs(dau->pdgId()) == 15)
					zFour += dau->p4();
			}
			z_p4_byDecay->push_back(zFour);
		}
		if (abs(p->pdgId()) == 15 && p->status() == 2) {
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4Invis,
					p4Vis;
			bool hadronic = true, muon = false, electron = false;
			unsigned int prongs = 0, piZeros = 0, KProngs = 0, KZeros = 0;
			bool verbose = false;
			checkDaughters(&*p, muon, electron, prongs, piZeros, KProngs,
					KZeros, verbose);
			// additional loop for fourvectors
			for (size_t d = 0; d < p->numberOfDaughters(); d++) {
				const reco::Candidate* dau = p->daughter(d);
				switch (abs(dau->pdgId())) {
				case 12:
				case 14:
				case 16:
					p4Invis += dau->p4();
					break;
				default:
					p4Vis += dau->p4();
				}
			}
			if (muon || electron)
				hadronic = false;

			double dR_TrueVis, dR_VisInvis, dR_InvisTrue;
			dR_TrueVis = deltaR(p->p4(), p4Vis);
			dR_VisInvis = deltaR(p4Vis, p4Invis);
			dR_InvisTrue = deltaR(p4Invis, p->p4());

			tau_p4->push_back(p->p4());
			tau_p4Invis->push_back(p4Invis);
			tau_p4Vis->push_back(p4Vis);
			tau_nDaughters->push_back(p->numberOfDaughters());
			tau_hadronic->push_back(hadronic);
			tau_electron->push_back(electron);
			tau_muon->push_back(muon);
			tau_prongs->push_back(prongs);
			tau_piZeros->push_back(piZeros);
			tau_KProngs->push_back(KProngs);
			tau_KZeros->push_back(KZeros);
			tau_mother->push_back(p->mother()->pdgId());
			tau_grandma->push_back(p->mother()->mother()->pdgId());
			deltaR_TrueVis->push_back(dR_TrueVis);
			deltaR_VisInvis->push_back(dR_VisInvis);
			deltaR_InvisTrue->push_back(dR_InvisTrue);

			if (verbose) {
				std::cout << "hadronic: " << hadronic << ", muon: " << muon
						<< ", electron: " << electron << std::endl;
				std::cout << "prongs: " << prongs << ", piZeros: " << piZeros
						<< std::endl;
			}

//			if (p->mother()->mother()->pdgId() != 23) {
//				const reco::Candidate* obj = &*p;
//				std::cout << "Found non-Z Tau, print history:" << obj->pdgId()
//						<< " <-- ";
//				while (obj->mother()) {
//					obj = obj->mother();
//					std::cout << obj->pdgId() << " <-- ";
//				}
//				std::cout << "BANG" << std::endl;
//			}
		}
	}

	if (tau_grandma->at(0) == 23 && tau_grandma->at(1) == 23) {
		deltaR_taus = deltaR(tau_p4->at(0), tau_p4->at(1));
		double phi0 = tau_p4->at(0).Phi();
		double phi1 = tau_p4->at(1).Phi();
		deltaPhi_taus = deltaPhi(phi0, phi1);
		deltaEta_taus = tau_p4->at(0).Eta() - tau_p4->at(1).Eta();
	} else {
		deltaR_taus = 99.;
		deltaPhi_taus = 99.;
		deltaEta_taus = 99.;
	}

	tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void GenTauChecker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void GenTauChecker::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void GenTauChecker::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void GenTauChecker::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void GenTauChecker::beginLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void GenTauChecker::endLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenTauChecker::fillDescriptions(
		edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

// method to search daughters for interesting stuff
void GenTauChecker::checkDaughters(const reco::Candidate* object, bool& muon,
		bool& electron, unsigned int& prongs, unsigned int& piZeros,
		unsigned int& KProngs, unsigned int& KZeros, bool& verbose) {
	// loop over daughters
	for (size_t d = 0; d < object->numberOfDaughters(); d++) {
		const reco::Candidate* dau = object->daughter(d);
		if (verbose) {
			std::cout << "Daughter No " << d + 1 << ": " << dau->pdgId()
					<< " | " << dau->status() << std::endl;
		}
		switch (abs(dau->pdgId())) {
		case 11:
			if (dau->status() != 1)
				std::cout << "ERROR: electron with wrong status" << std::endl;
			electron = true;
			break;
		case 13:
			if (dau->status() != 1)
				std::cout << "ERROR: muon with wrong status" << std::endl;
			muon = true;
			break;
		case 12:
		case 14:
		case 16:
			if (dau->status() != 1)
				std::cout << "ERROR: neutrino with wrong status" << std::endl;
			break;
		case 211:
			if (dau->status() != 1)
				std::cout << "ERROR: prong with wrong status" << std::endl;
			prongs++;
			break;
		case 111:
			piZeros++;
			break;
		case 321:
			KProngs++;
			break;
		case 311:
			KZeros++;
			break;
		default:
			if (verbose)
				std::cout << "recursive call:" << std::endl;
			checkDaughters(dau, muon, electron, prongs, piZeros, KProngs,
					KZeros, verbose);
		}
	}
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenTauChecker)
;
