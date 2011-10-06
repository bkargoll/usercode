// -*- C++ -*-
//
// Package:    PrepMatchCollectionsChargedMatching
// Class:      PrepMatchCollectionsChargedMatching
// 
/**\class PrepMatchCollectionsChargedMatching PrepMatchCollectionsChargedMatching.cc GenLevelTools/PrepMatchCollectionsChargedMatching/src/PrepMatchCollectionsChargedMatching.cc

 Description: prepare genParticles and pfTaus for matching using charged constituents

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Bastian Kargoll
//         Created:  Tue Aug 30 20:02:31 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/TauReco/interface/PFTau.h"

#include "DataFormats/Common/interface/Association.h"

#include "Math/Vector4D.h"

//
// class declaration
//

class PrepMatchCollectionsChargedMatching: public edm::EDProducer {
public:
	explicit PrepMatchCollectionsChargedMatching(const edm::ParameterSet&);
	~PrepMatchCollectionsChargedMatching();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
	virtual void beginJob();
	virtual void produce(edm::Event&, const edm::EventSetup&);
	virtual void endJob();

	virtual void beginRun(edm::Run&, edm::EventSetup const&);
	virtual void endRun(edm::Run&, edm::EventSetup const&);
	virtual void beginLuminosityBlock(edm::LuminosityBlock&,
			edm::EventSetup const&);
	virtual void endLuminosityBlock(edm::LuminosityBlock&,
			edm::EventSetup const&);

	void genChargedMomentum(const reco::Candidate*, ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> >&, bool);
	void pfChargedMomentum(reco::PFTauRef, ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> >&, bool);

	// ----------member data ---------------------------
	bool verbose_;

	edm::InputTag genParticles_;
	edm::InputTag pfTaus_;
};

//
// constants, enums and typedefs
//
typedef edm::Association<reco::GenParticleCollection> MCOriginalMap;
typedef edm::Association<reco::PFTauCollection> pfTauOriginalMap;

//
// static data member definitions
//

//
// constructors and destructor
//
PrepMatchCollectionsChargedMatching::PrepMatchCollectionsChargedMatching(
		const edm::ParameterSet& iConfig) {
	verbose_ = iConfig.getParameter<bool> ("verbose");
	genParticles_ = iConfig.getParameter<edm::InputTag> ("genParticles");
	pfTaus_ = iConfig.getParameter<edm::InputTag> ("pfTaus");

	produces<reco::GenParticleCollection> ("genTausChargedP4").setBranchAlias(
			"genTausChargedP4");
	produces<MCOriginalMap> ("originalGenParticleMap").setBranchAlias(
			"originalGenParticleMap");
	produces<reco::PFTauCollection> ("hpsTausChargedP4").setBranchAlias(
			"hpsTausChargedP4");
	produces<pfTauOriginalMap> ("originalHPSTauMap").setBranchAlias(
			"originalHPSTauMap");
}

PrepMatchCollectionsChargedMatching::~PrepMatchCollectionsChargedMatching() {

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called to produce the data  ------------
void PrepMatchCollectionsChargedMatching::produce(edm::Event& iEvent,
		const edm::EventSetup& iSetup) {
	// read generator particles
	edm::Handle<reco::GenParticleCollection> genCollection;
	iEvent.getByLabel(genParticles_, genCollection);
	// read HPS taus
	edm::Handle<reco::PFTauCollection> hpsCollection;
	iEvent.getByLabel(pfTaus_, hpsCollection);

	// create the vectors to put in the event stream
	std::auto_ptr<reco::GenParticleCollection> genTausChargedP4(
			new reco::GenParticleCollection);
	std::auto_ptr<MCOriginalMap> originalGenParticleMap(new MCOriginalMap(
			genCollection));
	std::auto_ptr<reco::PFTauCollection> hpsTausChargedP4(
			new reco::PFTauCollection);
	std::auto_ptr<pfTauOriginalMap> originalHPSTauMap(new pfTauOriginalMap(
			hpsCollection));

	// initialize the edm::Association
	MCOriginalMap::Filler genFiller(*originalGenParticleMap);
	std::vector<int> genIndices;

	pfTauOriginalMap::Filler hpsFiller(*originalHPSTauMap);
	std::vector<int> hpsIndices;

	// loop over gen particles
	for (unsigned i = 0; i != genCollection->size(); i++) {
		reco::GenParticleRef p(genCollection, i);
		// select taus
		if (abs(p->pdgId()) == 15 && p->status() == 2) {
			// calculate p4 of charged daughters
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4Charged;
			genChargedMomentum(&(*p), p4Charged, verbose_);
			reco::GenParticle genTau = *p;
			genTau.setP4(p4Charged);
			// safe tau as modified genParticle
			genTausChargedP4->push_back(genTau);
			// safe index of associated original genParticle
			genIndices.push_back(i);
		}
	}
	// put the new genTaus in the event stream
	edm::OrphanHandle<reco::GenParticleCollection> genTauHandle = iEvent.put(
			genTausChargedP4, "genTausChargedP4");
	// create the association
	genFiller.insert(genTauHandle, genIndices.begin(), genIndices.end());
	genFiller.fill();
	iEvent.put(originalGenParticleMap, "originalGenParticleMap");

	// loop over hps taus
	for (unsigned i = 0; i != hpsCollection->size(); i++) {
		reco::PFTauRef p(hpsCollection, i);
		// select taus
		if (p->decayMode() != -1 && p->signalPFChargedHadrCands().size() > 0) {
			// calculate p4 of charged substituents
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4Charged;
			pfChargedMomentum(p, p4Charged, verbose_);
			reco::PFTau hpsTau = *p;
			hpsTau.setP4(p4Charged);
			hpsTausChargedP4->push_back(hpsTau);
			hpsIndices.push_back(i);
		}
	}
	// put the new hpsTaus in the event stream
	edm::OrphanHandle<reco::PFTauCollection> hpsTauHandle = iEvent.put(
			hpsTausChargedP4, "hpsTausChargedP4");
	// create the association
	hpsFiller.insert(hpsTauHandle, hpsIndices.begin(), hpsIndices.end());
	hpsFiller.fill();
	iEvent.put(originalHPSTauMap, "originalHPSTauMap");
}

// ------------ method called once each job just before starting event loop  ------------
void PrepMatchCollectionsChargedMatching::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void PrepMatchCollectionsChargedMatching::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void PrepMatchCollectionsChargedMatching::beginRun(edm::Run&,
		edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void PrepMatchCollectionsChargedMatching::endRun(edm::Run&,
		edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void PrepMatchCollectionsChargedMatching::beginLuminosityBlock(
		edm::LuminosityBlock&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void PrepMatchCollectionsChargedMatching::endLuminosityBlock(
		edm::LuminosityBlock&, edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PrepMatchCollectionsChargedMatching::fillDescriptions(
		edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

// calculate charged part of 4momentum of given particle
void PrepMatchCollectionsChargedMatching::genChargedMomentum(
		const reco::Candidate* object, ROOT::Math::LorentzVector<
				ROOT::Math::PxPyPzE4D<double> >& charP4, bool verbose = false) {
	// loop over daughters
	for (size_t d = 0; d < object->numberOfDaughters(); d++) {
		const reco::Candidate* dau = object->daughter(d);
		if (verbose) {
			std::cout << "Daughter No " << d + 1 << ": " << dau->pdgId()
					<< " | " << dau->status() << std::endl;
		}
		switch (dau->status()) {
		case 1:
			if (dau->charge() != 0) {
				charP4 += dau->p4();
				if (verbose)
					std::cout << "charged particle, pdgID " << dau->pdgId()
							<< std::endl;
			} else if (verbose) {
				std::cout << "neutral particle, pdgID " << dau->pdgId()
						<< std::endl;
			}
			break;
		case 2:
		case 3:
			if (verbose)
				std::cout << "recursive call " << std::endl;
			genChargedMomentum(dau, charP4, verbose);
			break;
		default:
			std::cout << "ERROR: genParticle with status != 1,2 or 3"
					<< std::endl;
		}
	}
}

// calculate p4 sum of charged tracks
void PrepMatchCollectionsChargedMatching::pfChargedMomentum(reco::PFTauRef tau,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& charP4,
		bool verbose = false) {
	if (verbose)
		std::cout << "Adding up " << tau->signalPFChargedHadrCands().size()
				<< " constituents." << std::endl;
	for (size_t d = 0; d < tau->signalPFChargedHadrCands().size(); d++) {
		charP4 += tau->signalPFChargedHadrCands().at(d)->p4();
	}
}
//define this as a plug-in
DEFINE_FWK_MODULE(PrepMatchCollectionsChargedMatching)
;
