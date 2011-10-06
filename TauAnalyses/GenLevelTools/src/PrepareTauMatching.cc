// -*- C++ -*-
//
// Package:    PrepareTauMatching
// Class:      PrepareTauMatching
// 
/**\class PrepareTauMatching PrepareTauMatching.cc GenLevelTools/PrepareTauMatching/src/PrepareTauMatching.cc

 Description: prepare genParticle collection containing information about visible momentum and link to original genParticle

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Bastian Kargoll
//         Created:  Wed Jul 13 15:31:01 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Common/interface/Association.h"
//#include "DataFormats/Common/interface/AssociationMap.h"

#include "Math/Vector4D.h"

//
// class declaration
//

class PrepareTauMatching: public edm::EDProducer {
public:
	explicit PrepareTauMatching(const edm::ParameterSet&);
	~PrepareTauMatching();

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

	void visibleMomentum(const reco::Candidate*, ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> >&, bool);

	// ----------member data ---------------------------
	bool verbose_;

	edm::InputTag genParticles_;
};

//
// constants, enums and typedefs
//
typedef edm::Association<reco::GenParticleCollection> MCOriginalMap;
//typedef edm::AssociationMap<edm::OneToValue<reco::GenParticleCollection,
//		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >
//		GenTauP4Map;

//
// static data member definitions
//

//
// constructors and destructor
//
PrepareTauMatching::PrepareTauMatching(const edm::ParameterSet& iConfig) {

	verbose_ = iConfig.getParameter<bool> ("verbose");
	genParticles_ = iConfig.getParameter<edm::InputTag> ("genParticles");

	produces<reco::GenParticleCollection> ("genTausVisibleP4").setBranchAlias(
			"genTausVisibleP4");
	produces<MCOriginalMap> ("originalGenParticleMap").setBranchAlias(
			"originalGenParticleMap");
	//	produces<GenTauP4Map> ("genTauTrueP4Map").setBranchAlias("genTauTrueP4Map");
	//	produces<std::vector<unsigned int> > ("genParticleIndex").setBranchAlias(
	//			"genParticleIndex");

}

PrepareTauMatching::~PrepareTauMatching() {

}

//
// member functions
//

// ------------ method called to produce the data  ------------
void PrepareTauMatching::produce(edm::Event& iEvent,
		const edm::EventSetup& iSetup) {

	// read generator particles
	edm::Handle<reco::GenParticleCollection> particles;
	iEvent.getByLabel(genParticles_, particles);

	// create the vectors to put in the event stream
	std::auto_ptr<reco::GenParticleCollection> genTausVisibleP4(
			new reco::GenParticleCollection);
	std::auto_ptr<MCOriginalMap> originalGenParticleMap(new MCOriginalMap(
			particles));
	//	std::auto_ptr<GenTauP4Map> genTauTrueP4Map(new GenTauP4Map);
	//	std::auto_ptr<std::vector<unsigned int> > genParticleIndex(new std::vector<
	//			unsigned int>);

	// initialize the edm::Association
	MCOriginalMap::Filler filler(*originalGenParticleMap);
	std::vector<int> indices;

	// loop over particles
	//	unsigned int cntTaus = 0;
	for (unsigned i = 0; i != particles->size(); i++) {
		reco::GenParticleRef p(particles, i);
		// select taus
		if (abs(p->pdgId()) == 15 && p->status() == 2) {
			//			//safe true p4
			//			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4True;
			//			p4True = p->p4();
			if (verbose_)
				std::cout << p->numberOfDaughters() << "daughters, True pT = "
						<< p->pt() << std::endl;
			// calculate p4 of visible daughters
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4Vis;
			visibleMomentum(&(*p), p4Vis, verbose_);
			reco::GenParticle genTau = *p;
			genTau.setP4(p4Vis);
			if (verbose_)
				std::cout << ", Vis. pT = " << genTau.pt() << std::endl;
			// safe tau as modified genParticle
			genTausVisibleP4->push_back(genTau);
			// safe index of associated original genParticle
			//			genParticleIndex->push_back(i);
			indices.push_back(i);
			//			reco::GenParticleRef tau(genTausVisibleP4,cntTaus);
			//			genTauTrueP4Map->insert(tau,p4True);
			//			++cntTaus;
		}
	}

	edm::OrphanHandle<reco::GenParticleCollection> genTauHandle = iEvent.put(
			genTausVisibleP4, "genTausVisibleP4");
	//	iEvent.put(genTauTrueP4Map, "genTauTrueP4Map");
	//	iEvent.put(genParticleIndex, "genParticleIndex");

	// create the association
	filler.insert(genTauHandle, indices.begin(), indices.end());
	filler.fill();
	iEvent.put(originalGenParticleMap, "originalGenParticleMap");

	//	if (verbose_) {
	//		std::cout << "1" << std::endl;
	//		std::cout << "There are " << genTausVisibleP4->size() << "="
	//				<< genParticleIndex->size() << " taus in this event.";
	//		for (unsigned int i = 0; i < genParticleIndex->size(); ++i) {
	//			std::cout << "\nTau no. " << i << " can be found at position " << genParticleIndex->at(i);
	//			if (genTausVisibleP4->at(i).charge() != particles->at(genParticleIndex->at(i)).charge() ||
	//				genTausVisibleP4->at(i).status() != particles->at(genParticleIndex->at(i)).status()	){
	//				std::cout << "@@@@!!!!WARNING: NOT THE SAME PARTICLE!!!!@@@@";
	//			}
	//		}
	//		std::cout << std::endl;
	//	}

}

// ------------ method called once each job just before starting event loop  ------------
void PrepareTauMatching::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void PrepareTauMatching::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void PrepareTauMatching::beginRun(edm::Run&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void PrepareTauMatching::endRun(edm::Run&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void PrepareTauMatching::beginLuminosityBlock(edm::LuminosityBlock&,
		edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void PrepareTauMatching::endLuminosityBlock(edm::LuminosityBlock&,
		edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void PrepareTauMatching::fillDescriptions(
		edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

// calculate visible part of 4momentum of given particle
void PrepareTauMatching::visibleMomentum(const reco::Candidate* object,
		ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >& visP4,
		bool verbose = false) {
	// loop over daughters
	for (size_t d = 0; d < object->numberOfDaughters(); d++) {
		const reco::Candidate* dau = object->daughter(d);
		if (verbose) {
			std::cout << "Daughter No " << d + 1 << ": " << dau->pdgId()
					<< " | " << dau->status() << std::endl;
		}
		switch (dau->status()) {
		case 1:
			switch (abs(dau->pdgId())) {
			case 12:
			case 14:
			case 16:
				if (verbose)
					std::cout << "neutrino " << std::endl;
				break;
			default:
				if (verbose)
					std::cout << "visible particle " << std::endl;
				visP4 += dau->p4();
			}
			break;
		case 2:
		case 3:
			if (verbose)
				std::cout << "recursive call " << std::endl;
			visibleMomentum(dau, visP4, verbose);
			break;
		default:
			std::cout << "ERROR: genParticle with status != 1,2 or 3"
					<< std::endl;
		}

	}
}
//define this as a plug-in
DEFINE_FWK_MODULE(PrepareTauMatching)
;
