// -*- C++ -*-
//
// Package:    MatchingEfficiency
// Class:      MatchingEfficiency
// 
/**\class MatchingEfficiency MatchingEfficiency.cc TestAnalyses/MatchingEfficiency/src/MatchingEfficiency.cc

 Description: determine matching efficiencies

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Bastian Kargoll
//         Created:  Thu Jul  7 12:09:21 CEST 2011
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

#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TauAnalyses/GenLevelTools/src/GenTauTools.cc"

#include "Math/Vector4D.h"
#include "TTree.h"
#include "TMath.h"
//#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/Math/interface/deltaPhi.h"
//
// class declaration
//

class MatchingEfficiency: public edm::EDAnalyzer {
public:
	explicit MatchingEfficiency(const edm::ParameterSet&);
	~MatchingEfficiency();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

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

	bool addIfNew(const reco::GenParticle*, std::vector<
			const reco::GenParticle*> &);

	//	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > chargedP4(reco::PFTauRef);

	// ----------member data ---------------------------
	bool verbose_;
	double minGenPt_;
	std::vector<int> motherId_;

	std::vector<unsigned int>* cntGenTaus_;
	std::vector<unsigned int>* cntMatchedGenTaus_;
	unsigned int cntMatchedRecTaus_;

	edm::InputTag genParticles_;
	edm::InputTag genTausForMatching_;
	edm::InputTag genTausMapToGenParticles_;
	edm::InputTag tauCandidates_;
	edm::InputTag tauCandidatesForMatching_;
	edm::InputTag MatchingTausMapToTauCandidates_;
	edm::InputTag matching_;

	TTree* tree;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* genTau_p4;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* genTau_p4Matching;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* genTau_p4Vis;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_p4;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_p4Matching;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_genp4;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_genp4Matching;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_genp4Vis;
	std::vector<unsigned int>* genTau_decayMode;
	std::vector<unsigned int>* recoTau_genDecayMode;
	std::vector<int>* recoTau_HPSdecayMode;
	std::vector<unsigned int>* recoTau_nChargedHadrCand;
	std::vector<bool>* genTau_matched;
	std::vector<bool>* recoTau_matched;

//	std::vector<float>* recoTau_discr_DecayModeFinding;
	std::vector<std::string> tauDiscr_names;
	std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscr_handles;
	std::vector<std::vector<bool>* > tauDiscr_discrVec;
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
MatchingEfficiency::MatchingEfficiency(const edm::ParameterSet& iConfig)

{
	verbose_ = iConfig.getParameter<bool> ("verbose");
	genParticles_ = iConfig.getParameter<edm::InputTag> ("genParticles");
	genTausForMatching_ = iConfig.getParameter<edm::InputTag> (
			"genTausForMatching");
	genTausMapToGenParticles_ = iConfig.getParameter<edm::InputTag> (
			"genTausMapToGenParticles");
	tauCandidates_ = iConfig.getParameter<edm::InputTag> ("tauCandidates");
	tauCandidatesForMatching_ = iConfig.getParameter<edm::InputTag> (
			"tauCandidatesForMatching");
	MatchingTausMapToTauCandidates_ = iConfig.getParameter<edm::InputTag> (
			"MatchingTausMapToTauCandidates");
	matching_ = iConfig.getParameter<edm::InputTag> ("matching");
	minGenPt_ = iConfig.getParameter<double> ("minGenPt");
	motherId_ = iConfig.getParameter<std::vector<int> > ("motherId");

	// use the TFileService
	edm::Service<TFileService> fs;
	// create the Tree
	tree = fs->make<TTree> ("Tree", "MatchingEfficiency");

	cntGenTaus_ = new std::vector<unsigned int>(35);
	cntMatchedGenTaus_ = new std::vector<unsigned int>(35);

	genTau_p4 = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	genTau_p4Matching = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	genTau_p4Vis = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	recoTau_p4 = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	recoTau_p4Matching = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	recoTau_genp4 = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	recoTau_genp4Matching = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	recoTau_genp4Vis = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	genTau_decayMode = new std::vector<unsigned int>;
	recoTau_genDecayMode = new std::vector<unsigned int>;
	recoTau_nChargedHadrCand = new std::vector<unsigned int>;
	;
	recoTau_HPSdecayMode = new std::vector<int>;
	genTau_matched = new std::vector<bool>;
	recoTau_matched = new std::vector<bool>;
//	recoTau_discr_DecayModeFinding = new std::vector<float>;
	for (unsigned int iDis = 0; iDis < 10; ++iDis) {
		tauDiscr_discrVec.push_back(new std::vector<bool>);
	}

	// define what discriminators to use
	tauDiscr_names.push_back("hpsPFTauDiscriminationByDecayModeFinding");
	tauDiscr_names.push_back("hpsPFTauDiscriminationByLooseElectronRejection");
	tauDiscr_names.push_back("hpsPFTauDiscriminationByLooseIsolation");
	tauDiscr_names.push_back("hpsPFTauDiscriminationByLooseMuonRejection");
	tauDiscr_names.push_back("hpsPFTauDiscriminationByMediumElectronRejection");
	tauDiscr_names.push_back("hpsPFTauDiscriminationByMediumIsolation");
	tauDiscr_names.push_back("hpsPFTauDiscriminationByTightElectronRejection");
	tauDiscr_names.push_back("hpsPFTauDiscriminationByTightIsolation");
	tauDiscr_names.push_back("hpsPFTauDiscriminationByTightMuonRejection");
	tauDiscr_names.push_back("hpsPFTauDiscriminationByVLooseIsolation");

	tauDiscr_handles.resize(tauDiscr_names.size());

	// register the branches
	tree->Branch(
			"genTau_p4",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&genTau_p4);
	tree->Branch(
			"genTau_p4Matching",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&genTau_p4Matching);
	tree->Branch(
			"genTau_p4Vis",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&genTau_p4Vis);
	tree->Branch(
			"recoTau_p4",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&recoTau_p4);
	tree->Branch(
			"recoTau_p4Matching",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&recoTau_p4Matching);
	tree->Branch(
			"recoTau_genp4",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&recoTau_genp4);
	tree->Branch(
			"recoTau_genp4Matching",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&recoTau_genp4Matching);
	tree->Branch(
			"recoTau_genp4Vis",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&recoTau_genp4Vis);
	tree->Branch("genTau_decayMode", "std::vector<unsigned int>",
			&genTau_decayMode);
	tree->Branch("recoTau_nChargedHadrCand", "std::vector<unsigned int>",
			&recoTau_nChargedHadrCand);
	tree->Branch("recoTau_genDecayMode", "std::vector<unsigned int>",
			&recoTau_genDecayMode);
	tree->Branch("recoTau_HPSdecayMode", "std::vector<int>",
			&recoTau_HPSdecayMode);
	tree->Branch("genTau_matched", "std::vector<bool>", &genTau_matched);
	tree->Branch("recoTau_matched", "std::vector<bool>", &recoTau_matched);
//	tree->Branch("recoTau_discr_DecayModeFinding", "std::vector<float>",
//			&recoTau_discr_DecayModeFinding);
	for (unsigned int iDis = 0; iDis < 10; ++iDis) {
		tree->Branch((tauDiscr_names.at(iDis)).c_str(), "std::vector<bool>",	&(tauDiscr_discrVec.at(iDis)));
	}

}

MatchingEfficiency::~MatchingEfficiency() {

	delete cntGenTaus_;
	delete cntMatchedGenTaus_;
}

//
// member functions
//

// ------------ method called for each event  ------------
void MatchingEfficiency::analyze(const edm::Event& iEvent,
		const edm::EventSetup& iSetup) {

	// clear the vectors
	genTau_p4->clear();
	genTau_p4Matching->clear();
	genTau_p4Vis->clear();
	recoTau_p4->clear();
	recoTau_p4Matching->clear();
	recoTau_genp4->clear();
	recoTau_genp4Matching->clear();
	recoTau_genp4Vis->clear();
	genTau_decayMode->clear();
	recoTau_genDecayMode->clear();
	recoTau_HPSdecayMode->clear();
	recoTau_nChargedHadrCand->clear();
	genTau_matched->clear();
	recoTau_matched->clear();
//	recoTau_discr_DecayModeFinding->clear();
	for (unsigned int i = 0; i < tauDiscr_discrVec.size(); i++) {
		tauDiscr_discrVec.at(i)->clear();
	}

	// read genParticles
	edm::Handle<reco::GenParticleCollection> genCollection;
	iEvent.getByLabel(genParticles_, genCollection);
	// read genTaus
	edm::Handle<reco::GenParticleCollection> genTauCollection;
	iEvent.getByLabel(genTausForMatching_, genTauCollection);
	// read association from genTaus to original genParticles
	edm::Handle<MCOriginalMap> genTausMapToGenParticles;
	iEvent.getByLabel(genTausMapToGenParticles_, genTausMapToGenParticles);
	// read tau collection
	edm::Handle<reco::PFTauCollection> tauCollection;
	iEvent.getByLabel(tauCandidates_, tauCollection);
	// read tau collection from matching
	edm::Handle<reco::PFTauCollection> tauMatchingCollection;
	iEvent.getByLabel(tauCandidatesForMatching_, tauMatchingCollection);
	// read association from matchingTaus to original HPSTaus
	edm::Handle<pfTauOriginalMap> MatchingTausMapToTauCandidates;
	iEvent.getByLabel(MatchingTausMapToTauCandidates_,
			MatchingTausMapToTauCandidates);
	// read gen-reco matching
	edm::Handle<reco::GenParticleMatch> matching;
	iEvent.getByLabel(matching_, matching);
	// read discriminator
//	edm::Handle<reco::PFTauDiscriminator> decayModeFinding;
//	iEvent.getByLabel("hpsPFTauDiscriminationByDecayModeFinding",
//			decayModeFinding);
	// read PFTauDiscriminators
	if(tauDiscr_names.size() != 10 || tauDiscr_names.size() != tauDiscr_handles.size() )
		std::cout << "WARNING WRONG SIZES" << std::endl;
	for (unsigned int iDis = 0; iDis < tauDiscr_names.size(); ++iDis) {
		edm::Handle<reco::PFTauDiscriminator> tmpHandle;
		iEvent.getByLabel(tauDiscr_names.at(iDis), tmpHandle);
		tauDiscr_handles.at(iDis) = tmpHandle;
	}

	std::vector<const reco::GenParticle*> genTaus;
	std::vector<const reco::GenParticle*> matchedGenTaus;

	GenTauTools genHelper;

	// count generator level taus
	for (unsigned int iTau = 0; iTau < genTauCollection->size(); iTau++) {
		reco::GenParticleRef genTauRef(genTauCollection, iTau);
		if (genTauRef->pt() > minGenPt_ && abs(genTauRef->eta()) < 2.5) {
			genTau_p4Matching->push_back(genTauRef->p4());
			reco::GenParticleRef trueP4Ref =
					(*genTausMapToGenParticles)[genTauRef];
			if (trueP4Ref.isNull() || genTauRef->charge()
					!= trueP4Ref->charge() || genTauRef->mother()->p4()
					!= trueP4Ref->mother()->p4()) {
				std::cout << "ERROR: genTau association is broken" << std::endl;
				ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
						nullVector;
				genTau_p4->push_back(nullVector);
				genTau_p4Vis->push_back(nullVector);
			} else {
				genTau_p4->push_back(trueP4Ref->p4());
				ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > visP4;
				genHelper.visibleP4(&(*trueP4Ref), visP4, verbose_);
				genTau_p4Vis->push_back(visP4);
			}
			unsigned int mode = genHelper.genDecayMode(genTauRef.get());
			genTau_decayMode->push_back(mode);
			(cntGenTaus_->at(mode))++;
			genTaus.push_back(genTauRef.get());
			genTau_matched->push_back(false); // decision is made later
		}
	}

	// count reconstructed
	int matchedRecTaus = 0;
	for (unsigned iTau = 0; iTau != tauMatchingCollection->size(); ++iTau) {
		reco::PFTauRef tauMatchingCandidate(tauMatchingCollection, iTau);
		recoTau_p4Matching->push_back(tauMatchingCandidate->p4());
		recoTau_HPSdecayMode->push_back(tauMatchingCandidate->decayMode());
		recoTau_nChargedHadrCand->push_back(
				tauMatchingCandidate->signalPFChargedHadrCands().size());
		// access original HPS tau
		reco::PFTauRef trueP4Ref =
				(*MatchingTausMapToTauCandidates)[tauMatchingCandidate];
		if (trueP4Ref.isNull() || tauMatchingCandidate->decayMode()
							!= trueP4Ref->decayMode()) {
						std::cout << "ERROR: hpsTau association is broken" << std::endl;
						ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
								nullVector;
						recoTau_p4->push_back(nullVector);
//						recoTau_discr_DecayModeFinding->push_back(false);
						for (unsigned iDis = 0; iDis != tauDiscr_names.size(); ++iDis) {
							tauDiscr_discrVec.at(iDis)->push_back(false);
						}
		}
		else {
			recoTau_p4->push_back(trueP4Ref->p4());
//			recoTau_discr_DecayModeFinding->push_back(
//					(*decayModeFinding)[trueP4Ref]);
			//check discriminators
			for (unsigned iDis = 0; iDis != tauDiscr_names.size(); ++iDis) {
				edm::Handle<reco::PFTauDiscriminator> tmpHandle = tauDiscr_handles.at(iDis);
				bool discrValue = (*tmpHandle)[trueP4Ref];
				tauDiscr_discrVec.at(iDis)->push_back(discrValue);
			}
		}
		// check gen-reco-matching
		reco::GenParticleRef mcMatch = (*matching)[tauMatchingCandidate];
		if (mcMatch.isNonnull() && mcMatch->pt() > minGenPt_ && abs(
				mcMatch->eta()) < 2.5) {
			cntMatchedRecTaus_++;
			matchedRecTaus++;
			recoTau_matched->push_back(true);
			int refIndex = -1;
			for (unsigned i = 0; i < genTaus.size(); i++) {
				if (mcMatch->p4() == genTaus.at(i)->p4() && mcMatch->status()
						== genTaus.at(i)->status()) {
					bool isNew = addIfNew(genTaus.at(i), matchedGenTaus);
					if (isNew)
						refIndex = i;
					else
						std::cout << "WARNING: Duplicate matching" << std::endl;
				}
			}
			recoTau_genDecayMode->push_back(genTau_decayMode->at(refIndex));
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4 =
					genTau_p4->at(refIndex);
			recoTau_genp4->push_back(p4);
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4Matching =
					genTau_p4Matching->at(refIndex);
			recoTau_genp4Matching->push_back(p4Matching);
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > p4Vis =
					genTau_p4Vis->at(refIndex);
			recoTau_genp4Vis->push_back(p4Vis);
			genTau_matched->at(refIndex) = true;
		} else {
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
					nullVector;
			recoTau_genDecayMode->push_back(-1);
			recoTau_genp4->push_back(nullVector);
			recoTau_genp4Matching->push_back(nullVector);
			recoTau_genp4Vis->push_back(nullVector);
			recoTau_matched->push_back(false);
		}
	}
	for (unsigned int i = 0; i < matchedGenTaus.size(); i++) {
		unsigned int mode = genHelper.genDecayMode(matchedGenTaus.at(i));
		(cntMatchedGenTaus_->at(mode))++;
	}

	if (verbose_) {
		std::cout << "gen Taus:          " << genTaus.size()
				<< "\nmatched gen Taus:  " << matchedGenTaus.size() << " ("
				<< float(matchedGenTaus.size()) / genTaus.size() << ")"
				<< "\nmatched reco Taus: " << matchedRecTaus << " ("
				<< float(matchedRecTaus) / genTaus.size() << ")" << std::endl;
	}

	tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void MatchingEfficiency::beginJob() {
	cntMatchedRecTaus_ = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void MatchingEfficiency::endJob() {
	unsigned int gen[3], mGen[3];
	// other decays
	gen[2] = cntGenTaus_->at(19);
	mGen[2] = cntMatchedGenTaus_->at(19);
	// leptonic decays
	gen[1] = cntGenTaus_->at(16) + cntGenTaus_->at(17);
	mGen[1] = cntMatchedGenTaus_->at(16) + cntMatchedGenTaus_->at(17);
	// hadronic decays
	gen[0] = 0;
	mGen[0] = 0;
	for (unsigned int i = 0; i < cntGenTaus_->size(); i++) {
		if (i != 16 && i != 17 && i != 19) {
			gen[0] += cntGenTaus_->at(i);
			mGen[0] += cntMatchedGenTaus_->at(i);
		}
	}

	std::cout << "hadronic, leptonic, other" << std::endl;

	std::cout << "gen Taus:          " << gen[0] << ", " << gen[1] << ", "
			<< gen[2] << "\nmatched gen Taus:  " << mGen[0] << ", " << mGen[1]
			<< ", " << mGen[2] << "\nmatched reco Taus: " << cntMatchedRecTaus_
			<< std::endl;
}

// ------------ method called when starting to processes a run  ------------
void MatchingEfficiency::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void MatchingEfficiency::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void MatchingEfficiency::beginLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void MatchingEfficiency::endLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MatchingEfficiency::fillDescriptions(
		edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

// add element to vector if  it is not in there already
bool MatchingEfficiency::addIfNew(const reco::GenParticle* particle,
		std::vector<const reco::GenParticle*>& vector) {
	double pPt = particle->pt();
	int pStatus = particle->status();
	for (unsigned i = 0; i < vector.size(); i++) {
		if (vector.at(i)->pt() == pPt && vector.at(i)->status() == pStatus)
			return false;
	}
	vector.push_back(particle);
	return true;
}
//// calculate charged part of 4momentum of given recoTau
//ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > MatchingEfficiency::chargedP4(reco::PFTauRef tau) {
//	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > result;
//	// loop over charged daughters
//	const reco::PFCandidateRefVector chargedCands = tau->signalPFChargedHadrCands();
//	edm::RefVector<reco::PFCandidateCollection>::const_iterator iCand;
//	for (iCand = chargedCands.begin(); iCand != chargedCands.end(); iCand++) {
//		const reco::PFCandidate& pfCand = *(iCand->get());
//		result += pfCand.p4();
//	}
//	return result;
//}
//define this as a plug-in
DEFINE_FWK_MODULE(MatchingEfficiency)
;
