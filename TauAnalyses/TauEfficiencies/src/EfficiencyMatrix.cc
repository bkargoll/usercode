// -*- C++ -*-
//
// Package:    EfficiencyMatrix
// Class:      EfficiencyMatrix
// 
/**\class EfficiencyMatrix EfficiencyMatrix.cc TestAnalyses/EfficiencyMatrix/src/EfficiencyMatrix.cc

 Description: calculate efficiencies of individual decay channels

 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Bastian Kargoll
//         Created:  Mon Jun  6 11:37:26 CEST 2011
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

#include "Math/Vector4D.h"
#include "TTree.h"
#include "TString.h"

#include "TauAnalyses/GenLevelTools/src/GenTauTools.cc"

//
// class declaration
//

class EfficiencyMatrix: public edm::EDAnalyzer {
public:
	explicit EfficiencyMatrix(const edm::ParameterSet&);
	~EfficiencyMatrix();

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

	bool hasGoodMother(const reco::GenParticle*, int);

	// ----------member data ---------------------------
	bool verbose_;
	std::vector<int> motherId_;

	edm::InputTag genParticles_;
	edm::InputTag trueTaus_;
	edm::InputTag matchingTaus_;
	edm::InputTag MatchingTausMapToTauCandidates_;
	edm::InputTag matching_;

	TTree* tree;

	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* recoTau_p4;
	std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >
			* genTau_p4Vis;
	std::vector<int>* genTau_mode;
	std::vector<bool>* genTau_goodMother;
	std::vector<int>* recoTau_mode;
	std::vector<int>* recoTau_matchedGenMode;
	std::vector<bool>* recoTau_matchedGenGoodMother;

//	std::vector<float>* recoTau_discr_DecayModeFinding;
//	std::vector<float>* recoTau_discr_LooseElectronRejection;
//	std::vector<float>* recoTau_discr_LooseIsolation;
//	std::vector<float>* recoTau_discr_LooseMuonRejection;
//	std::vector<float>* recoTau_discr_MediumElectronRejection;
//	std::vector<float>* recoTau_discr_MediumIsolation;
//	std::vector<float>* recoTau_discr_TightElectronRejection;
//	std::vector<float>* recoTau_discr_TightIsolation;
//	std::vector<float>* recoTau_discr_TightMuonRejection;
//	std::vector<float>* recoTau_discr_VLooseIsolation;

	std::vector<std::string> tauDiscr_names;
	std::vector<edm::Handle<reco::PFTauDiscriminator> > tauDiscr_handles;
	std::vector<std::vector<bool>* > tauDiscr_discrVec;
//	std::vector<float>& tauDiscr_discrVec[] = {
//		*recoTau_discr_DecayModeFinding,
//		*recoTau_discr_LooseElectronRejection,
//		*recoTau_discr_LooseIsolation,
//		*recoTau_discr_LooseMuonRejection,
//		*recoTau_discr_MediumElectronRejection,
//		*recoTau_discr_MediumIsolation,
//		*recoTau_discr_TightElectronRejection,
//		*recoTau_discr_TightIsolation,
//		*recoTau_discr_TightMuonRejection,
//		*recoTau_discr_VLooseIsolation
//	};
};

//
// constants, enums and typedefs
//
typedef edm::Association<reco::PFTauCollection> pfTauOriginalMap;
//
// static data member definitions
//

//
// constructors and destructor
//
EfficiencyMatrix::EfficiencyMatrix(const edm::ParameterSet& iConfig)

{
	verbose_ = iConfig.getParameter<bool> ("verbose");
	genParticles_ = iConfig.getParameter<edm::InputTag> ("genParticles");
	trueTaus_ = iConfig.getParameter<edm::InputTag> ("originalPFTaus");
	matchingTaus_ = iConfig.getParameter<edm::InputTag> ("PFTausForMatching");
	MatchingTausMapToTauCandidates_ = iConfig.getParameter<edm::InputTag> ("PFTauMatchingOriginalMap");
	matching_ = iConfig.getParameter<edm::InputTag> ("matching");
	motherId_ = iConfig.getParameter<std::vector<int> > ("motherId");

	// use the TFileService
	edm::Service<TFileService> fs;
	// create the Tree
	tree = fs->make<TTree> ("Tree", "EfficiencyMatrix");

	// create the needed vectors
	recoTau_p4 = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	genTau_p4Vis = new std::vector<ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> > >;
	genTau_mode = new std::vector<int>;
	genTau_goodMother = new std::vector<bool>;
	recoTau_mode = new std::vector<int>;
	recoTau_matchedGenMode = new std::vector<int>;
	recoTau_matchedGenGoodMother = new std::vector<bool>;
	for (unsigned int iDis = 0; iDis < 10; ++iDis) {
		tauDiscr_discrVec.push_back(new std::vector<bool>);
	}
//	recoTau_discr_DecayModeFinding = new std::vector<float>;
//	recoTau_discr_LooseElectronRejection = new std::vector<float>;
//	recoTau_discr_LooseIsolation = new std::vector<float>;
//	recoTau_discr_LooseMuonRejection = new std::vector<float>;
//	recoTau_discr_MediumElectronRejection = new std::vector<float>;
//	recoTau_discr_MediumIsolation = new std::vector<float>;
//	recoTau_discr_TightElectronRejection = new std::vector<float>;
//	recoTau_discr_TightIsolation = new std::vector<float>;
//	recoTau_discr_TightMuonRejection = new std::vector<float>;
//	recoTau_discr_VLooseIsolation = new std::vector<float>;

	// define what discriminators to use
	tauDiscr_names.push_back("hpsPFTauDiscriminationByDecayModeFinding");
//	tauDiscr_discrVec.push_back(recoTau_discr_DecayModeFinding);
	tauDiscr_names.push_back("hpsPFTauDiscriminationByLooseElectronRejection");
//	tauDiscr_discrVec.push_back(recoTau_discr_LooseElectronRejection);
	tauDiscr_names.push_back("hpsPFTauDiscriminationByLooseIsolation");
//	tauDiscr_discrVec.push_back(recoTau_discr_LooseIsolation);
	tauDiscr_names.push_back("hpsPFTauDiscriminationByLooseMuonRejection");
//	tauDiscr_discrVec.push_back(recoTau_discr_LooseMuonRejection);
	tauDiscr_names.push_back("hpsPFTauDiscriminationByMediumElectronRejection");
//	tauDiscr_discrVec.push_back(recoTau_discr_MediumElectronRejection);
	tauDiscr_names.push_back("hpsPFTauDiscriminationByMediumIsolation");
//	tauDiscr_discrVec.push_back(recoTau_discr_MediumIsolation);
	tauDiscr_names.push_back("hpsPFTauDiscriminationByTightElectronRejection");
//	tauDiscr_discrVec.push_back(recoTau_discr_TightElectronRejection);
	tauDiscr_names.push_back("hpsPFTauDiscriminationByTightIsolation");
//	tauDiscr_discrVec.push_back(recoTau_discr_TightIsolation);
	tauDiscr_names.push_back("hpsPFTauDiscriminationByTightMuonRejection");
//	tauDiscr_discrVec.push_back(recoTau_discr_TightMuonRejection);
	tauDiscr_names.push_back("hpsPFTauDiscriminationByVLooseIsolation");
//	tauDiscr_discrVec.push_back(recoTau_discr_VLooseIsolation);

//	tauDiscr_discrVec[0] = *recoTau_discr_DecayModeFinding;
//	tauDiscr_discrVec[1] = *recoTau_discr_LooseElectronRejection;
//	tauDiscr_discrVec[2] = *recoTau_discr_LooseIsolation;
//	tauDiscr_discrVec[3] = *recoTau_discr_LooseMuonRejection;
//	tauDiscr_discrVec[4] = *recoTau_discr_MediumElectronRejection;
//	tauDiscr_discrVec[5] = *recoTau_discr_MediumIsolation;
//	tauDiscr_discrVec[6] = *recoTau_discr_TightElectronRejection;
//	tauDiscr_discrVec[7] = *recoTau_discr_TightIsolation;
//	tauDiscr_discrVec[8] = *recoTau_discr_TightMuonRejection;
//	tauDiscr_discrVec[9] = *recoTau_discr_VLooseIsolation;

	tauDiscr_handles.resize(tauDiscr_names.size());

	//register the branches
	tree->Branch(
			"recoTau_p4",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&recoTau_p4);
	tree->Branch(
			"genTau_p4Vis",
			"std::vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >",
			&genTau_p4Vis);
	tree->Branch("genTau_mode", "std::vector<int>", &genTau_mode);
	tree->Branch("genTau_goodMother", "std::vector<bool>", &genTau_goodMother);
	tree->Branch("recoTau_mode", "std::vector<int>", &recoTau_mode);
	tree->Branch("recoTau_matchedGenMode", "std::vector<int>",
			&recoTau_matchedGenMode);
	tree->Branch("recoTau_matchedGenGoodMother", "std::vector<bool>",
			&recoTau_matchedGenGoodMother);
	for (unsigned int iDis = 0; iDis < 10; ++iDis) {
		tree->Branch((tauDiscr_names.at(iDis)).c_str(), "std::vector<bool>",	&(tauDiscr_discrVec.at(iDis)));
	}

//	tree->Branch("recoTau_discr_DecayModeFinding", "std::vector<float>",
//			&recoTau_discr_DecayModeFinding);
//	tree->Branch("recoTau_discr_LooseElectronRejection", "std::vector<float>",
//			&recoTau_discr_LooseElectronRejection);
//	tree->Branch("recoTau_discr_LooseIsolation", "std::vector<float>",
//			&recoTau_discr_LooseIsolation);
//	tree->Branch("recoTau_discr_LooseMuonRejection", "std::vector<float>",
//			&recoTau_discr_LooseMuonRejection);
//	tree->Branch("recoTau_discr_MediumElectronRejection", "std::vector<float>",
//			&recoTau_discr_MediumElectronRejection);
//	tree->Branch("recoTau_discr_MediumIsolation", "std::vector<float>",
//			&recoTau_discr_MediumIsolation);
//	tree->Branch("recoTau_discr_TightElectronRejection", "std::vector<float>",
//			&recoTau_discr_TightElectronRejection);
//	tree->Branch("recoTau_discr_TightIsolation", "std::vector<float>",
//			&recoTau_discr_TightIsolation);
//	tree->Branch("recoTau_discr_TightMuonRejection", "std::vector<float>",
//			&recoTau_discr_TightMuonRejection);
//	tree->Branch("recoTau_discr_VLooseIsolation", "std::vector<float>",
//			&recoTau_discr_VLooseIsolation);
}

EfficiencyMatrix::~EfficiencyMatrix() {
}

//
// member functions
//

// ------------ method called for each event  ------------
void EfficiencyMatrix::analyze(const edm::Event& iEvent,
		const edm::EventSetup& iSetup) {

	// clear vectors
	recoTau_p4->clear();
	genTau_p4Vis->clear();
	genTau_mode->clear();
	genTau_goodMother->clear();
	recoTau_mode->clear();
	recoTau_matchedGenMode->clear();
	recoTau_matchedGenGoodMother->clear();

	//	recoTau_discr_DecayModeFinding->clear();
	//	recoTau_discr_LooseElectronRejection->clear();
	//	recoTau_discr_LooseIsolation->clear();
	//	recoTau_discr_LooseMuonRejection->clear();
	//	recoTau_discr_MediumElectronRejection->clear();
	//	recoTau_discr_MediumIsolation->clear();
	//	recoTau_discr_TightElectronRejection->clear();
	//	recoTau_discr_TightIsolation->clear();
	//	recoTau_discr_TightMuonRejection->clear();
	//	recoTau_discr_VLooseIsolation->clear();
	for (unsigned int i = 0; i < tauDiscr_discrVec.size(); i++) {
//		tauDiscr_discrVec[i].clear();
		tauDiscr_discrVec.at(i)->clear();
	}

	// read genParticles
	edm::Handle<reco::GenParticleCollection> genCollection;
	iEvent.getByLabel(genParticles_, genCollection);
	// read tau collection
	edm::Handle<reco::PFTauCollection> tauCollection;
	iEvent.getByLabel(matchingTaus_, tauCollection);
	// read original tau collection
	edm::Handle<reco::PFTauCollection> tauCollectionOriginalP4;
	iEvent.getByLabel(trueTaus_, tauCollectionOriginalP4);
	// read association from matchingTaus to original HPSTaus
	edm::Handle<pfTauOriginalMap> MatchingTausMapToTauCandidates;
	iEvent.getByLabel(MatchingTausMapToTauCandidates_,
			MatchingTausMapToTauCandidates);


	// read PFTauDiscriminators
	if(tauDiscr_names.size() != 10 || tauDiscr_names.size() != tauDiscr_handles.size() )
		std::cout << "WARNING WRONG SIZES" << std::endl;
	for (unsigned int iDis = 0; iDis < tauDiscr_names.size(); ++iDis) {
		edm::Handle<reco::PFTauDiscriminator> tmpHandle;
		iEvent.getByLabel(tauDiscr_names.at(iDis), tmpHandle);
		tauDiscr_handles.at(iDis) = tmpHandle;
	}

	// read matching
	edm::Handle<reco::GenParticleMatch> matching;
	iEvent.getByLabel(matching_, matching);

	GenTauTools genTool;

	// count generator level decays
	for (reco::GenParticleCollection::const_iterator itGen =
			genCollection->begin(); itGen != genCollection->end(); ++itGen) {
		genTau_p4Vis->push_back(itGen->p4());
		int mode = genTool.genDecayMode(&(*itGen));
		if (mode == 19)
			mode = -1;
		genTau_mode->push_back(mode);
		genTau_goodMother->push_back(hasGoodMother(&(*itGen), 2)); // status 2 taus have depth 2: Z->tau(3)->tau(2)
	}

	// determine reconstructed decay modes
	for (unsigned iTau = 0; iTau != tauCollection->size(); ++iTau) {
		reco::PFTauRef tauCandidate(tauCollection, iTau);
		recoTau_p4->push_back(tauCandidate->p4());
		recoTau_mode->push_back(tauCandidate->decayMode());
		reco::GenParticleRef mcMatch = (*matching)[tauCandidate];
		if (mcMatch.isNonnull()) {
			int mode = genTool.genDecayMode(&*mcMatch);
			if (mode == 19)
				mode = -1;
			recoTau_matchedGenMode->push_back(mode);
			recoTau_matchedGenGoodMother->push_back(hasGoodMother(
					mcMatch.get(), 2));
		} else {
			recoTau_matchedGenMode->push_back(-2);
			recoTau_matchedGenGoodMother->push_back(false);
		}
		// access original HPS tau
		reco::PFTauRef trueP4Ref =
				(*MatchingTausMapToTauCandidates)[tauCandidate];
		//check discriminators
		for (unsigned iDis = 0; iDis != tauDiscr_names.size(); ++iDis) {
			edm::Handle<reco::PFTauDiscriminator> tmpHandle = tauDiscr_handles.at(iDis);
			bool discrValue = (*tmpHandle)[trueP4Ref];
//			tauDiscr_discrVec[iDis].push_back(discrValue);
			tauDiscr_discrVec.at(iDis)->push_back(discrValue);
		}
	}

	tree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void EfficiencyMatrix::beginJob() {
	//std::cout << "GO!" << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void EfficiencyMatrix::endJob() {

}

// ------------ method called when starting to processes a run  ------------
void EfficiencyMatrix::beginRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a run  ------------
void EfficiencyMatrix::endRun(edm::Run const&, edm::EventSetup const&) {
}

// ------------ method called when starting to processes a luminosity block  ------------
void EfficiencyMatrix::beginLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method called when ending the processing of a luminosity block  ------------
void EfficiencyMatrix::endLuminosityBlock(edm::LuminosityBlock const&,
		edm::EventSetup const&) {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void EfficiencyMatrix::fillDescriptions(
		edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

// check if any ancestor of genParticle is in given list; check status of particle before use to avoid double counting!
bool EfficiencyMatrix::hasGoodMother(const reco::GenParticle* tau,
		int maxDepth = 50) {
	if (motherId_.empty() || motherId_.at(0) == 0) // pdgId=0: no check for mother
		return true;
	int depth = 0;
	int id = 0;
	const reco::Candidate* object = tau->mother();
	while (object->mother() && depth < maxDepth) {
		id = std::abs(object->pdgId());
		for (unsigned int i = 0; i < motherId_.size(); i++) {
			if (id == motherId_.at(i))
				return true;
		}
		++depth;
		object = object->mother();
	}
	return false;
}
//define this as a plug-in
DEFINE_FWK_MODULE(EfficiencyMatrix)
;
