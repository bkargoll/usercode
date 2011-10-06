#ifndef GenLevelTools_GenTauTools_cc
#define GenLevelTools_GenTauTools_cc

// -*- C++ -*-
//
// Package:    GenTauTools
// Class:      GenTauTools
//
/**\class GenTauTools GenTauTools.cc GenLevelTools/GenTauTools/src/GenTauTools.cc

 Description: tools to analyze generator information

 */
//
// Original Author:  Bastian Kargoll
//         Created:  Thu Jun 30 16:37:00 CEST 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/View.h"
#include "Math/Vector4D.h"
#include "TMath.h"

// class declaration
class GenTauTools {
public:
	GenTauTools() {
	}

	~GenTauTools() {
	}
	;

	/// determine decay mode of a given generator tau
	unsigned int genDecayMode(const reco::GenParticle* tau) {
		bool muon = false, electron = false;
		unsigned int prongs = 0, piZeros = 0, KProngs = 0, KZeros = 0;
		bool verbose = false;
		GenTauTools::checkDaughters(tau, muon, electron, prongs, piZeros,
				KProngs, KZeros, verbose);
		if (electron)
			return 16;
		if (muon)
			return 17;
		int Neutrals = piZeros + KZeros;
		int Prongs = prongs + KProngs;
		int result = (Prongs - 1) * 5 + Neutrals;
		if (result < 0)
			result = 19; // not determined = 19
		return result;
	}

	// calculate charged part of 4momentum of given particle
	void chargedP4(const reco::Candidate* object, ROOT::Math::LorentzVector<
			ROOT::Math::PxPyPzE4D<double> >& p4, bool verbose = false) {
		// loop over daughters
		for (size_t d = 0; d < object->numberOfDaughters(); d++) {
			const reco::Candidate* dau = object->daughter(d);
			if (verbose) {
				std::cout << "Daughter No " << d + 1 << ": " << dau->pdgId()
						<< " | " << dau->status() << std::endl;
			}
			switch (dau->status()) {
			case 1:
				if (dau->charge()) {
					p4 += dau->p4();
					if (verbose)
						std::cout << "charged particle, pdgId  "
								<< dau->pdgId() << std::endl;
				} else {
					if (verbose)
						std::cout << "neutral particle, pdgId  "
								<< dau->pdgId() << std::endl;
				}
				break;
			case 2:
			case 3:
				if (verbose)
					std::cout << "recursive call " << std::endl;
				chargedP4(dau, p4, verbose);
				break;
			default:
				std::cout << "ERROR: genParticle with status != 1,2 or 3"
						<< std::endl;
			}

		}
	}

	// calculate visible part of 4momentum of given particle
	void visibleP4(const reco::Candidate* object,
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
				visibleP4(dau, visP4, verbose);
				break;
			default:
				std::cout << "ERROR: genParticle with status != 1,2 or 3"
						<< std::endl;
			}

		}
	}

private:

	/// method to search daughters for interesting stuff
	void checkDaughters(const reco::Candidate* object, bool& muon,
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
					std::cout << "ERROR: electron with wrong status"
							<< std::endl;
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
					std::cout << "ERROR: neutrino with wrong status"
							<< std::endl;
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

	// ----------member data ---------------------------
};

#endif
