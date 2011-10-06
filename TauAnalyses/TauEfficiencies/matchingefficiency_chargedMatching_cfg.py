import FWCore.ParameterSet.Config as cms

process = cms.Process("MatchingEfficiencyChargedMatching")

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START42_V11::All' #off. DYToTauTau
#process.GlobalTag.globaltag = 'START42_V12::All' #newest
#process.GlobalTag.globaltag = 'START42_V6::All' #Vladimirs Sample
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

MaxEvents = 20000
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(MaxEvents))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# input files
if MaxEvents < 1001:
    process.load("GenLevelTools.GenTauChecker.DYToTauTau_local_cff") # local disk
else:
    process.load("GenLevelTools.GenTauChecker.DYToTauTau_cff") # dCache


# initialize TFileService and define output-file
process.TFileService = cms.Service("TFileService",
#    fileName=cms.string('test/MatchingEfficiency_DYToTauTau.root'),
    fileName=cms.string('/user/kargoll/results/MatchingEfficiency_DYToTauTau_chargedMatchingdR002dPt02Charge.root'),
    closeFileFast=cms.untracked.bool(True)
)

# rerun PF tau sequence
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# prepare genTau collection with visible p4 for matching
process.prepMatch = cms.EDProducer('PrepMatchCollectionsChargedMatching',
                                   verbose=cms.bool(False),
                                   genParticles=cms.InputTag("genParticles"),
                                   pfTaus=cms.InputTag("hpsPFTauProducer")
)

# match taus to genParticles with their charged constituents
process.mcMatcher = cms.EDProducer("MCMatcher",     # cut on deltaR, deltaPt/Pt; pick best by deltaR
    src     = cms.InputTag("prepMatch","hpsTausChargedP4"),     # RECO objects to match
    matched = cms.InputTag("prepMatch","genTausChargedP4"), # mc-truth particle collection
    mcPdgId     = cms.vint32(15),           # one or more PDG ID (15 = taus); absolute values (see below)
    checkCharge = cms.bool(True),           # True = require RECO and MC objects to have the same charge
    mcStatus = cms.vint32(2),               # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
    maxDeltaR = cms.double(0.02),            # Minimum deltaR for the match
    maxDPtRel = cms.double(0.2),           # Minimum deltaPt/Pt for the match
    resolveAmbiguities = cms.bool(True),    # Forbid two RECO objects to match to the same GEN object
    resolveByMatchQuality = cms.bool(True) # False = just match input in order; True = pick lowest deltaR pair first
)

# configure the efficiency matrix analyzer
process.matchEff = cms.EDAnalyzer('MatchingEfficiency',
                           verbose=cms.bool(False),
                           genParticles=cms.InputTag("genParticles"),
                           genTausForMatching=cms.InputTag("prepMatch","genTausChargedP4"),
                           genTausMapToGenParticles=cms.InputTag("prepMatch","originalGenParticleMap"),
                           tauCandidates=cms.InputTag("hpsPFTauProducer"),
                           tauCandidatesForMatching=cms.InputTag("prepMatch","hpsTausChargedP4"),
                           MatchingTausMapToTauCandidates=cms.InputTag("prepMatch","originalHPSTauMap"),
                           matching=cms.InputTag("mcMatcher"),
                           minGenPt=cms.double(0.0), # Minimum pt of generator taus
                           motherId=cms.vint32(0) # pdgID of tau mother, set to 0 if not needed
)

process.p = cms.Path(process.PFTau *
                     process.prepMatch *
                     process.mcMatcher *
                     process.matchEff)
