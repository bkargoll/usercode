import FWCore.ParameterSet.Config as cms

process = cms.Process("EfficiencyMatrix")

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START42_V11::All' #off. DYToTauTau
#process.GlobalTag.globaltag = 'START42_V12::All' #most recent GT
#process.GlobalTag.globaltag = 'START42_V6::All' #Vladimirs Sample
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")

MaxEvents = 30000
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(MaxEvents))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# initialize TFileService and define output-file
process.TFileService = cms.Service("TFileService",
    fileName=cms.string('testFile.root'),
    closeFileFast=cms.untracked.bool(True)
)

# input files
if MaxEvents < 1001:
    process.load("GenLevelTools.GenTauChecker.DYToTauTau_local_cff") # local disk
    process.TFileService.fileName = 'test/EfficiencyMatrix_DYToTauTau_test.root'
else:
    process.load("GenLevelTools.GenTauChecker.DYToTauTau_cff") # dCache
    process.TFileService.fileName = '/user/kargoll/results/EfficiencyMatrix_DYToTauTau_chargedMatchingdR002dPt02Charge.root'

# rerun PF tau sequence
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# prepare genTau collection with visible p4 for matching
#process.prepMatch = cms.EDProducer('PrepareTauMatching',
#                                   verbose=cms.bool(False),
#                                   genParticles=cms.InputTag("genParticles")
#)
process.prepMatch = cms.EDProducer('PrepMatchCollectionsChargedMatching',
                                   verbose=cms.bool(False),
                                   genParticles=cms.InputTag("genParticles"),
                                   pfTaus=cms.InputTag("hpsPFTauProducer")
)

# match taus to genParticles
#process.mcMatcher = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
#    src=cms.InputTag("hpsPFTauProducer"), # RECO objects to match
#    matched=cms.InputTag("prepMatch", "genTausVisibleP4"), # mc-truth particle collection
#    mcPdgId=cms.vint32(15), # one or more PDG ID (15 = taus); absolute values (see below)
#    checkCharge=cms.bool(False), # True = require RECO and MC objects to have the same charge
#    mcStatus=cms.vint32(2), # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
#    maxDeltaR=cms.double(0.3), # Minimum deltaR for the match
#    maxDPtRel=cms.double(1000.0), # Minimum deltaPt/Pt for the match
#    resolveAmbiguities=cms.bool(True), # Forbid two RECO objects to match to the same GEN object
#    resolveByMatchQuality=cms.bool(True) # False = just match input in order; True = pick lowest deltaR pair first
#)
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
process.effMatrix = cms.EDAnalyzer('EfficiencyMatrix',
                           verbose=cms.bool(True),
                           genParticles=cms.InputTag("prepMatch", "genTausChargedP4"),
                           originalPFTaus=cms.InputTag("hpsPFTauProducer"),
                           PFTausForMatching=cms.InputTag("prepMatch","hpsTausChargedP4"),
                           PFTauMatchingOriginalMap=cms.InputTag("prepMatch","originalHPSTauMap"),
                           matching=cms.InputTag("mcMatcher"),
                           motherId=cms.vint32(23) # pdgID of tau mother, set to 0 if not needed
)

process.p = cms.Path(process.PFTau * 
                     process.prepMatch * 
                     process.mcMatcher * 
                     process.effMatrix)
