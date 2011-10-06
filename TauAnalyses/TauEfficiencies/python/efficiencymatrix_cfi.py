import FWCore.ParameterSet.Config as cms

effMatrix = cms.EDAnalyzer('EfficiencyMatrix',
                           verbose=cms.bool(True),
                           genParticles=cms.InputTag("genParticles"),
                           tauCandidates=cms.InputTag(""),
                           tauDiscriminators=cms.vstring(""),
                           matching=cms.InputTag("")
)
