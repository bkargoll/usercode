import FWCore.ParameterSet.Config as cms

process = cms.Process("Test")

process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    'file:/home/home2/institut_3b/kargoll/CMSSW/CMSSW_2_2_13/src/PhysicsTools/PatExamples/input/FullSimTTBar-2_2_X_2008-11-03-STARTUP_V7-AODSIM.100.root'

  )
)

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('STARTUP_V7::All')
process.load("Configuration.StandardSequences.MagneticField_cff")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.load("PhysicsTools.PatAlgos.patSequences_cff")

process.p = cms.Path(
    process.patDefaultSequence  
)

process.MessageLogger = cms.Service("MessageLogger")

process.analyzeBasicPat = cms.EDFilter("PatBasicAnalyzer",
  photonSrc   = cms.untracked.InputTag("cleanLayer1Photons"),
  electronSrc = cms.untracked.InputTag("cleanLayer1Electrons"),
  muonSrc     = cms.untracked.InputTag("cleanLayer1Muons"),                                             
  tauSrc      = cms.untracked.InputTag("cleanLayer1Taus"),
  jetSrc      = cms.untracked.InputTag("cleanLayer1Jets"),
  metSrc      = cms.untracked.InputTag("layer1METs")
)

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('analyzePatBasics.root')
                                   )

process.q = cms.Path(process.analyzeBasicPat)

#process.Tracer = cms.Service('Tracer')


