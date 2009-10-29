######## config-file to test the cleaning tool ########

import FWCore.ParameterSet.Config as cms

process = cms.Process("TestClean")

# load message logger and define output report
process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

# define maximum number of processed events and how often the processed event is reported
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 50

# source: Heikos Test-PAT-Tupel
process.source = cms.Source("PoolSource",
       fileNames = cms.untracked.vstring( 'file:/user/kargoll/TTbar_Summer09_MC_31X_V3_7TeV-v1_GENSIMRECO.geenen.test.pat.root' )
)

# load TFileService
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('testCleaning.root')
                                   )

# load the cleaning tool
process.load("Tools.Cleaning.cleaning_cff")

# define analyzer to check the new selection
process.content = cms.EDAnalyzer("CleaningChecker",
  photonSrc   = cms.untracked.InputTag("cleanLayer1Photons"),
  electronSrc = cms.untracked.InputTag("cleanLayer1Electrons"),
  muonSrc     = cms.untracked.InputTag("cleanLayer1Muons"),                                             
  tauSrc      = cms.untracked.InputTag("cleanLayer1Taus"),
  jetSrc      = cms.untracked.InputTag("cleanLayer1Jets"),
  cleanphotonSrc   = cms.untracked.InputTag("cleanPhotons"),
  cleanelectronSrc = cms.untracked.InputTag("cleanElectrons"),
  cleanmuonSrc     = cms.untracked.InputTag("cleanMuons"),                                             
  cleantauSrc      = cms.untracked.InputTag("cleanTaus"),
  cleanjetSrc      = cms.untracked.InputTag("cleanJets")
)

process.p = cms.Path( process.cleanObjects * process.content)
