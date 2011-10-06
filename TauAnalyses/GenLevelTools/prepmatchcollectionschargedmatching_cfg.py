import FWCore.ParameterSet.Config as cms

process = cms.Process("chargedMatching")

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

MaxEvents = 1000
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(MaxEvents))
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# input files
if MaxEvents < 1001:
    process.load("GenLevelTools.GenTauChecker.DYToTauTau_local_cff") # local disk
else:
    process.load("GenLevelTools.GenTauChecker.DYToTauTau_cff") # dCache

process.prepMatch = cms.EDProducer('PrepMatchCollectionsChargedMatching',
                                   verbose=cms.bool(False),
                                   genParticles=cms.InputTag("genParticles"),
                                   pfTaus=cms.InputTag("hpsPFTauProducer")
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName=cms.untracked.string('/user/kargoll/results/prepareTauMatching_test.root'),
    outputCommands=cms.untracked.vstring('drop *',
      "keep *_genParticles_*_*",
      "keep *_prepMatch_*_*")
)

  
process.p = cms.Path(process.prepMatch)

process.e = cms.EndPath(process.out)