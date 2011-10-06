import FWCore.ParameterSet.Config as cms

process = cms.Process("GenTauChecker")

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'START42_V11::All' #off. DYToTauTau
process.GlobalTag.globaltag = 'START42_V6::All' #Vladimirs Sample
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(100000))
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

# input files
#inputDir = '/user/kargoll/ZtautauSM/3a40cf5c572f03c73d0aee82b8fe5cca/'
#filesDat = open("/user/kargoll/ZtautauSM/ZtautauSM.dat", "r")
##inputDir = '/user/kargoll/WJet_14_05_11_HLT/3a40cf5c572f03c73d0aee82b8fe5cca/'
##filesDat = open("/user/kargoll/WJet_14_05_11_HLT/WJet_14_05_11_HLT.dat","r")
#fileList = []
#for file in filesDat:
#    fileList.append('file:' + inputDir + file.strip())
#process.source = cms.Source("PoolSource",
#    fileNames=cms.untracked.vstring(fileList)
#)
process.load("GenLevelTools.GenTauChecker.DYToTauTau_cff")

# initialize TFileService and define output-file
process.TFileService = cms.Service("TFileService",
    fileName=cms.string('test/GenTauChecker_DYToTauTau.root'),
    closeFileFast=cms.untracked.bool(True)
)

process.check = cms.EDAnalyzer('GenTauChecker',
                               genParticles=cms.InputTag("genParticles")
)


process.p = cms.Path(process.check)
