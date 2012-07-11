import FWCore.ParameterSet.Config as cms


process = cms.Process("MUTAU")

process.load("FWCore.MessageService.MessageLogger_cfi")

# enable the TrigReport and TimeReport
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

#process.load("Configuration.StandardSequences.GeometryPilot2_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_52_V9::All'



process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring(
                                                               '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v1/000/195/165/0000/4428CB4B-74AC-E111-A125-20CF3019DF0F.root'
                                                               ),
                            secondaryFileNames = cms.untracked.vstring(
                                                                       ),
                            inputCommands = cms.untracked.vstring( 
                                                                  'keep *'
                                                                  )
                            )

#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'goodList.json').getVLuminosityBlockRange()

# SingleMu trigger
#HLT_IsoMu24_eta2p1_v
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.hltSingleMu = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltSingleMu.HLTPaths = ["HLT_IsoMu24_eta2p1_v*"]

# calculate Jet Energy Corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')
#process.path = cms.Path(process.ak5PFJetsL1L2L3Residual)

# rerun PFTau sequence
process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

# Medium MVA isolation
process.offlineSelectedTausXXMisomvaLelecrej = cms.EDFilter( "PFTauSelector",
                                            src = cms.InputTag( "hpsPFTauProducer" ),
                                            cut = cms.string("pt>15."),
                                            discriminators = cms.VPSet(
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByDecayModeFinding" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumIsolationMVA" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightMuonRejection" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet( discriminator = cms.InputTag( "hpsPFTauDiscriminationByLooseElectronRejection") ,
                                                                                 #discriminator = cms.InputTag( "hpsPFTauDiscriminationAgainstElectronLoose" ),
                                                                                 selectionCut = cms.double( 0.5 )
                                                                                 ),
                                                                       ),
                                           filter = cms.bool(True)
                                           )
process.offlineSelectedTausXXMisomvaMelecrej = cms.EDFilter( "PFTauSelector",
                                            src = cms.InputTag( "hpsPFTauProducer" ),
                                            cut = cms.string("pt>15."),
                                            discriminators = cms.VPSet(
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByDecayModeFinding" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumIsolationMVA" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightMuonRejection" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet( discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumElectronRejection") ,
                                                                                 #discriminator = cms.InputTag( "hpsPFTauDiscriminationAgainstElectronLoose" ),
                                                                                 selectionCut = cms.double( 0.5 )
                                                                                 ),
                                                                       ),
                                           filter = cms.bool(True)
                                           )                                           
process.offlineSelectedTausXXMisomvaTelecrej = cms.EDFilter( "PFTauSelector",
                                            src = cms.InputTag( "hpsPFTauProducer" ),
                                            cut = cms.string("pt>15."),
                                            discriminators = cms.VPSet(
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByDecayModeFinding" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumIsolationMVA" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightMuonRejection" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet( discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightElectronRejection") ,
                                                                                 #discriminator = cms.InputTag( "hpsPFTauDiscriminationAgainstElectronLoose" ),
                                                                                 selectionCut = cms.double( 0.5 )
                                                                                 ),
                                                                       ),
                                           filter = cms.bool(True)
                                           )                                             
process.offlineSelectedTausXXMisomvaMVAelecrej = cms.EDFilter( "PFTauSelector",
                                            src = cms.InputTag( "hpsPFTauProducer" ),
                                            cut = cms.string("pt>15."),
                                            discriminators = cms.VPSet(
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByDecayModeFinding" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumIsolationMVA" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightMuonRejection" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet( discriminator = cms.InputTag( "hpsPFTauDiscriminationByMVAElectronRejection") ,
                                                                                 #discriminator = cms.InputTag( "hpsPFTauDiscriminationAgainstElectronLoose" ),
                                                                                 selectionCut = cms.double( 0.5 )
                                                                                 ),
                                                                       ),
                                           filter = cms.bool(True)
                                           )                                              
        
# medium deltaBeta isolation                                 
process.offlineSelectedTausXXMisodbLelecrej = cms.EDFilter( "PFTauSelector",
                                            src = cms.InputTag( "hpsPFTauProducer" ),
                                            cut = cms.string("pt>15."),
                                            discriminators = cms.VPSet(
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByDecayModeFinding" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightMuonRejection" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet( discriminator = cms.InputTag( "hpsPFTauDiscriminationByLooseElectronRejection") ,
                                                                                 #discriminator = cms.InputTag( "hpsPFTauDiscriminationAgainstElectronLoose" ),
                                                                                 selectionCut = cms.double( 0.5 )
                                                                                 ),
                                                                       ),
                                           filter = cms.bool(True)
                                           )
process.offlineSelectedTausXXMisodbMelecrej = cms.EDFilter( "PFTauSelector",
                                            src = cms.InputTag( "hpsPFTauProducer" ),
                                            cut = cms.string("pt>15."),
                                            discriminators = cms.VPSet(
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByDecayModeFinding" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightMuonRejection" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet( discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumElectronRejection") ,
                                                                                 #discriminator = cms.InputTag( "hpsPFTauDiscriminationAgainstElectronLoose" ),
                                                                                 selectionCut = cms.double( 0.5 )
                                                                                 ),
                                                                       ),
                                           filter = cms.bool(True)
                                           )                                           
process.offlineSelectedTausXXMisodbTelecrej = cms.EDFilter( "PFTauSelector",
                                            src = cms.InputTag( "hpsPFTauProducer" ),
                                            cut = cms.string("pt>15."),
                                            discriminators = cms.VPSet(
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByDecayModeFinding" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightMuonRejection" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet( discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightElectronRejection") ,
                                                                                 #discriminator = cms.InputTag( "hpsPFTauDiscriminationAgainstElectronLoose" ),
                                                                                 selectionCut = cms.double( 0.5 )
                                                                                 ),
                                                                       ),
                                           filter = cms.bool(True)
                                           )                                             
process.offlineSelectedTausXXMisodbMVAelecrej = cms.EDFilter( "PFTauSelector",
                                            src = cms.InputTag( "hpsPFTauProducer" ),
                                            cut = cms.string("pt>15."),
                                            discriminators = cms.VPSet(
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByDecayModeFinding" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet(  discriminator = cms.InputTag( "hpsPFTauDiscriminationByTightMuonRejection" ),
                                                                                  selectionCut = cms.double( 0.5 )
                                                                                  ),
                                                                       cms.PSet( discriminator = cms.InputTag( "hpsPFTauDiscriminationByMVAElectronRejection") ,
                                                                                 #discriminator = cms.InputTag( "hpsPFTauDiscriminationAgainstElectronLoose" ),
                                                                                 selectionCut = cms.double( 0.5 )
                                                                                 ),
                                                                       ),
                                           filter = cms.bool(True)
                                           )                                           
 
process.offlineSelectedMuons = cms.EDFilter( "MuonSelector",
    src = cms.InputTag('muons'),
    cut = cms.string("pt > 15.0 && isGlobalMuon && isPFMuon &&  globalTrack().normalizedChi2<10. && globalTrack().hitPattern().numberOfValidMuonHits>0 && numberOfMatchedStations() > 1 && abs(innerTrack().dxy)<0.2 && innerTrack().hitPattern().numberOfValidPixelHits>0 && track().hitPattern().trackerLayersWithMeasurement > 5 && (pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - 0.5*pfIsolationR04().sumPUPt, 0)) < 0.1*pt"),
    filter = cms.bool(True) 
    )

#define paths
process.looseEleMVAiso = cms.Path(process.hltSingleMu +
                                 process.offlineSelectedMuons +
                                 process.PFTau +
                                 process.offlineSelectedTausXXMisomvaLelecrej +
                                 process.ak5PFJetsL1L2L3Residual)
process.mediumEleMVAiso = cms.Path(process.hltSingleMu +
                                process.offlineSelectedMuons +
                                process.PFTau +
                                process.offlineSelectedTausXXMisomvaMelecrej +
                                process.ak5PFJetsL1L2L3Residual)
process.tightEleMVAiso = cms.Path(process.hltSingleMu +
                                 process.offlineSelectedMuons +
                                 process.PFTau +
                                 process.offlineSelectedTausXXMisomvaTelecrej +
                                 process.ak5PFJetsL1L2L3Residual)
process.mvaEleMVAiso = cms.Path(process.hltSingleMu +
                                 process.offlineSelectedMuons +
                                 process.PFTau +
                                 process.offlineSelectedTausXXMisomvaMVAelecrej +
                                 process.ak5PFJetsL1L2L3Residual)
process.looseEleDBiso = cms.Path(process.hltSingleMu +
                                 process.offlineSelectedMuons +
                                 process.PFTau +
                                 process.offlineSelectedTausXXMisodbLelecrej +
                                 process.ak5PFJetsL1L2L3Residual)
process.mediumEleDBiso = cms.Path(process.hltSingleMu +
                                process.offlineSelectedMuons +
                                process.PFTau +
                                process.offlineSelectedTausXXMisodbMelecrej +
                                process.ak5PFJetsL1L2L3Residual)
process.tightEleDBiso = cms.Path(process.hltSingleMu +
                                 process.offlineSelectedMuons +
                                 process.PFTau +
                                 process.offlineSelectedTausXXMisodbTelecrej +
                                 process.ak5PFJetsL1L2L3Residual)
process.mvaEleDBiso = cms.Path(process.hltSingleMu +
                                 process.offlineSelectedMuons +
                                 process.PFTau +
                                 process.offlineSelectedTausXXMisodbMVAelecrej +
                                 process.ak5PFJetsL1L2L3Residual)

process.hltOutputA = cms.OutputModule( "PoolOutputModule",
                                       fileName = cms.untracked.string( "offlineMuTauSelectionv3_HLTSingleMu_CorrectedJets.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW-RECO" )
    ),
    outputCommands = cms.untracked.vstring( 'keep *'),
    SelectEvents = cms.untracked.PSet(
                                      SelectEvents = cms.vstring('looseEleMVAiso','mediumEleMVAiso','tightEleMVAiso','mvaEleMVAiso','looseEleDBiso','mediumEleDBiso','tightEleDBiso','mvaEleDBiso')
    )
                                       )
process.HLTDQMOutput = cms.EndPath( process.hltOutputA)
