import FWCore.ParameterSet.Config as cms


cleanMuons = cms.EDFilter("PATMuonSelector",
                                  src = cms.InputTag("cleanLayer1Muons"),
                                  cut = cms.string(' ')
                                  )

cleanElectrons = cms.EDFilter("PATElectronSelector",
                                 src = cms.InputTag("cleanLayer1Electrons"),
                                 cut = cms.string(' !hasOverlaps("muons") ')
                                      )

cleanPhotons = cms.EDFilter("PATPhotonSelector",    # this filter does nothing, as overlapping photons are discarded during PAT-production
                                    src = cms.InputTag("cleanLayer1Photons"),
                                    cut = cms.string(' !hasOverlaps("electrons") ')
                                    )

cleanTaus = cms.EDFilter("PATTauSelector",
                                    src = cms.InputTag("cleanLayer1Taus"),
                                    cut = cms.string(' !hasOverlaps("electrons") & !hasOverlaps("muons") ')
                                    )

cleanJets = cms.EDFilter("PATJetSelector",
                                    src = cms.InputTag("cleanLayer1Jets"),
                                    cut = cms.string(' !hasOverlaps("electrons") & !hasOverlaps("muons") & !hasOverlaps("taus") & !hasOverlaps("photons") & !hasOverlaps("tkIsoElectrons") ')  # tkIsoElectrons is a subselection of electrons, so is needless here
                                 )   

cleanJetsSC5 = cms.EDFilter("PATJetSelector",
                                    src = cms.InputTag("cleanLayer1JetsSC5"),
                                    cut = cms.string(' !hasOverlaps("electrons") & !hasOverlaps("muons") & !hasOverlaps("taus") & !hasOverlaps("photons") & !hasOverlaps("tkIsoElectrons") ')  # tkIsoElectrons is a subselection of electrons, so is needless here
                                 )

# define sequence to clean all collections
cleanObjects = cms.Sequence( cleanMuons * cleanElectrons * cleanPhotons * cleanTaus * cleanJets * cleanJetsSC5 )
