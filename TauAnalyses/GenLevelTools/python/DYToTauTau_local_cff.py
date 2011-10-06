import FWCore.ParameterSet.Config as cms

source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(
        'file:///user/kargoll/DYToTauTau/4E971DBC-B97C-E011-9CC4-003048D4DFA8.root'
        )
    )

