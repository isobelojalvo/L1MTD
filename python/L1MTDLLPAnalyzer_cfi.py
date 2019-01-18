import FWCore.ParameterSet.Config as cms

L1MTDLLPAnalyzer = cms.EDAnalyzer('L1MTDLLPAnalyzer',
                                  FTLBarrel         = cms.InputTag("", "", ""),
                                  FTLEndcap         = cms.InputTag("", "", ""),
                                  recHitBarrel      = cms.InputTag("", "", ""),
                                  recHitEndcap      = cms.InputTag("", "", ""),
                                  mtdClusterBarrel  = cms.InputTag("", "", ""),
                                  mtdClusterEndcap  = cms.InputTag("", "", ""),
                                  l1Taus            = cms.InputTag("", "", ""),
                                  time_cut          = cms.double(0),
                                  genParticles      = cms.InputTag("", "", "")
                                  )
