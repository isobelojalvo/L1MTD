import FWCore.ParameterSet.Config as cms

TrackAnalyzer = cms.EDAnalyzer('L1MTDAnalyzer',
                               #L1Clusters       = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                               L1TrackInputTag     = cms.InputTag("", "", ""), 
                               timingValuesNominal = cms.InputTag("", "", ""), 
                               timingValuesSmeared = cms.InputTag("", "", ""), 
                               genParticles_xyz    = cms.InputTag("", "", ""),
                               genParticles_t      = cms.InputTag("", "", "")
                               #packedCandidates = cms.InputTag("packedPFCandidates","","RECO"),
                               #ecalTPGsBarrel = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT")
)
