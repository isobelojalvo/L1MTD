# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --python_filename=rerun_step2_L1_onMCL1_FEVTHLTDEBUG.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 1 --era Phase2_timing --eventcontent FEVTDEBUGHLT --filein file:/afs/cern.ch/user/r/rekovic/release/CMSSW_9_3_2/src/step2_DIGI_PU200_10ev.root --conditions 93X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D17 --fileout file:step2_ZEE_PU200_1ev_rerun-L1-L1Ntuple.root --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleEMU
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('L1',eras.Phase2_trigger)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(4000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        $inputFileNames        
        ),
    inputCommands = cms.untracked.vstring("keep *", 
        "drop l1tHGCalTowerMapBXVector_hgcalTriggerPrimitiveDigiProducer_towerMap_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT",
        "drop l1tEMTFHit2016s_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016s_simEmtfDigis__HLT")
     #skipEvents = cms.untracked.uint32(80)
)

process.options = cms.untracked.PSet(
)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:test_reprocess.root'),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')


process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff')
process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives)

process.load('SimCalorimetry.EcalEBTrigPrimProducers.ecalEBTriggerPrimitiveDigis_cff')
process.EcalEBtp_step = cms.Path(process.simEcalEBTriggerPrimitiveDigis)

process.L1TrackTrigger_step = cms.Path(process.L1TrackletTracksWithAssociators)

process.VertexProducer.l1TracksInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

#track times
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)

#track times
process.load("SimTracker.TrackTriggerAssociation.ttTrackTimeValueMapProducer_cfi")

process.ttTrackTimeValueMapProducer50ps                     = process.ttTrackTimeValueMapProducer.clone()
process.ttTrackTimeValueMapProducer50ps.tkTriggerTrackTruth = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks")
process.ttTrackTimeValueMapProducer50ps.tkTriggerTrackSrc   = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")

process.ttTrackTimeValueMapProducer30ps                     = process.ttTrackTimeValueMapProducer.clone()
process.ttTrackTimeValueMapProducer30ps.tkTriggerTrackTruth = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks")
process.ttTrackTimeValueMapProducer30ps.tkTriggerTrackSrc   = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")
process.ttTrackTimeValueMapProducer30ps.resolutionModels    = cms.VPSet(cms.PSet( modelName = cms.string('L1ConfigurableFlatResolutionModel'),
                                                                                  resolutionInNs = cms.double(0.030) ),
                                                                        cms.PSet( modelName = cms.string('L1PerfectResolutionModel') ) )

process.ttTrackTimeValueMapProducer100ps                     = process.ttTrackTimeValueMapProducer.clone()
process.ttTrackTimeValueMapProducer100ps.tkTriggerTrackTruth = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks")
process.ttTrackTimeValueMapProducer100ps.tkTriggerTrackSrc   = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")
process.ttTrackTimeValueMapProducer100ps.resolutionModels    = cms.VPSet(cms.PSet( modelName = cms.string('L1ConfigurableFlatResolutionModel'),
                                                                                   resolutionInNs = cms.double(0.100) ),
                                                                         cms.PSet( modelName = cms.string('L1PerfectResolutionModel') ) )



process.timingtracks = cms.Path(process.ttTrackTimeValueMapProducer50ps*process.ttTrackTimeValueMapProducer30ps*process.ttTrackTimeValueMapProducer100ps)

############################################################
# L1 pf object
###########################################################
process.load("L1Trigger.Phase2L1ParticleFlow.pfTracksFromL1Tracks_cfi")
from L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff import *
process.l1pf = cms.Path(process.pfTracksFromL1Tracks+process.l1ParticleFlow)

############################################################
# L1 Tau object
############################################################
process.load("L1Trigger.Phase2L1Taus.L1PFTauProducer_cff")
process.L1PFTauProducer.min_pi0pt = cms.double(2.5);
process.L1PFTauProducer.L1PFObjects = cms.InputTag("l1pfProducer","PF")
process.L1PFTauProducer.L1Neutrals = cms.InputTag("l1pfProducer")
process.L1PFTauProducer.L1Clusters = cms.InputTag("l1pfProducer","PF")
process.L1PFTaus = cms.Path(process.L1PFTauProducer)

##the analyzer
process.L1MTDAnalyzer50ps = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                           timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer50ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                           L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                           l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                           l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                           l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           time_cut         = cms.double(0.150),
                                           isoConeDeltaR    = cms.double(0.4),
                                           isoConeDeltaZ    = cms.double(1.0) 
                                           )

process.L1MTDAnalyzer50ps1sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                                 timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer50ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                                 L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                                 l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                                 l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                                 l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                 l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                 time_cut         = cms.double(0.050),
                                                 isoConeDeltaR    = cms.double(0.4),
                                                 isoConeDeltaZ    = cms.double(1.0) 
                                                 )

process.L1MTDAnalyzer50ps2sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                                 timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer50ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                                 L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                                 l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                                 l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                                 l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                 l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                 time_cut         = cms.double(0.100),
                                                 isoConeDeltaR    = cms.double(0.4),
                                                 isoConeDeltaZ    = cms.double(1.0) 
                                                 )
##30ps
process.L1MTDAnalyzer30ps = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                           timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer30ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                           L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                           l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                           l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                           l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           time_cut         = cms.double(0.090),
                                           isoConeDeltaR    = cms.double(0.4),
                                           isoConeDeltaZ    = cms.double(1.0) 
                                           )

process.L1MTDAnalyzer30ps1sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                                 timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer30ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                                 L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                                 l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                                 l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                                 l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                 l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                 time_cut         = cms.double(0.030),
                                                 isoConeDeltaR    = cms.double(0.4),
                                                 isoConeDeltaZ    = cms.double(1.0) 
                                                 )

process.L1MTDAnalyzer30ps2sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                                 timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer30ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                                 L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                                 l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                                 l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                                 l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                 l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                 time_cut         = cms.double(0.060),
                                                 isoConeDeltaR    = cms.double(0.4),
                                                 isoConeDeltaZ    = cms.double(1.0) 
                                                 )

process.L1MTDAnalyzer100ps = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                           timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer100ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                           L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                           l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                           l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                           l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           time_cut         = cms.double(0.300),
                                           isoConeDeltaR    = cms.double(0.4),
                                           isoConeDeltaZ    = cms.double(1.0) 
                                           )

process.L1MTDAnalyzer100ps1sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                                  timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer100ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                                  L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                                  l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                                  l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                                  l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                  l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                  time_cut         = cms.double(0.100),
                                                  isoConeDeltaR    = cms.double(0.4),
                                                  isoConeDeltaZ    = cms.double(1.0) 
                                                  )

process.L1MTDAnalyzer100ps2sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                                  timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer100ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                                  L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                                  l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                                  l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                                  l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                  l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                                  time_cut         = cms.double(0.200),
                                                  isoConeDeltaR    = cms.double(0.4),
                                                  isoConeDeltaZ    = cms.double(1.0) 
                                                  )



process.analyzer = cms.Path(process.L1MTDAnalyzer50ps+process.L1MTDAnalyzer100ps+process.L1MTDAnalyzer30ps+process.L1MTDAnalyzer50ps1sigma+process.L1MTDAnalyzer100ps1sigma+process.L1MTDAnalyzer30ps1sigma+process.L1MTDAnalyzer50ps2sigma+process.L1MTDAnalyzer100ps2sigma+process.L1MTDAnalyzer30ps2sigma)


process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("analyzer.root"), 
   closeFileFast = cms.untracked.bool(True)
)

# Schedule definition

process.schedule = cms.Schedule(process.EcalEBtp_step,process.L1TrackTrigger_step,process.L1simulation_step,process.timingtracks,process.l1pf,process.L1PFTaus,process.analyzer,process.endjob_step) 

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )



process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("$outputFileName")
)
