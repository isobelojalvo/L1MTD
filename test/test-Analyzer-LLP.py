# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --python_filename=rerun_step2_L1_onMCL1_FEVTHLTDEBUG.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 1 --era Phase2_timing --eventcontent FEVTDEBUGHLT --filein file:/afs/cern.ch/user/r/rekovic/release/CMSSW_9_3_2/src/step2_DIGI_PU200_10ev.root --conditions 93X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D17 --fileout file:step2_ZEE_PU200_1ev_rerun-L1-L1Ntuple.root --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleEMU
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras
from Configuration.Eras.Modifier_phase2_hgcalV9_cff import phase2_hgcalV9 # BBT, 01-18-19

#process = cms.Process('ANALYSIS')
process = cms.Process('ANALYSIS',eras.Phase2_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
#process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff')
#process.load("Configuration.Geometry.GeometryExtended2023D24_cff")

process.load("Geometry.MTDNumberingBuilder.mtdNumberingGeometry_cfi")

process.load("Geometry.MTDNumberingBuilder.mtdTopology_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdGeometry_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdParameters_cfi")
process.mtdGeometry.applyAlignment = cms.bool(False)

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(4000)
    input = cms.untracked.int32(40)
)

# Isobel
#outFileName = "timing-StopToBL-M-300-CTau1000-noPU.root"
#inFileName = "file:/hdfs/store/user/ojalvo/Timing-Samples/DisplacedSUSY_StopToBL_M-300/DisplacedSUSY_StopToBL_M-300_CTau-1000-noPU-gen-sim-digi-raw-reco.root"

# Ben
outFileName = "timing-DisplacedSUSY_StopToBL_M-300_CTau-1-200PU.root"
inFileName = "root://cmseos.fnal.gov//store/user/benjtann/MTD_LLP/DisplacedSUSY_StopToBL_M-300/DisplacedSUSY_StopToBL_M-300_CTau-1-gen-sim-digi-raw-reco.root"


#outFileName = "file:dyll-llp-official.root"
#inFileName = "/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/GEN-SIM-RECO/PU200_pilot_103X_upgrade2023_realistic_v2_ext2-v2/00000/FA39381B-9BBC-8C4A-90B8-A28CF803D044.root"
#inFileName = "file:/hdfs/store/user/ojalvo/Timing-Samples/DYLL-200PU/DYll-gen-sim-raw-reco.root"
#outFileName = "timing-SmuonToMuNu_M-200_CTau-1-noPU.root"
#inFileName = "file:/hdfs/store/user/ojalvo/Timing-Samples/DisplacedSUSY_SmuonToNeutralino_M-300/displacedSUSY_SmuonToMuNeutralino_M-200_CTau-1-noPU-gen-sim-digi-raw-reco.root"
#outFileName = "timing-StopToBL-M-300-CTau100-noPU.root"
#inFileName = "file:/hdfs/store/user/ojalvo/Timing-Samples/DisplacedSUSY_StopToBL_M-300/DisplacedSUSY_StopToBL_M-300_CTau-100-noPU-gen-sim-digi-raw-reco.root"


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        inFileName
        #"file:/afs/hep.wisc.edu/cms/ojalvo/triggerPhaseII/2018_devel/MTD_devel/2018/test-cmsdriver-10_4_X/PPD-PhaseIIMTDTDRAutumn18DR-00002.root"
        #"root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_4_0_mtd3/RelValTTbar_Tauola_14TeV/GEN-SIM-RECO/PU25ns_103X_upgrade2023_realistic_v2_2023D35PU200_4-v2/20000/E75806D2-0A52-7240-912A-2BBBF690B16A.root"
        #"file:DisplacedSUSY_StopToBL_M-300_CTau-1000-gen-sim-digi-raw-reco-20.root"
        #"root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_4_0_mtd2_patch1/RelValMinBias_14TeV/MINIAODSIM/103X_upgrade2023_realistic_v2_2023D35noPU-v1/20000/0F0F8B8D-DE01-FB49-9CC4-F1641C26E063.root"
        #"root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/MINIAODSIM/PU25ns_93X_upgrade2023_realistic_v5_2023D17PU200-v1/10000/6CE39BE9-EA2D-E811-8FDA-0242AC130002.root"
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


#process.source.secondaryFileNames = cms.untracked.vstring(
# "root://cmsxrootd.fnal.gov///store/relval/CMSSW_10_4_0_mtd2_patch1/RelValMinBias_14TeV/GEN-SIM-DIGI-RAW/103X_upgrade2023_realistic_v2_2023D35noPU-v1/20000/75DB809A-A348-7645-8DCC-8DA7070D24D9.root"
# )
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:282")

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
#process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '102X_upgrade2023_realistic_v5', '') # BBT, 01-18-19


process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

process.load('L1Trigger.L1THGCal.hgcalTriggerPrimitives_cff') # BBT, 01-18-19, uncomment
process.hgcl1tpg_step = cms.Path(process.hgcalTriggerPrimitives) # BBT, 01-18-19, uncomment

process.load('SimCalorimetry.EcalEBTrigPrimProducers.ecalEBTriggerPrimitiveDigis_cff')
process.EcalEBtp_step = cms.Path(process.simEcalEBTriggerPrimitiveDigis)

#process.L1TrackTrigger_step = cms.Path(process.L1TrackletTracksWithAssociators)

#process.VertexProducer.l1TracksInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")

process.simCaloStage2Layer1Digis.ecalToken = cms.InputTag("simEcalTriggerPrimitiveDigis") # BBT, 01-18-19
process.simCaloStage2Layer1Digis.hcalToken = cms.InputTag("simHcalTriggerPrimitiveDigis") # BBT, 01-18-19
process.load('L1Trigger.L1TCalorimeter.simCaloStage2Digis_cfi') # BBT, 01-18-19

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)

#process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string("fullDump.root"),
#    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
#)
#process.endjob_step = cms.EndPath(process.out)

process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

from RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi import *

process.load('RecoLocalFastTime.FTLClusterizer.mtdClusters_cfi')
process.mtdCluster_step = cms.Path(process.mtdClusters)


#track times
#process.load("SimTracker.TrackTriggerAssociation.ttTrackTimeValueMapProducer_cfi")

#process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
#    ComponentName = cms.string('TransientTrackBuilder')
#)

#process.ttTrackTimeValueMapProducer.tkTriggerTrackTruth = cms.InputTag("TTTrackAssociatorFromPixelDigis", "Level1TTTracks")
#process.ttTrackTimeValueMapProducer.tkTriggerTrackSrc = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks")

#process.timingtracks = cms.Path(process.ttTrackTimeValueMapProducer)

#process.load("L1Trigger.L1MTD.L1TrackerEtMissProducertime_cfi")
#process.L1TrackerEtMisstime.timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer","TTTracksFromTrackletL1ConfigurableFlatResolutionModel")
#process.L1TrackerEtMisstime.time_cut = cms.double(0.150)
#process.mettime = cms.Path(process.L1TrackerEtMisstime)#

############################################################
# L1 pf object
###########################################################
#process.load("L1Trigger.Phase2L1ParticleFlow.pfTracksFromL1Tracks_cfi")
#from L1Trigger.Phase2L1ParticleFlow.l1ParticleFlow_cff import *
#process.l1pf = cms.Path(process.pfTracksFromL1Tracks+process.l1ParticleFlow)

############################################################
# L1 Tau object
############################################################
#process.load("L1Trigger.Phase2L1Taus.L1PFTauProducer_cff")
#process.L1PFTauProducer.min_pi0pt = cms.double(2.5);
#process.L1PFTauProducer.L1PFObjects = cms.InputTag("l1pfProducer","PF")
#process.L1PFTauProducer.L1Neutrals = cms.InputTag("l1pfProducer")
#process.L1PFTauProducer.L1Clusters = cms.InputTag("l1pfProducer","PF")
#process.L1PFTaus = cms.Path(process.L1PFTauProducer)

##the analyzer

process.L1MTDAnalyzer = cms.EDAnalyzer('L1MTDLLPAnalyzer',
                                       BTLMinimumEnergy = cms.double(2),
                                       FTLBarrel = cms.InputTag("mix","FTLBarrel","HLT"),
                                       FTLEndcap = cms.InputTag("mix","FTLEndcap","HLT"),
                                       recHitBarrel = cms.InputTag("mtdRecHits","FTLBarrel","RECO"),
                                       recHitEndcap = cms.InputTag("mtdRecHits","FTLEndcap","RECO"),
                                       mtdClusterBarrel = cms.InputTag("mtdClusters", "FTLBarrel", "ANALYSIS"),
                                       mtdClusterEndcap = cms.InputTag("mtdClusters", "FTLEndcap", "ANALYSIS"),
                                       l1Taus = cms.InputTag("l1extraParticles","Central","RECO"),
                                       time_cut  = cms.double(0.150),
                                       genParticles = cms.InputTag("genParticles","","HLT"),
                                       )

process.analyzer = cms.Path(process.L1MTDAnalyzer)


process.TFileService = cms.Service("TFileService", 
   fileName = cms.string(outFileName), 
   closeFileFast = cms.untracked.bool(True)
)

# Schedule definition

process.schedule = cms.Schedule(process.EcalEBtp_step,process.mtdCluster_step,process.analyzer)
#process.schedule = cms.Schedule(process.EcalEBtp_step,process.hgcl1tpg_step,process.mtdCluster_step,process.analyzer) # BBT, 01-18-19
#process.schedule = cms.Schedule(process.EcalEBtp_step,process.mtdCluster_step,process.analyzer) # BBT, 01-18-19
#process.schedule = cms.Schedule(process.L1simulation_step,process.EcalEBtp_step,process.mtdCluster_step,process.analyzer) # BBT, 01-18-19
#process.schedule = cms.Schedule(process.EcalEBtp_step,process.mtdCluster_step,process.analyzer,process.L1simulation_step) # BBT, 01-18-19
#process.L1simulation_step,process.timingtracks,process.l1pf,process.L1PFTaus,process.mettime,process.analyzer,process.endjob_step) 

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
