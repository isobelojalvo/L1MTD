# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --python_filename=rerun_step2_L1_onMCL1_FEVTHLTDEBUG.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 1 --era Phase2_timing --eventcontent FEVTDEBUGHLT --filein file:/afs/cern.ch/user/r/rekovic/release/CMSSW_9_3_2/src/step2_DIGI_PU200_10ev.root --conditions 93X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D17 --fileout file:step2_ZEE_PU200_1ev_rerun-L1-L1Ntuple.root --customise=L1Trigger/L1TNtuples/customiseL1Ntuple.L1NtupleEMU
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

#process = cms.Process('ANALYSIS')
process = cms.Process('ANALYSIS',eras.Phase2_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

### not tested
#process.load('L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff')
process.load("Geometry.MTDNumberingBuilder.mtdNumberingGeometry_cfi")

process.load("Geometry.MTDNumberingBuilder.mtdTopology_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdGeometry_cfi")
process.load("Geometry.MTDGeometryBuilder.mtdParameters_cfi")
process.mtdGeometry.applyAlignment = cms.bool(False)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)

outFileName = "timing-SmuonToMuNu_M-200_CTau-100-200PU.root"
inFileName  = "file:/scratch/ojalvo/DisplacedSUSY_Smuon-M-200_CTau-100-200PU/displacedSUSY_SmuonToMuNeutralino_M-200_CTau-100-gen-sim-digi-raw-reco-100.root"
#outFileName = "file:dyll-llp-official.root"
#inFileName = "/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/GEN-SIM-RECO/PU200_pilot_103X_upgrade2023_realistic_v2_ext2-v2/00000/FA39381B-9BBC-8C4A-90B8-A28CF803D044.root"
#inFileName = "file:/hdfs/store/user/ojalvo/Timing-Samples/DYLL-200PU/DYll-gen-sim-raw-reco.root"
#outFileName = "timing-SmuonToMuNu_M-200_CTau-100-noPU.root"
#inFileName = "file:/scratch/ojalvo/DisplacedSUSY_Smuon-M-200_CTau-100-0PU/displacedSUSY_SmuonToMuNeutralino_M-200_CTau-100-noPU-gen-sim-digi-raw-reco.root"
#inFileName = "file:pickevents.root" #
#outFileName = "timing-StopToBL-M-300-CTau100-noPU.root"
#inFileName = "file:/hdfs/store/user/ojalvo/Timing-Samples/DisplacedSUSY_StopToBL_M-300/DisplacedSUSY_StopToBL_M-300_CTau-100-noPU-gen-sim-digi-raw-reco.root"


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        inFileName
        ),
    inputCommands = cms.untracked.vstring("keep *", 
        "drop l1tHGCalTowerMapBXVector_hgcalTriggerPrimitiveDigiProducer_towerMap_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT",
        "drop l1tEMTFHit2016s_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016s_simEmtfDigis__HLT"),
     #skipEvents = cms.untracked.uint32(80)
    eventsToSkip = cms.untracked.VEventRange("1:144","1:519","1:507","1:515")
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
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')


process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

process.load('SimCalorimetry.EcalEBTrigPrimProducers.ecalEBTriggerPrimitiveDigis_cff')
process.EcalEBtp_step = cms.Path(process.simEcalEBTriggerPrimitiveDigis)

# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)

process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

from RecoLocalFastTime.FTLClusterizer.MTDCPEESProducer_cfi import *

process.load('RecoLocalFastTime.FTLClusterizer.mtdClusters_cfi')
process.mtdCluster_step = cms.Path(process.mtdClusters)


process.L1MTDAnalyzer = cms.EDAnalyzer('L1MTDLLPAnalyzer',
                                       BTLMinimumEnergy = cms.double(2),
                                       FTLBarrel = cms.InputTag("mix","FTLBarrel","HLT"),
                                       FTLEndcap = cms.InputTag("mix","FTLEndcap","HLT"),
                                       recHitBarrel = cms.InputTag("mtdRecHits","FTLBarrel","RECO"),
                                       recHitEndcap = cms.InputTag("mtdRecHits","FTLEndcap","RECO"),
                                       mtdClusterBarrel = cms.InputTag("mtdClusters", "FTLBarrel", "ANALYSIS"),
                                       mtdClusterEndcap = cms.InputTag("mtdClusters", "FTLEndcap", "ANALYSIS"),
                                       BTLSimHits   = cms.InputTag("g4SimHits","FastTimerHitsBarrel","SIM"),
                                       l1Taus    = cms.InputTag("l1extraParticles","Central","RECO"),
                                       HepMCProduct     = cms.InputTag("generatorSmeared","","SIM"),
                                       genParticles     = cms.InputTag("genParticles","","HLT"),
                                       time_cut  = cms.double(0.150)
                                       )

process.analyzer = cms.Path(process.L1MTDAnalyzer)


process.TFileService = cms.Service("TFileService", 
   fileName = cms.string(outFileName), 
   closeFileFast = cms.untracked.bool(True)
)

# Schedule definition

process.schedule = cms.Schedule(process.EcalEBtp_step,process.mtdCluster_step,process.analyzer)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#dump_file = open('dump.py','w')
#dump_file.write(process.dumpPython())
