# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: repr --processName=REPR --python_filename=reprocess_test_10_5_0_pre1.py --no_exec -s L1 --datatier GEN-SIM-DIGI-RAW -n 2 --era Phase2 --eventcontent FEVTDEBUGHLT --filein root://cms-xrd-global.cern.ch//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2_ext4-v1/280000/FF5C31D5-D96E-5E48-B97F-61A0E00DF5C4.root --conditions 103X_upgrade2023_realistic_v2 --beamspot HLLHC14TeV --geometry Extended2023D28 --fileout file:step2_2ev_reprocess_slim.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('L1',eras.Phase2C4_trigger)
#process = cms.Process('REPR',eras.Phase2C4_timing_layer_bar)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D35Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D35_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    #fileNames = cms.untracked.vstring('file:/hdfs/store/mc/PhaseIIMTDTDRAutumn18DR/TTbar_14TeV_TuneCP5_Pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/80000/B79DA3D0-B7E3-6246-864B-4F62E4621D15.root')
    fileNames = cms.untracked.vstring('/store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DYToLL_M-50_14TeV_TuneCP5_pythia8/MINIAODSIM/NoPU_103X_upgrade2023_realistic_v2-v2/90000/93BDDDAE-6FE2-834B-A94A-0ADE4FC6AAC3.root'),
    secondaryFileNames = cms.untracked.vstring('file:/hdfs/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/485E90DF-FEF5-1D49-9AD5-161B430BBC80.root',
                                               'file:/hdfs/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/4E30FBFB-E0FB-414C-96BC-FC3CAD3D3E1F.root',
                                               'file:/hdfs/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/4FB6298E-9700-A844-86A7-CAC802657F78.root',
                                               'file:/hdfs/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/CBF8AF2E-06DE-1844-9A67-704D678674EB.root',
                                               'file:/hdfs/store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_TuneCP5_pythia8/FEVT/NoPU_103X_upgrade2023_realistic_v2-v2/90000/E36263AC-70DE-FD44-826F-898852F4598D.root'
                                               )
    #secondaryFileNames = cms.untracked.vstring('/store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DYToLL_M-50_14TeV_TuneCP5_pythia8/MINIAODSIM/NoPU_103X_upgrade2023_realistic_v2-v2/90000/93BDDDAE-6FE2-834B-A94A-0ADE4FC6AAC3.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('repr nevts:2'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step2_2ev_reprocess_slim.root'),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

process.load('L1Trigger.mtdTPGenerator.mtdTPGenerator_cfi')
process.mtdTPGenerator.L1PFObjects = cms.InputTag("l1pfProducer","PF","L1")
process.MTD_step = cms.Path(process.mtdTPGenerator)


# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

process.GlobalTag = GlobalTag(process.GlobalTag, '103X_upgrade2023_realistic_v2', '') 

process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
process.load('CalibCalorimetry.CaloTPG.CaloTPGTranscoder_cfi')

############################################################
# L1 Tau object
############################################################

process.load("L1Trigger.Phase2L1Taus.L1PFTauProducer_cff")
process.L1PFTauProducer.debug = cms.untracked.bool(True)
process.L1PFTauProducer.L1PFObjects = cms.InputTag("mtdTPGenerator","Time") ##FIXME
process.L1PFTauProducer.L1Neutrals = cms.InputTag("mtdTPGenerator","Time")
process.L1PFTaus = cms.Path(process.L1PFTauProducer)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("fullDump.root"),
    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
)
#process.endjob_step = cms.EndPath(process.out)

process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("analyzer.root"), 
   closeFileFast = cms.untracked.bool(True)
)


##the analyzer
process.L1MTDAnalyzer = cms.EDAnalyzer('L1MTDPFAnalyzer',
                                       #timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer","TTTracksFromTrackletL1PerfectResolutionModel"),
                                       L1VertexInputTag = cms.InputTag("L1TkPrimaryVertex"),
                                       L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                       l1PFCands        = cms.InputTag("mtdTPGenerator","Time"),
                                       l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus","L1"),##check me
                                       l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                       l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                       recoElectrons    = cms.InputTag("slimmedElectrons", "", "PAT"),
                                       recoPhotons      = cms.InputTag("slimmedPhotons","","PAT"),
                                       recoMuons        = cms.InputTag("slimmedMuons", "", "PAT"),
                                       recoTaus         = cms.InputTag("slimmedTaus", "", "PAT"),
                                       recoMet          = cms.InputTag("slimmedMETs","","PAT"),
                                       genParticles     = cms.InputTag("genParticles", "", "HLT"),
                                       genJets          = cms.InputTag("ak4GenJetsNoNu","","HLT"),
                                       time_cut         = cms.double(0.90),
                                       isoConeDeltaR    = cms.double(0.5),
                                       isoConeDeltaZ    = cms.double(0.5) 
                                       )

process.analyzer = cms.Path(process.L1MTDAnalyzer)


# Path and EndPath definitions
process.L1simulation_step = cms.Path(process.SimL1Emulator)
#process.endjob_step = cms.EndPath(process.out*process.endOfProcess)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)


# Schedule definition
process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.MTD_step,process.L1PFTaus,process.analyzer)
#process.schedule = cms.Schedule(process.L1simulation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Customisation from command line

from L1Trigger.Configuration.customiseUtils import L1TrackTriggerTracklet
process = L1TrackTriggerTracklet(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

