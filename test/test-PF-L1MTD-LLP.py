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
    input = cms.untracked.int32(1000)
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
                                       isoConeDeltaZ    = cms.double(1.0) 
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

##delete me 
'''
process.source.secondaryFileNames = cms.untracked.vstring(
    "/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/110000/F9075F34-D392-6A43-AB9B-8475DF216C23.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/110000/E6713AE6-896D-744E-87B2-37842139EE74.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/110000/E00A3494-F1B1-DF40-886F-D1FFB45497D4.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/110000/48A2D642-EF40-EC47-97D7-3745DA496237.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/110000/3C633436-B60E-0245-9236-ADC04CE2C57B.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/110000/2CA0DA87-EC05-EE4E-AF8C-4E4DAFDF9B70.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/110000/20A2B7FA-D042-2D4D-8896-E0B019182520.root",
    "/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/110000/0DE10720-7C0F-324A-9684-8B7B89187380.root"
)
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:177")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

# Input source
process.source.fileNames = cms.untracked.vstring("file:/hdfs/store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm_TuneCP5_14TeV_pythia8/MINIAODSIM/PU200_103X_upgrade2023_realistic_v2-v1/110000/F3E6C49B-16EB-F240-9964-117F923BA060.root")


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DisplacedSUSY_stopToChi_Gravitino_M_1000_700_100mm.root")
)
'''

'''
process.source.secondaryFileNames = cms.untracked.vstring(
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/E788F43D-E0C0-3D43-B0DB-537663F5AB25.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/CCB975E4-2E01-A740-B684-CBC51230A4B5.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/B9B94B69-8E74-6543-980B-EFDEE94C7A98.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/930FFE41-A18B-E846-8E9E-A8DE837A7B32.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/6A43E37E-F02B-2048-8F85-12D4BC5AF38D.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/5E136D57-727A-9547-B998-09E42C02C062.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/58AEFE7B-B6B6-0344-B00D-94F3C493FE09.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/3C4536FA-F9A2-E445-AA18-84ACE50362EE.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/0CAEE7A1-F85A-5441-823C-4C96DF87FD25.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/09E54504-1F9E-7E43-904D-0ADE944A3174.root',
    '/store/mc/PhaseIIMTDTDRAutumn18DR/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/70000/06313F21-AFF8-7248-B06B-76D311DD0F0F.root'
    )
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:177")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

# Input source
process.source.fileNames = cms.untracked.vstring("file:/hdfs/store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm_TuneCP5_14TeV_pythia8/MINIAODSIM/PU200_103X_upgrade2023_realistic_v2-v1/70000/017EB303-130E-5B48-8B48-AB6B763602B6.root")


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DisplacedSUSY_stopToChi_Gravitino_M_1000_700_10mm.root")
)
'''

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:178")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.source.fileNames = cms.untracked.vstring("file:/hdfs/store/mc/PhaseIIMTDTDRAutumn18MiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8/MINIAODSIM/PU200_103X_upgrade2023_realistic_v2-v1/80000/9EDEBB17-1D14-8B4F-B506-0B3EA2D219AC.root")
process.source.secondaryFileNames = cms.untracked.vstring("root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/GluGluHToTauTau_M125_14TeV_powheg_pythia8/FEVT/PU200_103X_upgrade2023_realistic_v2-v1/80000/C69A0FB2-9A88-3B4C-916B-A943865A8856.root")
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("GluGluHToTauTau_M125_14TeV-200PU.root")
)
