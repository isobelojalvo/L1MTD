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
    input = cms.untracked.int32(40000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/00157B11-405C-E811-89CA-0CC47AFB81B4.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/003906EE-A95C-E811-9573-0025905C96A4.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/00AA00A5-645C-E811-8D0E-0025904C6564.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/00C881AB-EC5B-E811-95FB-0025905C3E38.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/00F464AE-025C-E811-BD82-0CC47AF9B2FE.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0234A560-735C-E811-AAF6-0025904C6214.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/02C1440B-FC5B-E811-BD1A-0025904CF75A.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/04306E5F-BD5B-E811-B29A-0CC47AF9B2E6.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/048260EE-2F5C-E811-BC93-0CC47AFB7D5C.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/048F32FC-B05C-E811-822F-0025904CF766.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/04DF9918-B05C-E811-B1BB-0025905C96E8.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/04EF8249-A95C-E811-832C-0025905C3DD6.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/060A7682-F75B-E811-8FCB-0CC47AF9B2EA.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0619A992-B95B-E811-8B12-0025904C51DA.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/067ABC7A-435C-E811-A99A-0025905C96A6.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/069EDA2A-9C5C-E811-B9EF-0025905C3D6A.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/06A6A880-AC5C-E811-8618-0025904C66A2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/06A8BCA5-685B-E811-BBF3-0025905C3D40.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/06C36403-605C-E811-AB7A-0025905C54FC.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/06CB4DED-895C-E811-83B2-0025905C3D40.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/06FE488B-315C-E811-95F4-0CC47AF9B2E6.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/083F71F8-2C5B-E811-99A9-0025905C5486.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/08844A57-9A5C-E811-908D-0025905C53D2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/08C91432-F75B-E811-8771-0CC47AF9B2CA.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0A039FF7-885C-E811-8508-0025905C5430.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0A0444C1-F15B-E811-8572-0025904C6566.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0A10EB5B-B25C-E811-B84B-0CC47AF9B2D2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0A356C7A-435C-E811-8BA4-0025905C94D2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0A3C1F22-2D5C-E811-86F0-0CC47AFB7D90.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0A63F333-845C-E811-99B2-0CC47AFB81BC.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0C2C8C83-765C-E811-983C-0025905D1CB4.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0C3D1B37-785C-E811-80AB-0025905C54C4.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0CFCDED3-EF5B-E811-B468-0025905C543A.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0CFEBE60-865B-E811-BA08-0025904C51DA.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/0E37E046-5D5B-E811-BAC9-0CC47AF9B2D2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/102D0B9C-2F5C-E811-8676-0025905C5500.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/1031849D-385C-E811-98E1-0025905C94D2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/10501ED7-A95C-E811-8D81-0CC47AF9B2E6.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/108E92DE-725C-E811-B4DA-0CC47AFB7DA8.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/10983F4E-795C-E811-B3FE-0025905C4300.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/10FC5070-915C-E811-B914-0025904C6566.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/12021EB3-8D5C-E811-813A-0025905C5488.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/125164D1-765B-E811-B9D6-0025905C2CE8.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/12585B8B-AA5B-E811-9D4D-0025904C5DE2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/12624FF2-B25C-E811-BF7D-0CC47AFB7D5C.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/12D3C2E2-AB5C-E811-93E9-0025905C53D8.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/12E0FCBE-385C-E811-9B89-0025904C66A0.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/146D5771-6A5C-E811-9A78-0025904C66A0.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/147EC2B6-2B5C-E811-AC83-0025905C2CD2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/1497C865-2B5C-E811-BD21-0025905C53D2.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/149CB995-785C-E811-8808-0025905C431A.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/14B72250-995C-E811-8AB0-0025905D1D50.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/14B8E52E-FC5B-E811-98E2-0025905C53F0.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/16003A48-A85C-E811-902E-0CC47AFB8104.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/16973893-3B5C-E811-A3BA-0025905C53D8.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/169CE0E6-B25C-E811-9178-0CC47AF9B2B6.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/16EAF7D2-FD5B-E811-8428-0025905D1E0A.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/16FF82F4-7C5C-E811-B422-0CC47AFB7CEC.root',
        '/store/mc/PhaseIIFall17D/SingleNeutrino/GEN-SIM-DIGI-RAW/L1TPU200_93X_upgrade2023_realistic_v5-v1/80000/1808B0B0-6A5C-E811-9CFF-0025904C66E4.root'
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



# z = 1
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
                                           isoConeDeltaZ    = cms.double(1.0),
                                           WP90ele          = cms.double(0.080),
                                           WPtime90ele      = cms.double(0.020),
                                           WP95ele          = cms.double(0.120),
                                           WPtime95ele      = cms.double(0.065),
                                           WP90mu           = cms.double(0.080),
                                           WPtime90mu       = cms.double(0.020),
                                           WP95mu           = cms.double(0.120),
                                           WPtime95mu       = cms.double(0.08)
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
                                           isoConeDeltaZ    = cms.double(1.0),
                                           WP90ele          = cms.double(0.080),
                                           WPtime90ele      = cms.double(0.020),
                                           WP95ele          = cms.double(0.120),
                                           WPtime95ele      = cms.double(0.065),
                                           WP90mu           = cms.double(0.080),
                                           WPtime90mu       = cms.double(0.020),
                                           WP95mu           = cms.double(0.120),
                                           WPtime95mu       = cms.double(0.065)
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
                                           isoConeDeltaZ    = cms.double(1.0),
                                           WP90ele          = cms.double(0.080),
                                           WPtime90ele      = cms.double(0.065),
                                           WP95ele          = cms.double(0.120),
                                           WPtime95ele      = cms.double(0.090),
                                           WP90mu           = cms.double(0.080),
                                           WPtime90mu       = cms.double(0.065),
                                           WP95mu           = cms.double(0.120),
                                           WPtime95mu       = cms.double(0.090)
                                           )

## 1sigma
process.L1MTDAnalyzer50ps1sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                           timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer50ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                           L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                           l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                           l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                           l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           time_cut         = cms.double(0.050),
                                           isoConeDeltaR    = cms.double(0.4),
                                           isoConeDeltaZ    = cms.double(1.0),
                                           WP90ele          = cms.double(0.080),
                                           WPtime90ele      = cms.double(0.010),
                                           WP95ele          = cms.double(0.120),
                                           WPtime95ele      = cms.double(0.050),
                                           WP90mu           = cms.double(0.080),
                                           WPtime90mu       = cms.double(0.010),
                                           WP95mu           = cms.double(0.120),
                                           WPtime95mu       = cms.double(0.050)
                                           )
##30ps
process.L1MTDAnalyzer30ps1sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                           timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer30ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                           L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                           l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                           l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                           l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           time_cut         = cms.double(0.030),
                                           isoConeDeltaR    = cms.double(0.4),
                                           isoConeDeltaZ    = cms.double(1.0),
                                           WP90ele          = cms.double(0.080),
                                           WPtime90ele      = cms.double(0.010),
                                           WP95ele          = cms.double(0.120),
                                           WPtime95ele      = cms.double(0.065),
                                           WP90mu           = cms.double(0.080),
                                           WPtime90mu       = cms.double(0.010),
                                           WP95mu           = cms.double(0.120),
                                           WPtime95mu       = cms.double(0.065)
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
                                           isoConeDeltaZ    = cms.double(1.0),
                                           WP90ele          = cms.double(0.080),
                                           WPtime90ele      = cms.double(0.010),
                                           WP95ele          = cms.double(0.120),
                                           WPtime95ele      = cms.double(0.055),
                                           WP90mu           = cms.double(0.080),
                                           WPtime90mu       = cms.double(0.010),
                                           WP95mu           = cms.double(0.120),
                                           WPtime95mu       = cms.double(0.055)
                                           )



## 2sigma
process.L1MTDAnalyzer50ps2sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                           timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer50ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                           L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                           l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                           l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                           l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           time_cut         = cms.double(0.100),
                                           isoConeDeltaR    = cms.double(0.4),
                                           isoConeDeltaZ    = cms.double(1.0),
                                           WP90ele          = cms.double(0.080),
                                           WPtime90ele      = cms.double(0.010),
                                           WP95ele          = cms.double(0.120),
                                           WPtime95ele      = cms.double(0.050),
                                           WP90mu           = cms.double(0.080),
                                           WPtime90mu       = cms.double(0.010),
                                           WP95mu           = cms.double(0.120),
                                           WPtime95mu       = cms.double(0.050)
                                           )
##30ps
process.L1MTDAnalyzer30ps2sigma = cms.EDAnalyzer('mtdIsoAnalyzerRates',
                                           timingValuesNominal = cms.InputTag("ttTrackTimeValueMapProducer30ps","TTTracksFromTrackletL1ConfigurableFlatResolutionModel"),
                                           L1TrackInputTag  = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),
                                           l1PFCands        = cms.InputTag("l1pfProducer","PF"),
                                           l1PFTaus         = cms.InputTag("L1PFTauProducer","L1PFTaus"),
                                           l1PFMET          = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           l1PFMETTime      = cms.InputTag("L1TrackerEtMiss","MET","L1"),
                                           time_cut         = cms.double(0.060),
                                           isoConeDeltaR    = cms.double(0.4),
                                           isoConeDeltaZ    = cms.double(1.0),
                                           WP90ele          = cms.double(0.080),
                                           WPtime90ele      = cms.double(0.010),
                                           WP95ele          = cms.double(0.120),
                                           WPtime95ele      = cms.double(0.060),
                                           WP90mu           = cms.double(0.080),
                                           WPtime90mu       = cms.double(0.010),
                                           WP95mu           = cms.double(0.120),
                                           WPtime95mu       = cms.double(0.060)
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
                                           isoConeDeltaZ    = cms.double(1.0),
                                           WP90ele          = cms.double(0.080),
                                           WPtime90ele      = cms.double(0.010),
                                           WP95ele          = cms.double(0.120),
                                           WPtime95ele      = cms.double(0.050),
                                           WP90mu           = cms.double(0.080),
                                           WPtime90mu       = cms.double(0.010),
                                           WP95mu           = cms.double(0.120),
                                           WPtime95mu       = cms.double(0.050)
                                           )
############# z = 2
process.L1MTDAnalyzer50psZ2 = process.L1MTDAnalyzer50ps.clone()
process.L1MTDAnalyzer50psZ2.isoConeDeltaZ = cms.double(2.0)

############# z = 5


############ z = 10


process.analyzer = cms.Path(process.L1MTDAnalyzer50ps
                            +process.L1MTDAnalyzer30ps
                            +process.L1MTDAnalyzer100ps
                            +process.L1MTDAnalyzer50ps1sigma
                            +process.L1MTDAnalyzer30ps1sigma
                            +process.L1MTDAnalyzer100ps1sigma
                            +process.L1MTDAnalyzer50ps2sigma
                            +process.L1MTDAnalyzer30ps2sigma
                            +process.L1MTDAnalyzer100ps2sigma)


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
