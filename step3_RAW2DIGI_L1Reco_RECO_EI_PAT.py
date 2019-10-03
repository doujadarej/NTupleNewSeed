# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions auto:phase1_2017_realistic -n 10 --era Run2_2017 --eventcontent RECOSIM --runUnscheduled -s RAW2DIGI,L1Reco,RECO,EI,PAT --datatier GEN-SIM-RECO --geometry DB:Extended --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('Configuration.StandardSequences.PATMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#####################ADDED 

from RecoLocalTracker.Configuration.RecoLocalTracker_cff import *
from SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi import *
from SimTracker.TrackerHitAssociation.tpClusterProducer_cfi import *
from SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi import *

from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEGeneric_cfi import *
from SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi import *

# Track Associators
from SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi import *

#process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v7'  


### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

### validation-specific includes
process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("Validation.RecoTrack.MultiTrackValidatorGenPs_cfi")
process.load("DQMServices.Components.EDMtoMEConverter_cff")
process.load("Validation.Configuration.postValidation_cff")
process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')

process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByChi2_cfi")
process.trackAssociatorByChi2.chi2cut = cms.double(500.0)


process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.TFileService = cms.Service("TFileService", fileName = cms.string("trackingNTuple.root") )

process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

from CommonTools.RecoAlgos.trackingParticleRefSelector_cfi import trackingParticleRefSelector as _trackingParticleRefSelector
process.trackingParticlesIntime = _trackingParticleRefSelector.clone(
    signalOnly = False,
    intimeOnly = False,
    chargedOnly = True,
    tip = 1e5,
    lip = 1e5,
    minRapidity = -2.4,
    maxRapidity =  2.4,
    ptMin = 1,
)

process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.MixingModule.trackingTruthProducerSelection_cfi")

process.trackingParticles.simHitCollections = cms.PSet( )
process.mix.playback = cms.untracked.bool(True)
process.mix.digitizers = cms.PSet(
     mergedtruth = cms.PSet(process.trackingParticles)
)

############################################# 




process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/opt/sbg/data/safe1/cms/jandrea/TrackingSeedDefinition/10024.0_TTbar_13+TTbar_13TeV_TuneCUETP8M1_2017_GenSimFull+DigiFull_2017+RecoFull_2017+ALCAFull_2017+HARVESTFull_2017/step2.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:step3.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2017_realistic', '')





############# CHANGED PLACE OF THIS  

#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('CommonTools.ParticleFlow.EITopPAG_cff')
#process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
#process.load('Configuration.StandardSequences.PATMC_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff').load('Configuration.StandardSequences.L1Reco_cff')


############ ADDED THIS 

process.trackingPerf = cms.EDAnalyzer('TrackingPerf',
     tracks               = cms.untracked.InputTag('generalTracks'),
     trackLabel           = cms.InputTag('generalTracks'),
      #trackAssociator      = cms.untracked.InputTag('trackingParticleRecoTrackAsssociation'),
      #trackAssociator      = cms.untracked.InputTag('trackAssociatorByHits'),
      trackAssociator      = cms.untracked.InputTag('quickTrackAssociatorByHits'),
      beamSpot             = cms.untracked.InputTag('offlineBeamSpot'),
      parametersDefiner    = cms.untracked.string('LhcParametersDefinerForTP'),
      #trackingParticles    = cms.untracked.InputTag('trackingParticlesIntime'),
      trackingParticles    = cms.untracked.InputTag('trackingParticlesIntime'),
      trackingParticlesRef = cms.untracked.bool(True),
      TTRHBuilder          = cms.string('WithTrackAngle'),
      useCluster           = cms.untracked.bool(False),
      vertices             = cms.untracked.InputTag('offlinePrimaryVertices'),
      jetInput             = cms.InputTag('ak4PFJets'), 
      gsftrack             = cms.InputTag('electronGsfTracks'),

)
############### 

# Path and EndPath definitions
#process.trackingPerformance = cms.Path(process.trackingPerf)
process.p = cms.Path( 
	#process.LHCTransport*
	#process.g4SimHits*
#	process.mix *
	#process.simHitTPAssocProducer*
	process.tpClusterProducer*
	process.quickTrackAssociatorByHits*
	process.trackingParticleRecoTrackAsssociation*
	process.trackingParticlesIntime*
	process.trackingPerf)

process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.p,process.endjob_step,process.RECOSIMoutput_step)

process.schedule.associate(process.patTask)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

#call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
