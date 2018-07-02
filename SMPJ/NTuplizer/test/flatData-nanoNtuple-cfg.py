import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTreeFileData-nano.root'))
process.load("Configuration.EventContent.EventContent_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.globaltag = '94X_dataRun2_v6'
#process.GlobalTag.globaltag = '94X_mc2017_realistic_v13'

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#-1
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
        '/store/data/Run2017B/JetHT/MINIAOD/PromptReco-v2/000/299/180/00000/4A96909D-F66B-E711-AFDB-02163E0144F4.root',
        ),
)
#############   Format MessageLogger #################
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

from RecoJets.JetProducers.PFJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

ak4PFJetsPuppi = ak4PFJets.clone(
    src = cms.InputTag("puppi")
)

process.load('CommonTools.PileupAlgos.Puppi_cff')

process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag('slimmedJets')
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets, ak8GenJets

genParticleCollection = 'prunedGenParticles'

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from RecoJets.JetProducers.GenJetParameters_cfi import *

process.ak4GenJetsCustom = ak4GenJets.clone(
    src = genParticleCollection,
    rParam = cms.double(0.4),
    jetAlgorithm = cms.string("AntiKt")
)

process.ak8GenJetsCustom = ak4GenJets.clone(
    src = genParticleCollection,
    rParam = cms.double(0.8),
    jetAlgorithm = cms.string("AntiKt")
)

genParticleCollection = 'prunedGenParticles'
genJetCollection = 'slimmedGenJets'

from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
    particles = genParticleCollection
)

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.ak4genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets = genJetCollection,
    rParam = cms.double(0.4),
)

genJetCollection = 'slimmedGenJetsAK8'

process.ak8genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets = genJetCollection,
    rParam = cms.double(0.8),
)

#only needed if information of the associated b hadrons are required
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
    genParticles = genParticleCollection,
    jetFlavourInfos = cms.InputTag("ak4genJetFlavourInfos"),
    flavour = cms.int32(5),
    onlyJetClusteredHadrons = cms.bool(True),
    noBBbarResonances = cms.bool(False),
)


# ---------------------------------------------------------
# set up TransientTrackBuilder
process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName=cms.string('TransientTrackBuilder')
)

##-------------------- User analyzer  --------------------------------
process.boostedAK4 = cms.EDAnalyzer('nanoNtupleProducer',
  jets             = cms.untracked.InputTag('slimmedJets'),
  met              = cms.InputTag('slimmedMETs'),
  vertices         = cms.InputTag('offlineSlimmedPrimaryVertices'),
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  candidates       = cms.InputTag('packedPFCandidates'),
  triggerPrescales = cms.InputTag('patTrigger'),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  nJetsMin         = cms.int32(6),
  nBJetsMin        = cms.int32(2),
  ptMin            = cms.double(50), 
  minMuPt          = cms.double(30),                              
  minElPt          = cms.double(30),                              
  ptMinLeading     = cms.double(50),   
  massMin          = cms.double(50),
  htMin            = cms.double(5),
  etaMax           = cms.double(5.0),
  kinfit           = cms.string('kinFitTtFullHadEvent'),
  btagMinThreshold = cms.double(0.814),
  btagMin          = cms.double(0.1),                     
  btagMaxThreshold = cms.double(1.1),
  btaggerSubjet   = cms.string('pfDeepCSVJetTags:probb'),
  btagger          = cms.string('pfDeepCSVJetTags:probbb'),
#  btagger          = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
  qgtagger         = cms.InputTag('QGTagger','qgLikelihood'),
  pu               = cms.untracked.string("addPileupInfo"),
  genparticles     = cms.untracked.InputTag('prunedGenParticles'),  
  triggerNames     = cms.vstring('HLT_PFJet40_v1','HLT_PFJet60_v1', 'HLT_PFJet80_v1', 'HLT_PFJet140_v1', 'HLT_PFJet200_v1', 'HLT_PFJet260_v1','HLT_PFJet320_v1', 'HLT_PFJet400_v1', 'HLT_PFJet450_v1','HLT_PFJet500_v1','HLT_PFJet550_v1','HLT_PFJet40_v','HLT_PFJet60_v', 'HLT_PFJet80_v', 'HLT_PFJet140_v', 'HLT_PFJet200_v', 'HLT_PFJet260_v','HLT_PFJet320_v', 'HLT_PFJet400_v', 'HLT_PFJet450_v','HLT_PFJet500_v','HLT_PFJet550_v','HLT_AK8PFJet40_v','HLT_AK8PFJet60_v', 'HLT_AK8PFJet80_v', 'HLT_AK8PFJet140_v', 'HLT_AK8PFJet200_v', 'HLT_AK8PFJet260_v','HLT_AK8PFJet320_v', 'HLT_AK8PFJet400_v', 'HLT_AK8PFJet450_v','HLT_AK8PFJet500_v','HLT_AK8PFJet550_v','HLT_AK8PFJet40_v1','HLT_AK8PFJet60_v1', 'HLT_AK8PFJet80_v1', 'HLT_AK8PFJet140_v1', 'HLT_AK8PFJet200_v1', 'HLT_AK8PFJet260_v1','HLT_AK8PFJet320_v1', 'HLT_AK8PFJet400_v1', 'HLT_AK8PFJet450_v1','HLT_AK8PFJet500_v1','HLT_AK8PFJet550_v1','HLT_AK8PFJetFwd40_v1','HLT_AK8PFJetFwd60_v1', 'HLT_AK8PFJetFwd80_v1', 'HLT_AK8PFJetFwd140_v1', 'HLT_AK8PFJetFwd200_v1', 'HLT_AK8PFJetFwd260_v1','HLT_AK8PFJetFwd320_v1', 'HLT_AK8PFJetFwd400_v1', 'HLT_AK8PFJetFwd450_v1','HLT_AK8PFJetFwd500_v1','HLT_AK8PFJetFwd40_v1','HLT_AK8PFJetFwd60_v1', 'HLT_AK8PFJetFwd80_v1', 'HLT_AK8PFJetFwd140_v1', 'HLT_AK8PFJetFwd200_v1', 'HLT_AK8PFJetFwd260_v1','HLT_AK8PFJetFwd320_v1', 'HLT_AK8PFJetFwd400_v1', 'HLT_AK8PFJetFwd450_v1','HLT_AK8PFJetFwd500_v1','HLT_PFJetFwd40_v1','HLT_PFJetFwd60_v1', 'HLT_PFJetFwd80_v1', 'HLT_PFJetFwd140_v1', 'HLT_PFJetFwd200_v1', 'HLT_PFJetFwd260_v1','HLT_PFJetFwd320_v1', 'HLT_PFJetFwd400_v1', 'HLT_PFJetFwd450_v1','HLT_PFJetFwd500_v1','HLT_PFJetFwd40_v','HLT_PFJetFwd60_v', 'HLT_PFJetFwd80_v', 'HLT_PFJetFwd140_v', 'HLT_PFJetFwd200_v', 'HLT_PFJetFwd260_v','HLT_PFJetFwd320_v', 'HLT_PFJetFwd400_v', 'HLT_PFJetFwd450_v','HLT_PFJetFwd500_v'),
  triggerResults   = cms.InputTag('TriggerResults','','HLT'),
  metNames         = cms.vstring('Flag_globalTightHalo2016Filter','Flag_goodVertices','Flag_HBHENoiseFilter','Flag_HBHENoiseIsoFilter','Flag_EcalDeadCellTriggerPrimitiveFilter','Flag_BadChargedCandidateFilter','Flag_BadPFMuonFilter','Flag_eeBadScFilter'),
  metResults       = cms.InputTag('TriggerResults','','RECO'),
  isMC             = cms.untracked.bool(False),
  isAK8            = cms.untracked.bool(False),
  genjets          = cms.untracked.InputTag('slimmedGenJets'),
  GenptMin         = cms.untracked.double(50),
  GenetaMax        = cms.untracked.double(2.5),
  jetFlavourInfos  = cms.InputTag("ak4genJetFlavourInfos"),                 
  isPrint          = cms.untracked.bool(False),
  saveWeights      = cms.untracked.bool(False),
)


process.boostedAK8 = process.boostedAK4.clone(
    jets             = cms.untracked.InputTag('slimmedJetsAK8'),
    genjets          = cms.untracked.InputTag('slimmedGenJetsAK8'),
    jetFlavourInfos = cms.InputTag('ak8genJetFlavourInfos'),
#    jetFlavourInfos = cms.InputTag('slimmedGenJetsFlavourInfos'),
    isAK8           = cms.untracked.bool(True),
    isMC             = cms.untracked.bool(False),
)


process.p = cms.Path(
   process.boostedAK4 *
   process.boostedAK8
)








