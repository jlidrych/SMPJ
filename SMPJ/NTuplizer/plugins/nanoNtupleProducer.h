#ifndef nanoNtupleProducer_h
#define nanoNtupleProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"//add
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"//add

#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include <algorithm>

using namespace reco;

class nanoNtupleProducer : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit nanoNtupleProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
    virtual void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~nanoNtupleProducer();


  private:  
    virtual bool isGoodJet(const pat::Jet &jet);

    void initialize();
    //---- configurable parameters --------  
    edm::EDGetTokenT<pat::JetCollection> jetsToken;
    edm::EDGetTokenT<GenJetCollection> genjetsToken;
    edm::EDGetTokenT<pat::METCollection> metToken;
    edm::EDGetTokenT<pat::PackedCandidateCollection> candsToken;
    edm::EDGetTokenT<double> rhoToken;
    edm::EDGetTokenT<reco::VertexCollection> recVtxsToken;
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
    edm::EDGetTokenT<edm::TriggerResults> metResultsToken;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken;
    edm::EDGetTokenT<edm::View<PileupSummaryInfo> > pupInfoToken;
    edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken;
    edm::EDGetTokenT<LHEEventProduct> lheEvtInfoToken;
    edm::EDGetTokenT<LHERunInfoProduct> runInfoToken;
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
 
    std::string srcBtag_;
    std::vector<std::string> triggerNames_;
    std::vector<std::string> metNames_;
    double etaMax_,ptMin_,ptMinLeading_,massMin_,btagMin_,minMuPt_,minElPt_,GenetaMax_,GenptMin_;
    bool   isMC_,isPrint_,saveWeights_,debug_,isAK8_; 
    //---------------------------------------------------
    edm::Service<TFileService> fs_;
    TTree *outTree_; 
    TH1F *cutFlowHisto_;
    //---- TRIGGER -------------------------  
    TH1F *triggerPassHisto_,*triggerNamesHisto_;
    //---- output TREE variables ------
    //---- global event variables -----
    int   run_,evt_,nVtx_,lumi_,nJets_,nBJets_,nLeptons_,nGenJets_;
    float rho_,met_,metSig_,ht_,pvRho_,pvz_,pvndof_,pvchi2_,metGenSig_,metGen_;

    int nTriggerObjects_;
    std::vector<float> *TrigObjpt_;
    std::vector<float> *TrigObjeta_;
    std::vector<float> *TrigObjphi_;

    std::vector<string> *TrigObjcollection_;

    std::vector<bool> *triggerBit_;
    std::vector<bool> *metBit_;
    std::vector<int>  *triggerPre_;
    //---- jet variables --------------
    std::vector<bool>  *isBtag_;
    std::vector<int>   *flavor_,*nSubJets_,*nSubGenJets_,*nBSubJets_,*flavorHadron_;
    std::vector<float> *cor_, *unc_,*pt_,*eta_,*phi_,*mass_,*massSoftDrop_,*energy_,*chf_,*nhf_,*phf_,*elf_,*muf_,*btag_,*tau1_,*tau2_,*tau3_, *massPruned_;
    std::vector<float> *btagSub0_,*btagSub1_,*massSub0_,*massSub1_,*ptSub0_,*ptSub1_,*etaSub0_,*etaSub1_,*phiSub0_,*phiSub1_;
    std::vector<float> *PuppiCHSpT_;
    std::vector<int>   *flavorSub0_,*flavorSub1_,*flavorHadronSub0_,*flavorHadronSub1_;
    //---- MC variables ---------------
    int npu_,decay_;
    float genEvtWeight_,lheOriginalXWGTUP_;
    std::vector<float> *scaleWeights_;
    std::vector<float> *pdfWeights_;

    
    //candidate information
    std::vector<float> *candPdgId_;
    std::vector<int> *candJetID_;
    std::vector<float> *candPhi_,*candEta_,*candPt_,*candMass_,*candPvAssociationQuality_,*candDXY_,*candDZ_,*candDZAssociatedPV_,*candPuppiWeight_;
    std::vector<int> *candFatherSubJet_;


    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;

    //gen jets
    std::vector<float> *GenJetpt_;
    std::vector<float> *GenJetphi_;
    std::vector<float> *GenJeteta_;
    std::vector<float> *GenJetenergy_;
    std::vector<float> *GenJetmass_;
    std::vector<bool> *isBJetGen_;

    edm::Handle<pat::JetCollection> jets;
    edm::Handle<GenJetCollection> genjets;
    edm::Handle<pat::METCollection> met;
    edm::Handle<pat::PackedCandidateCollection> cands;
    edm::Handle<double> rho;
    edm::Handle<reco::VertexCollection> recVtxs;
    edm::Handle<edm::TriggerResults> triggerResults;
    edm::Handle<edm::TriggerResults> metResults;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<edm::View<PileupSummaryInfo> > pupInfo;
    edm::Handle<edm::View<reco::GenParticle> > genParticles;
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    edm::Handle<LHEEventProduct> lheEvtInfo;
    edm::Handle<LHERunInfoProduct> runInfo;

    JetCorrectionUncertainty *mPFUncCHS;

};





#endif
