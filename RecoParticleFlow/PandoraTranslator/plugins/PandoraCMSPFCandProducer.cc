#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "PandoraCMSPFCandProducer.h"
#include "RecoParticleFlow/PandoraTranslator/interface/CMSBFieldPlugin.h"
#include "RecoParticleFlow/PandoraTranslator/interface/CMSPseudoLayerPlugin.h"
#include "LCContent.h"
#include "LCContentFast.h"
#include "RecoParticleFlow/PandoraTranslator/interface/CMSAlgorithms.h"
#include "RecoParticleFlow/PandoraTranslator/interface/CMSTemplateAlgorithm.h"
#include "RecoParticleFlow/PandoraTranslator/interface/CMSGlobalHadronCompensationPlugin.h"
#include "RecoParticleFlow/PandoraTranslator/interface/MuonCoilCorrection.h"

//#include "PandoraMonitoringApi.h"

#include "Managers/PluginManager.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"

#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

// Addition for HGC geometry
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/FlatTrd.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCEEDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCHEDetId.h"
#include "Geometry/FCalGeometry/interface/HGCalGeometry.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// for track propagation through HGC  
// N.B. we are only propogating to first layer, so check these later
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"

//We need the speed of light
#include "CLHEP/Units/PhysicalConstants.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/RootAutoLibraryLoader/interface/RootAutoLibraryLoader.h"

#include "TGClient.h"
#include "TVirtualX.h"
#include "TROOT.h"
#include "TRint.h"
#include "TGraphErrors.h"
 
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <iterator>
#include <unordered_set>
#include <numeric>

using namespace edm;
using namespace reco;
using namespace pandora;  
using namespace lc_content;
using namespace cms_content;


// SCZ utility stuffs
typedef std::pair<float, float> HitEnergyDistance;
typedef std::vector<HitEnergyDistance> HitEnergyDistanceVector;

bool SortHitsByDistance(const HitEnergyDistance &lhs, const HitEnergyDistance &rhs)
{
  return (lhs.second < rhs.second);
}

template<typename A, typename B>
std::pair<B,A> flip_pair(const std::pair<A,B> &p)
{
  return std::pair<B,A>(p.second, p.first);
}

template<typename A, typename B>
std::multimap<B,A> flip_map(const std::map<A,B> &src)
{
  std::multimap<B,A> dst;
  std::transform(src.begin(), src.end(), std::inserter(dst, dst.begin()),
		 flip_pair<A,B>);
  return dst;
}

namespace cms_content {
  pandora::StatusCode RegisterBasicPlugins(const pandora::Pandora &pandora)
  {
    LC_ENERGY_CORRECTION_LIST(PANDORA_REGISTER_ENERGY_CORRECTION);
    LC_PARTICLE_ID_LIST(PANDORA_REGISTER_PARTICLE_ID);

    PANDORA_REGISTER_ENERGY_CORRECTION("CMSGlobalHadronCompensation",          pandora::HADRONIC,     cms_content::GlobalHadronCompensation);
    PANDORA_REGISTER_ENERGY_CORRECTION("CMSMuonCoilCorrection",                pandora::HADRONIC,     cms_content::MuonCoilCorrection);

    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(pandora, new cms_content::CMSPseudoLayerPlugin));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetShowerProfilePlugin(pandora, new lc_content::LCShowerProfilePlugin));

    return pandora::STATUS_CODE_SUCCESS;
  }
}

//pandora::Pandora * PandoraCMSPFCandProducer::m_pPandora = NULL;
//
// constructors and destructor
//
PandoraCMSPFCandProducer::PandoraCMSPFCandProducer(const edm::ParameterSet& iConfig) : 
  debugPrint(iConfig.getParameter<bool>("debugPrint")), debugHisto(iConfig.getParameter<bool>("debugHisto")), useRecoTrackAsssociation(iConfig.getParameter<bool>("useRecoTrackAsssociation")),
  m_calibEE(ForwardSubdetector::HGCEE,"EE",debugPrint), m_calibHEF(ForwardSubdetector::HGCHEF,"HEF",debugPrint), m_calibHEB(ForwardSubdetector::HGCHEB,"HEB",debugPrint), calibInitialized(false),
  _deltaPtOverPtForPfo(iConfig.getParameter<double>("MaxDeltaPtOverPtForPfo")),
  _deltaPtOverPtForClusterlessPfo(iConfig.getParameter<double>("MaxDeltaPtOverPtForClusterlessPfo"))
{  
  produces<reco::PFClusterCollection>();
  produces<reco::PFBlockCollection>();
  produces<reco::PFCandidateCollection>();

  // Produce some of the extra collections from PFProducer
  electronOutputCol_
    = iConfig.getParameter<std::string>("pf_electron_output_col");
  produces<reco::PFCandidateCollection>(electronOutputCol_);
  //    produces<reco::PFCandidateElectronExtraCollection>(electronExtraOutputCol_);
  //    produces<reco::PFCandidatePhotonExtraCollection>(photonExtraOutputCol_);
  
  inputTagSimTracks_ = iConfig.getUntrackedParameter<InputTag>("simTracks",InputTag("g4SimHits"));
  inputTagSimHGCHitsEE_ = iConfig.getUntrackedParameter<InputTag>("SimHGCHitsEE",InputTag("g4SimHits:HGCHitsEE"));
  inputTagSimHGCHitsHEback_ = iConfig.getUntrackedParameter<InputTag>("SimHGCHitsHEback",InputTag("g4SimHits:HGCHitsHEback"));
  inputTagSimHGCHitsHEfront_ = iConfig.getUntrackedParameter<InputTag>("SimHGCHitsHEfront",InputTag("g4SimHits:HGCHitsHEfront"));

  inputTagHGCrechit_ = iConfig.getParameter<InputTag>("HGCrechitCollection");
  inputTagGeneralTracks_ = iConfig.getParameter<InputTag>("generaltracks");
  inputTagtPRecoTrackAsssociation_ = iConfig.getParameter<InputTag>("tPRecoTrackAsssociation");
  inputTagGenParticles_ = iConfig.getParameter<InputTag>("genParticles");
  inputTagSimTracks_ = iConfig.getParameter<InputTag>("simTracks");
  inputTagSimVertices_ = iConfig.getParameter<InputTag>("simVertices");
  inputSimHits_ = iConfig.getParameter<std::vector<InputTag> >("HGCSimHits");
  m_pandoraSettingsXmlFile = iConfig.getParameter<edm::FileInPath>("inputconfigfile");

  m_calibrationParameterFile = iConfig.getParameter<edm::FileInPath>("calibrParFile");
  m_energyCorrMethod = iConfig.getParameter<std::string>("energyCorrMethod");
  m_energyWeightingFilename  = iConfig.getParameter<edm::FileInPath>("energyWeightFile");
  m_layerDepthFilename = iConfig.getParameter<edm::FileInPath>("layerDepthFile");
  m_overburdenDepthFilename = iConfig.getParameter<edm::FileInPath>("overburdenDepthFile");
  m_useOverburdenCorrection = iConfig.getParameter<bool>("useOverburdenCorrection");
  _outputFileName = iConfig.getParameter<std::string>("outputFile");

  stm = new steerManager(m_energyWeightingFilename.fullPath().c_str());
  
  speedoflight = (CLHEP::c_light/CLHEP::cm)/CLHEP::ns;
}

PandoraCMSPFCandProducer::~PandoraCMSPFCandProducer()
{
  delete stm;

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   
  // SCZ NEIN
  /*
  if(debugHisto){
    delete m_hitEperLayer_EM[ForwardSubdetector::HGCEE];
    delete m_hitEperLayer_EM[ForwardSubdetector::HGCHEF];
    delete m_hitEperLayer_EM[ForwardSubdetector::HGCHEB];
    delete m_hitEperLayer_HAD[ForwardSubdetector::HGCEE];
    delete m_hitEperLayer_HAD[ForwardSubdetector::HGCHEF];
    delete m_hitEperLayer_HAD[ForwardSubdetector::HGCHEB]; 
  }
  */
}


//
// member functions
//

// ------------ method called for each event  ------------
void PandoraCMSPFCandProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // SCZ NEIN
  //  if(debugHisto) resetVariables();  

  if(debugPrint) std::cout << "Analyzing events" << std::endl ; 

  // std::cout << "Analyzing events 1 " << std::endl ;
  
  prepareMCParticleAssociations(iEvent);

  prepareTrack(iEvent);
  
  prepareHits(iEvent);

  //preparePandoraMCParticles(iEvent); // this should happen after pandora is loaded with the various reconstructed things in the event


  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,PandoraApi::ProcessEvent(*m_pPandora));
  //  if(debugPrint || debugHisto) preparePFO(iEvent); //not necessary for production // SCZ NEIN
  
  edm::Handle<reco::PFRecHitCollection> HGCRecHitHandle;
  iEvent.getByLabel(inputTagHGCrechit_, HGCRecHitHandle);
  
  edm::Handle<reco::PFRecTrackCollection> tkRefCollection;
  iEvent.getByLabel(inputTagGeneralTracks_, tkRefCollection);

  // now we do all the hard work?
  convertPandoraToCMSSW(tkRefCollection,HGCRecHitHandle,iEvent);

   PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,PandoraApi::Reset(*m_pPandora));

   std::unordered_map<const void*,unsigned int>().swap(recHitMap);
   std::unordered_map<const void*,unsigned int>().swap(recTrackMap);

   // clear out the MC particle maps
   RecoDetIdToSimHitsMap().swap( m_recoDetIdToSimHit );
   SimHitsToSimTrackMap().swap( m_simHitsToSimTracks );
   SimTrackToVertexMap().swap( m_simTrackToStartSimVertex );
   SimTrackToVertexMap().swap( m_simTrackToStartSimVertexEnd );
   SimVertexToTracksMap().swap( m_simVertexToSimTracks );
   SimVertexToTracksMap().swap( m_simTrackToSimVertex );
   SimVertexToTracksMap().swap( m_simVertexToSimTrackParent );
   BarCodeMap().swap(m_barCodesToSimTrack);
   BarCodeMap().swap(m_barCodesToSimVertices);
   BarCodeMap().swap(m_barCodesToSimHits);
   BarCodeMap().swap(m_barCodesToGenParticle);
}

void PandoraCMSPFCandProducer::initPandoraCalibrParameters()
{
   m_firstMCpartEta = -99999.;
   m_firstMCpartPhi = -99999.;
   m_secondMCpartEta = -99999.;
   m_secondMCpartPhi = -99999.;

   m_calibEE.m_Calibr_ADC2GeV  = 0.000012 ; //w/o absorber thickness correction
   m_calibHEF.m_Calibr_ADC2GeV = 0.0000176; //w/o absorber thickness correction
   m_calibHEB.m_Calibr_ADC2GeV = 0.0003108; //w/o absorber thickness correction

   m_calibEE.m_EM_addCalibr   = 10.;
   m_calibHEF.m_EM_addCalibr  = 10.;
   m_calibHEB.m_EM_addCalibr  = 10.;
   m_calibEE.m_HAD_addCalibr  = 10.;
   m_calibHEF.m_HAD_addCalibr = 10.;
   m_calibHEB.m_HAD_addCalibr = 10.;

   m_calibEE.m_CalThresh  = 27.55e-6; //EE
   m_calibHEF.m_CalThresh = 42.50e-6;
   m_calibHEB.m_CalThresh = 742.2e-6;

   m_calibEE.m_CalMipThresh  = 0.5;
   m_calibHEF.m_CalMipThresh = 0.5;
   m_calibHEB.m_CalMipThresh = 0.5;

   m_calibEE.m_CalToMip     = 18149.;
   m_calibHEF.m_CalToMip    = 11765.;
   m_calibHEB.m_CalToMip    = 667.4;

   m_calibEE.m_CalToEMGeV   = 1.;
   m_calibHEF.m_CalToEMGeV  = 1.;
   m_calibHEB.m_CalToEMGeV  = 1.;

   m_calibEE.m_CalToHADGeV  = 1.;
   m_calibHEF.m_CalToHADGeV = 1.;
   m_calibHEB.m_CalToHADGeV = 1.;
   m_muonToMip             = 1.;

   return;
}

void PandoraCMSPFCandProducer::readCalibrParameterFile()
{
  std::ifstream calibrParFile(m_calibrationParameterFile.fullPath().c_str() , std::ifstream::in );

  if (!calibrParFile.is_open()) {
    if(debugPrint) std::cout << "PandoraCMSPFCandProducer::readCalibrParameterFile: calibrParFile does not exist ("
        << m_calibrationParameterFile << ")" << std::endl;
    return;
  }

   while ( !calibrParFile.eof() ) {
      std::string linebuf;
      getline( calibrParFile, linebuf );
      if (linebuf.substr(0,1) == "#") continue;
      if (linebuf.substr(0,2) == "//") continue;

      if (linebuf.empty()) continue;

      std::string paraName;
      double paraValue;
      std::stringstream ss(linebuf);
      ss >> paraName >> paraValue;

           if (paraName=="Calibr_ADC2GeV_EE"     ) {m_calibEE.m_Calibr_ADC2GeV  = paraValue;}
      else if (paraName=="Calibr_ADC2GeV_HEF"    ) {m_calibHEF.m_Calibr_ADC2GeV = paraValue;}
      else if (paraName=="Calibr_ADC2GeV_HEB"    ) {m_calibHEB.m_Calibr_ADC2GeV = paraValue;}

      else if (paraName=="EMaddCalibrEE"         ) {m_calibEE.m_EM_addCalibr    = paraValue;}
      else if (paraName=="EMaddCalibrHEF"        ) {m_calibHEF.m_EM_addCalibr   = paraValue;}
      else if (paraName=="EMaddCalibrHEB"        ) {m_calibHEB.m_EM_addCalibr   = paraValue;}
      else if (paraName=="HADaddCalibrEE"        ) {m_calibEE.m_HAD_addCalibr   = paraValue;}
      else if (paraName=="HADaddCalibrHEF"       ) {m_calibHEF.m_HAD_addCalibr  = paraValue;}
      else if (paraName=="HADaddCalibrHEB"       ) {m_calibHEB.m_HAD_addCalibr  = paraValue;}
 
      else if (paraName=="ECalThresEndCap"       ) {m_calibEE.m_CalThresh       = paraValue;}
      else if (paraName=="HCalThresEndCapHEF"    ) {m_calibHEF.m_CalThresh      = paraValue;}
      else if (paraName=="HCalThresEndCapHEB"    ) {m_calibHEB.m_CalThresh      = paraValue;}
 
      else if (paraName=="ECalMipThresEndCap"    ) {m_calibEE.m_CalMipThresh    = paraValue;}
      else if (paraName=="HCalMipThresEndCapHEF" ) {m_calibHEF.m_CalMipThresh   = paraValue;}
      else if (paraName=="HCalMipThresEndCapHEB" ) {m_calibHEB.m_CalMipThresh   = paraValue;}

      else if (paraName=="ECalToMipEndCap"       ) {m_calibEE.m_CalToMip        = paraValue;}
      else if (paraName=="HCalToMipEndCapHEF"    ) {m_calibHEF.m_CalToMip       = paraValue;}
      else if (paraName=="HCalToMipEndCapHEB"    ) {m_calibHEB.m_CalToMip       = paraValue;}

      else if (paraName=="ECalToEMGeVEndCap"     ) {m_calibEE.m_CalToEMGeV      = paraValue;}
      else if (paraName=="HCalToEMGeVEndCapHEF"  ) {m_calibHEF.m_CalToEMGeV     = paraValue;}
      else if (paraName=="HCalToEMGeVEndCapHEB"  ) {m_calibHEB.m_CalToEMGeV     = paraValue;}

      else if (paraName=="ECalToHadGeVEndCap"    ) {m_calibEE.m_CalToHADGeV     = paraValue;}
      else if (paraName=="HCalToHadGeVEndCapHEF" ) {m_calibHEF.m_CalToHADGeV    = paraValue;}
      else if (paraName=="HCalToHadGeVEndCapHEB" ) {m_calibHEB.m_CalToHADGeV    = paraValue;}

      else if (paraName=="MuonToMip"             ) {m_muonToMip                 = paraValue;}
      else continue;
      
      if(debugPrint) std::cout << "reading calibr parameter " << paraName << " " << paraValue << std::endl;
   }

   calibrParFile.close();

   return;
}


void PandoraCMSPFCandProducer::readEnergyWeight()
{
   if(debugPrint) std::cout << "PandoraCMSPFCandProducer::readEnergyWeight" << std::endl;

   //FIXME : for the moment, everything is given in unit of MIP

   stm->addArrayParameter("layerSet_EE");
   stm->addArrayParameter("layerSet_HEF");
   stm->addArrayParameter("layerSet_HEB");
   stm->addArrayParameter("energyWeight_EM_EE");
   stm->addArrayParameter("energyWeight_EM_HEF");
   stm->addArrayParameter("energyWeight_EM_HEB");
   stm->addArrayParameter("energyWeight_Had_EE");
   stm->addArrayParameter("energyWeight_Had_HEF");
   stm->addArrayParameter("energyWeight_Had_HEB");
   
   stm->addSingleParameter("offSet_EM");
   stm->addSingleParameter("offSet_Had");

   stm->read();
   stm->printPars();

   offSet_EM  = stm->getSinglePara("offSet_EM");
   offSet_Had = stm->getSinglePara("offSet_Had");

   return;
}



void PandoraCMSPFCandProducer::prepareGeometry(){ // function to setup a geometry for pandora
  
  // Get the ecal/hcal barrel, endcap geometry
  const EcalBarrelGeometry* ecalBarrelGeometry = dynamic_cast< const EcalBarrelGeometry* > (geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel));
  const HcalGeometry* hcalBarrelGeometry = dynamic_cast< const HcalGeometry* > (geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalBarrel));
  assert( ecalBarrelGeometry );
  assert( hcalBarrelGeometry );
  
  // Additions
  const HGCalGeometry* HGCEEGeometry = hgceeGeoHandle.product() ; 
  const HGCalGeometry* HGCHEFGeometry = hgchefGeoHandle.product() ; 
  const HGCalGeometry* HGCHEBGeometry = hgchebGeoHandle.product() ; 
  
  std::vector<DetId> ecalBarrelCells = geoHandle->getValidDetIds(DetId::Ecal, EcalBarrel);
  std::vector<DetId> ecalEndcapCells = hgceeGeoHandle->getValidDetIds(DetId::Forward, HGCEE);
  std::vector<DetId> hcalBarrelCells = geoHandle->getValidDetIds(DetId::Hcal, HcalBarrel);
  std::vector<DetId> hcalEndcapCellsFront = hgchefGeoHandle->getValidDetIds(DetId::Forward, HGCHEF);
  std::vector<DetId> hcalEndcapCellsBack  = hgchebGeoHandle->getValidDetIds(DetId::Forward, HGCHEB);

  PandoraApi::Geometry::SubDetector::Parameters *ebParameters = new PandoraApi::Geometry::SubDetector::Parameters();
  PandoraApi::Geometry::SubDetector::Parameters *eeParameters = new PandoraApi::Geometry::SubDetector::Parameters();
  PandoraApi::Geometry::SubDetector::Parameters *hbParameters = new PandoraApi::Geometry::SubDetector::Parameters();
  PandoraApi::Geometry::SubDetector::Parameters *heParameters = new PandoraApi::Geometry::SubDetector::Parameters();

  PandoraApi::Geometry::LayerParameters *ebLayerParameters    = new PandoraApi::Geometry::LayerParameters();
  PandoraApi::Geometry::LayerParameters *hbLayerParameters    = new PandoraApi::Geometry::LayerParameters();

  std::vector<PandoraApi::Geometry::LayerParameters*> hgcEELayerParameters;
  std::vector<PandoraApi::Geometry::LayerParameters*> hgcHEFLayerParameters;
  std::vector<PandoraApi::Geometry::LayerParameters*> hgcHEBLayerParameters;

  //initialize HGC layer calibrations (only happens once)
  if(!calibInitialized){
    calibInitialized = true;
    
    //get # of layers from geometry
    const HGCalGeometry* hgcGeometries[3] = {HGCEEGeometry, HGCHEFGeometry, HGCHEBGeometry};
    unsigned int nHGClayers[3];
    for(unsigned int h = 0; h < 3; h++){
      const HGCalDDDConstants &dddCons = hgcGeometries[h]->topology().dddConstants();
      auto firstLayerIt = dddCons.getFirstTrForm();
      auto lastLayerIt = dddCons.getLastTrForm();
      std::unordered_set<float> uniqueZvals;
      for(auto layerIt=firstLayerIt; layerIt !=lastLayerIt; layerIt++) {
        if(layerIt->h3v.z()>0) uniqueZvals.insert(layerIt->h3v.z());
      }
      nHGClayers[h] = uniqueZvals.size();
    }
    
    //divide # of HGCHEB layers by 2 to account for offset in each sector
    nHGCeeLayers = nHGClayers[0]; nHGChefLayers = nHGClayers[1]; nHGChebLayers = nHGClayers[2]/2;
    if(debugPrint) std::cout << "HGCEE layers = " << nHGCeeLayers << ", HGCHEF layers = " << nHGChefLayers << ", HGCHEB layers = " << nHGChebLayers << std::endl;

    //set layer max vals in calib objects
    m_calibEE.m_TotalLayers = nHGCeeLayers;
    m_calibHEF.m_TotalLayers = nHGChefLayers;
    m_calibHEB.m_TotalLayers = nHGChebLayers;

    //open ROOT file with histograms containing layer depths
    TFile* file = TFile::Open(m_layerDepthFilename.fullPath().c_str(),"READ");
    TH1F* h_x0 = (TH1F*)file->Get("x0");
    TH1F* h_lambda = (TH1F*)file->Get("lambda");
    unsigned h_max = h_x0->GetNbinsX();
    for(unsigned ih = 1; ih <= h_max; ih++){
        if(ih <= nHGCeeLayers) {
            m_calibEE.nCellRadiationLengths.push_back(h_x0->GetBinContent(ih));
            m_calibEE.nCellInteractionLengths.push_back(h_lambda->GetBinContent(ih));
        }
        else if(ih <= nHGCeeLayers + nHGChefLayers){
            m_calibHEF.nCellRadiationLengths.push_back(h_x0->GetBinContent(ih));
            m_calibHEF.nCellInteractionLengths.push_back(h_lambda->GetBinContent(ih));
        }
        else if(ih <= nHGCeeLayers + nHGChefLayers + nHGChebLayers){
            m_calibHEB.nCellRadiationLengths.push_back(h_x0->GetBinContent(ih));
            m_calibHEB.nCellInteractionLengths.push_back(h_lambda->GetBinContent(ih));
        }
    }
    //close file
    file->Close();
    //set calibration layers: EE*2* for EE, HEF1 for HEF and HEB
    m_calibEE.calibrationRadiationLength = m_calibEE.nCellRadiationLengths[2];
    m_calibEE.calibrationInteractionLength = m_calibEE.nCellInteractionLengths[2];
    m_calibHEF.calibrationRadiationLength = m_calibHEF.nCellRadiationLengths[1];
    m_calibHEF.calibrationInteractionLength = m_calibHEF.nCellInteractionLengths[1];
    m_calibHEB.calibrationRadiationLength = m_calibHEF.nCellRadiationLengths[1];
    m_calibHEB.calibrationInteractionLength = m_calibHEF.nCellInteractionLengths[1];

    //toggle use of overburden corrections
    m_calibEE.useOverburdenCorrection = m_useOverburdenCorrection;

    //open ROOT file with graphs containing overburden depths vs eta
    //(for HGCEE layer 1 eta-dependent depth correction)
    if(m_calibEE.useOverburdenCorrection){
      file = TFile::Open(m_overburdenDepthFilename.fullPath().c_str(),"READ");
      TGraphErrors* g_x0 = (TGraphErrors*)file->Get("x0Overburden");
      TGraphErrors* g_lambda = (TGraphErrors*)file->Get("lambdaOverburden");
      double *x_x0, *xe_x0, *y_x0, *x_lambda, *xe_lambda, *y_lambda;
      x_x0 = g_x0->GetX();
      xe_x0 = g_x0->GetEX();
      y_x0 = g_x0->GetY();
      x_lambda = g_lambda->GetX();
      xe_lambda = g_lambda->GetEX();
      y_lambda = g_lambda->GetY();
      int nbins = g_x0->GetN();
      //fill map with low x edges and y values
      for(int i = 0; i < nbins; i++){
          m_calibEE.nOverburdenRadiationLengths[x_x0[i] - xe_x0[i]] = y_x0[i];
          m_calibEE.nOverburdenInteractionLengths[x_lambda[i] - xe_lambda[i]] = y_lambda[i];
      }
      //close file
      file->Close();
    }
    
    //initialize corrections after getting all calibrations (in beginJob) and layer depths
    m_calibEE.initialize();
    m_calibHEF.initialize();
    m_calibHEB.initialize();
    
  }
  
  std::vector<double> min_innerR_depth_ee, min_innerZ_depth_ee ; 
  std::vector<double> min_innerR_depth_hef, min_innerZ_depth_hef ; 
  std::vector<double> min_innerR_depth_heb, min_innerZ_depth_heb ;
  for (unsigned int i=0; i<=nHGCeeLayers; i++) { 
    PandoraApi::Geometry::LayerParameters *eeLayerParameters;
    eeLayerParameters = new PandoraApi::Geometry::LayerParameters();
    hgcEELayerParameters.push_back( eeLayerParameters ) ; 
    min_innerR_depth_ee.push_back( 99999.0 ) ; 
    min_innerZ_depth_ee.push_back( 99999.0 ) ;
  }
  for (unsigned int i=0; i<=nHGChefLayers; i++) { 
    PandoraApi::Geometry::LayerParameters *hefLayerParameters;
    hefLayerParameters = new PandoraApi::Geometry::LayerParameters();
    hgcHEFLayerParameters.push_back( hefLayerParameters ) ; 
    min_innerR_depth_hef.push_back( 99999.0 ) ; 
    min_innerZ_depth_hef.push_back( 99999.0 ) ;
  }
  for (unsigned int i=0; i<=nHGChebLayers; i++) { 
      PandoraApi::Geometry::LayerParameters *hebLayerParameters;
      hebLayerParameters = new PandoraApi::Geometry::LayerParameters();
      hgcHEBLayerParameters.push_back( hebLayerParameters ) ; 
      min_innerR_depth_heb.push_back( 99999.0 ) ; 
      min_innerZ_depth_heb.push_back( 99999.0 ) ; 
  }

  // dummy vectors for CalculateCornerSubDetectorParameters function
  std::vector<double> min_innerR_depth_eb, min_innerZ_depth_eb ; 
  std::vector<double> min_innerR_depth_hb, min_innerZ_depth_hb ; 
  
  SetDefaultSubDetectorParameters("EcalBarrel", pandora::ECAL_BARREL, *ebParameters);
  SetDefaultSubDetectorParameters("EcalEndcap", pandora::ECAL_ENDCAP, *eeParameters);
  SetDefaultSubDetectorParameters("HcalBarrel", pandora::HCAL_BARREL, *hbParameters);
  SetDefaultSubDetectorParameters("HcalEndcap", pandora::HCAL_ENDCAP, *heParameters);

  //corner parameters for EB
  double min_innerRadius = 99999.0 ; double max_outerRadius = 0.0 ;
  double min_innerZ = 99999.0 ; double max_outerZ = 0.0 ;
  CalculateCornerSubDetectorParameters(ecalBarrelGeometry, ecalBarrelCells, pandora::ECAL_BARREL, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       false, min_innerR_depth_eb, min_innerZ_depth_eb);
  SetCornerSubDetectorParameters(*ebParameters, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ); // One ECAL layer
  SetSingleLayerParameters(*ebParameters,*ebLayerParameters);
  ebParameters->m_nLayers = 1 ; // One ECAL layer
  
  //corner parameters for HB
  min_innerRadius = 99999.0 ; max_outerRadius = 0.0 ;
  min_innerZ = 99999.0 ; max_outerZ = 0.0 ;
  CalculateCornerSubDetectorParameters(hcalBarrelGeometry, hcalBarrelCells, pandora::HCAL_BARREL, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       false, min_innerR_depth_hb, min_innerZ_depth_hb);
  SetCornerSubDetectorParameters(*hbParameters, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ);
  SetSingleLayerParameters(*hbParameters, *hbLayerParameters);
  hbParameters->m_nLayers = 1 ; //todo: include HCAL layers
  
  //corner & layer parameters for EE
  min_innerRadius = 99999.0 ; max_outerRadius = 0.0 ;
  min_innerZ = 99999.0 ; max_outerZ = 0.0 ;
  CalculateCornerSubDetectorParameters(HGCEEGeometry, ecalEndcapCells, pandora::ECAL_ENDCAP, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       true, min_innerR_depth_ee, min_innerZ_depth_ee);
  SetCornerSubDetectorParameters(*eeParameters, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ);
  SetMultiLayerParameters(*eeParameters, hgcEELayerParameters, min_innerR_depth_ee, min_innerZ_depth_ee, nHGCeeLayers, m_calibEE);
  eeParameters->m_nLayers = nHGCeeLayers;

  //corner & layer parameters for HE
  //consider both HEF and HEB together
  min_innerRadius = 99999.0 ; max_outerRadius = 0.0 ;
  min_innerZ = 99999.0 ; max_outerZ = 0.0 ;
  CalculateCornerSubDetectorParameters(HGCHEFGeometry, hcalEndcapCellsFront, pandora::HCAL_ENDCAP, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       true, min_innerR_depth_hef, min_innerZ_depth_hef);
  CalculateCornerSubDetectorParameters(HGCHEBGeometry, hcalEndcapCellsBack, pandora::HCAL_ENDCAP, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ,
                                       true, min_innerR_depth_heb, min_innerZ_depth_heb);
  SetCornerSubDetectorParameters(*heParameters, min_innerRadius, max_outerRadius, min_innerZ, max_outerZ);
  SetMultiLayerParameters(*heParameters, hgcHEFLayerParameters, min_innerR_depth_hef, min_innerZ_depth_hef, nHGChefLayers, m_calibHEF);
  SetMultiLayerParameters(*heParameters, hgcHEBLayerParameters, min_innerR_depth_heb, min_innerZ_depth_heb, nHGChebLayers, m_calibHEB);
  heParameters->m_nLayers = nHGChefLayers + nHGChebLayers;

  PandoraApi::Geometry::SubDetector::Create(*m_pPandora, *ebParameters);
  PandoraApi::Geometry::SubDetector::Create(*m_pPandora, *hbParameters);
  PandoraApi::Geometry::SubDetector::Create(*m_pPandora, *eeParameters);
  PandoraApi::Geometry::SubDetector::Create(*m_pPandora, *heParameters);

  // make propagator
  constexpr float m_pion = 0.1396;
  _mat_prop.reset( new PropagatorWithMaterial(alongMomentum, m_pion, magneticField.product()) );
  // setup HGC layers for track propagation
  Surface::RotationType rot; //unit rotation matrix

  _minusSurface.clear();
  _plusSurface.clear();
  const HGCalDDDConstants &dddCons = HGCEEGeometry->topology().dddConstants();
  std::map<float,float> zrhoCoord;
  auto firstLayerIt = dddCons.getFirstTrForm();
  auto lastLayerIt = dddCons.getLastTrForm();
  for(auto layerIt=firstLayerIt; layerIt !=lastLayerIt; layerIt++) {
    float Z(fabs(layerIt->h3v.z()));
    auto lastmod = std::reverse_iterator<std::vector<HGCalDDDConstants::hgtrap>::const_iterator>(dddCons.getLastModule(true));
    float Radius(lastmod->tl+layerIt->h3v.perp());
    zrhoCoord[Z]=Radius;
  }
  //take only innermost layer
  auto it=zrhoCoord.begin();
  float Z(it->first);
  float Radius(it->second);
  _minusSurface.push_back(new BoundDisk( Surface::PositionType(0,0,-Z), rot, new SimpleDiskBounds( 0, Radius, -0.001, 0.001)));
  _plusSurface.push_back(new BoundDisk( Surface::PositionType(0,0,+Z), rot, new SimpleDiskBounds( 0, Radius, -0.001, 0.001)));

}

void PandoraCMSPFCandProducer::SetDefaultSubDetectorParameters(const std::string &subDetectorName, const pandora::SubDetectorType subDetectorType, PandoraApi::Geometry::SubDetector::Parameters &parameters) const {
  //identification
  parameters.m_subDetectorName = subDetectorName;
  parameters.m_subDetectorType = subDetectorType;
  
  // Phi Coordinate is when start drawing the detector, wrt x-axis.  
  // Assuming this is 0 since CMS ranges from -pi to pi
  parameters.m_innerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  parameters.m_outerPhiCoordinate = 0.0 ; // -1.0 * CLHEP::pi ; 
  
  // Symmetry order is how you draw the "polygon" detector.  
  // Circle approximation for now (0), but can be configured to match N(cells)
  parameters.m_innerSymmetryOrder = 0 ; 
  parameters.m_outerSymmetryOrder = 0 ; 
  
  parameters.m_isMirroredInZ = true ; // Duplicate detector +/- z
}

void PandoraCMSPFCandProducer::CalculateCornerSubDetectorParameters(const CaloSubdetectorGeometry* geom,  const std::vector<DetId>& cells, const pandora::SubDetectorType subDetectorType, 
                                                      double& min_innerRadius, double& max_outerRadius, double& min_innerZ, double& max_outerZ,
                                                      bool doLayers, std::vector<double>& min_innerR_depth, std::vector<double>& min_innerZ_depth) const
{
  //barrel:
  // Inner radius taken as average magnitude of (x,y) for corners 0-3
  // Outer radius taken as average magnitude of (x,y) for corners 4-7  
  double ci_barrel[4] = {0,1,2,3};
  double co_barrel[4] = {4,5,6,7};
  
  //endcap:
  // Inner radius taken as average magnitude of (x,y) for corners 0,3,4,7 
  // Outer radius taken as average magnitude of (x,y) for corners 1,2,5,6  
  double ci_endcap[4] = {0,3,4,7};
  double co_endcap[4] = {1,2,5,6};
  
  double *ci, *co;
  if(subDetectorType==pandora::ECAL_BARREL || subDetectorType==pandora::HCAL_BARREL){
    ci = ci_barrel;
    co = co_barrel;
  }
  else { //if(subDetectorType==pandora::ECAL_ENDCAP || subDetectorType==pandora::HCAL_ENDCAP){
    ci = ci_endcap;
    co = co_endcap;
  }
  
  // Determine: inner/outer detector radius
  for (std::vector<DetId>::const_iterator ib=cells.begin(); ib!=cells.end(); ib++) {
    const CaloCellGeometry *thisCell = geom->getGeometry(*ib);
    const CaloCellGeometry::CornersVec& corners = thisCell->getCorners();
    
    //kind of hacky
    unsigned int layer = 0;
    if(doLayers && subDetectorType==pandora::ECAL_ENDCAP) layer = (unsigned int) ((HGCEEDetId)(*ib)).layer() ;
    else if(doLayers && subDetectorType==pandora::HCAL_ENDCAP) layer = (unsigned int) ((HGCHEDetId)(*ib)).layer() ;    
    
    //inner radius calculation
    double avgX_inner = 0.25 * (corners[ci[0]].x() + corners[ci[1]].x() + corners[ci[2]].x() + corners[ci[3]].x()) ;
    double avgY_inner = 0.25 * (corners[ci[0]].y() + corners[ci[1]].y() + corners[ci[2]].y() + corners[ci[3]].y()) ;
    double innerRadius = sqrt( avgX_inner * avgX_inner + avgY_inner * avgY_inner ) ;
    if ( innerRadius < min_innerRadius ) min_innerRadius = innerRadius ;
    if ( doLayers && innerRadius < min_innerR_depth.at(layer) ) min_innerR_depth.at(layer) = innerRadius ; 
    
    //outer radius calculation
    double avgX_outer = 0.25 * (corners[co[0]].x() + corners[co[1]].x() + corners[co[2]].x() + corners[co[3]].x()) ;
    double avgY_outer = 0.25 * (corners[co[0]].y() + corners[co[1]].y() + corners[co[2]].y() + corners[co[3]].y()) ;
    double outerRadius = sqrt( avgX_outer * avgX_outer + avgY_outer * avgY_outer ) ;
    if ( outerRadius > max_outerRadius ) max_outerRadius = outerRadius ;
    
    //z calculation
    for( unsigned int isubcell = 0; isubcell<8; isubcell++){
        if( fabs(corners[isubcell].z()) < min_innerZ ) min_innerZ = fabs(corners[isubcell].z()) ;
        if ( corners[isubcell].z() > max_outerZ ) max_outerZ = corners[isubcell].z();
        if ( doLayers && fabs(corners[isubcell].z()) < min_innerZ_depth.at(layer) ) min_innerZ_depth.at(layer) = fabs(corners[isubcell].z()) ;
    }
  }

}

void PandoraCMSPFCandProducer::SetCornerSubDetectorParameters(PandoraApi::Geometry::SubDetector::Parameters &parameters, 
                                                const double& min_innerRadius, const double& max_outerRadius, const double& min_innerZ, const double& max_outerZ) const 
{
  parameters.m_innerRCoordinate = min_innerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  parameters.m_outerRCoordinate = max_outerRadius * 10.0 ; // CMS units cm, Pandora expects mm
  parameters.m_innerZCoordinate = min_innerZ * 10.0 ; // CMS units cm, Pandora expects mm
  parameters.m_outerZCoordinate = max_outerZ * 10.0 ; // CMS units cm, Pandora expects mm
}

void PandoraCMSPFCandProducer::SetSingleLayerParameters(PandoraApi::Geometry::SubDetector::Parameters &parameters, PandoraApi::Geometry::LayerParameters &layerParameters) const {
  layerParameters.m_closestDistanceToIp = parameters.m_innerRCoordinate; 
  layerParameters.m_nInteractionLengths = 0.0 ;
  layerParameters.m_nRadiationLengths = 0.0 ;
  parameters.m_layerParametersList.push_back(layerParameters) ; 
}

void PandoraCMSPFCandProducer::SetMultiLayerParameters(PandoraApi::Geometry::SubDetector::Parameters &parameters, std::vector<PandoraApi::Geometry::LayerParameters*> &layerParameters,
                                         std::vector<double>& min_innerR_depth, std::vector<double>& min_innerZ_depth, const unsigned int& nTotalLayers, CalibHGC& calib) const 
{
  for (unsigned int i=1; i<=nTotalLayers; i++) { //skip nonexistent layer 0
    double distToIP = 10.0 * sqrt(min_innerR_depth.at(i)*min_innerR_depth.at(i) + min_innerZ_depth.at(i)*min_innerZ_depth.at(i)) ; 
    layerParameters.at(i)->m_closestDistanceToIp = distToIP ; 
    layerParameters.at(i)->m_nInteractionLengths = calib.nCellInteractionLengths[i] ;
    layerParameters.at(i)->m_nRadiationLengths = calib.nCellRadiationLengths[i] ;
    parameters.m_layerParametersList.push_back(*(layerParameters.at(i))) ; 
  }
}

void PandoraCMSPFCandProducer::prepareTrack(edm::Event& iEvent){ // function to setup tracks in an event for pandora
  //Why PF uses global point (0,0,0) for all events?
  math::XYZVector B_(math::XYZVector(magneticField->inTesla(GlobalPoint(0,0,0))));

  edm::Handle<reco::RecoToSimCollection > rectosimCollection;
  if(useRecoTrackAsssociation) iEvent.getByLabel(inputTagtPRecoTrackAsssociation_, rectosimCollection);
  
  PandoraApi::Track::Parameters trackParameters;
  //We need the speed of light
  if(debugPrint) std::cout << speedoflight << " mm/ns" << std::endl;

  edm::Handle<reco::PFRecTrackCollection> tkRefCollection;
  bool found = iEvent.getByLabel(inputTagGeneralTracks_, tkRefCollection);
  if(!found ) {
    std::ostringstream err;
    err<<"cannot find generalTracks: "<< inputTagGeneralTracks_;
    LogError("PandoraCMSPFCandProducer")<<err.str()<<std::endl;
    throw cms::Exception( "MissingProduct", err.str());
  }

  for(unsigned int i=0; i<tkRefCollection->size(); i++) {
    const reco::PFRecTrack * pftrack = &(*tkRefCollection)[i];
    const reco::TrackRef track = pftrack->trackRef();
    
    //std::cout << "got track with algo == " << track->algo() << std::endl;

    //For the d0 = -dxy
    trackParameters.m_d0 = track->d0() * 10. ; //in mm
    //For the z0
    trackParameters.m_z0 = track->dz() * 10. ; //in mm
    //For the Track momentum at the 2D distance of closest approach
    //For tracks reconstructed in the CMS Tracker, the reference position is the point of closest approach to the centre of CMS. (math::XYZPoint posClosest = track->referencePoint();)
    // According to TrackBase.h the momentum() method returns momentum vector at innermost (reference) point on track
    const pandora::CartesianVector momentumAtClosestApproach(track->momentum().x(),track->momentum().y(),track->momentum().z()); //in GeV
    trackParameters.m_momentumAtDca = momentumAtClosestApproach;
 
    //For the track of the state at the start in mm and GeV
    const pandora::CartesianVector positionAtStart(track->innerPosition().x()* 10.,track->innerPosition().y()* 10., track->innerPosition().z() * 10. );
    const pandora::CartesianVector momentumAtStart(track->innerMomentum().x(),track->innerMomentum().y(), track->innerMomentum().z() );
    trackParameters.m_trackStateAtStart = pandora::TrackState(positionAtStart,momentumAtStart);
    //For the track of the state at the end in mm and GeV
    const pandora::CartesianVector positionAtEnd(track->outerPosition().x() * 10.,track->outerPosition().y() * 10., track->outerPosition().z() * 10.);
    const pandora::CartesianVector momentumAtEnd(track->outerMomentum().x(),track->outerMomentum().y(), track->outerMomentum().z() );
    trackParameters.m_trackStateAtEnd = pandora::TrackState(positionAtEnd,momentumAtEnd);
    //For the charge
    double charge = track->charge();
    trackParameters.m_charge = charge;
    //Associate the reconstructed Track (in the Tracker) with the corresponding MC true simulated particle
    trackParameters.m_particleId = 211; // INITIALIZATION // NS

    std::vector<std::pair<TrackingParticleRef, double> > tp;
    TrackingParticleRef tpr; 
    edm::RefToBase<reco::Track> tr(track);
    if(useRecoTrackAsssociation){
      const reco::RecoToSimCollection pRecoToSim = *(rectosimCollection.product());		
      if(pRecoToSim.find(tr) != pRecoToSim.end()){
        tp = pRecoToSim[tr];
        if(debugPrint) std::cout << "Reco Track pT: "  << track->pt() <<  " matched to " << tp.size() << " MC Tracks" << " associated with quality: " << tp.begin()->second << std::endl;
        tpr = tp.begin()->first;
        
        trackParameters.m_particleId = (tpr->pdgId());
        if(debugPrint) std::cout << "the pdg id of this track is " << (tpr->pdgId()) << std::endl;
        //The parent vertex (from which this track was produced) has daughter particles.
        //These are the desire siblings of this track which we need. 
        TrackingParticleRefVector simSiblings = getTpSiblings(tpr);
        
        const TrackingParticle * sib; 
        int numofsibs = 0;
        std::vector<int> pdgidofsibs; pdgidofsibs.clear();
        
        if (simSiblings.isNonnull()) {
          for(TrackingParticleRefVector::iterator si = simSiblings.begin(); si != simSiblings.end(); si++){
            //Check if the current sibling is the track under study
            if ( (*si) ==  tpr  ) {continue;}
            sib = &(**si);
            pdgidofsibs.push_back(sib->pdgId());
            ++numofsibs;
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(*m_pPandora, pftrack, sib)); 
          }
          if(debugPrint) {
            std::cout << "This track has " << numofsibs << " sibling tracks with pdgids:" << std::endl; 
            for (std::vector<int>::iterator sib_pdg_it = pdgidofsibs.begin(); sib_pdg_it != pdgidofsibs.end(); sib_pdg_it++){
              std::cout << (*sib_pdg_it) << std::endl;
            }
          }
        }
        else {
          if(debugPrint) std::cout << "Particle pdgId = "<< (tpr->pdgId()) << " produced at rho = " << (tpr->vertex().Rho()) << ", z = " << (tpr->vertex().Z()) << ", has NO siblings!" << std::endl;
        }
      
        //Now the track under study has daughter particles. To find them we study the decay vertices of the track
        TrackingParticleRefVector simDaughters = getTpDaughters(tpr);
        const TrackingParticle * dau; 
        int numofdaus = 0;
        std::vector<int> pdgidofdaus; pdgidofdaus.clear();
      
        if (simDaughters.isNonnull()) {
          for(TrackingParticleRefVector::iterator di = simDaughters.begin(); di != simDaughters.end(); di++){
            //We have already checked that simDaughters don't contain the track under study
            dau = &(**di);
            pdgidofdaus.push_back(dau->pdgId());
            ++numofdaus;
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackParentDaughterRelationship(*m_pPandora, pftrack, dau)); 
          }
          if(debugPrint) {
            std::cout << "This track has " << numofdaus << " daughter tracks with pdgids:" << std::endl; 
            for (std::vector<int>::iterator dau_pdg_it = pdgidofdaus.begin(); dau_pdg_it != pdgidofdaus.end(); dau_pdg_it++){
              std::cout << (*dau_pdg_it) << std::endl;
            }
          }
        }        
        else {
          if(debugPrint) std::cout << "Particle pdgId = "<< (tpr->pdgId()) << " produced at rho = " << (tpr->vertex().Rho()) << ", z = " << (tpr->vertex().Z()) << ", has NO daughters!" << std::endl;
        }
      
      }
    }

    //The mass 
    trackParameters.m_mass = pandora::PdgTable::GetParticleMass(trackParameters.m_particleId.Get());

    //For the ECAL entrance
    // Starting from outermost hit position of the track and propagating to ECAL Entrance

    const TrajectoryStateOnSurface myTSOS = trajectoryStateTransform::outerStateOnSurface(*track, *(tkGeom.product()), magneticField.product());
    if(debugPrint) {
      std::cout << "magnetic field z " << B_.z() << std::endl;
      std::cout << "theOutParticle x position before propagation in cm "<< myTSOS.globalPosition().x()<< std::endl;
      std::cout << "theOutParticle x momentum before propagation in cm "<< myTSOS.globalMomentum().x()<< std::endl;
    }

    const auto& layer( (myTSOS.globalPosition().z() > 0 ? _plusSurface[0] : _minusSurface[0]) );
    //std::cout << "BoundDisk inner radius = " << layer->innerRadius() << ", outer radius = " << layer->outerRadius() << std::endl;
    TrajectoryStateOnSurface piStateAtSurface = _mat_prop->propagate(myTSOS, *layer);
    
    bool reachesCalorimeter = false;
    bool isonendcap = false;
    double mom = 0;
    if(piStateAtSurface.isValid()){
      // std::cout<< "!!!Reached ECAL!!! "<< std::endl;
      reachesCalorimeter = true;
      // std::cout<< "It is on the endcaps "<< std::endl;
      isonendcap = true;
      if(debugPrint) {
        std::cout << "theOutParticle x position after propagation to ECAL "<< piStateAtSurface.globalPosition().x()<< std::endl;
        std::cout << "theOutParticle x momentum after propagation to ECAL "<< piStateAtSurface.globalMomentum().x()<< std::endl;
      }
      mom = sqrt((piStateAtSurface.globalMomentum().x()*piStateAtSurface.globalMomentum().x()
                 +piStateAtSurface.globalMomentum().y()*piStateAtSurface.globalMomentum().y()
                 +piStateAtSurface.globalMomentum().z()*piStateAtSurface.globalMomentum().z()));
    }
    
    trackParameters.m_reachesCalorimeter = reachesCalorimeter;
    if (reachesCalorimeter){
      const pandora::CartesianVector positionAtCalorimeter(piStateAtSurface.globalPosition().x() * 10.,piStateAtSurface.globalPosition().y() * 10.,piStateAtSurface.globalPosition().z() * 10.);//in mm
      const pandora::CartesianVector momentumAtCalorimeter(piStateAtSurface.globalMomentum().x(),piStateAtSurface.globalMomentum().y(),piStateAtSurface.globalMomentum().z());
      trackParameters.m_trackStateAtCalorimeter = pandora::TrackState(positionAtCalorimeter, momentumAtCalorimeter);
      // For the time at calorimeter we need the speed of light
      //This is in BaseParticlePropagator c_light() method in mm/ns but is protected (299.792458 mm/ns)
      //So we take it from CLHEP
      trackParameters.m_timeAtCalorimeter = positionAtCalorimeter.GetMagnitude() / speedoflight; // in ns
    }
    else { 
      trackParameters.m_trackStateAtCalorimeter = trackParameters.m_trackStateAtEnd.Get();
      //trackParameters.m_timeAtCalorimeter = std::numeric_limits<double>::max();
      trackParameters.m_timeAtCalorimeter = 999999999.;
    }

    trackParameters.m_isProjectedToEndCap = isonendcap;

    bool canFormPfo = false;
    bool canFormClusterlessPfo = false;

    if(trackParameters.m_reachesCalorimeter.Get()>0 && mom>1. && track->quality(reco::TrackBase::tight) ){
      canFormPfo = true;
      canFormClusterlessPfo = true;
    }
    const float deltaPtRel = track->ptError()/track->pt();
    canFormClusterlessPfo *= (deltaPtRel < _deltaPtOverPtForClusterlessPfo); //tighter requirements for tracks with no cluster

    trackParameters.m_canFormPfo = canFormPfo;
    trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;

    //The parent address
    trackParameters.m_pParentAddress =  (void *) pftrack;
    recTrackMap.emplace((void*)pftrack,i); //associate parent address with collection index    
 
    if( deltaPtRel < _deltaPtOverPtForPfo ) {
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(*m_pPandora, trackParameters));    
    }
  }

}

void PandoraCMSPFCandProducer::prepareHits( edm::Event& iEvent)
{
  // Get the primary vertex
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByLabel("offlinePrimaryVertices", pvHandle);
  reco::Vertex pv = pvHandle->at(0) ; 

  // get the HGC PFRecHit collection
  edm::Handle<reco::PFRecHitCollection> HGCRecHitHandle;
  
  bool found = iEvent.getByLabel(inputTagHGCrechit_, HGCRecHitHandle);
  
  if(!found ) {
    std::ostringstream err;
    err<<"cannot find rechits: "<< inputTagHGCrechit_;
    LogError("PandoraCMSPFCandProducer")<<err.str()<<std::endl;
    throw cms::Exception( "MissingProduct", err.str());
  }

  // get the sim hits (make sure they are same order in the python!
  std::vector<edm::Handle<std::vector<PCaloHit> > > simHits;  
  for( const auto& tag : inputSimHits_ ) {
    simHits.emplace_back( edm::Handle<std::vector<PCaloHit> >() );
    iEvent.getByLabel(tag,simHits.back());
  }

  // get sim tracks to associate to
  edm::Handle<std::vector<SimTrack> > simTracks;
  iEvent.getByLabel(inputTagSimTracks_,simTracks);
  

  PandoraApi::RectangularCaloHitParameters caloHitParameters;

  if(debugPrint) std::cout<< speedoflight << " cm/ns" << std::endl;

  sumCaloEnergy = 0.;
  sumCaloEnergyEM = 0.;
  sumCaloEnergyHAD = 0.;

  simDir_sumCaloEnergyEM = 0.;
  simDir_sumCaloEnergyHAD = 0.;

  sumCaloECALEnergyEM = 0.;
  sumCaloHCALEnergyEM  = 0.;
  sumCaloECALEnergyHAD= 0.;
  sumCaloECALEnergyHAD_unc= 0.;
  sumCaloHCALEnergyHAD = 0.;
  
  // Get the HGC geometry
  const HGCalGeometry* HGCEEGeometry = hgceeGeoHandle.product() ; 
  const HGCalGeometry* HGCHEFGeometry = hgchefGeoHandle.product() ; 
  const HGCalGeometry* HGCHEBGeometry = hgchebGeoHandle.product() ; 

  //
  // Process HGC rec hits 
  // 
  int nNotFoundEE = 0, nCaloHitsEE = 0 ; 
  int nNotFoundHEF = 0, nCaloHitsHEF = 0 ; 
  int nNotFoundHEB = 0, nCaloHitsHEB = 0 ; 
  for(unsigned i=0; i<HGCRecHitHandle->size(); i++) {
    const reco::PFRecHit* rh = &(*HGCRecHitHandle)[i];
    const DetId detid(rh->detId());
    if (detid.subdetId() == 3) ProcessRecHits(simTracks,simHits,rh, i, HGCEEGeometry, m_calibEE, nCaloHitsEE, nNotFoundEE, pv, pandora::ECAL, pandora::ENDCAP, caloHitParameters);
    else if (detid.subdetId() == 4) ProcessRecHits(simTracks,simHits,rh, i, HGCHEFGeometry, m_calibHEF, nCaloHitsHEF, nNotFoundHEF, pv, pandora::HCAL, pandora::ENDCAP, caloHitParameters);
    else if (detid.subdetId() == 5) ProcessRecHits(simTracks,simHits,rh, i, HGCHEBGeometry, m_calibHEB, nCaloHitsHEB, nNotFoundHEB, pv, pandora::HCAL, pandora::ENDCAP, caloHitParameters);
    else continue;
  }

  // SCZ nein
  if (false) {
    //  if(debugHisto){
    h_sumCaloE->Fill(sumCaloEnergy);
    h_sumCaloEM->Fill(sumCaloEnergyEM);
    h_sumCaloHad->Fill(sumCaloEnergyHAD);
    h_simDir_sumCaloEM ->Fill(simDir_sumCaloEnergyEM );
    h_simDir_sumCaloHad->Fill(simDir_sumCaloEnergyHAD);
    
    h2_Calo_EM_hcalEecalE->Fill(sumCaloECALEnergyEM, sumCaloHCALEnergyEM);
    h2_Calo_Had_hcalEecalE->Fill(sumCaloECALEnergyHAD, sumCaloHCALEnergyHAD);
    h_sumEcalEEM->Fill(sumCaloECALEnergyEM);
    h_sumHcalEEM->Fill(sumCaloHCALEnergyEM);
    
    h_sumEcalEHad->Fill(sumCaloECALEnergyHAD);
    h_sumEcalEHad_unc->Fill(sumCaloECALEnergyHAD_unc);
    if(sumCaloECALEnergyHAD>=0.5) h_sumEcalEHadc->Fill(sumCaloECALEnergyHAD);
    if(sumCaloECALEnergyHAD<0.5) h_sumHcalEHadc->Fill(sumCaloHCALEnergyHAD);
    h_sumHcalEHad->Fill(sumCaloHCALEnergyHAD);
    h_sumEHad->Fill(sumCaloHCALEnergyHAD+sumCaloECALEnergyHAD);
    
    for (int ilay=0; ilay<100; ilay++) {
       int ibin = ilay+1;
       
       h_hitEperLayer_EM[ForwardSubdetector::HGCEE] ->SetBinContent(ibin,m_hitEperLayer_EM[ForwardSubdetector::HGCEE] [ilay]+ h_hitEperLayer_EM[ForwardSubdetector::HGCEE] ->GetBinContent(ibin));
       h_hitEperLayer_EM[ForwardSubdetector::HGCHEF]->SetBinContent(ibin,m_hitEperLayer_EM[ForwardSubdetector::HGCHEF][ilay]+ h_hitEperLayer_EM[ForwardSubdetector::HGCHEF]->GetBinContent(ibin));
       h_hitEperLayer_EM[ForwardSubdetector::HGCHEB]->SetBinContent(ibin,m_hitEperLayer_EM[ForwardSubdetector::HGCHEB][ilay]+ h_hitEperLayer_EM[ForwardSubdetector::HGCHEB]->GetBinContent(ibin));
       
       h_hitEperLayer_HAD[ForwardSubdetector::HGCEE] ->SetBinContent(ibin,m_hitEperLayer_HAD[ForwardSubdetector::HGCEE] [ilay]+ h_hitEperLayer_HAD[ForwardSubdetector::HGCEE] ->GetBinContent(ibin));
       h_hitEperLayer_HAD[ForwardSubdetector::HGCHEF]->SetBinContent(ibin,m_hitEperLayer_HAD[ForwardSubdetector::HGCHEF][ilay]+ h_hitEperLayer_HAD[ForwardSubdetector::HGCHEF]->GetBinContent(ibin));
       h_hitEperLayer_HAD[ForwardSubdetector::HGCHEB]->SetBinContent(ibin,m_hitEperLayer_HAD[ForwardSubdetector::HGCHEB][ilay]+ h_hitEperLayer_HAD[ForwardSubdetector::HGCHEB]->GetBinContent(ibin));
    }
  }
  if(debugPrint) {
    std::cout << "sumCaloEnergy = " << sumCaloEnergy << std::endl;
    std::cout << "sumCaloEnergyEM  = " << sumCaloEnergyEM  << std::endl;
    std::cout << "sumCaloEnergyHAD = " << sumCaloEnergyHAD << std::endl;
    
    std::cout << "prepareHits HGC summary: " << std::endl ; 
    std::cout << "HGC Calo Hits               : " << nCaloHitsEE << " (HGC EE) " 
          << nCaloHitsHEF << " (HGC HEF) " << nCaloHitsHEB << " (HGC HEB) " << std::endl ;
    std::cout << "DetIDs not found in geometry: " << nNotFoundEE << " (HGC EE) " 
          << nNotFoundHEF << " (HGC HEF) " << nNotFoundHEB << " (HGC HEB) " << std::endl ;
  }
}

void PandoraCMSPFCandProducer::ProcessRecHits(const edm::Handle<std::vector<SimTrack> >& simTracks,
                                              const std::vector<edm::Handle<std::vector<PCaloHit> > >& simHits,
                                              const reco::PFRecHit* rh, unsigned int index, const HGCalGeometry* geom, CalibHGC& calib, int& nCaloHits, int& nNotFound, reco::Vertex& pv,
                 const pandora::HitType hitType, const pandora::HitRegion hitRegion, PandoraApi::RectangularCaloHitParameters& caloHitParameters)
{
    const DetId detid(rh->detId());
    double eta = fabs(rh->position().Eta());
    //    double cos_theta = std::tanh(rh->position().Eta()); 
    double cos_theta = 1.; // We do not need this correction because increasing thickness of absorber material compensates
                           // for increasing path length at non-normal incidence
    double energy = rh->energy() * cos_theta * calib.GetADC2GeV(); 
        
    if (energy < calib.m_CalThresh) return;
    
    double time = rh->time();
    // std::cout << "energy " << energy <<  " time " << time <<std::endl;
    
    const CaloCellGeometry *thisCell = geom->getGeometry(detid);
    if(!thisCell) {
        LogError("PandoraCMSPFCandProducerPrepareHits") << "warning detid " << detid.rawId() << " not found in geometry" << std::endl;
        nNotFound++;
        return;
    }
    
    unsigned int layer = 0;
    if(calib.m_id==ForwardSubdetector::HGCEE) layer = (unsigned int) ((HGCEEDetId)(detid)).layer() ;
    else if(calib.m_id==ForwardSubdetector::HGCHEF || calib.m_id==ForwardSubdetector::HGCHEB) layer = (unsigned int) ((HGCHEDetId)(detid)).layer() ;
    
    //hack because calo and HGC CornersVec are different formats
    const HGCalGeometry::CornersVec corners = ( std::move( geom->getCorners( detid ) ) );
    assert( corners.size() == 8 );
    
    // Various thickness measurements: 
    // m_cellSizeU --> Should be along beam for barrel, so along z...take as 0 <--> 1
    // m_cellSizeV --> Perpendicular to U and to thickness, but what is thickness?...take as 0 <--> 3
    // m_cellThickness --> Equivalent to depth?...take as 0 <--> 4
    const pandora::CartesianVector corner0( corners[0].x(), corners[0].y(), corners[0].z() );
    const pandora::CartesianVector corner1( corners[1].x(), corners[1].y(), corners[1].z() );
    const pandora::CartesianVector corner3( corners[3].x(), corners[3].y(), corners[3].z() );
    const pandora::CartesianVector corner4( corners[4].x(), corners[4].y(), corners[4].z() );
    caloHitParameters.m_cellSizeU     = 10.0 * (corner0 - corner1).GetMagnitude() ; 
    caloHitParameters.m_cellSizeV     = 10.0 * (corner0 - corner3).GetMagnitude() ; 
    caloHitParameters.m_cellThickness = 10.0 * (corner0 - corner4).GetMagnitude() ; 
    
//   for (unsigned int i=0; i<8; i++) { 
//     std::cout << "Corners " << i << ": x " << corners[i].x() << " y " << corners[i].y() << " z " << corners[i].z() << std::endl ; 
//   }
    
    // Position is average of all eight corners, convert from cm to mm
    double x = 0.0, y = 0.0, z = 0.0 ; 
    double xf = 0.0, yf = 0.0, zf = 0.0 ; 
    double xb = 0.0, yb = 0.0, zb = 0.0 ; 
    for (unsigned int i=0; i<8; i++) {
      if ( i < 4 ) { xf += corners[i].x() ; yf += corners[i].y() ; zf += corners[i].z() ; }
      else { xb += corners[i].x() ; yb += corners[i].y() ; zb += corners[i].z() ; }
      x += corners[i].x() ; y += corners[i].y() ; z += corners[i].z() ; 
    }
    // Average x,y,z position 
    x = x / 8.0 ; y = y / 8.0 ; z = z / 8.0 ; 
    xf = xf / 8.0 ; yf = yf / 8.0 ; zf = zf / 8.0 ; 
    xb = xb / 8.0 ; yb = yb / 8.0 ; zb = zb / 8.0 ; 
    const pandora::CartesianVector positionVector(10.0*x,10.0*y,10.0*z);
    caloHitParameters.m_positionVector = positionVector;
    
    // Expected direction (currently) drawn from primary vertex to front face of calorimeter cell
    const pandora::CartesianVector axisVector(10.0*(xf-pv.x()),10.0*(yf-pv.y()),10.0*(zf-pv.z())) ; 
    caloHitParameters.m_expectedDirection = axisVector.GetUnitVector();
    
    // Cell normal vector runs from front face to back of cell
    const pandora::CartesianVector normalVector(10.0*(xb-xf),10.0*(yb-yf),10.0*(zb-zf)) ; 
    caloHitParameters.m_cellNormalVector = normalVector.GetUnitVector();
    
    double distToFrontFace = sqrt( xf*xf + yf*yf + zf*zf ) ;
    // dist = cm, c = cm/nsec, rechit t in psec
    caloHitParameters.m_time = (distToFrontFace / speedoflight) + (time/1000.0) ; 
    
    //set hit and energy values
    caloHitParameters.m_hitType = hitType;
    caloHitParameters.m_hitRegion = hitRegion;
    caloHitParameters.m_inputEnergy = energy;
    caloHitParameters.m_electromagneticEnergy = calib.GetEMCalib(layer,eta) * energy;
    caloHitParameters.m_hadronicEnergy = calib.GetHADCalib(layer,eta) * energy; // = energy; 
    caloHitParameters.m_mipEquivalentEnergy = calib.m_CalToMip * energy;

    double angleCorrectionMIP(1.); 
    double hitR = distToFrontFace; 
    double hitZ = zf;
    angleCorrectionMIP = hitR/hitZ; 
    
    //---choose hits with eta-phi close to init. mc particle
    TVector3 hit3v(xf,yf,zf);
    double hitEta = hit3v.PseudoRapidity();
    double hitPhi = hit3v.Phi();
    //    if(debugHisto){ // SCZ NEIN
    if (false) {
      h_hit_Eta -> Fill(hitEta);
      h_hit_Phi -> Fill(hitPhi);
    }
    //    if (debugHisto && std::fabs(hitEta-m_firstMCpartEta) < 0.05 && std::fabs(hitPhi-m_firstMCpartPhi) < 0.05) {
    if (false) {
       h_MIP[calib.m_id] -> Fill(caloHitParameters.m_mipEquivalentEnergy.Get());
       h_MIP_Corr[calib.m_id] -> Fill(caloHitParameters.m_mipEquivalentEnergy.Get()*angleCorrectionMIP);
    }
    if (caloHitParameters.m_mipEquivalentEnergy.Get() < calib.m_CalMipThresh) {
       //std::cout << "EE MIP threshold rejected" << std::endl;
       return;
    }
    
    if(debugPrint || debugHisto){
      sumCaloEnergy += energy;
      sumCaloEnergyEM += caloHitParameters.m_electromagneticEnergy.Get();
      sumCaloEnergyHAD += caloHitParameters.m_hadronicEnergy.Get();
    }
    // SCZ NEIN
    if (false ){ 
      //      if(debugHisto){
      m_hitEperLayer_EM[calib.m_id][layer] += caloHitParameters.m_electromagneticEnergy.Get();
      m_hitEperLayer_HAD[calib.m_id][layer] += caloHitParameters.m_hadronicEnergy.Get();

      if( ( std::fabs(hitEta-m_firstMCpartEta) < 0.2 && std::fabs(hitPhi-m_firstMCpartPhi) < 0.2) || (std::fabs(hitEta-m_secondMCpartEta) < 0.2  && std::fabs(hitPhi-m_secondMCpartPhi) < 0.2) ){
        simDir_sumCaloEnergyEM  += caloHitParameters.m_electromagneticEnergy.Get();
        simDir_sumCaloEnergyHAD += caloHitParameters.m_hadronicEnergy.Get();
      }
      
      if(hitType==pandora::ECAL){
        sumCaloECALEnergyEM  += energy  ;
        sumCaloECALEnergyHAD += energy  ;
        sumCaloECALEnergyHAD_unc += energy ;
      }
      else if(debugHisto && hitType==pandora::HCAL){
        sumCaloHCALEnergyEM  += energy  ;
        sumCaloHCALEnergyHAD += energy  ;      
      }
    }
    
    caloHitParameters.m_layer = layer + (calib.m_id==ForwardSubdetector::HGCHEB ? nHGChefLayers : 0); //offset for HEB because combined with HEF in pandora
    caloHitParameters.m_nCellRadiationLengths = calib.nCellRadiationLengths[layer];
    caloHitParameters.m_nCellInteractionLengths = calib.nCellInteractionLengths[layer];
    caloHitParameters.m_isDigital = false;
    caloHitParameters.m_isInOuterSamplingLayer = ( calib.m_id==ForwardSubdetector::HGCHEB && layer == nHGChebLayers ) ;
    caloHitParameters.m_pParentAddress = (void *) rh;
    recHitMap.emplace((void*)rh,index); //associate parent address with collection index
    
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::CaloHit::Create(*m_pPandora, caloHitParameters));
    // also make the association to the paren sim particle
    unsigned detId = rh->detId();
    auto shrange = m_recoDetIdToSimHit.equal_range(detId);
    for( auto iter = shrange.first; iter != shrange.second; ++iter ) {
      const unsigned det = std::get<0>(iter->second);
      const unsigned hit_idx = std::get<1>(iter->second);
      const float fraction = std::get<2>(iter->second);
      const PCaloHit& simhit = simHits[det]->at(hit_idx);
      int tkid = simhit.geantTrackId();
      if( tkid > 0 ) { // <= 0 is a NULL track ID, i.e. no parent for the hit
        auto i_tk = m_barCodesToSimTrack.find(tkid);
        if( i_tk == m_barCodesToSimTrack.end() ) {
          throw cms::Exception("BadPCaloHitToTrackAssc")
            << "Track (" << tkid << ") associated to PCaloHit not found in list of barcodes!";
        }      
        const SimTrack& simtrack = simTracks->at(i_tk->second);
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(*m_pPandora, (const void *)rh,(const void*)&simtrack,fraction));
      }
    }
    
    nCaloHits++;
}


void PandoraCMSPFCandProducer::preparePandoraMCParticles(edm::Event& iEvent){ // function to setup a mcParticle for pandora
  
}

void PandoraCMSPFCandProducer::preparePFO(edm::Event& iEvent){
    // PandoraPFANew/v00-09/include/Pandora/PandoraInternal.h
    // typedef std::set<ParticleFlowObject *> PfoList;  
    //     PandoraPFANew/v00-09/include/Api/PandoraContentApi.h
        //    class ParticleFlowObject
        //    {
        //    public:
        //        /**
        //         *  @brief  Parameters class
        //         */
        //        class Parameters
        //        {
        //        public:
        //            pandora::InputInt               m_particleId;       ///< The particle flow object id (PDG code)
        //            pandora::InputInt               m_charge;           ///< The particle flow object charge
        //            pandora::InputFloat             m_mass;             ///< The particle flow object mass
        //            pandora::InputFloat             m_energy;           ///< The particle flow object energy
        //            pandora::InputCartesianVector   m_momentum;         ///< The particle flow object momentum
        //            pandora::ClusterList            m_clusterList;      ///< The clusters in the particle flow object
        //            pandora::TrackList              m_trackList;        ///< The tracks in the particle flow object
        //        };
        //        /**
        //         *  @brief  Create a particle flow object
        //         * 
        //         *  @param  algorithm the algorithm creating the particle flow object
        //         *  @param  particleFlowObjectParameters the particle flow object parameters
        //         */
        //        static pandora::StatusCode Create(const pandora::Algorithm &algorithm, const Parameters &parameters);
        //    };
        //    typedef ParticleFlowObject::Parameters ParticleFlowObjectParameters;

    // const pandora::CartesianVector momentum(1., 2., 3.);
    // for (pandora::PfoList::const_iterator itPFO = pPfoList->begin(), itPFOEnd = pPfoList->end(); itPFO != itPFOEnd; ++itPFO){
    //   (*itPFO)->SetParticleId();
    //   (*itPFO)->SetCharge();
    //   (*itPFO)->SetMass();
    //   (*itPFO)->SetEnergy();
    //   (*itPFO)->SetMomentum();
    // }
  
  const pandora::PfoList *pPfoList = NULL;
  // PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoList));

  const pandora::StatusCode statusCode(PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoList));  
  if (pandora::STATUS_CODE_SUCCESS != statusCode){throw pandora::StatusCodeException(statusCode);}

  edm::Handle<std::vector<reco::GenParticle> > genpart;
  iEvent.getByLabel(inputTagGenParticles_,genpart); //NS ADD
  if(debugPrint) std::cout << " GENPART SIZE IS " << genpart->size() << std::endl;

  // get the Calorimeter PFRecHit collections
  edm::Handle<reco::PFRecHitCollection> HGCRecHitHandle;
  iEvent.getByLabel(inputTagHGCrechit_, HGCRecHitHandle);
  
  Double_t found_energy = -1;
  Double_t ene_all_true = 0;

  //vars for histograms (only filled for i==0 in genparticle loop)
  double _sumPFOEnergy(0.);
  double sumClustEMEcalE(0.); //PFO cluster energy in Ecal
  double sumClustEMHcalE(0.); //PFO cluster energy in Hcal
  double sumClustHADEcalE(0.); //PFO cluster energy in Ecal
  double sumClustHADHcalE(0.); //PFO cluster energy in Hcal
  int sumPFOs(0);
  
  for(size_t i = 0; i < genpart->size(); ++ i) { // 1 
    const GenParticle * pa = &(*genpart)[i];
    ene_true    = pa->energy();

    if (pa->numberOfDaughters()==0)  ene_all_true = ene_all_true+ene_true;

    charge_true = 0;
    pT_true     = sqrt(pa->px()*pa->px()+pa->py()*pa->py());
    pid_true    = pa->pdgId();
    mass_true   = pa->mass();

    isDecBefCal = isDecayedBeforeCalo[i];

    if(pa->pdgId()>0) charge_true = 1;
    if(pa->pdgId()<0) charge_true = -1;
    if(pa->pdgId()==22 || pa->pdgId() == 130) charge_true = 0;
    //std::cout << " IN GENPART LOOP charge min " << std::endl;
    double diff_min  = 1e9;
    ene_match       = 0;
    mass_match      = 0;
    pid_match       = 0;
    pT_match        = 0;
    charge_match    = 0;
    ene_match_em    = 0;
    ene_match_had   = 0;
    ene_match_track = 0;
 
    nbPFOs = 0;
    for (pandora::PfoList::const_iterator itPFO = pPfoList->begin(), itPFOEnd = pPfoList->end(); itPFO != itPFOEnd; ++itPFO){ // 4
      nbPFOs++;
      double charge = (*itPFO)->GetCharge() ;
      double energy = (*itPFO)->GetEnergy();
      double pid    = (*itPFO)->GetParticleId();

      double mass   = (*itPFO)->GetMass() ;
      double pT     = sqrt(((*itPFO)->GetMomentum()).GetX()*((*itPFO)->GetMomentum()).GetX()+((*itPFO)->GetMomentum()).GetY()*((*itPFO)->GetMomentum()).GetY()) ; // EDW GetX

      const ClusterList &clusterList((*itPFO)->GetClusterList());
      //std::cout << " size of cluster list " << clusterList.size() << std::endl;
      ClusterVector clusterVector(clusterList.begin(), clusterList.end());
      ene_em  =0;
      ene_had =0;
 
      unsigned int firstLayer = std::numeric_limits<unsigned int>::max();
      unsigned int lastLayer = 0;

      for (ClusterVector::const_iterator clusterIter = clusterVector.begin(), clusterIterEnd = clusterVector.end(); clusterIter != clusterIterEnd; ++clusterIter) {
        const Cluster *pCluster = (*clusterIter);
        //ene_em  = pCluster->GetElectromagneticEnergy();
        //ene_had = pCluster->GetHadronicEnergy();
        // hits

	const float totalElectromagneticEnergy(pCluster->GetElectromagneticEnergy() - pCluster->GetIsolatedElectromagneticEnergy());
	mipFraction = pCluster->GetMipFraction();
	dCosR = 0.;
	clusterRms = 0.;

	const ClusterFitResult &clusterFitResult(pCluster->GetFitToAllHitsResult());

	if (clusterFitResult.IsFitSuccessful()) {
	  dCosR = clusterFitResult.GetRadialDirectionCosine();
	  clusterRms = clusterFitResult.GetRms();
	}

	const CartesianVector &clusterDirection(pCluster->GetFitToAllHitsResult().IsFitSuccessful() ?
						pCluster->GetFitToAllHitsResult().GetDirection() : pCluster->GetInitialDirection());

	const CartesianVector &clusterIntercept(pCluster->GetFitToAllHitsResult().IsFitSuccessful() ?
						pCluster->GetFitToAllHitsResult().GetIntercept() : CartesianVector(0.f, 0.f, 0.f));

	// Calculate properties of longitudinal shower profile: layer90 and shower max layer                                                                                           
	bool foundLayer90(false);
	float layer90EnergySum(0.f), /*showerMaxRadLengths(0.f),*/ energyAboveHighRadLengths(0.f);
	float nRadiationLengths(0.f), /*nRadiationLengths90(0.f),*/ nRadiationLengthsInLastLayer(0.f), maxEnergyInlayer(0.f);
	innerLayerRadLengths = 0.;
	nRadiationLengths90  = 0.;
	HitEnergyDistanceVector hitEnergyDistanceVector;

	const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
	const unsigned int innerPseudoLayer(pCluster->GetInnerPseudoLayer());
	const unsigned int firstPseudoLayer(m_pPandora->GetPlugins()->GetPseudoLayerPlugin()->GetPseudoLayerAtIp());

	for (unsigned int iLayer = innerPseudoLayer, outerPseudoLayer = pCluster->GetOuterPseudoLayer(); iLayer <= outerPseudoLayer; ++iLayer)
	  {
	    OrderedCaloHitList::const_iterator iter = orderedCaloHitList.find(iLayer);

	    if ((orderedCaloHitList.end() == iter) || (iter->second->empty()))
	      {
		nRadiationLengths += nRadiationLengthsInLastLayer;
		continue;
	      }

	    // Extract information from the calo hits                                                                                                                                  
	    float energyInLayer(0.f);
	    float nRadiationLengthsInLayer(0.f);

	    for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
	      {
		float cosOpeningAngle(std::fabs((*hitIter)->GetCellNormalVector().GetCosOpeningAngle(clusterDirection)));
		cosOpeningAngle = std::max(cosOpeningAngle, 0.3f); // ATTN Hard-coded default value of configurable parameter in particle id plugin                                    

		const float hitEnergy((*hitIter)->GetElectromagneticEnergy());
		energyInLayer += hitEnergy;
		nRadiationLengthsInLayer += (*hitIter)->GetNCellRadiationLengths() / cosOpeningAngle;

		const float radialDistance(((*hitIter)->GetPositionVector() - clusterIntercept).GetCrossProduct(clusterDirection).GetMagnitude());
		hitEnergyDistanceVector.push_back(HitEnergyDistance(hitEnergy, radialDistance));
	      }

	    layer90EnergySum += energyInLayer;
	    nRadiationLengthsInLayer /= static_cast<float>(iter->second->size());
	    nRadiationLengthsInLastLayer = nRadiationLengthsInLayer;
	    nRadiationLengths += nRadiationLengthsInLayer;

	    // Number of radiation lengths before cluster start                                                                                                                          
	    if (innerPseudoLayer == iLayer)
	      {
		nRadiationLengths *= static_cast<float>(innerPseudoLayer + 1 - firstPseudoLayer);
		innerLayerRadLengths = nRadiationLengths;
	      }

	    // Number of radiation lengths before longitudinal layer90                                                                                                                 
	    if (!foundLayer90 && (layer90EnergySum > 0.9f * totalElectromagneticEnergy))
	      {
		foundLayer90 = true;
		nRadiationLengths90 = nRadiationLengths;
	      }

	    // Number of radiation lengths before shower max layer                                                                                                                     
	    if (energyInLayer > maxEnergyInlayer)
	      {
		showerMaxRadLengths = nRadiationLengths;
		maxEnergyInlayer = energyInLayer;
	      }

	    // Energy above specified "high" number of radiation lengths                                                                                                               
	    if (nRadiationLengths > 40.f) // ATTN Hard-coded default value of configurable parameter in particle id plugin                                                             
	      {
		energyAboveHighRadLengths += energyInLayer;
	      }
	  }

	// Longitudinal shower profile properties                                                                                                                                      
	energyAboveHighRadLengthsFrac = (totalElectromagneticEnergy > 0.f) ? energyAboveHighRadLengths / totalElectromagneticEnergy : 0.f;

	// Transverse shower profile properties                                                                                                                                        
	std::sort(hitEnergyDistanceVector.begin(), hitEnergyDistanceVector.end(), SortHitsByDistance);
	float radial90EnergySum(0.f);
	radial90 = std::numeric_limits<float>::max();

	for (HitEnergyDistanceVector::const_iterator iter = hitEnergyDistanceVector.begin(), iterEnd = hitEnergyDistanceVector.end(); iter != iterEnd; ++iter)
	  {
	    radial90EnergySum += iter->first;
	    //	    std::cout << " GIGADEBUG radial90EnergySum " << radial90EnergySum << " totalElectromagneticEnergy " << std::endl;

	    if (radial90EnergySum > 0.9f * totalElectromagneticEnergy)
	      {
		radial90 = iter->second;
		break;
	      }
	  }

	if (radial90 == std::numeric_limits<float>::max()) {
	  std::cout << " SCZ GIGADEBUG radial90 not set radial90EnergySum " << radial90EnergySum << " totalElectromagneticEnergy " << totalElectromagneticEnergy
		    << " pCluster->GetElectromagneticEnergy() " << pCluster->GetElectromagneticEnergy() << std::endl;
	}
	
	// Explicit shower profile properties                                                                                                                                          
	showerProfileStart = pCluster->GetShowerProfileStart(*m_pPandora);
	showerProfileChi2 = pCluster->GetShowerProfileDiscrepancy(*m_pPandora);



	//        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        CaloHitList pCaloHitList;
        orderedCaloHitList.GetCaloHitList(pCaloHitList);
    
        for (CaloHitList::const_iterator hitIter = pCaloHitList.begin(), hitIterEnd = pCaloHitList.end(); hitIter != hitIterEnd; ++hitIter) {
          const CaloHit *pCaloHit = (*hitIter);
          // Determing extremal pseudolayers
          const unsigned int pseudoLayer(pCaloHit->GetPseudoLayer());
          if (pseudoLayer > lastLayer)
              lastLayer = pseudoLayer;
          if (pseudoLayer < firstLayer)
              firstLayer = pseudoLayer;
        }
      }

      double clusterEMenergyECAL  = 0.;
      double clusterEMenergyHCAL  = 0.;
      double clusterHADenergyECAL = 0.;
      double clusterHADenergyHCAL = 0.;
             
      const pandora::ClusterAddressList clusterAddressList((*itPFO)->GetClusterAddressList());
      for (pandora::ClusterAddressList::const_iterator itCluster = clusterAddressList.begin(), itClusterEnd = clusterAddressList.end(); itCluster != itClusterEnd; ++itCluster) {
        const unsigned int nHitsInCluster((*itCluster).size());
          
        for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit) {
          const reco::PFRecHitRef hgcHit(HGCRecHitHandle,recHitMap[(*itCluster)[iHit]]);               
          const DetId& detid(hgcHit->detId());
          if (!detid)
             continue;
          double eta = fabs(hgcHit->position().Eta());
          //double cos_theta = std::tanh(hgcHit->position().Eta());
          double cos_theta = 1.; // We do not need this correction because increasing thickness of absorber material compensates
                                 // for increasing path length at non-normal incidence
    
          ForwardSubdetector thesubdet = (ForwardSubdetector)detid.subdetId();
          if (thesubdet == 3) {
            int layer = (int) ((HGCEEDetId)(detid)).layer() ;
            clusterEMenergyECAL += hgcHit->energy() * cos_theta * m_calibEE.GetADC2GeV() * m_calibEE.GetEMCalib(layer,eta);
            clusterHADenergyECAL += hgcHit->energy() * cos_theta * m_calibEE.GetADC2GeV() * m_calibEE.GetHADCalib(layer,eta);
          }
          else if (thesubdet == 4) {
            int layer = (int) ((HGCHEDetId)(detid)).layer() ;
            clusterEMenergyHCAL += hgcHit->energy() * cos_theta * m_calibHEF.GetADC2GeV() * m_calibHEF.GetEMCalib(layer,eta);
            clusterHADenergyHCAL += hgcHit->energy() * cos_theta * m_calibHEF.GetADC2GeV() * m_calibHEF.GetHADCalib(layer,eta);
          }
          else if (thesubdet == 5) {
            int layer = (int) ((HGCHEDetId)(detid)).layer() ;
            clusterEMenergyHCAL += hgcHit->energy() * cos_theta * m_calibHEB.GetADC2GeV() * m_calibHEB.GetEMCalib(layer,eta);
            clusterHADenergyHCAL += hgcHit->energy() * cos_theta * m_calibHEB.GetADC2GeV() * m_calibHEB.GetHADCalib(layer,eta);
          }
          else {
          }
        }
      }
    
      ene_em = clusterEMenergyECAL + clusterEMenergyHCAL;
      ene_had = clusterHADenergyECAL + clusterHADenergyHCAL;
      
      if(i==0){
        if(debugHisto){
          sumPFOs++;
          sumClustEMEcalE += clusterEMenergyECAL;
          sumClustEMHcalE += clusterEMenergyHCAL;
          sumClustHADEcalE += clusterHADenergyECAL;
          sumClustHADHcalE += clusterHADenergyECAL;
        }
        _sumPFOEnergy += energy;
        
        if(debugPrint) std::cout << "Particle Id: " << pid << std::endl;
        if(debugPrint) std::cout << "Energy: " << energy << std::endl;
        
        const pandora::TrackAddressList trackAddressList((*itPFO)->GetTrackAddressList());
//        PANDORA_MONITORING_API(VisualizeTracks( trackAddressList  , "currentTrackList", AUTO, false, true  ) );    
//        PANDORA_MONITORING_API(ViewEvent() );
        for (pandora::TrackAddressList::const_iterator itTrack = trackAddressList.begin(), itTrackEnd = trackAddressList.end();itTrack != itTrackEnd; ++itTrack){
          reco::PFRecTrack * pftrack = (reco::PFRecTrack*)(*itTrack);
          const reco::TrackRef track =  pftrack->trackRef();
          if(debugPrint) std::cout << "Track from pfo charge " << track->charge() << std::endl;
          if(debugPrint) std::cout << "Track from pfo transverse momentum " << track->pt() << std::endl;
        }
        
        if(debugPrint) std::cout << " ENERGY  is  " << _sumPFOEnergy << std::endl;
      }
      
      const TrackList &trackList((*itPFO)->GetTrackList());
      TrackVector trackVector(trackList.begin(), trackList.end());
        
      ene_track=0; 
      
      for (TrackVector::const_iterator trackIter = trackVector.begin(), trackIterEnd = trackVector.end(); trackIter != trackIterEnd; ++trackIter) {
        const pandora::Track *pPandoraTrack = (*trackIter);
        // Extract pandora track states
        const TrackState &trackState(pPandoraTrack->GetTrackStateAtStart());
        const CartesianVector &momentum(trackState.GetMomentum());
        ene_track = momentum.GetMagnitude();
      }
      
      double diff = 1e10;  
      if(debugHisto && charge_true == charge && energy != found_energy){
        diff = abs(energy-ene_true);
        if(diff<diff_min){ // 3
         diff_min      = diff;
         ene_match     = energy;
         mass_match    = mass;
         charge_match  = charge;
         pT_match      = pT;
         pid_match     = pid;
         ene_match_em  = ene_em;
         ene_match_had = ene_had;
         ene_match_track = ene_track;
         first_layer_match = double(firstLayer*1.);
         last_layer_match  = double(lastLayer*1.);
         if(pid_match==11  && last_layer_match > 30) pid_match = -211;
         if(pid_match==-11 && last_layer_match > 30) pid_match =  211;
        } // 3
      }
    } // 4  

    if(debugHisto && ene_match>0){ // 2
     runno = iEvent.id().run();
     eventno = iEvent.id().event();
     lumi = iEvent.luminosityBlock();
     found_energy = ene_match;
     mytree->Fill();
    } // 2
  } // 1

  if(debugPrint) std::cout << " ENERGY ALL TRUE " << ene_all_true << std::endl;

  if(debugHisto){
    h2_EM_hcalEecalE->Fill(sumClustEMEcalE,sumClustEMHcalE);
    h2_Had_hcalEecalE->Fill(sumClustHADEcalE,sumClustHADHcalE);
    h_sumPfoE->Fill(_sumPFOEnergy);
    h_nbPFOs->Fill(sumPFOs);
  }
}

//Get the track siblings
TrackingParticleRefVector PandoraCMSPFCandProducer::getTpSiblings(TrackingParticleRef tp){

  if (tp.isNonnull() && tp->parentVertex().isNonnull() && !tp->parentVertex()->daughterTracks().empty()) {
    return tp->parentVertex()->daughterTracks();
  } else {
    return TrackingParticleRefVector();
  }

}
//Get the track daughters
TrackingParticleRefVector PandoraCMSPFCandProducer::getTpDaughters(TrackingParticleRef tp){

  TrackingVertexRefVector trvertexes;
  TrackingParticleRefVector trdaughter;

  if (tp.isNonnull() && tp->decayVertices().isNonnull() ) {
    trvertexes = tp->decayVertices();
    //Loop on vector of TrackingVertex objects where the TrackingParticle decays. 
    for(TrackingVertexRefVector::iterator vi = trvertexes.begin(); vi != trvertexes.end(); vi++){
      //loop on all daughter tracks 
      for(TrackingParticleRefVector::iterator di = (**vi).daughterTracks().begin(); di != (**vi).daughterTracks().end(); di++){
        //Check if the daughter is the same as our mother tp particle
        if ( (*di) == tp  ) {continue;}
        trdaughter.push_back( (*di) );
      }//end on loop over daughter
    }//end on loop over vertices
    return trdaughter;
  } else {
    return TrackingParticleRefVector();
  }

}

//make CMSSW PF objects from Pandora output
void PandoraCMSPFCandProducer::convertPandoraToCMSSW(const edm::Handle<reco::PFRecTrackCollection>& trackh,
                                                     const edm::Handle<reco::PFRecHitCollection>& rechith, 
                                                     edm::Event& iEvent)
{
  const auto& pfrectracks = *trackh;
  const auto& pfrechits = *rechith;

  const pandora::PfoList *pPfoList = NULL;
  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pPandora, pPfoList));


  // SCZ needed matching collections
  edm::Handle<std::vector<PCaloHit> > HGCHitsEE;
  edm::Handle<std::vector<PCaloHit> > HGCHitsHEback;
  edm::Handle<std::vector<PCaloHit> > HGCHitsHEfront;
  iEvent.getByLabel(inputTagSimHGCHitsEE_,HGCHitsEE);
  iEvent.getByLabel(inputTagSimHGCHitsHEback_,HGCHitsHEback);
  iEvent.getByLabel(inputTagSimHGCHitsHEfront_,HGCHitsHEfront);
  std::vector<std::vector<PCaloHit> > theHGCHits;
  theHGCHits.resize(3);
  theHGCHits[0].insert(theHGCHits[0].end(),HGCHitsEE->begin(),HGCHitsEE->end());
  theHGCHits[1].insert(theHGCHits[1].end(),HGCHitsHEfront->begin(),HGCHitsHEfront->end());
  theHGCHits[2].insert(theHGCHits[2].end(),HGCHitsHEback->begin(),HGCHitsHEback->end());
  edm::Handle<std::vector<SimTrack> > simTk;
  edm::Handle<std::vector<SimVertex> > simVtx;
  iEvent.getByLabel(inputTagSimTracks_,simTk);
  iEvent.getByLabel(inputTagSimTracks_,simVtx);
  std::vector<SimTrack> theSimTracks;
  std::vector<SimVertex> theSimVertexes;
  theSimTracks.insert(theSimTracks.end(),simTk->begin(),simTk->end());
  theSimVertexes.insert(theSimVertexes.end(),simVtx->begin(),simVtx->end());


  edm::Handle<std::vector<reco::GenParticle> > genpart;
  iEvent.getByLabel(inputTagGenParticles_,genpart);

  edm::Handle<std::vector<int> > genBarcodes;
  iEvent.getByLabel(inputTagGenParticles_,genBarcodes);


  // SCZ needed geometry for matching
  std::vector<const HGCalGeometry*> geom;
  geom.resize(3);
  geom[0] = hgceeGeoHandle.product();
  geom[1] = hgchefGeoHandle.product();
  geom[2] = hgchebGeoHandle.product();

  std::unordered_map<int,int> geantTrackIdToPdgId(theSimTracks.size());
  std::unordered_map<int,float> geantTrackIdToSimTrackEnergy(theSimTracks.size());
  std::unordered_map<int,float> geantTrackIdToStartingRadius(theSimTracks.size());

  for (unsigned int isimtk = 0; isimtk < theSimTracks.size(); isimtk++){
    int gid = (int)theSimTracks[isimtk].trackId();
    geantTrackIdToSimTrackEnergy[gid] = theSimTracks[isimtk].trackerSurfaceMomentum().E();
    geantTrackIdToPdgId[gid] = theSimTracks[isimtk].type();
    //    unsigned int vi = theSimTracks[isimtk].vertIndex();
    //    const math::XYZTLorentzVectorD start_pos = theSimVertexes[vi].position();
    //    geantTrackIdToStartingRadius[gid] = start_pos.R();
    geantTrackIdToStartingRadius[gid] = -1.;
    for (unsigned int igen = 0 ; igen < genBarcodes->size() ; igen++) {
      if (genBarcodes->at(igen) == theSimTracks[isimtk].genpartIndex()) {
	geantTrackIdToStartingRadius[gid] = genpart->at(igen).vertex().rho();
      }
    }

    //    std::cout << " SCZ GIGADEBUG gid " << gid << " r " << geantTrackIdToStartingRadius[gid] << std::endl;
  }

  // size is just an upper bound here
  std::unordered_multimap<uint32_t,std::pair<unsigned,unsigned> > recoDetIdToSimHitIndex(theHGCHits[0].size() + theHGCHits[1].size() + theHGCHits[2].size()); 

  //  for(std::vector<PCaloHit>::iterator hit_it = theHGCHits[i].begin(); hit_it != theHGCHits[i].end(); ++hit_it)
  for (unsigned int i = 0 ; i < 3 ; i++) {

    uint32_t mySubDet(ForwardSubdetector::HGCEE);                                                                                                                          
    if(i==1) mySubDet=ForwardSubdetector::HGCHEF;                                                                                                                          
    else if(i==2) mySubDet=ForwardSubdetector::HGCHEB;                                                                                                                     
    const HGCalTopology &topo=geom[i]->topology();                                                                                                                         
    const HGCalDDDConstants &dddConst=topo.dddConstants();                                                                                                                 

    for (unsigned int j = 0 ; j < theHGCHits[i].size() ; j++) {
      //gang SIM->RECO cells to get final layer assignment
      HGCalDetId simId(theHGCHits[i][j].id());
      int layer(simId.layer()),cell(simId.cell());
      std::pair<int,int> recoLayerCell=dddConst.simToReco(cell,layer,topo.detectorType());
      cell  = recoLayerCell.first;
      layer = recoLayerCell.second;
      if(layer < 0) continue;
      
      uint32_t recoDetId = ( (i==0) ?
			     (uint32_t)HGCEEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell) :
			     (uint32_t)HGCHEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell)
			     );
      //      std::cout << " SCZ GIGADEBUG insert recoDetId det simhit trackid " << recoDetId << " " << i << " " << j << " " << theHGCHits[i][j].geantTrackId() <<  std::endl;
      //      if (theHGCHits[i][j].geantTrackId() < 1) std::cout << " SCZ MEGADEBUG THE ABOVE HAS trackid 0" << std::endl;
      recoDetIdToSimHitIndex.insert(std::make_pair(recoDetId,std::make_pair(i,j)));


    }
  }


  //first loop: make clusters
  std::auto_ptr<reco::PFClusterCollection> clusters(new reco::PFClusterCollection);
  std::unordered_multimap<unsigned,unsigned> pfos_to_clusters; // since pfos can be many to one pfos
  std::unordered_multimap<unsigned,unsigned> pfos_to_tracks; // ditto for tracks
  auto pfoBegin = pPfoList->cbegin();
  unsigned cluster_idx = 0;
  for (auto itPFO = pfoBegin; itPFO != pPfoList->cend(); ++itPFO){
    size_t pfo_idx = std::distance(pfoBegin,itPFO);
    //loop over pandora clusters
    const ClusterList &clusterList((*itPFO)->GetClusterList());
    //ClusterVector clusterVector(clusterList.begin(), clusterList.end());
    
    // keep tabs on if this PFO was EM or HAD
    const bool pfoIsEM = ( 22 == (*itPFO)->GetParticleId() || 11 == std::abs((*itPFO)->GetParticleId()) ) ;

    std::cout << " SCZ MEGADEBUG PFO ID " << (*itPFO)->GetParticleId() << " " << " clusterList.size() " << clusterList.size()
              << " (*itPFO)->GetTrackList().size() " << (*itPFO)->GetTrackList().size() << std::endl;

    assert (clusterList.size()<=1); // SCZ
    assert ((*itPFO)->GetTrackList().size()<=1); // SCZ
    
    for (auto clusterIter = clusterList.cbegin(); clusterIter != clusterList.cend(); ++clusterIter){
      // keep track of clusters used by the PFOs
      
      pfos_to_clusters.emplace(pfo_idx,cluster_idx++);
      const Cluster *pCluster = (*clusterIter);
      reco::PFCluster temp;
      
      if( pCluster->GetOrderedCaloHitList().size() == 0 ) {
	throw cms::Exception("EmptyClusterOnPFO")
	  << " cluster is empty, wtf. " ;
      }

      // setup basic energy determination
      temp.setEmEnergy(pCluster->GetElectromagneticEnergy());
      temp.setHadEnergy(pCluster->GetHadronicEnergy());
      if( pfoIsEM ) {
        temp.setEnergy(temp.emEnergy());
      } else {
        temp.setEnergy(temp.hadEnergy());
      }

      pfPdgId = (*itPFO)->GetParticleId(); // SCZ
      energy = temp.energy(); // SCZ
      trackSeeded = (int)pCluster->IsTrackSeeded(); // SCZ

      const TrackList &trackList((*itPFO)->GetTrackList());
      TrackVector trackVector(trackList.begin(), trackList.end());
      assert(trackVector.size() <= 1);

      trackP = -1.;
      for (TrackVector::const_iterator trackIter = trackVector.begin(), trackIterEnd = trackVector.end(); trackIter != trackIterEnd; ++trackIter) {
        const pandora::Track *pPandoraTrack = (*trackIter);
        // Extract pandora track states
	const TrackState &trackState(pPandoraTrack->GetTrackStateAtStart());
        const CartesianVector &momentum(trackState.GetMomentum());
        trackP = momentum.GetMagnitude(); // SCZ
      }

      // set cluster axis and position information
      const auto& pandoraAxis = pCluster->GetInitialDirection();
      math::XYZVector axis(pandoraAxis.GetX(),pandoraAxis.GetY(),pandoraAxis.GetZ());
      temp.setAxis(axis);
      const auto& clusterFit = pCluster->GetFitToAllHitsResult();
      if( clusterFit.IsFitSuccessful() ) {
        const auto& pandoraPos = clusterFit.GetIntercept();
        math::XYZPoint pos( pandoraPos.GetX(), pandoraPos.GetY(), pandoraPos.GetZ() );
        temp.setPosition(pos);
      } else {
        const auto& pandoraPos = pCluster->GetCentroid(pCluster->GetInnerPseudoLayer());
        math::XYZPoint pos( pandoraPos.GetX(), pandoraPos.GetY(), pandoraPos.GetZ() );
        temp.setPosition(pos);
      }

      eta = temp.eta(); // SCZ
      phi = temp.phi(); // SCZ

      if( pCluster->IsTrackSeeded() ) {
        const auto* pandoraTrack = pCluster->GetTrackSeed();
        auto iter = recTrackMap.find(pandoraTrack->GetParentTrackAddress());
        if( iter != recTrackMap.end() ) {
          temp.setTrack(pfrectracks[iter->second].trackRef());
          pfos_to_tracks.emplace(pfo_idx,iter->second);
        } else {
          throw cms::Exception("TrackUsedButNotFound")
            << "Track used in PandoraPFA was not found in the original input track list!";
        }
      }
      
      ngp = 0;
      sumEneRecHit = 0.;
      sumEneGp = 0.;

      std::map<int,float> geantTrackIdToSimHitEnergy;

      //loop over calo hits in cluster
      const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());      
      if( orderedCaloHitList.size() ) {
        const auto& firstlayer = *(orderedCaloHitList.begin());
        const auto* firsthit = *(firstlayer.second->begin());
        auto iter = recHitMap.find(firsthit->GetParentCaloHitAddress());
        if( iter != recHitMap.end() ) {
          temp.setLayer(pfrechits.at(iter->second).layer());
	}
      } else {
	throw cms::Exception("TrackUsedButNotFound")
	  << "Hit used in PandoraPFA was not found in the original input hit list!";
      }

      for (auto layerIter = orderedCaloHitList.begin(), layerIterEnd = orderedCaloHitList.end(); layerIter != layerIterEnd; ++layerIter) {
        const CaloHitList& hits_in_layer = *(layerIter->second);
        for( auto hitIter = hits_in_layer.cbegin(), hitIterEnd = hits_in_layer.cend(); hitIter != hitIterEnd; ++hitIter ) {
          auto hit_index = recHitMap.find((*hitIter)->GetParentCaloHitAddress());
          if( hit_index != recHitMap.end() ) {
	    reco::PFRecHitRef ref(rechith,hit_index->second);
            temp.addRecHitFraction(reco::PFRecHitFraction(ref,(*hitIter)->GetWeight()));
	    auto currentRecoDetId = pfrechits.at(hit_index->second).detId();
	    int layer = pfrechits.at(hit_index->second).layer();
	    float currentRecoDetEnergy = pfrechits.at(hit_index->second).energy();
	    const DetId& detid(currentRecoDetId);

	    if (detid.subdetId() == 3) currentRecoDetEnergy *= (m_calibEE.GetEMCalib(layer,eta) * m_calibEE.GetADC2GeV());
	    if (detid.subdetId() == 4) currentRecoDetEnergy *= (m_calibHEF.GetEMCalib(layer,eta) * m_calibHEF.GetADC2GeV());
	    if (detid.subdetId() == 5) currentRecoDetEnergy *= (m_calibHEB.GetEMCalib(layer,eta) * m_calibHEB.GetADC2GeV());

	    sumEneRecHit += currentRecoDetEnergy;

	    //	    std::cout << " SCZ GIGADEBUG rechit layer " << layer << " x " << pfrechits.at(hit_index->second).position().x()
	    //		      << " y " << pfrechits.at(hit_index->second).position().y()
	    //		      << " z " << pfrechits.at(hit_index->second).position().z() << std::endl;


	    //	    std::pair<unsigned,unsigned> simHitIndices = recoDetIdToSimHitIndex.
	    //	      auto range = my_map.equal_range(reco det id);
	    auto range = recoDetIdToSimHitIndex.equal_range(currentRecoDetId);
	    //	    std::cout << "   SCZ GIGADEBUG just called equal_range on " << currentRecoDetId << std::endl;

	    for (auto index_iter = range.first ; index_iter != range.second ; index_iter++) {
	      //	      std::cout << " SCZ GIGADEBUG " << (index_iter->first) << " " << (index_iter->second.first) << " " << (index_iter->second.first) << std::endl;
	      auto my_hit = theHGCHits[index_iter->second.first][index_iter->second.second];
	      float currentSimHitEnergy = my_hit.energy();                                                                                                                    
	      

	      if (index_iter->second.first == 0) currentSimHitEnergy *= (m_calibEE.GetEMCalib(layer,eta) * m_calibEE.GetADC2GeV());
	      if (index_iter->second.first == 1) currentSimHitEnergy *= (m_calibHEF.GetEMCalib(layer,eta) * m_calibHEF.GetADC2GeV());
	      if (index_iter->second.first == 2) currentSimHitEnergy *= (m_calibHEB.GetEMCalib(layer,eta) * m_calibHEB.GetADC2GeV());                                                                 
	      sumEneRecHit += currentRecoDetEnergy;                                                                                                                            
	      
	      int gid = my_hit.geantTrackId();                                                                                                                                
	      sumEneGp += currentSimHitEnergy;                                                                                                                                 
	      if (geantTrackIdToSimHitEnergy.count(gid)) {                                                                                                                     
		geantTrackIdToSimHitEnergy[gid] += currentSimHitEnergy;                                                                                                        
	      } else {                                                                                                                                                         
		ngp++;
		geantTrackIdToSimHitEnergy[gid] = currentSimHitEnergy;                                                                                                         
	      }
	    }

	    // SCZ big loop
	    /*
	    for (unsigned int i = 0 ; i < theHGCHits.size() ; i++) {
	      uint32_t mySubDet(ForwardSubdetector::HGCEE);
	      if(i==1) mySubDet=ForwardSubdetector::HGCHEF;
	      else if(i==2) mySubDet=ForwardSubdetector::HGCHEB;
	      
	      const HGCalTopology &topo=geom[i]->topology();
	      const HGCalDDDConstants &dddConst=topo.dddConstants();
	      
	      for(std::vector<PCaloHit>::iterator hit_it = theHGCHits[i].begin(); hit_it != theHGCHits[i].end(); ++hit_it)
		{
		  //gang SIM->RECO cells to get final layer assignment                                                                                                                         
		  HGCalDetId simId(hit_it->id());
		  int layer(simId.layer()),cell(simId.cell());
		  std::pair<int,int> recoLayerCell=dddConst.simToReco(cell,layer,topo.detectorType());
		  cell  = recoLayerCell.first;
		  layer = recoLayerCell.second;
		  if(layer < 0) continue;
		  

		  uint32_t recoDetId = ( (i==0) ?
					 (uint32_t)HGCEEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell) :
					 (uint32_t)HGCHEDetId(ForwardSubdetector(mySubDet),simId.zside(),layer,simId.sector(),simId.subsector(),cell)
					 );
		  //		std::cout << "   SCZ GIGADEBUG has recoDetId " << recoDetId << " and energy " << hit_it->energy() << std::endl;
		  if (recoDetId == currentRecoDetId) {
		    float currentSimHitEnergy = hit_it->energy();

		    if (i == 0) currentSimHitEnergy *= (m_calibEE.GetEMCalib(layer,eta) * m_calibEE.GetADC2GeV());
                    if (i == 1) currentSimHitEnergy *= (m_calibHEF.GetEMCalib(layer,eta) * m_calibHEF.GetADC2GeV());
                    if (i == 2) currentSimHitEnergy *= (m_calibHEB.GetEMCalib(layer,eta) * m_calibHEB.GetADC2GeV());
		    sumEneRecHit += currentRecoDetEnergy;

		    int gid = hit_it->geantTrackId();
		    sumEneGp += currentSimHitEnergy;
		    if (geantTrackIdToSimHitEnergy.count(gid)) {
		      geantTrackIdToSimHitEnergy[gid] += currentSimHitEnergy;
		    } else {
		      ngp++;
		      geantTrackIdToSimHitEnergy[gid] = currentSimHitEnergy;
		      geantTrackIdToSimTrackEnergy[gid] = -1.;
		      geantTrackIdToPdgId[gid] = 0;
		      for (unsigned int isimtk = 0; isimtk < theSimTracks.size(); isimtk++){
			if ((int)theSimTracks[isimtk].trackId() == gid) {
			  geantTrackIdToSimTrackEnergy[gid] = theSimTracks[isimtk].trackerSurfaceMomentum().E();
			  geantTrackIdToPdgId[gid] = theSimTracks[isimtk].type();
			}
		      }
		    }
		  } 
		} // loop over sim hits
	    } // loop over HGC detectors
	    */

	  } else {
	    throw cms::Exception("TrackUsedButNotFound")
	      << "Hit used in PandoraPFA was not found in the original input hit list!";
	  }
        } // loop over hits in layer
      } // loop over layers
      clusters->push_back(temp);

      // SCZ: we have these now                                                                                                                                                
      //    map<int,int> geantTrackIdToPdgId;                                                                                                                                  
      //    map<int,float> geantTrackIdToSimHitEnergy;                                                                                                                         
      //    map<int,float> geantTrackIdToSimTrackEnergy;                                                                                                                       

      std::multimap<float, int> simHitEnergyToGeantTrackEnergy  = flip_map(geantTrackIdToSimHitEnergy);

      std::cout << " SCZ MEGADEBUG ngp=" << ngp << " energy=" << energy << " sumEneGp=" << sumEneGp << " sumEneRecHit=" << sumEneRecHit << std::endl;
      int counter = 0;

      for(std::multimap<float, int>::reverse_iterator pos = simHitEnergyToGeantTrackEnergy.rbegin(); pos != simHitEnergyToGeantTrackEnergy.rend(); ++pos) {
	gpEneFrac[counter] = pos->first / sumEneGp; // SCZ changed here denominator from energy
	gpEneTotal[counter] = geantTrackIdToSimTrackEnergy[pos->second];
	gpPdgId[counter] = geantTrackIdToPdgId[pos->second];
	gpStartRadius[counter] = geantTrackIdToStartingRadius[pos->second];
	std::cout << "     SCZ GIGADEBUG " << counter++ << " " << pos->first << " " << pos->second << " " << geantTrackIdToSimHitEnergy[pos->second]
		  << " " << geantTrackIdToSimTrackEnergy[pos->second] << " " << geantTrackIdToPdgId[pos->second] << std::endl;
      }
      
      std::cout << " SCZ MEGADEBUG Fill clusterTree" << std::endl;
      if (debugHisto) clusterTree->Fill();
      std::cout << " SCZ MEGADEBUG Filled clusterTree" << std::endl;

    }//end clusters
    // tail catcher for the tracks on the pfo
    for( const auto* pandoraTrack : (*itPFO)->GetTrackList() ) {
      auto iter = recTrackMap.find(pandoraTrack->GetParentTrackAddress());
      if( iter != recTrackMap.end() ) {
	pfos_to_tracks.emplace(pfo_idx,iter->second);
      } else {
	throw cms::Exception("TrackUsedButNotFound")
	  << "Track used in PandoraPFA was not found in the original input track list!";
      }
    }
  }//end pfos
  // put the clusters in the event and get the handle back
  edm::OrphanHandle<PFClusterCollection> clusterHandle = iEvent.put(clusters);
  
  //make blocks (no need for maps here since blocks are 1:1 to pandora candidates
  std::auto_ptr<reco::PFBlockCollection> pfblocks(new reco::PFBlockCollection);
  for (auto itPFO = pPfoList->cbegin(); itPFO != pPfoList->cend(); ++itPFO){
    const size_t pfo_idx = std::distance(pfoBegin,itPFO);
    auto tk_range =  pfos_to_tracks.equal_range(pfo_idx);
    auto clus_range =  pfos_to_clusters.equal_range(pfo_idx);
    // setup block to add to
    pfblocks->emplace_back( reco::PFBlock() );
    reco::PFBlock& block = pfblocks->back();
    // process tracks
    for( auto tk_itr = tk_range.first; tk_itr != tk_range.second; ++tk_itr ) {
      reco::PFRecTrackRef pftrackref = reco::PFRecTrackRef(trackh,tk_itr->second);  
      std::unique_ptr<reco::PFBlockElementTrack> tk_elem( new reco::PFBlockElementTrack( pftrackref ) ); 
      block.addElement(tk_elem.get());
    }
    // process clusters
    for( auto clus_itr = clus_range.first; clus_itr != clus_range.second; ++clus_itr ) {
      reco::PFClusterRef pfclusterref = reco::PFClusterRef(clusterHandle,clus_itr->second);  
      std::unique_ptr<reco::PFBlockElementCluster> clus_elem( nullptr );
      switch( pfclusterref->layer() ) {
      case PFLayer::HGC_ECAL:
        clus_elem.reset( new reco::PFBlockElementCluster( pfclusterref, reco::PFBlockElement::HGC_ECAL ) );
        break;
      case PFLayer::HGC_HCALF:
        clus_elem.reset( new reco::PFBlockElementCluster( pfclusterref, reco::PFBlockElement::HGC_HCALF ) );
        break;
      case PFLayer::HGC_HCALB:
        clus_elem.reset( new reco::PFBlockElementCluster( pfclusterref, reco::PFBlockElement::HGC_HCALB ) );
        break;
      default:
        throw cms::Exception("PoorlyDefinedCluster")
          << pfclusterref->layer() << " PFCluster from PandoraPFA does not have an assigned layer in HGC!";
      }
      block.addElement(clus_elem.get());
    }
  }
  edm::OrphanHandle<PFBlockCollection> blockHandle = iEvent.put(pfblocks);
  
  //make candidates
  std::auto_ptr<reco::PFCandidateCollection> pandoraCands(new reco::PFCandidateCollection);
  std::auto_ptr<reco::PFCandidateCollection> pandoraElectronCands(new reco::PFCandidateCollection); // SZ Feb 25 not filled yet
  for (auto itPFO = pPfoList->cbegin(); itPFO != pPfoList->cend(); ++itPFO){ 
    const unsigned pfo_idx = std::distance(pfoBegin,itPFO);
    const pandora::ParticleFlowObject& the_pfo = *(*itPFO);
    // translate pandora pdgids to pf pi+/-, photon, k_long notation
    const int pandoraPID = the_pfo.GetParticleId();
    const int pandoraCharge = the_pfo.GetCharge();
    reco::PFCandidate::ParticleType cmspfPID = reco::PFCandidate::X;
    if( 22 == pandoraPID ) {
      cmspfPID = reco::PFCandidate::gamma;
    } else if( 0 == pandoraCharge ) { // anything that has no charge and isn't a photon is a neutral hadron
      cmspfPID = reco::PFCandidate::h0;
    } else if( 0 != pandoraCharge ) {
      // LG disregarding electrons for now, ID is funky 20 Feb, 2015
      if( 13 == std::abs(pandoraPID) ) {
        cmspfPID = reco::PFCandidate::mu;
      } else { // anything that's not a muon and has charge is a h+
        cmspfPID = reco::PFCandidate::h;
      }
    } else {
      throw cms::Exception("StrangeParticleID")
        << "PandoraPID " << pandoraPID << " appears to not be in the set of real numbers!";
    }    
    // get the p4
    const float pandoraE = the_pfo.GetEnergy();
    const pandora::CartesianVector& pandoraP3 = the_pfo.GetMomentum();
    math::XYZTLorentzVector thep4(pandoraP3.GetX(),pandoraP3.GetY(),pandoraP3.GetZ(),pandoraE);
    
    pandoraCands->emplace_back( pandoraCharge, thep4, cmspfPID  );
    reco::PFCandidate& cand = pandoraCands->back();
    
    // set track as well as raw and corrected calorimeter energies    
    
    // setup any remaining information about what's in this PF Candidate
    reco::PFBlockRef blockref(blockHandle,pfo_idx);   
    reco::TrackRef hardest_track;
    const edm::OwnVector<reco::PFBlockElement>& elements = blockref->elements();
    for( unsigned elem_idx = 0; elem_idx < elements.size(); ++elem_idx ) {
      if( reco::PFBlockElement::TRACK == elements[elem_idx].type() ) {
	if( hardest_track.isNull() || hardest_track->pt() < elements[elem_idx].trackRef()->pt() ) {
	  hardest_track = elements[elem_idx].trackRef();
	}
      }
      cand.addElementInBlock(blockref,elem_idx);
    }
    if( hardest_track.isNonnull() && pandoraCharge != 0) cand.setTrackRef(hardest_track);   
  }
  iEvent.put(pandoraCands);
  iEvent.put(pandoraElectronCands,electronOutputCol_); // SZ Feb 25 not filled yet
}



// ------------ method called once each job just before starting event loop  ------------
void PandoraCMSPFCandProducer::beginJob()
{   
  // setup our maps for processing
  recHitMap.clear();
  recTrackMap.clear();

  const char *pDisplay(::getenv("DISPLAY"));
  if (NULL == pDisplay) {
    if(debugPrint) std::cout << "DISPLAY environment not set" << std::endl;
  }  else {
    if(debugPrint) std::cout << "DISPLAY environment set to " << pDisplay << std::endl;
  }
  int argc = 0;
  char* argv = (char *)"";
  TApplication *m_pApplication;
  m_pApplication = gROOT->GetApplication();
  if(debugPrint) std::cout << "In PandoraCMSPFCandProducer::beginJob gVirtualX->GetDisplay()" << gVirtualX->GetDisplay() << std::endl;
  if(!m_pApplication){
    if(debugPrint) std::cout << "In if of m_pApplication in PandoraCMSPFCandProducer::beginJob " << std::endl;
    m_pApplication = new TApplication("PandoraMonitoring", &argc, &argv);
  } 
// END AP

  if(debugHisto){
    file = new TFile(_outputFileName.c_str(),"recreate");
    
    const bool oldAddDir = TH1::AddDirectoryStatus(); 
    TH1::AddDirectory(true); 
    
    /*
    h_sumCaloE = new TH1F("sumCaloE","sum hit E in Calos",1000,0,400);
    h_sumCaloEM = new TH1F("sumCaloEM","sum hit E in Calos",1000,0,400);
    h_sumCaloHad = new TH1F("sumCaloHad","sum hit E in Calos",1000,0,400);
    h_simDir_sumCaloEM  = new TH1F("sumCaloEMsimDir","sum hit E in Calos",1000,0,400);
    h_simDir_sumCaloHad = new TH1F("sumCaloHadsimDir","sum hit E in Calos",1000,0,400);
    
    h_sumEcalEEM   = new TH1F("sumEcalEEM","sum hit EM E in Ecal",1000,0,350);
    h_sumHcalEEM   = new TH1F("sumHcalEEM","sum hit EM E in Hcal",1000,0,350);
    h_sumEcalEHad  = new TH1F("sumEcalEHad","sum hit Had E in Ecal",1000,0,400);
    h_sumEcalEHad_unc = new TH1F("sumEcalEHad_unc","sum hit Had E in Ecal",1000,0,400);
    h_sumHcalEHad  = new TH1F("sumHcalEHad","sum hit Had E in Hcal",1000,0,400);
    h_sumHcalEHadc = new TH1F("sumHcalEHadc","sum hit Had E in Hcal",1000,0,400);
    h_sumEcalEHadc = new TH1F("sumEcalEHadc","sum hit Had E in Hcal",1000,0,400);
    h_sumEHad      = new TH1F("sumEHad","sum hit Had E ",1000,0,400);
    
    h_sumPfoE = new TH1F("hsumPfoE","sumPFOenergy",1000,0.,1000.);
    h_nbPFOs = new TH1F("hnbPfos","nb of rec PFOs",30,0.,30.);
    
    h2_Calo_EM_hcalEecalE = new TH2F("CalohcalEecalEem","",1000,0,100.0,1000,0,100.0);
    h2_Calo_Had_hcalEecalE = new TH2F("CalohcalEecalEhad","",1000,0,100.0,1000,0,100.0);
    
    h2_EM_hcalEecalE = new TH2F("hcalEecalEem","",1000,0,400,1000,0,400);
    h2_Had_hcalEecalE = new TH2F("hcalEecalEhad","",1000,0,400,400,0,400);
    
    h_MCp_Eta = new TH1F("MCp_Eta","MCp_Eta",300,-3.5,3.5);
    h_MCp_Phi = new TH1F("MCp_Phi","MCp_Phi",300,-3.5,3.5);
    h_hit_Eta = new TH1F("hit_Eta","hit_Eta",300,-3.5,3.5);
    h_hit_Phi = new TH1F("hit_Phi","hit_Phi",300,-3.5,3.5);
    
    h_hitEperLayer_EM[ForwardSubdetector::HGCEE]  = new TH1F("hitEperLayer_EM_EE"  ,"sum hit E per layer EE",100,-0.5,99.5);
    h_hitEperLayer_EM[ForwardSubdetector::HGCHEF] = new TH1F("hitEperLayer_EM_HEF" ,"sum hit E per layer HEF",100,-0.5,99.5);
    h_hitEperLayer_EM[ForwardSubdetector::HGCHEB] = new TH1F("hitEperLayer_EM_HEB" ,"sum hit E per layer HEB",100,-0.5,99.5);
    h_hitEperLayer_HAD[ForwardSubdetector::HGCEE]  = new TH1F("hitEperLayer_HAD_EE"  ,"sum hit E per layer EE",100,-0.5,99.5);
    h_hitEperLayer_HAD[ForwardSubdetector::HGCHEF] = new TH1F("hitEperLayer_HAD_HEF" ,"sum hit E per layer HEF",100,-0.5,99.5);
    h_hitEperLayer_HAD[ForwardSubdetector::HGCHEB] = new TH1F("hitEperLayer_HAD_HEB" ,"sum hit E per layer HEB",100,-0.5,99.5);
    
    h_MIP[ForwardSubdetector::HGCEE]  = new TH1F("MIPEE" ,"Mip in EE ",1000,0,10);
    h_MIP[ForwardSubdetector::HGCHEF] = new TH1F("MIPHEF","Mip in HEF",1000,0,10);
    h_MIP[ForwardSubdetector::HGCHEB] = new TH1F("MIPHEB","Mip in HEB",1000,0,10);
    
    h_MIP_Corr[ForwardSubdetector::HGCEE]  = new TH1F("MIPCorrEE" ,"Mip corrected in EE ",1000,0,10);
    h_MIP_Corr[ForwardSubdetector::HGCHEF] = new TH1F("MIPCorrHEF","Mip corrected in HEF",1000,0,10);
    h_MIP_Corr[ForwardSubdetector::HGCHEB] = new TH1F("MIPCorrHEB","Mip corrected in HEB",1000,0,10);
    
    mytree = new TTree("mytree","mytree");
    mytree->Branch("ene_true",&ene_true);
    mytree->Branch("mass_true",&mass_true);
    mytree->Branch("pT_true",&pT_true);
    mytree->Branch("charge_true",&charge_true);
    mytree->Branch("pid_true",&pid_true);
    mytree->Branch("ene_match",&ene_match);
    mytree->Branch("ene_match_em",&ene_match_em);
    mytree->Branch("ene_match_had",&ene_match_had);
    mytree->Branch("ene_match_track",&ene_match_track);
    mytree->Branch("mass_match",&mass_match);
    mytree->Branch("pT_match",&pT_match);
    mytree->Branch("charge_match",&charge_match);
    mytree->Branch("pid_match",&pid_match);
    mytree->Branch("isDecBefCal",&isDecBefCal);
    mytree->Branch("first_layer_match",&first_layer_match);
    mytree->Branch("last_layer_match",&last_layer_match);
    mytree->Branch("runno",&runno);
    mytree->Branch("eventno",&eventno);
    mytree->Branch("lumi",&lumi);
    mytree->Branch("nbPFOs",&nbPFOs);
    */

    /*
  # define CLMAX 100
  TTree *gpTree;
  float gpTree_energy, gpTree_eta, gpTree_phi, gpTree_startRadius;
  int gpTree_pdgId;
  int gpTree_ncl;
  int gpTree_clTrackSeeded[CLMAX];
  float gpTree_clTotalEnergy[CLMAX];
  float gpTree_clDeposit[CLMAX];
  float gpTree_mipFraction[CLMAX];
  float gpTree_dCosR[CLMAX];
  float gpTree_clusterRms[CLMAX];
  float gpTree_innerLayerRadLengths[CLMAX];
  float gpTree_nRadiationLengths90[CLMAX];
  float gpTree_showerMaxRadLengths[CLMAX];
  float gpTree_energyAboveHighRadLengths[CLMAX];
  float gpTree_energyAboveHighRadLengthsFrac[CLMAX];
  float gpTree_radial90[CLMAX];
  float gpTree_showerProfileStart[CLMAX];
  float gpTree_showerProfileChi2[CLMAX];
    */


    gpTree = new TTree("gpTree","gpTree");
    gpTree->Branch("energy",&gpTree_energy);
    gpTree->Branch("eta",&gpTree_eta);
    gpTree->Branch("phi",&gpTree_phi);
    gpTree->Branch("pdgId",&gpTree_pdgId);
    gpTree->Branch("startRadius",&gpTree_startRadius);
    gpTree->Branch("ncl",&gpTree_ncl);
    gpTree->Branch("clTrackSeeded",gpTree_clTrackSeeded,"clTrackSeeded[ncl]/I");
    gpTree->Branch("clDeposit",gpTree_clDeposit,"clDeposit[ncl]/F");
    gpTree->Branch("clTotalEnergy",gpTree_clTotalEnergy,"clTotalEnergy[ncl]/F");
    gpTree->Branch("mipFraction",gpTree_mipFraction,"mipFraction[ncl]/F");
    gpTree->Branch("dCosR",gpTree_dCosR,"dCosR[ncl]/F");
    gpTree->Branch("clusterRms",gpTree_clusterRms,"clusterRms[ncl]/F");
    gpTree->Branch("innerLayerRadLengths",gpTree_innerLayerRadLengths,"innerLayerRadLengths[ncl]/F");
    gpTree->Branch("showerMaxRadLengths",gpTree_showerMaxRadLengths,"showerMaxRadLengths[ncl]/F");
    gpTree->Branch("energyAboveHighRadLengths",gpTree_energyAboveHighRadLengths,"energyAboveHighRadLengths[ncl]/F");
    gpTree->Branch("energyAboveHighRadLengthsFrac",gpTree_energyAboveHighRadLengthsFrac,"energyAboveHighRadLengthsFrac[ncl]/F");
    gpTree->Branch("radial90",gpTree_radial90,"radial90[ncl]/F");
    gpTree->Branch("showerProfileStart",gpTree_showerProfileStart,"showerProfileStart[ncl]/F");
    gpTree->Branch("showerProfileChi2",gpTree_showerProfileChi2,"showerProfileChi2[ncl]/F");


    /*
    TTree *clusterTree;
    float energy, trackP, eta, phi;
    int pfPdgId, trackSeeded;
    int ngp;
    floatgpEneFrac[GPMAX];
    floatgpPdgId[GPMAX];
    */

    std::cout << " SCZ MEGADEBUG new clusterTree" << std::endl;
    clusterTree = new TTree("clusterTree","clusterTree");
    clusterTree->Branch("energy",&energy);
    clusterTree->Branch("sumEneGp",&sumEneGp);
    clusterTree->Branch("sumEneRecHit",&sumEneRecHit);
    clusterTree->Branch("trackP",&trackP);
    clusterTree->Branch("eta",&eta);
    clusterTree->Branch("phi",&phi);
    clusterTree->Branch("pfPdgId",&pfPdgId);
    clusterTree->Branch("trackSeeded",&trackSeeded);
    clusterTree->Branch("ngp",&ngp);
    clusterTree->Branch("gpEneFrac",gpEneFrac,"gpEneFrac[ngp]/F");
    clusterTree->Branch("gpEneTotal",gpEneTotal,"gpEneTotal[ngp]/F");
    clusterTree->Branch("gpPdgId",gpPdgId,"gpPdgId[ngp]/F");
    clusterTree->Branch("gpStartRadius",gpStartRadius,"gpStartRadius[ngp]/F");
    clusterTree->Branch("mipFraction",&mipFraction);
    clusterTree->Branch("dCosR",&dCosR);
    clusterTree->Branch("clusterRms",&clusterRms);
    clusterTree->Branch("innerLayerRadLengths",&innerLayerRadLengths);
    clusterTree->Branch("showerMaxRadLengths",&showerMaxRadLengths);
    clusterTree->Branch("energyAboveHighRadLengths",&energyAboveHighRadLengths);
    clusterTree->Branch("energyAboveHighRadLengthsFrac",&energyAboveHighRadLengthsFrac);
    clusterTree->Branch("radial90",&radial90);
    clusterTree->Branch("showerProfileStart",&showerProfileStart);
    clusterTree->Branch("showerProfileChi2",&showerProfileChi2);

    
    TH1::AddDirectory(oldAddDir); 
    

    // SCZ NEIN
    /*
    m_hitEperLayer_EM[ForwardSubdetector::HGCEE]  = new double[100];
    m_hitEperLayer_EM[ForwardSubdetector::HGCHEF] = new double[100];
    m_hitEperLayer_EM[ForwardSubdetector::HGCHEB] = new double[100];
    m_hitEperLayer_HAD[ForwardSubdetector::HGCEE]  = new double[100];
    m_hitEperLayer_HAD[ForwardSubdetector::HGCHEF] = new double[100];
    m_hitEperLayer_HAD[ForwardSubdetector::HGCHEB] = new double[100];
    */
  }
  
  // read in calibration parameters
  m_calibEE.m_energyCorrMethod = m_calibHEF.m_energyCorrMethod = m_calibHEB.m_energyCorrMethod = m_energyCorrMethod;
  m_calibEE.m_stm = m_calibHEF.m_stm = m_calibHEB.m_stm = stm;
  initPandoraCalibrParameters();
  readCalibrParameterFile();
  if (m_energyCorrMethod == "WEIGHTING")
     readEnergyWeight();

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PandoraCMSPFCandProducer::endJob() 
{
  if(debugHisto){
    file->Write();
    file->Close();
  }
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PandoraCMSPFCandProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PandoraCMSPFCandProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------

void 
PandoraCMSPFCandProducer::beginLuminosityBlock(edm::LuminosityBlock const& iLumiBlock, const edm::EventSetup& iSetup)
{
  //refresh geometry handles
  // Get the ecal/hcal barrel geometry  
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  
  // Get the HGC geometry
  iSetup.get<IdealGeometryRecord>().get("HGCalEESensitive",hgceeGeoHandle) ; 
  iSetup.get<IdealGeometryRecord>().get("HGCalHESiliconSensitive",hgchefGeoHandle) ; 
  iSetup.get<IdealGeometryRecord>().get("HGCalHEScintillatorSensitive",hgchebGeoHandle) ; 
  
  //get tracker geometry
  iSetup.get<TrackerDigiGeometryRecord>().get(tkGeom);

  //get the magnetic field
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  
  // reset the pandora instance
  m_pPandora.reset( new pandora::Pandora() );

  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContent::RegisterAlgorithms(*m_pPandora));

  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, LCContentFast::RegisterAlgorithms(*m_pPandora));

  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, CMSAlgorithms::RegisterAlgorithms(*m_pPandora));

  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, cms_content::RegisterBasicPlugins(*m_pPandora));

  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetBFieldPlugin(*m_pPandora, new CMSBFieldPlugin()));

  PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::RegisterAlgorithmFactory(*m_pPandora, "Template", new CMSTemplateAlgorithm::Factory));

  // reset all the calibration info since it depends on the geometry
  m_calibEE  = CalibHGC(ForwardSubdetector::HGCEE,"EE",debugPrint);
  m_calibHEF = CalibHGC(ForwardSubdetector::HGCHEF,"HEF",debugPrint);
  m_calibHEB = CalibHGC(ForwardSubdetector::HGCHEB,"HEB",debugPrint);
  // read in calibration parameters
  m_calibEE.m_energyCorrMethod = m_calibHEF.m_energyCorrMethod = m_calibHEB.m_energyCorrMethod = m_energyCorrMethod;
  m_calibEE.m_stm = m_calibHEF.m_stm = m_calibHEB.m_stm = stm;
  initPandoraCalibrParameters();
  readCalibrParameterFile();
  if (m_energyCorrMethod == "WEIGHTING")
    readEnergyWeight();
  calibInitialized = false;
  
  // prepare the geometry
  prepareGeometry();

  //rebuild pandora
  if( pandora::STATUS_CODE_SUCCESS != PandoraApi::ReadSettings(*m_pPandora, m_pandoraSettingsXmlFile.fullPath()) ) {
    throw cms::Exception("InvalidXMLConfig")
      << "Unable to parse pandora configuration file";
  }    
  //PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*m_pPandora, m_pandoraSettingsXmlFile.fullPath()));
}


// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PandoraCMSPFCandProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

void 
PandoraCMSPFCandProducer::prepareMCParticleAssociations(edm::Event& iEvent) {
  // ensure we are in a good starting state
  m_recoDetIdToSimHit.clear();
  m_simHitsToSimTracks.clear();
  m_simTrackToStartSimVertex.clear();
  m_simTrackToStartSimVertexEnd.clear();
  m_simTrackToSimVertex.clear();
  m_simVertexToSimTracks.clear();
  m_barCodesToSimTrack.clear();
  m_barCodesToSimVertices.clear();
  m_barCodesToSimHits.clear();
  m_barCodesToGenParticle.clear();
  m_simVertexToSimTrackParent.clear();

  // get all the necessary bits from the event
  edm::Handle<reco::GenParticleCollection> genParticles;
  edm::Handle<std::vector<int> > genBarcodes;
  
  edm::Handle<std::vector<SimTrack> > simTracks;
  edm::Handle<std::vector<SimVertex> > simVertices;
  std::vector<edm::Handle<std::vector<PCaloHit> > > simHits;

  iEvent.getByLabel(inputTagGenParticles_,genParticles);
  iEvent.getByLabel(inputTagGenParticles_,genBarcodes);  
  iEvent.getByLabel(inputTagSimTracks_,simTracks);
  iEvent.getByLabel(inputTagSimVertices_,simVertices);
  for( const auto& tag : inputSimHits_ ) {
    simHits.emplace_back( edm::Handle<std::vector<PCaloHit> >() );
    iEvent.getByLabel(tag,simHits.back());
  }

  // setup the reco det id to sim-hit match
  std::unordered_set<unsigned> reco_detIds;
  for( unsigned i = 0; i < simHits.size(); ++i ) {
    const auto& hits = *(simHits[i]);
    for( unsigned j = 0; j < hits.size(); ++j ) {
      // get topo
      HGCalDetId simId(hits[j].id());
      ForwardSubdetector mysubdet = (ForwardSubdetector)(i+3);
      const HGCalGeometry* geom = nullptr;
      switch(mysubdet) {
      case HGCEE:
        geom = hgceeGeoHandle.product();
        break;
      case HGCHEF:
        geom = hgchefGeoHandle.product();
        break;
      case HGCHEB:
        geom = hgchebGeoHandle.product();
        break;
      default:
        throw cms::Exception("InvalidDetector")
          << "Got invalid HGC subdet: " << mysubdet;
      }
      const HGCalTopology& topo = geom->topology();
      const HGCalDDDConstants& dddConst = topo.dddConstants();
      
      int layer(simId.layer()), cell(simId.cell());
      std::pair<int,int> recoLayerCell = dddConst.simToReco(cell,layer,topo.detectorType());
      cell  = recoLayerCell.first;
      layer = recoLayerCell.second;
      if(layer < 0) continue;
      
      uint32_t recoDetId = ( ( geom == hgceeGeoHandle.product() ) ?
			     (uint32_t)HGCEEDetId(ForwardSubdetector(mysubdet),simId.zside(),layer,simId.sector(),simId.subsector(),cell) :
			     (uint32_t)HGCHEDetId(ForwardSubdetector(mysubdet),simId.zside(),layer,simId.sector(),simId.subsector(),cell)
			     );
      reco_detIds.insert(recoDetId);
      m_recoDetIdToSimHit.emplace(recoDetId,std::make_tuple(i,j,0.0f));
    }
  }
  
  // calculate and store the weights for particles associated to 
  // pcalohits
  for( const unsigned detid : reco_detIds ) {
    auto range = m_recoDetIdToSimHit.equal_range(detid);
    double e_tot = 0.0;
    for( auto iter = range.first; iter != range.second; ++iter ) {
      const PCaloHit& hit = simHits[std::get<0>(iter->second)]->at(std::get<1>(iter->second));
      if( hit.geantTrackId() > 0 ) e_tot += hit.energy();
    }
    for( auto iter = range.first; iter != range.second; ++iter ) {
      const PCaloHit& hit = simHits[std::get<0>(iter->second)]->at(std::get<1>(iter->second));
      if( hit.geantTrackId() > 0 ) {
        float fraction = hit.energy()/e_tot;
        std::make_tuple(std::get<0>(iter->second),std::get<1>(iter->second),fraction).swap(iter->second);
      }
    }
  }

  // stash barcodes for genparticles/simtracks/simvertices/simhits
  // gen particles
  const reco::GenParticleCollection& genParts = *genParticles;
  const std::vector<int> genBCs = *genBarcodes;
  for( unsigned i = 0; i < genParts.size(); ++i ) {
    m_barCodesToGenParticle.emplace(genBCs[i], i);
  }
  // simtracks (indexed by sim vertex or by gen particle barcode)
  const std::vector<SimTrack>& simTks = *simTracks;
  for( unsigned i = 0; i < simTks.size(); ++i ) {
    m_barCodesToSimTrack.emplace( simTks[i].trackId(), i );
    if( !simTks[i].noVertex() ) 
      m_simVertexToSimTracks.emplace( simTks[i].vertIndex(), i);
    if( !simTks[i].noGenpart() ) 
      m_barCodesToSimTrack.emplace( simTks[i].genpartIndex(), i);
  }
  // sim vertices (indexed by sim track indices) 
  const std::vector<SimVertex>& simVtxs = *simVertices;
  for( unsigned i = 0; i < simVtxs.size(); ++i ) {
    m_barCodesToSimVertices.emplace( simVtxs[i].vertexId(), i);
    if( !simVtxs[i].noParent() ) {
      m_simTrackToSimVertex.emplace( simVtxs[i].parentIndex(), i );
      m_simVertexToSimTrackParent.emplace( simVtxs[i].vertexId(), simVtxs[i].parentIndex() );
    }
  }
  // now that we have these maps we can fill all the associations in the pandora event
  // later when we fill tracks and hits we can make the final associations there
  // first pass: import all status one gen particles
  for( unsigned i = 0; i < genParts.size(); ++i ) {
    const reco::GenParticle& genp = genParts[i];
    if( genp.status() == 1 ) {
      PandoraApi::MCParticle::Parameters parameters;
      // basic info
      parameters.m_energy = genp.energy(); 
      parameters.m_momentum = pandora::CartesianVector(genp.px() , genp.py(),  genp.pz() );
      // gen particles don't really have a concept of physical extent... these are here for bookkeeping
      parameters.m_vertex = pandora::CartesianVector( genp.vx() * 10. , genp.vy() * 10., genp.vz() * 10. );
      parameters.m_endpoint = pandora::CartesianVector( genp.vx() * 10. , genp.vy() * 10., genp.vz() * 10. );
      parameters.m_particleId = genp.pdgId();
      parameters.m_mcParticleType = pandora::MCParticleType::MC_3D;
      parameters.m_pParentAddress = (void*)&genp;  
      //std::cout << "Adding genparticle : " << i << ' ' << genp.pdgId() << " to pandora!" << std::endl;
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, parameters));  
    }
  }
  // second pass: import all of the sim tracks
  for( unsigned i = 0; i < simTks.size(); ++i ) {
    // loop over the sim tracks adding each sim track as a MC particle in the event
    const SimTrack& simtrack = simTks[i];
    //std::cout << simtrack.noVertex() << ' ' << simtrack.noGenpart() << std::endl;
    PandoraApi::MCParticle::Parameters parameters;
    // first case,  this track comes from a gen particle
    if( !simtrack.noVertex() && !simtrack.noGenpart() ) {
      // get the parent gen particle of this  sim track
      auto i_start = m_barCodesToGenParticle.find(simtrack.genpartIndex());
      if( i_start != m_barCodesToGenParticle.end() ) {
        const reco::GenParticle& start_vtx = genParts[i_start->second];
        parameters.m_vertex = pandora::CartesianVector( start_vtx.vx() * 10. , start_vtx.vy() * 10., start_vtx.vz() * 10. );
      } else {
        throw cms::Exception("SimTrackHasNoGenParticle")
          << "SimTrack: " << simtrack.trackId() << ' ' << simtrack.type() << " has no associated gen particle!";
      }
    }
    // second case this track comes from a sim vertex
    if( simtrack.noGenpart() ) {
      // get the parent sim vertex of this sim track
      auto i_start = m_simVertexToSimTracks.find(simtrack.vertIndex());
      if( i_start != m_simVertexToSimTracks.end() ) {
        const SimVertex& start_vtx = simVtxs[i_start->first];
        const auto& pos = start_vtx.position();
        parameters.m_vertex = pandora::CartesianVector( pos.x() * 10. , pos.y() * 10., pos.z() * 10. );
      } else {
        throw cms::Exception("SimTrackHasNoSimVertex")
          << "SimTrack: " << simtrack.trackId() << ' ' << simtrack.type() << " has no associated sim vertex!";
      }        
    }      
    
    // find the vertex that has this track as a parent:
    auto i_end = m_simTrackToSimVertex.find(simtrack.trackId());
    if( i_end != m_simTrackToSimVertex.end() ) {
      const SimVertex& end_vtx = simVtxs[i_end->second];
      const auto& pos = end_vtx.position();
      parameters.m_endpoint = pandora::CartesianVector( pos.x() * 10. , pos.y() * 10., pos.z() * 10. );
    } else {
      constexpr float floatmax = 9999.f;
      parameters.m_endpoint = pandora::CartesianVector( floatmax, floatmax, floatmax );
    }
    // basic info
    const auto& sim_p = simtrack.momentum();
    parameters.m_energy = sim_p.T(); 
    parameters.m_momentum = pandora::CartesianVector(sim_p.X() , sim_p.Y(),  sim_p.Z() );
    parameters.m_particleId = simtrack.type();
    parameters.m_mcParticleType = pandora::MCParticleType::MC_3D;
    parameters.m_pParentAddress = (void*)&simtrack;
    
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, parameters));  
  }
  // third pass: set gen particle -> sim track associations
  for( unsigned i = 0; i < genParts.size(); ++i ) {
    const reco::GenParticle& genp = genParts[i];
    const void* gen_ptr = static_cast<const void*>(&genp);
    if( genp.status() == 1 ) {
      auto i_tk = m_barCodesToSimTrack.find(genBCs[i]);
      if( i_tk != m_barCodesToSimTrack.end() ) {
        const SimTrack& simtrack = simTks[i_tk->second];
        const void* sim_ptr = static_cast<const void*>(&simtrack);
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*m_pPandora,gen_ptr,sim_ptr));
      }
    }
  }
  // fourth pass: set sim track -> sim track associations (skipping the gen -> sim associations)
  for( unsigned i = 0; i < simTks.size(); ++i ) {
    // only find parents!
    if( !simTks[i].noVertex() && !simTks[i].noGenpart() ) continue; // skip sim tracks whose parents are gen particles (we already covered those)
    if( !simTks[i].noGenpart() ) {
      const SimTrack& daughter = simTks[i];
      const void* daughter_ptr = static_cast<const void*>(&daughter);    
    // get the sim vertex that this daughter belongs to    
      const void* parent_ptr = nullptr;
      auto i_vtxparent = m_simVertexToSimTrackParent.find(daughter.vertIndex());
      if( i_vtxparent != m_simVertexToSimTrackParent.end() ) {
        auto i_parent = m_barCodesToSimTrack.find(i_vtxparent->second);
        if( i_parent != m_barCodesToSimTrack.end() ) {
          parent_ptr = &simTks[i_parent->second];
        } else {
          throw cms::Exception("SimTrackHasNoBarcode")
            << "SimTrack Barcode: " << i_vtxparent->second << " has no associated sim track!";
        }
      } else {
        throw cms::Exception("SimVertexHasNoParent")
          << "SimTrack: " << i << " has no associated parent simtrack!";
      } 
      PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*m_pPandora,parent_ptr,daughter_ptr));
    }
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PandoraCMSPFCandProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PandoraCMSPFCandProducer);

