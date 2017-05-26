#ifndef CalibratedElectronProducerRun2_h
#define CalibratedElectronProducerRun2_h

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "CondFormats/DataRecord/interface/GBRDWrapperRcd.h"
#include "CondFormats/EgammaObjects/interface/GBRForestD.h"
#include "EgammaAnalysis/ElectronTools/interface/EpCombinationToolSemi.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEnergyCalibratorRun2.h"

#include <vector>
#include <random>
#include <TRandom2.h>

template<typename T>
class CalibratedElectronProducerRun2T: public edm::stream::EDProducer<>
{
    public:
        explicit CalibratedElectronProducerRun2T( const edm::ParameterSet & ) ;
        ~CalibratedElectronProducerRun2T() override;
        void produce( edm::Event &, const edm::EventSetup & ) override ;

    private:
        edm::EDGetTokenT<edm::View<T> >         theElectronToken;
        std::vector<std::string>                theGBRForestName;
        std::vector<const GBRForestD* > theGBRForestHandle;

        EpCombinationToolSemi        theEpCombinationTool;
        ElectronEnergyCalibratorRun2 theEnCorrectorRun2;
        std::unique_ptr<TRandom> theSemiDeterministicRng;
        edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionEBToken_;
        edm::EDGetTokenT<EcalRecHitCollection> recHitCollectionEEToken_;
};

template<typename T>
CalibratedElectronProducerRun2T<T>::CalibratedElectronProducerRun2T( const edm::ParameterSet & conf ) :
  theElectronToken(consumes<edm::View<T> >(conf.getParameter<edm::InputTag>("electrons"))),
  theGBRForestName(conf.getParameter< std::vector<std::string> >("gbrForestName")),
  theEpCombinationTool(),
  theEnCorrectorRun2(theEpCombinationTool, conf.getParameter<bool>("isMC"), conf.getParameter<bool>("isSynchronization"), conf.getParameter<std::string>("correctionFile")),
  recHitCollectionEBToken_(consumes<EcalRecHitCollection>(conf.getParameter<edm::InputTag>( "recHitCollectionEB" ))),
  recHitCollectionEEToken_(consumes<EcalRecHitCollection>(conf.getParameter<edm::InputTag>( "recHitCollectionEE" )))
{
  if (conf.existsAs<bool>("semiDeterministic") && conf.getParameter<bool>("semiDeterministic")) {
    theSemiDeterministicRng.reset(new TRandom2());
    theEnCorrectorRun2.initPrivateRng(theSemiDeterministicRng.get());
  }
  produces<std::vector<T> >();
}

template<typename T>
CalibratedElectronProducerRun2T<T>::~CalibratedElectronProducerRun2T()
{
}

template<typename T>
void
CalibratedElectronProducerRun2T<T>::produce( edm::Event & iEvent, const edm::EventSetup & iSetup ) 
{

    for (auto&& forestName : theGBRForestName) {
      edm::ESHandle<GBRForestD> forestHandle;
      iSetup.get<GBRDWrapperRcd>().get(forestName, forestHandle);
      theGBRForestHandle.emplace_back(forestHandle.product());      
    }

    theEpCombinationTool.init(theGBRForestHandle);

    edm::Handle<edm::View<T> > in;
    iEvent.getByToken(theElectronToken, in);

    std::unique_ptr<std::vector<T> > out(new std::vector<T>());

	edm::Handle<EcalRecHitCollection> recHitCollectionEBHandle;
	edm::Handle<EcalRecHitCollection> recHitCollectionEEHandle;

	iEvent.getByToken(recHitCollectionEBToken_, recHitCollectionEBHandle);
	iEvent.getByToken(recHitCollectionEEToken_, recHitCollectionEEHandle);

    out->reserve(in->size());   

    if (theSemiDeterministicRng && !in->empty()) { // no need to set a seed if in is empty
        const auto & first = in->front();
        std::seed_seq seeder = {int(iEvent.id().event()), int(iEvent.id().luminosityBlock()), int(iEvent.id().run()),
          int(in->size()), int(std::numeric_limits<int>::max()*first.phi()/M_PI) & 0xFFF, int(first.pdgId())};
        uint32_t seed = 0, tries = 10;
        do {
            seeder.generate(&seed,&seed+1); tries++;
        } while (seed == 0 && tries < 10);
        theSemiDeterministicRng->SetSeed(seed ? seed : iEvent.id().event());
    }

    for (const T &ele : *in) {
        out->push_back(ele);
		const EcalRecHitCollection* recHits = (ele.isEB()) ? recHitCollectionEBHandle.product() : recHitCollectionEEHandle.product();
        theEnCorrectorRun2.calibrate(out->back(), iEvent.id().run(), recHits, iEvent.streamID());
    }
    
    iEvent.put(std::move(out));
}

typedef CalibratedElectronProducerRun2T<reco::GsfElectron> CalibratedElectronProducerRun2;
typedef CalibratedElectronProducerRun2T<pat::Electron> CalibratedPatElectronProducerRun2;

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(CalibratedElectronProducerRun2);
DEFINE_FWK_MODULE(CalibratedPatElectronProducerRun2);

#endif
