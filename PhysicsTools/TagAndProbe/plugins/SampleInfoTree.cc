#include "PhysicsTools/TagAndProbe/plugins/SampleInfoTree.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

tnp::SampleInfoTree::SampleInfoTree(const edm::ParameterSet& iConfig) : 
  weightSrcToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genInfo")))  {
  
  // make trees as requested
  edm::Service<TFileService> fs;
  addTree_ = fs->make<TTree>("sampleInfo", "sampleInfo");
  //hNvtx = fs->make<TH1F>("nvtx", "nvtx", 100, 0., 100.);
    
  addTree_->Branch("sumWeight", &totSumWeight_, "sumWeight/D");
  addTree_->Branch("nEvents", &totNEvents_, "nEvents/D");

  totSumWeight_ = 0.0;
  totNEvents_ = 0.0;

  produces<edm::MergeableDouble, edm::InLumi> ("totalGenWeight"); 
  produces<edm::MergeableDouble, edm::InLumi> ("totalEvent");
}
    
void tnp::SampleInfoTree::beginLuminosityBlock(const edm::LuminosityBlock & theLuminosityBlock, const edm::EventSetup & theSetup) {
  sumWeight_ = 0.0;
  nEvents_ = 0.0;
}

void tnp::SampleInfoTree::endLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup) 
{}
//
//  edm::Handle<edm::MergeableDouble> totWeight;
//  iLumi.getByLabel(edm::InputTag("weightsCount","totalWeight"), totWeight);
//  if (totWeight.isValid())
//    totalGenWeight_ += (double)totWeight->value;
//  
//}

void tnp::SampleInfoTree::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<GenEventInfoProduct> genInfo;
  //edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
  
  if(!iEvent.isRealData()) {
    iEvent.getByToken(weightSrcToken_, genInfo);
    
    const auto & weights = genInfo->weights(); 
    if(!weights.empty()) {
      sumWeight_ += weights[0];
    }
    //if( doTruePileup_ || doObsPileup_ ) {
    //	iEvent.getByToken(puInfoToken_,puInfo);
    //	hysto_type::value_type truePu=0., obsPu=0.;
    //	for( auto & frame : *puInfo ) {
    //	  /// std::cout << frame.getBunchCrossing() << std::endl;
    //	  if( frame.getBunchCrossing() == 0 ) {
    //	    truePu = frame.getTrueNumInteractions();
    //	    obsPu = frame.getPU_NumInteractions();
    //	    break;
    //	  }
    //	}
    //	/// cout << truePu << " " << obsPu << endl; 
    //
    //	if( doTruePileup_ ) {
    //	  size_t bin = 0;
    //	  if( truePu >= maxTruePileup_ ) { bin = truePileup_.size() - 1; }
    //	  else if( truePu >= minTruePileup_ ) { bin = (size_t)std::floor( (truePu-minTruePileup_) / widthTruePileup_) + 1; }
    //	  truePileup_[bin] += weights[0];
    //	}
    //	
    //	if( doObsPileup_ ) {
    //	  size_t bin = 0;
    //	  if( obsPu >= maxObsPileup_ ) { bin = obsPileup_.size() - 1; }
    //	  else if( obsPu >= minObsPileup_ ) { bin = (size_t)std::floor( (obsPu-minObsPileup_) / widthObsPileup_) + 1; }
    //	  /// cout << bin << " " << std::floor( (obsPu-minObsPileup_) / widthObsPileup_) << endl;
    //	  obsPileup_[bin] += weights[0];
    //	}
    //}     
  } else {
    sumWeight_ += 1.;
    totSumWeight_ += 1.;
  }

  totNEvents_ += 1.0;
  nEvents_ += 1.;
}  

void tnp::SampleInfoTree::endJob() {
  addTree_->Fill();
}

void tnp::SampleInfoTree::endLuminosityBlockProduce(edm::LuminosityBlock & theLuminosityBlock, const edm::EventSetup & theSetup) {
  //LogTrace("WeightsCounting") << "endLumi: adding " << weightProcessedInLumi_ << " events" << endl;
  
  std::auto_ptr<edm::MergeableDouble> numWeightssPtr(new edm::MergeableDouble);
  numWeightssPtr->value = sumWeight_;
  theLuminosityBlock.put(numWeightssPtr, "totalGenWeight");
  
  std::auto_ptr<edm::MergeableDouble> numEventsPtr(new edm::MergeableDouble);
  numEventsPtr->value = nEvents_;
  theLuminosityBlock.put(numEventsPtr, "totalEvent");
  //return;
  //addTree_->Fill();
  
}

DEFINE_FWK_MODULE(tnp::SampleInfoTree);
