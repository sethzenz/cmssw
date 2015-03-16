/**
 *  @file   RecoParticleFlow/PandoraTranslator/interace/CMSAlgorithms.h
 * 
 *  @brief  Header file detailing faster versions of algorithms in the LCContent library, using e.g. KD-trees and *relying on c++11*
 * 
 *  $Log: $
 */
#ifndef CMS_ALGORITHMS_H
#define CMS_ALGORITHMS_H 1

#include "RecoParticleFlow/PandoraTranslator/interface/ExitingTrackAlg.h"
#include "RecoParticleFlow/PandoraTranslator/interface/LoopingTrackAssociationAlgorithm.h"
#include "RecoParticleFlow/PandoraTranslator/interface/MainFragmentRemovalAlgorithmFast.h"
#include "RecoParticleFlow/PandoraTranslator/interface/MipPhotonSeparationAlgorithm.h"
#include "RecoParticleFlow/PandoraTranslator/interface/MuonClusterAssociationAlgorithm.h"
#include "RecoParticleFlow/PandoraTranslator/interface/ResolveTrackAssociationsAlg.h"
#include "RecoParticleFlow/PandoraTranslator/interface/SplitMergedClustersAlg.h"
#include "RecoParticleFlow/PandoraTranslator/interface/SplitTrackAssociationsAlg.h"
#include "RecoParticleFlow/PandoraTranslator/interface/TrackDrivenAssociationAlg.h"
#include "RecoParticleFlow/PandoraTranslator/interface/TrackDrivenMergingAlg.h"
#include "RecoParticleFlow/PandoraTranslator/interface/TrackRecoveryHelixAlgorithm.h"
#include "RecoParticleFlow/PandoraTranslator/interface/TrackRecoveryInteractionsAlgorithm.h"

#include "RecoParticleFlow/PandoraTranslator/interface/MuonCoilCorrection.h"

/**
 *  @brief  CMSAlgorithms class
 */
class CMSAlgorithms
{
public:
    #define CMS_ALGORITHM_LIST(d)                                                                            \
        d("CMSMainFragmentRemovalFast",     cms_content_fast::MainFragmentRemovalAlgorithm::Factory)         \
        d("CMSExitingTrack",               cms_content::ExitingTrackAlg::Factory)                            \
        d("CMSLoopingTrackAssociation",    cms_content::LoopingTrackAssociationAlgorithm::Factory)           \
        d("CMSMipPhotonSeparation",        cms_content::MipPhotonSeparationAlgorithm::Factory)               \
        d("CMSMuonClusterAssociation",     cms_content::MuonClusterAssociationAlgorithm::Factory)            \
        d("CMSResolveTrackAssociations",   cms_content::ResolveTrackAssociationsAlg::Factory)                \
        d("CMSSplitMergedClusters",        cms_content::SplitMergedClustersAlg::Factory)                     \
        d("CMSSplitTrackAssociations",     cms_content::SplitTrackAssociationsAlg::Factory)                  \
        d("CMSTrackDrivenAssociation",     cms_content::TrackDrivenAssociationAlg::Factory)                  \
        d("CMSTrackDrivenMerging",         cms_content::TrackDrivenMergingAlg::Factory)                      \
        d("CMSTrackRecoveryHelix",         cms_content::TrackRecoveryHelixAlgorithm::Factory)                \
        d("CMSTrackRecoveryInteractions",  cms_content::TrackRecoveryInteractionsAlgorithm::Factory)

    /**
     *  @brief  Register all the linear collider algorithms with pandora
     * 
     *  @param  pandora the pandora instance with which to register content
     */
    static pandora::StatusCode RegisterAlgorithms(const pandora::Pandora &pandora);
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::StatusCode CMSAlgorithms::RegisterAlgorithms(const pandora::Pandora &pandora)
{
    CMS_ALGORITHM_LIST(PANDORA_REGISTER_ALGORITHM);

    return pandora::STATUS_CODE_SUCCESS;
}

#endif // #ifndef CMS_ALGORITHMS_H
