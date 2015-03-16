/**
 *  @file   PandoraSDK/src/LCHelpers/ReclusterHelper.cc
 * 
 *  @brief  Implementation of the recluster helper class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "RecoParticleFlow/PandoraTranslator/interface/ReclusterHelper.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"

using namespace pandora;

namespace cms_content
{

float ReclusterHelper::GetTrackClusterCompatibility(const Pandora &pandora, const Cluster *const pCluster, const TrackList &trackList)
{
    float trackEnergySum(0.);
    float trackEnergyErrorSum2(0.);

    for (TrackList::const_iterator trackIter = trackList.begin(), trackIterEnd = trackList.end(); trackIter != trackIterEnd; ++trackIter)
    {
        trackEnergySum += (*trackIter)->GetEnergyAtDca();
        const reco::PFRecTrack* pftrack = static_cast<const reco::PFRecTrack*>((*trackIter)->GetParentTrackAddress());
        const reco::TrackRef& trackref = pftrack->trackRef();
        const float deltaPtRel = trackref->ptError()/trackref->pt();
        const float deltaE = deltaPtRel*trackref->p();
        trackEnergyErrorSum2 += std::pow(deltaE,2.0f);
    }

    const float hadronicEnergyResolution(pandora.GetSettings()->GetHadronicEnergyResolution());

    if ((trackEnergySum < std::numeric_limits<float>::epsilon()) || (hadronicEnergyResolution < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float sigmaE_exp_calo(hadronicEnergyResolution * std::sqrt(trackEnergySum));
    const float sigmaE = std::sqrt( sigmaE_exp_calo*sigmaE_exp_calo + trackEnergyErrorSum2 );
    const float chi((pCluster->GetTrackComparisonEnergy(pandora) - trackEnergySum) / sigmaE);

    return chi;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float ReclusterHelper::GetTrackClusterCompatibility(const Pandora &pandora, const float clusterEnergy, const float trackEnergy,
                                                    const float trackEnergyError2)
{
    const float hadronicEnergyResolution(pandora.GetSettings()->GetHadronicEnergyResolution());

    if ((trackEnergy < std::numeric_limits<float>::epsilon()) || (hadronicEnergyResolution < std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const float sigmaE_exp_calo(hadronicEnergyResolution * std::sqrt(trackEnergy));
    const float sigmaE = std::sqrt( sigmaE_exp_calo*sigmaE_exp_calo + trackEnergyError2 );
    const float chi((clusterEnergy - trackEnergy) / sigmaE);

    return chi;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ReclusterHelper::ExtractReclusterResults(const Pandora &pandora, const ClusterList *const pReclusterCandidatesList, ReclusterResult &reclusterResult)
{
    unsigned int nExcessTrackAssociations(0);
    float chi(0.), chi2(0.), dof(0.), unassociatedEnergy(0.), minTrackAssociationEnergy(std::numeric_limits<float>::max());

    for (ClusterList::const_iterator iter = pReclusterCandidatesList->begin(), iterEnd = pReclusterCandidatesList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        const TrackList &trackList(pCluster->GetAssociatedTrackList());
        const unsigned int nTrackAssociations(trackList.size());
        const float clusterEnergy(pCluster->GetHadronicEnergy());

        if (0 == nTrackAssociations)
        {
            unassociatedEnergy += clusterEnergy;
            continue;
        }

        if (clusterEnergy < minTrackAssociationEnergy)
            minTrackAssociationEnergy = clusterEnergy;

        nExcessTrackAssociations += nTrackAssociations - 1;

        const float newChi(ReclusterHelper::GetTrackClusterCompatibility(pandora, pCluster, trackList));

        chi2 += newChi * newChi;
        chi += newChi;
        dof += 1.f;
    }

    if (dof < std::numeric_limits<float>::epsilon())
        return STATUS_CODE_FAILURE;

    reclusterResult.SetChi(chi);
    reclusterResult.SetChi2(chi2);
    reclusterResult.SetChiPerDof(chi /= dof);
    reclusterResult.SetChi2PerDof(chi2 /= dof);
    reclusterResult.SetUnassociatedEnergy(unassociatedEnergy);
    reclusterResult.SetMinTrackAssociationEnergy(minTrackAssociationEnergy);
    reclusterResult.SetNExcessTrackAssociations(nExcessTrackAssociations);

    return STATUS_CODE_SUCCESS;
}

} // namespace cms_content
