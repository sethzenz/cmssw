/**
 *  @file   RecoParticleFlow/PandoraTranslator/src/MuonCoilCorrection.cc
 * 
 *  @brief  Implementation of the cms muon correction plugins class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "RecoParticleFlow/PandoraTranslator/interface/ReclusterHelper.h"

#include "RecoParticleFlow/PandoraTranslator/interface/MuonCoilCorrection.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"

using namespace pandora;

namespace cms_content {

MuonCoilCorrection::MuonCoilCorrection() :
    m_muonHitEnergy(0.5f),
    m_coilEnergyLossCorrection(10.f),
    m_minMuonHitsInInnerLayer(3),
    m_coilEnergyCorrectionChi(3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonCoilCorrection::MakeEnergyCorrections(const Cluster *const pCluster, float &correctedHadronicEnergy) const
{
    bool containsMuonHit(false);
    unsigned int nMuonHitsInInnerLayer(0);
    unsigned int muonInnerLayer(std::numeric_limits<unsigned int>::max());

    // Extract muon-based properties from the cluster
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            if ((*hitIter)->GetHitType() == MUON)
            {
                containsMuonHit = true;
                ++nMuonHitsInInnerLayer;
            }
        }

        if (containsMuonHit)
        {
            muonInnerLayer = iter->first;
            break;
        }
    }

    if (!containsMuonHit)
        return STATUS_CODE_SUCCESS;;

    // Check whether energy deposits are likely to have been lost in coil region
    const CartesianVector muonInnerLayerCentroid(pCluster->GetCentroid(muonInnerLayer));
    const float centroidX(muonInnerLayerCentroid.GetX()), centroidY(muonInnerLayerCentroid.GetY());

    const float muonInnerLayerRadius(std::sqrt(centroidX * centroidX + centroidY * centroidY));
    const float coilInnerRadius(this->GetPandora().GetGeometry()->GetSubDetector(COIL).GetInnerRCoordinate());

    if (muonInnerLayerRadius < coilInnerRadius)
        return STATUS_CODE_SUCCESS;;

    const TrackList &trackList(pCluster->GetAssociatedTrackList());

    if (pCluster->GetInnerPseudoLayer() == muonInnerLayer)
    {
        // Energy correction for standalone muon cluster
        correctedHadronicEnergy += m_muonHitEnergy * static_cast<float>(nMuonHitsInInnerLayer);
    }
    else if (trackList.empty())
    {
        // Energy correction for neutral hadron cluster spilling into coil and muon detectors
        correctedHadronicEnergy += m_coilEnergyLossCorrection;
    }
    else
    {
        // Energy correction for charged hadron cluster spilling into coil and muon detectors
        if (nMuonHitsInInnerLayer < m_minMuonHitsInInnerLayer)
            return STATUS_CODE_SUCCESS;;

        float trackEnergySum(0.f);
        float trackEnergyErrorSum2(0.f);

        for (TrackList::const_iterator iter = trackList.begin(), iterEnd = trackList.end(); iter != iterEnd; ++iter)
        {
            trackEnergySum += (*iter)->GetEnergyAtDca();
            const reco::PFRecTrack* pftrack = static_cast<const reco::PFRecTrack*>((*iter)->GetParentTrackAddress());
            const reco::TrackRef& trackref = pftrack->trackRef();
            const float deltaPtRel = trackref->ptError()/trackref->pt();
            const float deltaE = deltaPtRel*trackref->p();
            trackEnergyErrorSum2 += std::pow(deltaE,2.0f);
        }

        const float oldChi(ReclusterHelper::GetTrackClusterCompatibility(this->GetPandora(), correctedHadronicEnergy, trackEnergySum, trackEnergyErrorSum2));
        const float newChi(ReclusterHelper::GetTrackClusterCompatibility(this->GetPandora(), correctedHadronicEnergy + m_coilEnergyLossCorrection, trackEnergySum, trackEnergyErrorSum2));

        if ((oldChi < m_coilEnergyCorrectionChi) && (std::fabs(newChi) < std::fabs(oldChi)))
        {
            correctedHadronicEnergy += m_coilEnergyLossCorrection;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonCoilCorrection::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MuonHitEnergy", m_muonHitEnergy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CoilEnergyLossCorrection", m_coilEnergyLossCorrection));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMuonHitsInInnerLayer", m_minMuonHitsInInnerLayer));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "CoilEnergyCorrectionChi", m_coilEnergyCorrectionChi));

    return STATUS_CODE_SUCCESS;
}

} // cms_content
