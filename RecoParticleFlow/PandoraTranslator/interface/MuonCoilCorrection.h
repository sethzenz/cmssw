/**
 *  @file   RecoParticleFlow/PandoraTranslator/interface/MuonCoilCorrection.h
 * 
 *  @brief  Header file for the lc energy correction plugins class.
 * 
 *  $Log: $
 */
#ifndef CMS_MUON_COIL_CORRECTION_H
#define CMS_MUON_COIL_CORRECTION_H 1

#include "Plugins/EnergyCorrectionsPlugin.h"

namespace cms_content {

/**
 *   @brief  MuonCoilCorrection class. Addresses issue of energy loss in uninstrumented coil region.
 */
class MuonCoilCorrection : public pandora::EnergyCorrectionPlugin
{
 public:
        /**
         *  @brief  Default constructor
         */
  MuonCoilCorrection();
  
  pandora::StatusCode MakeEnergyCorrections(const pandora::Cluster *const pCluster, float &correctedEnergy) const;
  
 private:
  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
  
  float           m_muonHitEnergy;                    ///< The energy for a digital muon calorimeter hit, units GeV
  float           m_coilEnergyLossCorrection;         ///< Energy correction due to missing energy deposited in coil, units GeV
  unsigned int    m_minMuonHitsInInnerLayer;          ///< Min muon hits in muon inner layer to correct charged cluster energy
  float           m_coilEnergyCorrectionChi;          ///< Track-cluster chi value used to assess need for coil energy correction
};

} // cms_content

#endif // #ifndef CMS_MUON_COIL_CORRECTION
