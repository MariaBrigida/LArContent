/**
 *  @file   larpandoracontent/LArUtility/TierPruningAlgorithm.h
 *
 *  @brief  Header file for the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_TIER_PRUNING_ALGORITHM_H
#define LAR_TIER_PRUNING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  TierPruningAlgorithm class
 */
class TierPruningAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;       ///< The pfo list names
    pandora::StringVector m_cluster3dListNames; ///< The 3D cluster list names
    pandora::StringVector m_cluster2dListNames; ///< The 2D cluster list names
};

} // namespace lar_content

#endif // #ifndef LAR_TIER_PRUNING_ALGORITHM_H
