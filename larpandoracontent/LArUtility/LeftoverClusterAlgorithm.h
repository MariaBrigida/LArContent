/**
 *  @file   larpandoracontent/LArUtility/LeftoverClusterAlgorithm.h
 *
 *  @brief  Header file for the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_LEFTOVER_CLUSTER_ALGORITHM_H
#define LAR_LEFTOVER_CLUSTER_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  LeftoverClusterAlgorithm class
 */
class LeftoverClusterAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_clusterListNames;       ///< The input cluster list names
    pandora::StringVector m_outputClusterListNames; ///< The output cluster list names
};

} // namespace lar_content

#endif // #ifndef LAR_LEFTOVER_CLUSTER_ALGORITHM_H
