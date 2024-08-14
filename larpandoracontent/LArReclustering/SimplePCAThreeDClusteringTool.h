/**
 *  @file   larpandoracontent/LArReclustering/SimplePCAThreeDClusteringTool.h
 *
 *  @brief  Header file for the reclustering algorithm that uses transverse calorimetric profiles.
 *
 *  $Log: $
 */

#ifndef LAR_SIMPLE_PCA_THREE_D_CLUSTERING_TOOL_H
#define LAR_SIMPLE_PCA_THREE_D_CLUSTERING_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h"

namespace lar_content
{
/**
 *  @brief  RecursivePfoMopUpAlgorithm class
 */
class SimplePCAThreeDClusteringTool : public ClusteringTool
{
public:

    SimplePCAThreeDClusteringTool();

private:

    bool Run(const pandora::Algorithm *const /*pAlgorithm*/, std::vector<pandora::CaloHitList*> &newCaloHitListsVector) override;
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle /*xmlHandle*/);
};

} // namespace lar_content

#endif // #endif LAR_SIMPLE_PCA_THREE_D_CLUSTERING_TOOL_H
