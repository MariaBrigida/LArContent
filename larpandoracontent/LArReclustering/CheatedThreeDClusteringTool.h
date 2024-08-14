/**
 *  @file   larpandoracontent/LArReclustering/CheatedThreeDClusteringTool.h
 *
 *  @brief  Header file for the reclustering algorithm that uses transverse calorimetric profiles.
 *
 *  $Log: $
 */

#ifndef LAR_CHEATED_THREE_D_CLUSTERING_TOOL_H
#define LAR_CHEATED_THREE_D_CLUSTERING_TOOL_H 1

#include "Pandora/Algorithm.h"
#include "larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h"

namespace lar_content
{

class ClusteringTool;

class CheatedThreeDClusteringTool : public ClusteringTool
{
public:

    CheatedThreeDClusteringTool();

private:

    //pandora::StatusCode Run();

    //bool Run(ThreeDReclusteringAlgorithm *const pAlgorithm, std::vector<pandora::CaloHitList> &newCaloHitListsVector);
    //virtual bool Run(std::vector<pandora::CaloHitList> &newCaloHitListsVector) = 0;
    bool Run(const pandora::Algorithm *const pAlgorithm, std::vector<pandora::CaloHitList*> &newCaloHitListsVector) override;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    int GetMainMcParticleIndex(const pandora::Algorithm *const pAlgorithm, const pandora::CaloHit *const pCaloHit);

    std::string m_mcParticleListName; ///< The mc particle list name 
};

} // namespace lar_content

#endif // #endif LAR_CHEATED_THREE_D_CLUSTERING_TOOL_H
