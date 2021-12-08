/**
 *  @file   larpandoracontent/LArUtility/TierPruningAlgorithm.cc
 *
 *  @brief  Implementation of the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArUtility/TierPruningAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode TierPruningAlgorithm::Run()
{
    std::map<HitType, std::string> viewToNameMap;
    for (const std::string listName : m_cluster2dListNames)
    {
        const ClusterList *pClusterList{nullptr};
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pClusterList));

        if (pClusterList && !pClusterList->empty())
            viewToNameMap[LArClusterHelper::GetClusterHitType(pClusterList->front())] = listName;
    }

    for (unsigned int i = 0; i < m_pfoListNames.size(); ++i)
    {
        // ATTN - one-to-one correspondance required between PFO and 3D cluster lists
        const std::string &pfoListName{m_pfoListNames.at(i)};
        const std::string &clusterListName{m_cluster3dListNames.at(i)};
        const PfoList *pPfoList(nullptr);
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, pfoListName, pPfoList));

        if (pPfoList && !pPfoList->empty())
        {
            PfoList nonPrimaryPfoList;
            for (const ParticleFlowObject *pPfo : *pPfoList)
            {
                if (LArPfoHelper::GetHierarchyTier(pPfo) != 1)
                    nonPrimaryPfoList.emplace_back(pPfo);
            }

            // Delete tier 2+ PFOs
            for (const ParticleFlowObject *pPfo : nonPrimaryPfoList)
            {
                // Remove 3D clusters
                ClusterList clustersToRemove;
                LArPfoHelper::GetClusters(pPfo, TPC_3D, clustersToRemove);

                for (const Cluster *pCluster : clustersToRemove)
                {
                    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfo, pCluster));
                    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                        PandoraContentApi::Delete<Cluster>(*this, pCluster, clusterListName));
                }
                clustersToRemove.clear();

                // Remove and delete 2D clusters and hits
                for (const HitType view : { TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W })
                {
                    LArPfoHelper::GetClusters(pPfo, view, clustersToRemove);
                    const std::string listName{viewToNameMap[view]};

                    for (const Cluster *pCluster : clustersToRemove)
                    {
                        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfo, pCluster));
                        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                            PandoraContentApi::Delete<Cluster>(*this, pCluster, listName));
                    }
                    clustersToRemove.clear();
                }

                PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                    PandoraContentApi::Delete<ParticleFlowObject>(*this, pPfo, pfoListName));
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TierPruningAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "Cluster3DListNames", m_cluster3dListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "Cluster2DListNames", m_cluster2dListNames));

    if (m_pfoListNames.size() != m_cluster3dListNames.size())
    {
        std::cout << "TierPruningAlgorithm: Mismatch between PFO and 3D Cluster list sizes" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
