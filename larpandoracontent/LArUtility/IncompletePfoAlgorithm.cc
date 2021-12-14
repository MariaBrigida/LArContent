/**
 *  @file   larpandoracontent/LArUtility/IncompletePfoAlgorithm.cc
 *
 *  @brief  Implementation of the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArUtility/IncompletePfoAlgorithm.h"

using namespace pandora;

namespace lar_content
{

IncompletePfoAlgorithm::IncompletePfoAlgorithm()
{
}

StatusCode IncompletePfoAlgorithm::Run()
{
    std::map<std::string, CaloHitList> nameToCaloHitListMap;
    // ATTN - You can't actually delete CaloHits, so you we need to use some list trickery to remove specific 3D hits from our lists
    // by renaming the list and then saving (moving) retained 3D hits back to the original list. We rename those lists here and create
    // mappings from the original names to the associated lists.
    for (std::string listName : m_caloHitListNames)
    {
        std::string tempListName{"_Temp_" + listName};
        const CaloHitList *pCaloHitList{nullptr};
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
            PandoraContentApi::GetList(*this, listName, pCaloHitList));
        nameToCaloHitListMap[listName] = CaloHitList(*pCaloHitList);
        PandoraContentApi::RenameList<CaloHitList>(*this, listName, tempListName);
    }

    for (std::string listName : m_pfoListNames)
    {
        const PfoList *pPfoList{nullptr};
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
            PandoraContentApi::GetList(*this, listName, pPfoList));
        PfoList allPfos;
        LArPfoHelper::GetAllDownstreamPfos(*pPfoList, allPfos);

        for (const ParticleFlowObject *pPfo : allPfos)
        {
            // Neutrino PFO doesn't contain hits, so skip it
            if (LArPfoHelper::IsNeutrino(pPfo))
                continue;

            // Check for the presence of 3D hits
            ClusterList clusters3d;
            LArPfoHelper::GetClusters(pPfo, HitType::TPC_3D, clusters3d);
            const bool no3dHits{clusters3d.empty()};
            int missingViews{0};
            if (!no3dHits)
            {
                // If we have 3D hits, check how many views produced those hits
                for (HitType view : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
                {
                    ClusterList clusters;
                    LArPfoHelper::GetClusters(pPfo, view, clusters);
                    if (clusters.empty())
                        ++missingViews;
                }
                if (missingViews)
                {
                    // If it was less than three views to make the 3D hits, dissolve the PFO and remove the 3D hits from the appropriate list
                    for (const Cluster *pCluster : clusters3d)
                    {
                        CaloHitList hitsToRemove(pCluster->GetIsolatedCaloHitList());
                        pCluster->GetOrderedCaloHitList().FillCaloHitList(hitsToRemove);
                        for (std::string caloHitListName : m_caloHitListNames)
                        {
                            CaloHitList &caloHitList{nameToCaloHitListMap[caloHitListName]};
                            for (const CaloHit *pCaloHit : hitsToRemove)
                            {
                                auto iter{std::find(caloHitList.begin(), caloHitList.end(), pCaloHit)};
                                if (iter != caloHitList.end())
                                    caloHitList.erase(iter);
                            }
                        }

                        // Remove the cluster from the PFO and delete the 3D cluster
                        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfo, pCluster));
                        for (std::string clusterListName : m_clusterListNames)
                        {
                            PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                                PandoraContentApi::Delete<Cluster>(*this, pCluster, clusterListName));
                        }
                    }
                }
            }
            if (no3dHits || missingViews)
            {
                // This PFO has some sub-optimal aspects, remove its links to parent and daughter PFOs, remove, but don't delete its
                // clusters (if it had a 3D cluster, it's alredy been deleted) and we'd like to try to better match the 2D clusters
                // in future algorithm runs)
                const PfoList childList(pPfo->GetDaughterPfoList());
                for (const ParticleFlowObject *pChild : childList)
                    PandoraContentApi::RemovePfoParentDaughterRelationship(*this, pPfo, pChild);
                const PfoList parentList(pPfo->GetParentPfoList());
                for (const ParticleFlowObject *pParent : parentList)
                    PandoraContentApi::RemovePfoParentDaughterRelationship(*this, pParent, pPfo);

                for (const HitType view : { TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W })
                {
                    ClusterList clustersToRemove;
                    LArPfoHelper::GetClusters(pPfo, view, clustersToRemove);

                    for (const Cluster *pCluster : clustersToRemove)
                    {
                        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfo, pCluster));
                    }
                }

                // Finally, delete the PFO
                PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
                    PandoraContentApi::Delete<ParticleFlowObject>(*this, pPfo, listName));
            }
        }
    }

    // ATTN - Whether we made edits or not, the earlier renaming means we need to save what's left back to the original list names here
    for (std::string caloHitListName : m_caloHitListNames)
    {
        CaloHitList &caloHitList{nameToCaloHitListMap[caloHitListName]};
        PandoraContentApi::SaveList<CaloHitList>(*this, caloHitList, caloHitListName);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode IncompletePfoAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
