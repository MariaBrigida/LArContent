/**
 *  @file   larpandoracontent/LArUtility/LeftoverClusterAlgorithm.cc
 *
 *  @brief  Implementation of the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArUtility/LeftoverClusterAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode LeftoverClusterAlgorithm::Run()
{
    for (size_t i = 0; i < m_clusterListNames.size(); ++i)
    {
        const std::string inputClusterListName{m_clusterListNames.at(i)};
        const ClusterList *pClusterList{nullptr};
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=,
            PandoraContentApi::GetList(*this, inputClusterListName, pClusterList));
        const std::string outputClusterListName{m_outputClusterListNames.at(i)};
        ClusterList outputClusterList;

        if (pClusterList && !pClusterList->empty())
        {
            for (const Cluster *const pCluster : *pClusterList)
            {
                if (pCluster->IsAvailable())
                    outputClusterList.emplace_back(pCluster);
            }
        }
        PandoraContentApi::SaveList<ClusterList>(*this, inputClusterListName, outputClusterListName, outputClusterList);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LeftoverClusterAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "ClusterListNames", m_clusterListNames));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "OutputClusterListNames", m_outputClusterListNames));

    if (m_clusterListNames.size() != m_outputClusterListNames.size())
    {
        std::cout << "LeftoverClusterAlgorithm: Mismatch between input and output cluster list sizes" << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
