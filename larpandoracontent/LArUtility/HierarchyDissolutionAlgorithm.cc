/**
 *  @file   larpandoracontent/LArUtility/HierarchyDissolutionAlgorithm.cc
 *
 *  @brief  Implementation of the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArUtility/HierarchyDissolutionAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode HierarchyDissolutionAlgorithm::Run()
{

    for (std::string listName : m_pfoListNames)
    {
        const PfoList *pPfoList{nullptr};
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, listName, pPfoList));

        if (pPfoList && !pPfoList->empty())
        {
            for (const ParticleFlowObject *pParent : *pPfoList)
            {
                // ATTN - RemovePfoParentDaughterRelationship modifies the daughter list, so you can't use a const reference of the daughter
                // PFO list, you must use a copy
                const PfoList childList(pParent->GetDaughterPfoList());
                for (const ParticleFlowObject *pChild : childList)
                    PandoraContentApi::RemovePfoParentDaughterRelationship(*this, pParent, pChild);
            }
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyDissolutionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
