/**
 *  @file   larpandoracontent/LArCheating/MyCheatingVertexSelectionAlgorithm.cc
 *
 *  @brief  Implementation of the cheating vertex selection algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/MyCheatingVertexSelectionAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"


using namespace pandora;

namespace lar_content
{

MyCheatingVertexSelectionAlgorithm::MyCheatingVertexSelectionAlgorithm() : m_replaceCurrentVertexList(true)
{
}

//MyCheatingVertexSelectionAlgorithm::~MyCheatingVertexSelectionAlgorithm(){};

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyCheatingVertexSelectionAlgorithm::Run()
{

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));


    const VertexList *pCandidateVertexList = NULL;
    if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, m_candidateVertexListName, pCandidateVertexList))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "MyCheatingVertexSelectionAlgorithm: candidate vertex list " << m_candidateVertexListName << " unavailable." << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }


    LArMCParticleHelper::MCContributionMap mcToHitsMap; 
    MCParticleVector primaries; 
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaries); 
    const MCParticle *pTrueNeutrino{nullptr}; 
    if (!primaries.empty())
    {
        for (const MCParticle *primary : primaries)
        {
            const MCParticleList &parents{primary->GetParentList()};
            if (parents.size() == 1 && LArMCParticleHelper::IsNeutrino(parents.front()))
            {
                pTrueNeutrino = parents.front();
                break;
            }
        }
    }
    const CartesianVector &trueVertexPosition{pTrueNeutrino->GetVertex()};

    const VertexList *pOutputVertexList{nullptr};
    std::string temporaryListName;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pOutputVertexList, temporaryListName));



    double dist(999999);
    CartesianVector closestVertexPosition(0.,0.,0.);

    for (VertexList::const_iterator iter = pCandidateVertexList->begin(), iterEnd = pCandidateVertexList->end(); iter != iterEnd; ++iter)
    {
        const Vertex *pCandidateVertex = *iter;

        const CartesianVector &candidateVertexPosition{pCandidateVertex->GetPosition()};

        double newdist = (trueVertexPosition-candidateVertexPosition).GetMagnitude(); 
        if(newdist<dist)
        {
          dist=newdist;
          closestVertexPosition=candidateVertexPosition;
        }

    }
    PandoraContentApi::Vertex::Parameters parameters;
    parameters.m_position = closestVertexPosition;
    parameters.m_vertexLabel = VERTEX_INTERACTION;
    parameters.m_vertexType = VERTEX_3D;

    const Vertex *pVertex(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, parameters, pVertex));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, m_outputVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, m_outputVertexListName));
   
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MyCheatingVertexSelectionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ReplaceCurrentVertexList", m_replaceCurrentVertexList));


    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CandidateVertexListNames", m_candidateVertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
