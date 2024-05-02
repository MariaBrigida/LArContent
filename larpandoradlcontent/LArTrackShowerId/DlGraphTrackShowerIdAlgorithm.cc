/**
 *  @file   larpandoradlcontent/LArReclustering/DlGraphTrackShowerIdAlgorithm.cc
 *
 *  @brief  Implementation of the GNN-based track-shower id algorithmw
 *
 *  $Log: $
 */

#include <chrono>
#include <cmath>

#include <torch/script.h>
#include <torch/torch.h>

#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoradlcontent/LArTrackShowerId/DlGraphTrackShowerIdAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlGraphTrackShowerIdAlgorithm::DlGraphTrackShowerIdAlgorithm() :
    m_trainingMode{false},
    m_trainingOutputFile{""},
    m_event{-1},
    m_rng(static_cast<std::mt19937::result_type>(std::chrono::high_resolution_clock::now().time_since_epoch().count()))
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlGraphTrackShowerIdAlgorithm::~DlGraphTrackShowerIdAlgorithm()
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlGraphTrackShowerIdAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    //else
    //    return this->Infer();

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlGraphTrackShowerIdAlgorithm::PrepareTrainingSample()
{

   	const std::string targetListName = "AllParticles3D";
    const PfoList *pPfoList(nullptr);

	const StatusCode statusCode(PandoraContentApi::SaveList<Pfo>(*this, "ShowerParticles3D", targetListName));
	const StatusCode statusCode2(PandoraContentApi::SaveList<Pfo>(*this, "TrackParticles3D", targetListName));

    if(STATUS_CODE_SUCCESS != statusCode) std::cout << "WARNING: the shower pfo list is empty" << std::endl;
    if(STATUS_CODE_SUCCESS != statusCode2) std::cout << "WARNING: the track pfo list is empty" << std::endl;

    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, "AllParticles3D", pPfoList)) && pPfoList)
    {
         if (pPfoList->size()==0) return STATUS_CODE_NOT_FOUND;
 
         for (const Pfo *const pPfo : *pPfoList)
         {
             ClusterList clusterList3D;
             LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);
 
             if (clusterList3D.empty())
                 continue;
 
             CaloHitList caloHitList3D;
             clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

		     if (caloHitList3D.empty()) 
			     continue;
        	 
             LArMvaHelper::MvaFeatureVector featureVector;

		     unsigned long nHits{0};
		     std::vector<int> mainMcParticleIndices;
		     for (const CaloHit *const pCaloHit : caloHitList3D)
		     {

		         const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
		         int mainMcParticleIndex = this->GetMainMcParticleIndex(pParentCaloHit);
		         int mainMcParticlePdgCode = this->GetMainMcParticlePdgCode(pParentCaloHit);
		         double caloHitX = pCaloHit->GetPositionVector().GetX();
		         double caloHitY = pCaloHit->GetPositionVector().GetY();
		         double caloHitZ = pCaloHit->GetPositionVector().GetZ();
		         double caloHitAdc = pCaloHit->GetMipEquivalentEnergy();
		
		         mainMcParticleIndices.emplace_back(mainMcParticleIndex);
		         featureVector.emplace_back(mainMcParticleIndex);
		         featureVector.emplace_back(mainMcParticlePdgCode);
		         featureVector.emplace_back(caloHitX);
		         featureVector.emplace_back(caloHitY);
		         featureVector.emplace_back(caloHitZ);
		         featureVector.emplace_back(caloHitAdc);
		
		         nHits++;
     		 }

		     const std::string trainingFilename{m_trainingOutputFile + ".csv"};
		     featureVector.insert(featureVector.begin(), static_cast<double>(nHits));
		
		     LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);
		}
	}	
    return STATUS_CODE_SUCCESS;
}


//-----------------------------------------------------------------------------------------------------------------------------------------
//Get vector of MC particles, sort them, loop into the vector and find in map (see "SortMCParticle")
int DlGraphTrackShowerIdAlgorithm::GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit)
{
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
    MCParticleVector mcParticleVector(pMCParticleList->begin(),pMCParticleList->end());
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);
    int iMcPart(0); 
    for (const auto &weightMapEntry : pCaloHit->GetMCParticleWeightMap()) 
    { 
      if(weightMapEntry.second>0.5) 
      { 
        iMcPart=0;  
        for(const MCParticle *const pMCParticle: mcParticleVector) 
        { 
            if(pMCParticle==weightMapEntry.first) { break;} 
            iMcPart++; 
        }
      } 
    }
    return iMcPart;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

int DlGraphTrackShowerIdAlgorithm::GetMainMcParticlePdgCode(const pandora::CaloHit *const pCaloHit)
{
    const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
    const int pdg{std::abs(pMCParticle->GetParticleId())};
    return pdg;
   
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlGraphTrackShowerIdAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        std::string modelName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", modelName));
        modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelName, m_model);
    }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    return STATUS_CODE_SUCCESS;
}
} // namespace lar_dl_content
