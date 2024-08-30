/**
 *  @file   larpandoradlcontent/LArReclustering/DlGraphClusteringAlgorithm.cc
 *
 *  @brief  Implementation of the GNN-based reclustering algorithm
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
//#include "larpandoracontent/LArHelpers/LArVertexHelper.h"

#include "larpandoradlcontent/LArReclustering/DlGraphClusteringAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlGraphClusteringAlgorithm::DlGraphClusteringAlgorithm() :
    m_trainingMode{false},
    m_trainingOutputFile{""},
    m_event{-1},
    m_writeTree{false},
    m_rng(static_cast<std::mt19937::result_type>(std::chrono::high_resolution_clock::now().time_since_epoch().count()))
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

DlGraphClusteringAlgorithm::~DlGraphClusteringAlgorithm()
{
    if (m_writeTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_rootTreeName, m_rootFileName, "RECREATE"));
        }
        catch (StatusCodeException e)
        {
            std::cout << "DlGraphClusteringAlgorithm: Unable to write to ROOT tree" << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlGraphClusteringAlgorithm::Run()
{
    if (m_trainingMode)
        return this->PrepareTrainingSample();
    else
        return this->Infer();

    return STATUS_CODE_SUCCESS;

//    return this->Test();
}

StatusCode DlGraphClusteringAlgorithm::PrepareTrainingSample()
{

     //InitializeReclustering in has set the hits to be reclustered as current 
     const CaloHitList *pCaloHitList = nullptr; 
     PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList)); 

     if (pCaloHitList->empty()) 
         return STATUS_CODE_SUCCESS;

     LArMvaHelper::MvaFeatureVector featureVector;

     unsigned long nHits{0};
     std::vector<int> mainMcParticleIndices;
     for (const CaloHit *const pCaloHit : *pCaloHitList)
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

	 //Count fraction of occurrences of second most contributing MC particle IDs
  	 std::map<int, int> mainMcParticleOccurrencesMap;
     std::vector<int> occurrences;
     for(int id : mainMcParticleIndices) {mainMcParticleOccurrencesMap[id]++;}
     for(const auto& pair : mainMcParticleOccurrencesMap) {occurrences.push_back(pair.second);}
     std::sort(occurrences.begin(), occurrences.end(), std::greater<int>());
     int secondMostFrequentOccurrence(0);
     if(occurrences.size() > 1) {secondMostFrequentOccurrence = occurrences[1];} // Assuming counts has at least 2 different elements
     else return STATUS_CODE_SUCCESS;

     double secondMostFrequentParticleHitFraction = (double)secondMostFrequentOccurrence/mainMcParticleIndices.size();
     if(secondMostFrequentParticleHitFraction>0.1) return STATUS_CODE_SUCCESS; //replace magic number with a variable!

     const std::string trainingFilename{m_trainingOutputFile + ".csv"};
     featureVector.insert(featureVector.begin(), static_cast<double>(nHits));

     LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlGraphClusteringAlgorithm::Infer()
{

     //InitializeReclustering in has set the hits to be reclustered as current 
     const CaloHitList *pCaloHitList = nullptr; 
     PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList)); 

     if (pCaloHitList->empty()) 
         return STATUS_CODE_SUCCESS;

    LArDLHelper::TorchInput inputNodeFeatures;
    LArDLHelper::TorchInput inputEdgeIndices;
    int nHits(pCaloHitList->size());
    int nEdges = nHits*nHits;

    LArDLHelper::InitialiseInput({nHits,4}, inputNodeFeatures);
    LArDLHelper::InitialiseInput({2,nEdges}, inputEdgeIndices, at::kLong);

    // Run the input through the trained model
    LArDLHelper::TorchInputVector inputs;

    auto node_accessor = inputNodeFeatures.accessor<float, 2>();
    auto edge_accessor = inputEdgeIndices.accessor<long, 2>();
    
    std::vector<int> mainMcParticleVector, mainMcParticlePdgCodeVector;
    unsigned long iHit{0};
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {

        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
        int mainMcParticleIndex = this->GetMainMcParticleIndex(pParentCaloHit);
        int mainMcParticlePdgCode = this->GetMainMcParticlePdgCode(pParentCaloHit);
        double caloHitX = pCaloHit->GetPositionVector().GetX();  //check the actual functions
        double caloHitY = pCaloHit->GetPositionVector().GetY();
        double caloHitZ = pCaloHit->GetPositionVector().GetZ();
        double caloHitAdc = pCaloHit->GetMipEquivalentEnergy();
        node_accessor[iHit][0] += caloHitX;
        node_accessor[iHit][1] += caloHitY;
        node_accessor[iHit][2] += caloHitZ;
        node_accessor[iHit][3] += caloHitAdc;
        mainMcParticleVector.push_back(mainMcParticleIndex);
        mainMcParticlePdgCodeVector.push_back(mainMcParticlePdgCode);

        iHit++;
     }

     //Fill edges tensor
     for(int kHit=0; kHit<nHits; kHit++)
     {
         for(int jHit=0; jHit<nHits; jHit++)
         {
             edge_accessor[0][kHit*nHits+jHit] += kHit;
             edge_accessor[1][kHit*nHits+jHit] += jHit;
         }

      }

    inputs.push_back(inputNodeFeatures);
    inputs.push_back(inputEdgeIndices);

    // Run the input through the trained model
    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_model, inputs, output);

    for(int kHit=0; kHit<nHits; kHit++)
    {
        torch::Tensor kHitNodeEmbeddings = output[kHit];
        for(int jHit=0; jHit<nHits; jHit++)
        {
            
            torch::Tensor jHitNodeEmbeddings = output[jHit];
            //Normalising the embeddings vectors prior to dot product
            kHitNodeEmbeddings = kHitNodeEmbeddings/kHitNodeEmbeddings.norm();
            jHitNodeEmbeddings = kHitNodeEmbeddings/jHitNodeEmbeddings.norm();
            //double kjEdgeScore=torch::dot(kHitNodeEmbeddings, jHitNodeEmbeddings).item<double>();
            //bool trueEdgeScore = (mainMcParticleVector.at(kHit)==mainMcParticleVector.at(jHit))?1:0;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

int DlGraphClusteringAlgorithm::GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit)
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

int DlGraphClusteringAlgorithm::GetMainMcParticlePdgCode(const pandora::CaloHit *const pCaloHit)
{
    const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
    const int pdg{std::abs(pMCParticle->GetParticleId())};
	return pdg;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlGraphClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
       if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        if (m_writeTree)
        {
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootTreeName", m_rootTreeName));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "RootFileName", m_rootFileName));
        }
        else
        {
            std::string modelName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileName", modelName));
            modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
            LArDLHelper::LoadModel(modelName, m_model);
        }
    }

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
