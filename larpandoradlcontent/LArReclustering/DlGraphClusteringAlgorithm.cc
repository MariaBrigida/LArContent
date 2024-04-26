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
    //m_rootTreeName{""},
    //m_rootFileName{""},
    //m_caloHitListNames{""},
    m_rng(static_cast<std::mt19937::result_type>(std::chrono::high_resolution_clock::now().time_since_epoch().count()))
{
}

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

     std::cout << "Hello" << std::endl;
  
     if (pCaloHitList->empty()) 
         return STATUS_CODE_SUCCESS;

     //I add a cut on the ratio between the number of hits contributed by the second most contributing MC particle
     //if (!this->PassesCuts(pCaloHitList))
     //    return STATUS_CODE_SUCCESS;

     LArMvaHelper::MvaFeatureVector featureVector;

     unsigned long nHits{0};
     std::vector<int> mainMcParticleIndices;
     for (const CaloHit *const pCaloHit : *pCaloHitList)
     {

         const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
         //if (TPC_VIEW_W != pParentCaloHit->GetHitType()) 
         //    continue; 
         //Find main contributing MC particle Id 
         int mainMcParticleIndex = this->GetMainMcParticleIndex(pParentCaloHit);
         //In reality I don't need the mainMcParticle index but a label: 1, 2...
         //Do I have to select showers that are only composed of 2 particles? In developing this algo I will assume so.
         double caloHitX = pCaloHit->GetPositionVector().GetX();  //check the actual functions
         double caloHitY = pCaloHit->GetPositionVector().GetY();
         double caloHitZ = pCaloHit->GetPositionVector().GetZ();
         double caloHitAdc = pCaloHit->GetMipEquivalentEnergy();
         std::cout << "calo hit mc particle index = " << mainMcParticleIndex << " coord x = " << caloHitX << " y = " << caloHitY << " z = " << caloHitZ << " adc = " << caloHitAdc << " parent hit type = " << pParentCaloHit->GetHitType() <<  std::endl;

		 mainMcParticleIndices.emplace_back(mainMcParticleIndex); 
         featureVector.emplace_back(mainMcParticleIndex);
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
     //Sort counts in descending order
     std::sort(occurrences.begin(), occurrences.end(), std::greater<int>());

     // Find the second most frequent count, if available
     int secondMostFrequentOccurrence(0);
     //for(auto occ: occurrences) std::cout << "occurrences = " << occ << std::endl;
     if(occurrences.size() > 1) {secondMostFrequentOccurrence = occurrences[1];} // Assuming counts has at least 2 different elements
     else return STATUS_CODE_SUCCESS;

	// std::cout << "secondMostFrequentOccurrence/mainMcParticleIndices.size() = "<< secondMostFrequentOccurrence/mainMcParticleIndices.size() << std::endl;
    // std::cout << "secondMostFrequentOccurrence = " << secondMostFrequentOccurrence << " mainMcParticleIndices.size() = " << mainMcParticleIndices.size() << std::endl;
     double secondMostFrequentParticleHitFraction = (double)secondMostFrequentOccurrence/mainMcParticleIndices.size();
     std::cout << "secondMostFrequentParticleHitFraction = " << secondMostFrequentParticleHitFraction << std::endl;
     if(secondMostFrequentParticleHitFraction<0.3) return STATUS_CODE_SUCCESS; //replace 0.3 with a variable!!

     const std::string trainingFilename{m_trainingOutputFile + ".csv"};
     //featureVector.insert(featureVector.begin() + 8, static_cast<double>(nHits));
     featureVector.insert(featureVector.begin(), static_cast<double>(nHits));
     std::cout << "nHits = " << nHits << std::endl;
     //if (nHits > 10)


     LArMvaHelper::ProduceTrainingExample(trainingFilename, true, featureVector);

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlGraphClusteringAlgorithm::Test(){
    // Create a tensor to scatter values into
    //auto tensor = torch::zeros({3, 5}, torch::dtype(torch::kFloat32).device(device));
    auto tensor = torch::zeros({3, 5}, torch::dtype(torch::kFloat32));

    // Create a tensor with values to scatter
    //auto src = torch::tensor({1, 2, 3, 4, 5}, torch::dtype(torch::kFloat32).device(device));
    auto src = torch::tensor({1, 2, 3, 4, 5}, torch::dtype(torch::kFloat32));

    // Create an index tensor indicating where to scatter the src values
    //auto index = torch::tensor({0, 1, 2, 0, 1}, torch::dtype(torch::kInt64).device(device));
    auto index = torch::tensor({0, 1, 2, 0, 1}, torch::dtype(torch::kInt64));

    // Perform the scatter operation (dim, index, src)
    // This will scatter the src values into the tensor according to the index
    tensor.scatter_(0, index.unsqueeze(0), src.unsqueeze(0));

    std::cout << tensor << std::endl;
    return STATUS_CODE_SUCCESS;
}

StatusCode DlGraphClusteringAlgorithm::Infer()
{

     //InitializeReclustering in has set the hits to be reclustered as current 
     const CaloHitList *pCaloHitList = nullptr; 
     PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList)); 

  
     if (pCaloHitList->empty()) 
         return STATUS_CODE_SUCCESS;

    //I don't know what shape my input needs to have...
    LArDLHelper::TorchInput inputNodeFeatures;
    LArDLHelper::TorchInput inputEdgeIndices;
    //tensor dimensionality.train_data.x shape =  torch.Size([67, 4])  train_data.edge_index shape =  torch.Size([2, 7722])
    //so [nhits,4] and [2,nedges]
    int nHits(pCaloHitList->size());
    int nEdges = nHits*nHits;
//    int nHits = 3;
//    int nEdges = 9;
    std::cout << "nHits = " << nHits << std::endl;
    std::cout << "nEdges = " << nEdges << std::endl;

    LArDLHelper::InitialiseInput({nHits,4}, inputNodeFeatures);
    LArDLHelper::InitialiseInput({2,nEdges}, inputEdgeIndices, at::kLong);

    // Run the input through the trained model
    LArDLHelper::TorchInputVector inputs;

    auto node_accessor = inputNodeFeatures.accessor<float, 2>();
    auto edge_accessor = inputEdgeIndices.accessor<long, 2>();
    
    //std::cout << " inputNodeFeatures = " << inputNodeFeatures << std::endl;
    //std::cout << "inputEdgeIndices = " << inputEdgeIndices << std::endl;


     /*std::cout << "Hello 3B" << std::endl;
    try{
    	auto edge_accessor = inputEdgeIndices.accessor<int, 2>();
        edge_accessor[0][0] += 1;
//   	} catch (const std::exception& e) {
//    	std::cerr << "Standard exception: " << e.what() << std::endl;
	} catch (const c10::Error& e) {
    	std::cerr << "C10/Torch exception: " << e.msg() << std::endl;
	} catch (...) {
    	std::cerr << "An unknown exception occurred" << std::endl;
	}
    
    	auto edge_accessor = inputEdgeIndices.accessor<int, 2>();*/
    std::vector<int> mainMcParticleVector;
    unsigned long iHit{0};
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {

        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
        //if (TPC_VIEW_W != pParentCaloHit->GetHitType()) 
     //    continue; 
        //Find main contributing MC particle Id 
        int mainMcParticleIndex = this->GetMainMcParticleIndex(pParentCaloHit);
        //In reality I don't need the mainMcParticle index but a label: 1, 2...
        //Do I have to select showers that are only composed of 2 particles? In developing this algo I will assume so.
        double caloHitX = pCaloHit->GetPositionVector().GetX();  //check the actual functions
        double caloHitY = pCaloHit->GetPositionVector().GetY();
        double caloHitZ = pCaloHit->GetPositionVector().GetZ();
        double caloHitAdc = pCaloHit->GetMipEquivalentEnergy();
        node_accessor[iHit][0] += caloHitX;
        node_accessor[iHit][1] += caloHitY;
        node_accessor[iHit][2] += caloHitZ;
        node_accessor[iHit][3] += caloHitAdc;
        mainMcParticleVector.push_back(mainMcParticleIndex);
        //std::cout << "calo hit mc particle index = " << mainMcParticleIndex << " coord x = " << caloHitX << " y = " << caloHitY << " z = " << caloHitZ << " adc = " << caloHitAdc <<  std::endl;

         iHit++;
     }


/*     for(int mHit=0; mHit<nHits; mHit++){
		node_accessor[mHit][0] += mHit;
		node_accessor[mHit][1] += 10*mHit;
		node_accessor[mHit][2] += 100*mHit;
		node_accessor[mHit][3] += 1000*mHit;
     }
*/

     //Fill edges tensor
     for(int kHit=0; kHit<nHits; kHit++)
     {
         for(int jHit=0; jHit<nHits; jHit++)
         {
             edge_accessor[0][kHit*nHits+jHit] += kHit;
             edge_accessor[1][kHit*nHits+jHit] += jHit;
         }

      }

//whatever this is for me
//        this->MakeNetworkInputFromHits(*pCaloHitList, view, driftMin, driftMax, wireMin[view], wireMax[view], input, pixelVector);

    //LArDLHelper::InitialiseInput({1, 1, m_height, m_width}, input);
    //auto accessor = networkInput.accessor<float, 4>();
    inputs.push_back(inputNodeFeatures);
    inputs.push_back(inputEdgeIndices);
    // Run the input through the trained model
    LArDLHelper::TorchOutput output;
    LArDLHelper::Forward(m_model, inputs, output);

    
	//CALCULATE EDGE SCORES, and print true score vs edge score!
    //This method uses a dot product. Alternatively, look into other methods e.d. concatenate embeddings and pass them to a MLP
    for(int kHit=0; kHit<nHits; kHit++)
    {
        torch::Tensor kHitNodeEmbeddings = output[kHit];
        for(int jHit=0; jHit<nHits; jHit++)
        {
            
            torch::Tensor jHitNodeEmbeddings = output[jHit];
            //Normalising the embeddings vectors prior to dot product
            kHitNodeEmbeddings = kHitNodeEmbeddings/kHitNodeEmbeddings.norm();
            jHitNodeEmbeddings = kHitNodeEmbeddings/jHitNodeEmbeddings.norm();
            double kjEdgeScore=torch::dot(kHitNodeEmbeddings, jHitNodeEmbeddings).item<double>();
            std::cout << "----------------- k = " << kHit << ", j = " << jHit << "----------------" << std::endl;
            std::cout << "kHitNodeEmbeddings = " << kHitNodeEmbeddings.view({1, -1}) << std::endl;
            std::cout << "jHitNodeEmbeddings = " << jHitNodeEmbeddings.view({1, -1}) << std::endl;
            std::cout << "mainMcParticleVector.at(kHit) = " << mainMcParticleVector.at(kHit) << " mainMcParticleVector.at(jHit) = " << mainMcParticleVector.at(jHit) << std::endl;
            bool trueEdgeScore = (mainMcParticleVector.at(kHit)==mainMcParticleVector.at(jHit))?1:0;
            std::cout << "true kjEdgeScore = " << trueEdgeScore << " predicted kjEdgeScore = " << kjEdgeScore << std::endl;
        }
    }



	//auto shape = output.sizes();
    //std::cout << "output = " << output << std::endl;
    //std::cout << "shape = " << shape << std::endl;
    /*auto outputAccessor{output.accessor<float, 2>()};
    for(int iHit=0; iHit<nHits; iHit++)
    {
        for(int jHit=0; jHit<nHits; jHit++){
           auto outputAccessor{output.accessor<float, 1>()};
           const float value{outputAccessor[i]};

        }

    }*/
    //auto classes{torch::argmax(output,1)};
    //auto classesAccessor{classes.accessor<long,3>()};

//What dimensionality does the torch tensor output of this model have?

//This is how I should deal with the output of this model
//h_src = h[batch.edge_label_index[0]]
//h_dst = h[batch.edge_label_index[1]]
//pred = (h_src * h_dst).sum(dim=-1)

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

//Get vector of MC particles, sort them, loop into the vector and find in map (see "SortMCParticle")
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



/*StatusCode DlGraphClusteringAlgorithm::MakeNetworkInputFromHits(const CaloHitList &caloHits, const HitType view, const float xMin,
    const float xMax, const float zMin, const float zMax, LArDLHelper::TorchInput &networkInput, PixelVector &pixelVector) const
{
    return STATUS_CODE_SUCCESS;
}
*/
//-----------------------------------------------------------------------------------------------------------------------------------------

/*StatusCode DlGraphClusteringAlgorithm::GetMCToHitsMap(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const
{
    const CaloHitList *pCaloHitList2D(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList2D));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_maxPhotonPropagation = std::numeric_limits<float>::max();
    LArMCParticleHelper::SelectReconstructableMCParticles(
        pMCParticleList, pCaloHitList2D, parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToHitsMap);

    return STATUS_CODE_SUCCESS;
}
*/

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlGraphClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));
    //PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));
       if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TrainingOutputFileName", m_trainingOutputFile));
    }
    else
    {
        //std::string modelName;
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelFileNameU", modelName));
        //modelName = LArFileHelper::FindFileInPath(modelName, "FW_SEARCH_PATH");
        //PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
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
 //       PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputVertexListName", m_outputVertexListName));
    }

    //think I don't need this
    //PANDORA_RETURN_RESULT_IF_AND_IF(
       // STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "CaloHitListNames", m_caloHitListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    return STATUS_CODE_SUCCESS;
}

//-----------------------------------------------------------------------------------------------------------------------------------------

//StatusCode DlGraphClusteringAlgorithm::PassesCuts(const CaloHitList *pCaloHitList)
//{
//}

//-----------------------------------------------------------------------------------------------------------------------------------------

} // namespace lar_dl_content
