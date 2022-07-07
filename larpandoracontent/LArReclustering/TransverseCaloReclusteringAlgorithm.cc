/**
 *  @file   larpandoracontent/LArReclustering/TransverseCaloReclusteringAlgorithm.cc
 *
 *  @brief  Implementation file for the reclustering algorithm that uses transverse calorimetric profiles.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArReclustering/TransverseCaloReclusteringAlgorithm.h"


using namespace pandora;

namespace lar_content
{

TransverseCaloReclusteringAlgorithm::TransverseCaloReclusteringAlgorithm()
{
}

StatusCode TransverseCaloReclusteringAlgorithm::Run()
{

    //Access the clusters
    //Split the clusters using truth info
    //Return the new cluster list

    //InitializeReclustering in has set the hits to be reclustered as current
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    if (pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

//    OrderedCaloHitList orderedCaloHitList;
//    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, orderedCaloHitList.Add(*pCaloHitList));

/*    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    MCParticleVector mcParticleVector; 
    for (const MCParticle *pMCParticle : *pMCParticleList) mcParticleVector.push_back(pMCParticle);*/


    ClusterVector clusterVector;

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));   
    PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), pCaloHitList, "CurrentClusterHits", BLUE));
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {

        //std::cout << "debug calo hit loop. clusterVector.size = " << clusterVector.size() << std::endl;
        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
        //if (TPC_VIEW_W != pParentCaloHit->GetHitType()) 
        //    continue; 
        //Find main contributing MC particle Id 
        //MCParticleVector mcParticleVector; 
        int mainMcParticleIndex = this->GetMainMcParticleIndex(pParentCaloHit);
        //std::cout <<"mainMcParticleIndex = " << mainMcParticleIndex << std::endl;


        //Loop over clusterVector entries
        const Cluster *pMainMCParticleCluster(nullptr);
        for (ClusterVector::const_iterator cluIter = clusterVector.begin(), cluIterEnd = clusterVector.end(); cluIter != cluIterEnd; ++cluIter)
        {
             const Cluster *const pCluster = *cluIter;
             //I only need to check the parent MC particle of one hit in the cluster. I'll take the front 
             //const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
             CaloHitList clusterCaloHitList;;
             pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterCaloHitList);
             //std::cout << "mainMcParticleIndex = " << mainMcParticleIndex << " this->GetMainMcParticleIndex(clusterCaloHitList.front()) = " << this->GetMainMcParticleIndex(clusterCaloHitList.front()) << std::endl;
             const CaloHit *const pFrontHitParentCaloHit = static_cast<const CaloHit *>(clusterCaloHitList.front()->GetParentAddress());
             if(this->GetMainMcParticleIndex(pFrontHitParentCaloHit)==mainMcParticleIndex) pMainMCParticleCluster=pCluster;
             
        }

        //If pMainMCParticle is null cluster (a cluster with matching MC particle hits is not found), create one
        if(!pMainMCParticleCluster)      
        {
            //std::cout << "debug calo hit loop: cluster with this mc hits index not found. clusterVector.size = " << clusterVector.size() << std::endl;
            const Cluster *pCluster = nullptr;
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList.push_back(pCaloHit);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
            //PandoraContentApi::AddToCluster(const pandora::Algorithm &algorithm, const pandora::Cluster *const pCluster, const T *const pT)
            clusterVector.push_back(pCluster);
        }
        else
        {
            //Attach calo hit to pMainMCParticleCluster
            //std::cout << "debug calo hit loop: cluster with this mc hits index was found and I will attach hits. clusterVector.size = " << clusterVector.size() << std::endl;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pMainMCParticleCluster, pCaloHit));
            //PandoraContentApi::AddToCluster(const pandora::Algorithm &algorithm, const pandora::Cluster *const pCluster, const T *const pT)
        }

    }

    std::cout << "At the end of the reclustering algorithm, clusterVector has size = " << clusterVector.size() << std::endl;
    //Now I want to display all of these new clusters!

    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));   
    for(const Cluster* pNewCluster: clusterVector)
    { 
      ClusterList newClusters;
      newClusters.push_back(pNewCluster);  
      PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &newClusters, "newClusters", AUTOITER));
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    return STATUS_CODE_SUCCESS;
}

int TransverseCaloReclusteringAlgorithm::GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit)
{
    //const CaloHit *const pParentCaloHit(static_cast<const CaloHit *>(pCaloHit->GetParentAddress()));

    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
//    std::cout << "debug pMCParticleList size = " << pMCParticleList->size() << std::endl;
//    std::cout << "debug pParentCaloHit->GetMCParticleWeightMap() size = " << pParentCaloHit->GetMCParticleWeightMap().size() << std::endl;
    int /*mcParticleIndex(-999),*/iMcPart(0); 
    for (const auto &weightMapEntry : pCaloHit->GetMCParticleWeightMap()) 
    { 
      if(weightMapEntry.second>0.5) 
      { 
//        std::cout << "particle ID = " << weightMapEntry.first << std::endl; 
        iMcPart=0;  
        for(const MCParticle *const pMCParticle: *pMCParticleList) 
        { 
            if(pMCParticle==weightMapEntry.first) {/*mcParticleIndex=iMcPart;*/ break;} 
            iMcPart++; 
        }
      } 
    }
    return iMcPart;
}


StatusCode TransverseCaloReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
