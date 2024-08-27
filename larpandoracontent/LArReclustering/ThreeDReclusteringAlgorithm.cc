/**
 *  @file   larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.cc
 *
 *  @brief  Implementation file for the reclustering algorithm that runs other algs.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Managers/ClusterManager.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h"

bool sortByCaloHits (pandora::CaloHitList a, pandora::CaloHitList b) { return (a.size()>b.size()); }

using namespace pandora;

namespace lar_content
{

ThreeDReclusteringAlgorithm::ThreeDReclusteringAlgorithm():
    m_visualDisplaysOn(false),
    m_hitThresholdForNewPfo(0)
{
}

ThreeDReclusteringAlgorithm::~ThreeDReclusteringAlgorithm()
{
}

StatusCode ThreeDReclusteringAlgorithm::Run()
{
    std::cout << "--------------------------------------------- NEW EVENT ------------------------------------------------" << std::endl;
    // Get shower pfos and then find the 3D cluster in the shower pfo.
    const PfoList *pShowerPfoList(nullptr);
    std::string initialPfosListName;
    m_newPfosListNameAllAfterReclustering = "newShowerParticles3D";
    std::string newPfosListNameUnchanged = "unchangedShowerParticles3D";
    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, "ShowerParticles3D", pShowerPfoList)) && pShowerPfoList->size())
    {
        std::cout << "In this event there are " << pShowerPfoList->size() << " shower pfos." << std::endl;
        if (pShowerPfoList->size()==0) return STATUS_CODE_NOT_FOUND;
        PfoList unchangedPfoList;   

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Pfo>(*this, initialPfosListName));

        //Some pfos are shower-like and yet include track-like 3D clusters. For the moment I don't want to deal with these.
        const ClusterList *pShowerClusters(nullptr);
        PandoraContentApi::GetList(*this, "ShowerClusters3D", pShowerClusters);
        if(!pShowerClusters) return STATUS_CODE_NOT_FOUND;


        int iPfo(0); //for debug
        //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
        for (const Pfo *const pShowerPfo : *pShowerPfoList)
        {
            //std::cout << "iPfo = " << iPfo << std::endl;

            ClusterList clusterList3D;
            LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);
            //PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterList3D, "ThisPfo3DClusters", RED));
            //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
           
            if(!this->PassesCutsForReclustering(pShowerPfo))
            {
                unchangedPfoList.push_back(pShowerPfo);
				std::cout << "This pfo does not pass cuts for reclustering" << std::endl;
                continue; // this function just checks it's a shower at the moment
            }
            if (clusterList3D.empty())
                continue;

            //Some pfos are shower-like and yet include track-like 3D clusters. For the moment I don't want to deal with these.
            if(pShowerClusters->end() == std::find(pShowerClusters->begin(), pShowerClusters->end(), clusterList3D.front())) continue;

            CaloHitList caloHitList3D;
            clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

            //Quality cuts
            if (caloHitList3D.size() < 2)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);


			//Calculate initial figure of merit 
            float initialFigureOfMerit=this->GetFigureOfMerit(caloHitList3D);
            //std::cout << "initialFigureOfMerit = " << initialFigureOfMerit << std::endl;
			
			//Create a variable for the minimum figure of merit and initialize to initial FOM
            float minimumFigureOfMerit(initialFigureOfMerit);

            //Free the hits in this cluster, so that they are available for reclustering!
            //Ask to remove the 3D cluster from the parent pfo, so that it's not owned any more
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pShowerPfo, clusterList3D.front()));
            //Pop this cluster in a local clusterlist
            //const ClusterList reclusterClusterList(1, clusterList3D.front());
            //const TrackList reclusterTrackList; //dummy track list

            // Initialize reclustering with these local lists
            std::string currentClustersListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentClustersListName));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, "ShowerClusters3D"));

            // Specify clusters and tracks to be used in reclustering
            //std::string originalClustersListName;
            //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, reclusterTrackList, reclusterClusterList, originalClustersListName));

			//I want to split the calo hit list into a new set of calo hit lists, taking the best outcome out of different algorithms
            //std::cout << "m_algorithmToolVector.size() = " << m_algorithmToolVector.size() << std::endl;
            std::vector<CaloHitList*> newCaloHitListsVector, minimumFigureOfMeritCaloHitListsVector, initialCaloHitListsVector;
            //newCaloHitListsVector.emplace_back(&caloHitList3D);
            initialCaloHitListsVector.emplace_back(&caloHitList3D);
            minimumFigureOfMeritCaloHitListsVector = initialCaloHitListsVector;

            for (auto toolIter = m_algorithmToolVector.begin(); toolIter != m_algorithmToolVector.end(); ++toolIter)
            {
				//std::cout << "ThreeDReclusteringAlgoritm : running a new clustering tool" << std::endl;
                // Split cluster into a set of new CaloHit lists
                //newCaloHitListsVector.emplace_back(&caloHitList3D);
               // minimumFigureOfMeritCaloHitListsVector = newCaloHitListsVector;
                //std::cout << "start of new tool: minimumFigureOfMeritCaloHitListsVector.size() = " << minimumFigureOfMeritCaloHitListsVector.size() << std::endl; 
                newCaloHitListsVector = initialCaloHitListsVector;

                try {
                    (*toolIter)->Run(this, newCaloHitListsVector);
                } catch (const std::exception &e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                } catch (...) {
                    std::cerr << "Unknown exception caught!" << std::endl;
                }

                //Calculate FOM for this set of new CaloHitLists
				float newFigureOfMerit = this->GetFigureOfMerit(newCaloHitListsVector);
                //std::cout << "newFigureOfMerit = " << newFigureOfMerit << " current newCaloHitListsVector.size() = " << newCaloHitListsVector.size() << " minimumFigureOfMeritCaloHitListsVector.size() = " << minimumFigureOfMeritCaloHitListsVector.size() << std::endl;                

                //Is this FOM smaller?
                if(newFigureOfMerit <= minimumFigureOfMerit)
                {
                     minimumFigureOfMerit = newFigureOfMerit;
                     //minimumFigureOfMeritCaloHitListsVector.clear();
                     minimumFigureOfMeritCaloHitListsVector = newCaloHitListsVector;
                     //std::cout << "after having found new minimum: minimumFigureOfMeritCaloHitListsVector.size() = " << minimumFigureOfMeritCaloHitListsVector.size() << std::endl; 
                     //std::cout << "after having found new minimum: newCaloHitListsVector.size() = " << newCaloHitListsVector.size() << std::endl; 
                }
                else //If not, clear the lists and vector!!
                {
                     //Delete the CaloHitList objects after processing
                     for (CaloHitList* caloHitList : newCaloHitListsVector)
                     {
                         delete caloHitList;  // Clean up the memory
                     }
                 }
                 newCaloHitListsVector.clear();
                 //std::cout << "end of tool: minimumFigureOfMeritCaloHitListsVector.size() = " << minimumFigureOfMeritCaloHitListsVector.size() << std::endl; 
                 //std::cout << "end of tool: newCaloHitListsVector.size() = " << newCaloHitListsVector.size() << std::endl; 
            }
            //std::cout << "minimumFigureOfMerit = " << minimumFigureOfMerit << std::endl;
            //std::cout << "minimumFigureOfMeritCaloHitListsVector.size() = " << minimumFigureOfMeritCaloHitListsVector.size() << std::endl; 


			//Now, based on this outcome, I want to build new 3D clusters - but only if the cluster list changed.


           if(minimumFigureOfMeritCaloHitListsVector==newCaloHitListsVector)
           {
               PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pShowerPfo, clusterList3D.front()));
               unchangedPfoList.push_back(pShowerPfo);
           }
           else
           {
               const ClusterList reclusterClusterList(1, clusterList3D.front());
               std::string clusterListToSaveName, clusterListToDeleteName;
 
               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
               PandoraContentApi::InitializeFragmentation(*this, reclusterClusterList, clusterListToDeleteName, clusterListToSaveName));



               //ClusterVector clusterVector;
               //ClusterList newClustersList;
               for(auto list: minimumFigureOfMeritCaloHitListsVector)
               {
                   //std::cout << "list size = " << list->size() << std::endl;
                   const Cluster *pCluster = nullptr;
                   PandoraContentApi::Cluster::Parameters parameters;
                   parameters.m_caloHitList = *list;
                   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));
                   m_newClustersList.push_back(pCluster);
               }

               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));
               //std::cout << "newClustersList size = " << newClustersList.size() << std::endl;


               //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RebuildPfo(pShowerPfo, m_newClustersList));
               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RebuildPfo(pShowerPfo));


        }
        iPfo++;
    }
        if(unchangedPfoList.size()>0) PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<PfoList>(*this, "ShowerParticles3D", m_newPfosListNameAllAfterReclustering,  unchangedPfoList));
        const PfoList *pNewPfosListAllAfterReclustering;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList<PfoList>(*this, m_newPfosListNameAllAfterReclustering, pNewPfosListAllAfterReclustering));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_newPfosListNameAllAfterReclustering, "ShowerParticles3D"));
        //const PfoList *pShowerParticles3DDebug;
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList<PfoList>(*this, "ShowerParticles3D", pShowerParticles3DDebug));
       PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, initialPfosListName));

    }

    return STATUS_CODE_SUCCESS;
}

bool ThreeDReclusteringAlgorithm::PassesCutsForReclustering(const pandora::ParticleFlowObject *const pShowerPfo)
{
    if (!LArPfoHelper::IsShower(pShowerPfo)) return false;
	std::cout << "It's a shower" << std::endl;
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);

    if (clusterList3D.empty()) return false;
	std::cout << "Has 3D clusters" << std::endl;
    CaloHitList caloHitList3D;
    clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

    //Quality cuts
    if (caloHitList3D.size() < 2) return false;
	std::cout << "Has more than 2 hits" << std::endl;

    //Some pfos are shower-like and yet include track-like 3D clusters. For the moment I don't want to deal with these.
    const ClusterList *pShowerClusters(nullptr);
    PandoraContentApi::GetList(*this, "ShowerClusters3D", pShowerClusters);
    if(!pShowerClusters) return false;
    if(pShowerClusters->end() == std::find(pShowerClusters->begin(), pShowerClusters->end(), clusterList3D.front())) return false;

	std::cout << "this->GetCheatedFigureOfMerit(caloHitList3D) = " << this->GetCheatedFigureOfMerit(caloHitList3D) << std::endl;
    if(this->GetCheatedFigureOfMerit(caloHitList3D)<0.3) return false; //30% impurity at least
	std::cout << "At least 30% impurity" << std::endl;
	
    return true;
}

float ThreeDReclusteringAlgorithm::GetCheatedFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{

    std::map<int,int> mainMcParticleMap;
    for (const CaloHit *const pCaloHit : mergedClusterCaloHitList3D)
    {
        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
        //Find main contributing MC particle Id 
        int mainMcParticleIndex = this->GetMainMcParticleIndex(pParentCaloHit);
        std::map<int, int>::iterator it = mainMcParticleMap.find(mainMcParticleIndex);
 
        if (it != mainMcParticleMap.end()) {
            it->second++;
        }
        else {
            mainMcParticleMap.insert(std::make_pair(mainMcParticleIndex, 1));
        }
    }
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
    MCParticleVector mcParticleVector(pMCParticleList->begin(),pMCParticleList->end());
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);
    
    auto maxSharedHits = std::max_element(mainMcParticleMap.begin(), mainMcParticleMap.end(), [](const auto &x, const auto &y) {
                    return x.second < y.second;
                });
    float mainMcParticleFraction = (float)maxSharedHits->second/mergedClusterCaloHitList3D.size();
    return (1-mainMcParticleFraction);
}

int ThreeDReclusteringAlgorithm::GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit)
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

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(std::string figureOfMeritName, CaloHitList mergedClusterCaloHitList3D)
{
    float figureOfMerit(-999);
    if(figureOfMeritName=="cheated") figureOfMerit=this->GetCheatedFigureOfMerit(mergedClusterCaloHitList3D);
    return figureOfMerit;
}

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(std::string figureOfMeritName, std::vector<CaloHitList*> newClustersCaloHitLists3D)
{
////    if(figureOfMeritName=="cheated")figureOfMerit=this->GetCheatedFigureOfMerit(mergedClusterCaloHitList3D, newClustersCaloHitLists3D);
////    return figureOfMerit;
      std::vector<float> newClustersFigureOfMeritVector;
      for(auto clusterCaloHitLists3D: newClustersCaloHitLists3D)
      {
        if(figureOfMeritName=="cheated")newClustersFigureOfMeritVector.push_back(this->GetCheatedFigureOfMerit(*clusterCaloHitLists3D));
//        newClustersFigureOfMeritVector.push_back(this->GetFigureOfMerit(*iter,mergedClusterCaloHitList3D));
      }
      float figureOfMerit=*(std::min_element(newClustersFigureOfMeritVector.begin(), newClustersFigureOfMeritVector.end()));
      return figureOfMerit;
}

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(std::vector<CaloHitList*> newClustersCaloHitLists3D)
{
    std::vector<float> figureOfMeritVector;
    for (StringVector::const_iterator iter = m_figureOfMeritNames.begin(), iterEnd = m_figureOfMeritNames.end(); iter != iterEnd; ++iter)
    {
        figureOfMeritVector.push_back(this->GetFigureOfMerit(*iter,newClustersCaloHitLists3D));
    }
    
    float figureOfMerit=*(std::min_element(figureOfMeritVector.begin(), figureOfMeritVector.end()));
    return figureOfMerit;
}

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    std::vector<float> figureOfMeritVector;
    for (StringVector::const_iterator iter = m_figureOfMeritNames.begin(), iterEnd = m_figureOfMeritNames.end(); iter != iterEnd; ++iter)
    {
        figureOfMeritVector.push_back(this->GetFigureOfMerit(*iter,mergedClusterCaloHitList3D)); 
    }
    
    float figureOfMerit=*(std::min_element(figureOfMeritVector.begin(), figureOfMeritVector.end()));
    return figureOfMerit;
}


StatusCode ThreeDReclusteringAlgorithm::RebuildPfo(const Pfo *pPfoToRebuild)
{
   m_newClustersUMap.clear();  
   m_newClustersVMap.clear();  
   m_newClustersWMap.clear();  
 
   this->BuildNewTwoDClusters(pPfoToRebuild);
   this->BuildNewPfos(pPfoToRebuild);

   return STATUS_CODE_SUCCESS;
}

StatusCode ThreeDReclusteringAlgorithm::BuildNewTwoDClusters(const Pfo *pPfoToRebuild){

    //std::cout << "Now I want to make 2D clusters" << std::endl;
    ClusterList clusterList2D;
    LArPfoHelper::GetTwoDClusterList(pPfoToRebuild, clusterList2D);
    //std::map<int,const Cluster*> newClustersUMap, newClustersVMap,newClustersWMap;

    for(const Cluster *const pTwoDCluster : clusterList2D)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfoToRebuild, pTwoDCluster));
        HitType hitType = LArClusterHelper::GetClusterHitType(pTwoDCluster);
        std::string clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

        std::string initialListName="";//, debugListName="";
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, initialListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

        // Fragmentation initialisation
        std::string originalListName, fragmentListName;
        ClusterList originalClusterList;
        originalClusterList.push_back(pTwoDCluster);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));
 
        const OrderedCaloHitList &twoDClusterOrderedCaloHitList(pTwoDCluster->GetOrderedCaloHitList());
        OrderedCaloHitList leftoverCaloHitList = twoDClusterOrderedCaloHitList;
 
        int iCluster(0);
        for(const Cluster *pNewCluster : m_newClustersList)
        {
     		if (!pNewCluster)
   			{
         		  std::cerr << "DEBUG, Error: pNewCluster is null!" << std::endl;
       			  continue;
  			}
            PandoraContentApi::Cluster::Parameters parameters;
            CaloHitList newClusterCaloHitList3D;
            pNewCluster->GetOrderedCaloHitList().FillCaloHitList(newClusterCaloHitList3D); 
            //std::cout << "pNewCluster->GetOrderedCaloHitList() size = " << pNewCluster->GetOrderedCaloHitList().size() << std::endl;
 
            for(const CaloHit *const p3DCaloHit : newClusterCaloHitList3D)
            {
                for (const OrderedCaloHitList::value_type &mapEntry : twoDClusterOrderedCaloHitList)
                {
                    for (const CaloHit *const pCaloHit : *mapEntry.second)
                    {
                        if(pCaloHit==static_cast<const CaloHit *>(p3DCaloHit->GetParentAddress())) 
                        {
                          parameters.m_caloHitList.push_back(static_cast<const CaloHit *>(pCaloHit));
                          leftoverCaloHitList.Remove(pCaloHit);
                        }
                    }
                }
            }
            //std::cout << "debug 0" << std::endl;
            const Cluster *pNewTwoDCluster(nullptr);
            if (!parameters.m_caloHitList.empty())
            {
                //std::cout << "Making a new 2D cluster with hitType = " << hitType << std::endl;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewTwoDCluster));
            }
            //PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &parameters.m_caloHitList, "_" + hitType, RED));
            //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            //std::cout << "debug 1" << std::endl;
            if(pNewTwoDCluster!=nullptr && !parameters.m_caloHitList.empty() && hitType==TPC_VIEW_U) {m_newClustersUMap.insert(std::make_pair(iCluster,pNewTwoDCluster));}
            else if(pNewTwoDCluster!=nullptr && !parameters.m_caloHitList.empty() && hitType==TPC_VIEW_V) {m_newClustersVMap.insert(std::make_pair(iCluster,pNewTwoDCluster));}
            else if(pNewTwoDCluster!=nullptr && !parameters.m_caloHitList.empty() && hitType==TPC_VIEW_W) {m_newClustersWMap.insert(std::make_pair(iCluster,pNewTwoDCluster));}
            //std::cout << "debug 2" << std::endl;
            
            iCluster++;
        }
        //PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &clusterListVDebug, "2DClustersBeforeMakingPfo_V", RED));
        //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

        //std::cout << "debug 3" << std::endl;
        
        //Check the leftover caloHits. Attach to the nearest cluster in the new cluster list (newClustersUVect, newClustersVVect or newClustersWVect?
        std::map<int,const Cluster*> clustersForLeftoverHitsMap;
        if(hitType==TPC_VIEW_U) clustersForLeftoverHitsMap = m_newClustersUMap;
        else if(hitType==TPC_VIEW_V) clustersForLeftoverHitsMap = m_newClustersVMap;
        else if(hitType==TPC_VIEW_W) clustersForLeftoverHitsMap = m_newClustersWMap;
        if(clustersForLeftoverHitsMap.size())  //THE QUESTION REMAINS OF WHAT TO DO WITH LEFTOVER HITS IF TEHRE IS NO CLUSTER TO ATTACH THEM TO (THIS CONDITION FAILS)!!!
        {
        for(const OrderedCaloHitList::value_type &mapEntry : leftoverCaloHitList)
            {
                for (const CaloHit *const pCaloHit : *mapEntry.second)
                {
                    const Cluster* pNearestCluster = nullptr;
                    double minimumDistance(std::numeric_limits<float>::max());
                    for(const auto & [clusterIndex, pNewTwoDCluster] : clustersForLeftoverHitsMap)
                    {
                        double dist = LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), pNewTwoDCluster);  
                        if (dist<minimumDistance)
                        {
                            minimumDistance=dist;
                            pNearestCluster=pNewTwoDCluster;
                        }
                    }
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this,pNearestCluster,pCaloHit));
                }
            }
        }


       //std::cout << "pre frag isAvailableU = " << isAvailableU << " isAvailableV = " << isAvailableV << " isAvailableW = " << isAvailableW << std::endl;
        //std::cout << "debug 5 " << std::endl;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
        //std::cout << "debug 6 " << std::endl;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, initialListName));
        //std::cout << "debug 7 " << std::endl;


    //Now I want to clear the lists and vector to avoid memory leaks
        // Delete the CaloHitList objects after processing
    /*for (CaloHitList* caloHitList : newCaloHitListsVector)
    {
        delete caloHitList;  // Clean up the memory
    }
    */
       //newCaloHitListsVector.clear();  // Clear the vector to avoid dangling pointers
    }
   return STATUS_CODE_SUCCESS;

}

StatusCode ThreeDReclusteringAlgorithm::BuildNewPfos(const Pfo *pPfoToRebuild){
    //Making new Pfos
    const PfoList *pNewPfoList(nullptr);
    std::string newPfoListName = "changedShowerParticles3D";
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewPfoList, newPfoListName));
    //std::cout << "debug 8 " << std::endl;
 
    std::string originalClusterListName="InitialCluster";
    int iCluster(0);

    //std::cout << "About to create new pfos. newClustersList.size() = " << newClustersList.size() << std::endl; 
    //ClusterList twoDClusterList;
    for(const Cluster *pNewThreeDCluster : m_newClustersList)
    {
            PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
            const bool isAvailableU((m_newClustersUMap.count(iCluster)) && m_newClustersUMap.at(iCluster)->IsAvailable());
            const bool isAvailableV((m_newClustersVMap.count(iCluster)) && m_newClustersVMap.at(iCluster)->IsAvailable());
            const bool isAvailableW((m_newClustersWMap.count(iCluster)) && m_newClustersWMap.at(iCluster)->IsAvailable());
            //std::cout << "isAvailableU = " << isAvailableU << " isAvailableV = " << isAvailableV << " isAvailableW = " << isAvailableW << std::endl;
            CaloHitList clusterUHits, clusterVHits, clusterWHits;
            if(isAvailableU) m_newClustersUMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterUHits);
            if(isAvailableV) m_newClustersVMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterVHits);
            if(isAvailableW) m_newClustersWMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterWHits);
            if(isAvailableU) pfoParameters.m_clusterList.push_back(m_newClustersUMap.at(iCluster));
            if(isAvailableV) pfoParameters.m_clusterList.push_back(m_newClustersVMap.at(iCluster));
            if(isAvailableW) pfoParameters.m_clusterList.push_back(m_newClustersWMap.at(iCluster));


            pfoParameters.m_clusterList.push_back(pNewThreeDCluster);
 
            pfoParameters.m_particleId = pPfoToRebuild->GetParticleId(); // SHOWER, placeholder for now... Are the new clusters all showers???
            pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
            pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
            pfoParameters.m_energy = 0.f;
            pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
 
            const ParticleFlowObject *pNewPfo(nullptr);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pNewPfo));
 
            //ClusterList newClusterList2D;
            //LArPfoHelper::GetTwoDClusterList(pNewPfo, newClusterList2D);
            //twoDClusterList.insert(twoDClusterList.end(),newClusterList2D.begin(),newClusterList2D.end());

            iCluster++;
            
    }
    //PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &twoDClusterList, "twoDClusterList", RED));
    //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    //std::cout << "debug 9 " << std::endl;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, newPfoListName, m_newPfosListNameAllAfterReclustering));
    //std::cout << "debug 10 " << std::endl;

        /////////////////////////////// END OF REBUILD PFO

   return STATUS_CODE_SUCCESS;
}

StatusCode ThreeDReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "FigureOfMeritNames", m_figureOfMeritNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualDisplaysOn", m_visualDisplaysOn));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "HitThresholdForNewPfo", m_hitThresholdForNewPfo));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "ClusteringTools", algorithmToolVector));

    for (auto algorithmTool : algorithmToolVector)
    {
        ClusteringTool *const pClusteringTool(dynamic_cast<ClusteringTool *>(algorithmTool));

        if (!pClusteringTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pClusteringTool);
     }
 


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
