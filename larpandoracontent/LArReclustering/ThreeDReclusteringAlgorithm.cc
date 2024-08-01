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

const float convertADCToMeV = 0.0075;
const int atomicNumberArgon = 18;
const float atomicMassArgon = 39.948;
const float criticalEnergyArgon = 32.84;
const float moliereRadiusCmArgon = 9.043; //cm

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
    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, "ShowerParticles3D", pShowerPfoList)) && pShowerPfoList)
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
        for (const Pfo *const pShowerPfo : *pShowerPfoList)
        {
            std::cout << "iPfo = " << iPfo << std::endl;

            ClusterList clusterList3D;
            LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);
            
            if(!this->PassesCutsForReclustering(pShowerPfo))
            {
                unchangedPfoList.push_back(pShowerPfo);
				std::cout << "This pfo does not pass cuts for reclustering" << std::endl;
                continue; // this function just checks it's a shower at the moment
            }
            std::cout << "DEBUG passed cuts for reclustering" << std::endl;
            if (clusterList3D.empty())
                continue;

            //Some pfos are shower-like and yet include track-like 3D clusters. For the moment I don't want to deal with these.
            if(pShowerClusters->end() == std::find(pShowerClusters->begin(), pShowerClusters->end(), clusterList3D.front())) continue;

            CaloHitList caloHitList3D;
            clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

            //Quality cuts
            if (caloHitList3D.size() < 2)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            //Free the hits in this cluster, so that they are available for reclustering!
            //Ask to remove the 3D cluster from the parent pfo, so that it's not owned any more
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pShowerPfo, clusterList3D.front()));
            //Pop this cluster in a local clusterlist
            const ClusterList reclusterClusterList(1, clusterList3D.front());
            const TrackList reclusterTrackList; //dummy track list

            // Initialize reclustering with these local lists
            std::string currentClustersListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentClustersListName));
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, "ShowerClusters3D"));

            // Specify clusters and tracks to be used in reclustering
            std::string originalClustersListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, reclusterTrackList, reclusterClusterList, originalClustersListName));

            //Call clustering loop
           std::string reclusterListName;
           ClusterList minimumFigureOfMeritClusterList;
           PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->FindBestThreeDClusters(clusterList3D, reclusterListName, minimumFigureOfMeritClusterList));
          //minimumFigureOfMeritClusterList=clusterList3D;

           if(minimumFigureOfMeritClusterList==clusterList3D)
           {
               std::cout << "NO CHANGE! Keep the original pfo" << std::endl;
               PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pShowerPfo, clusterList3D.front()));
               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, originalClustersListName));
               unchangedPfoList.push_back(pShowerPfo);
               
           }
           else
           {
               std::cout << "THE CLUSTER LIST CHANGED! I need to make " << minimumFigureOfMeritClusterList.size() << " new pfos." << std::endl;
               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, reclusterListName));
               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->RebuildPfo(pShowerPfo, minimumFigureOfMeritClusterList));
           }
            iPfo++;
        }
        if(unchangedPfoList.size()>0) PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<PfoList>(*this, "ShowerParticles3D", m_newPfosListNameAllAfterReclustering,  unchangedPfoList));
        const PfoList *pNewPfosListAllAfterReclustering;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList<PfoList>(*this, m_newPfosListNameAllAfterReclustering, pNewPfosListAllAfterReclustering));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_newPfosListNameAllAfterReclustering, "ShowerParticles3D"));
        const PfoList *pShowerParticles3DDebug;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList<PfoList>(*this, "ShowerParticles3D", pShowerParticles3DDebug));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, initialPfosListName));

    }

    return STATUS_CODE_SUCCESS;
}

//This can then go in the SDK with the name/format PandoraContentApi::RebuildLArTPCPfo(Algorithm &, Pfo  *pPfoToRebuild, Cluster *pTemplate3DCluster, MapOfListNameInfo &)
StatusCode ThreeDReclusteringAlgorithm::RebuildPfo(const Pfo *pPfoToRebuild, ClusterList &newThreeDClusters)
{
   ClusterList clusterList2D;
   LArPfoHelper::GetTwoDClusterList(pPfoToRebuild, clusterList2D);
   
   std::map<int,const Cluster*> newClustersUMap, newClustersVMap,newClustersWMap;

   for(const Cluster *const pTwoDCluster : clusterList2D)
   {
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pPfoToRebuild, pTwoDCluster));

       HitType hitType = LArClusterHelper::GetClusterHitType(pTwoDCluster);
       std::string clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

       std::string initialListName="", debugListName="";
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, initialListName));
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, debugListName));

       // Fragmentation initialisation
       std::string originalListName, fragmentListName;
       ClusterList originalClusterList;
       originalClusterList.push_back(pTwoDCluster);


       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

       const OrderedCaloHitList &twoDClusterOrderedCaloHitList(pTwoDCluster->GetOrderedCaloHitList());
       OrderedCaloHitList leftoverCaloHitList = twoDClusterOrderedCaloHitList;

       int iCluster(0);
       for(const Cluster *const pNewCluster : newThreeDClusters)
       {
           PandoraContentApi::Cluster::Parameters parameters;
           CaloHitList newClusterCaloHitList3D;
           pNewCluster->GetOrderedCaloHitList().FillCaloHitList(newClusterCaloHitList3D);                          
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
           const Cluster *pNewTwoDCluster(nullptr);
           if (!parameters.m_caloHitList.empty())
           {
               PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewTwoDCluster));
           }
           if(pNewTwoDCluster!=nullptr && !parameters.m_caloHitList.empty() && hitType==TPC_VIEW_U) {newClustersUMap.insert(std::make_pair(iCluster,pNewTwoDCluster));}
           else if(pNewTwoDCluster!=nullptr && !parameters.m_caloHitList.empty() && hitType==TPC_VIEW_V) {newClustersVMap.insert(std::make_pair(iCluster,pNewTwoDCluster));}
           else if(pNewTwoDCluster!=nullptr && !parameters.m_caloHitList.empty() && hitType==TPC_VIEW_W) {newClustersWMap.insert(std::make_pair(iCluster,pNewTwoDCluster));}

           iCluster++;
       }

       //Check the leftover caloHits. Attach to the nearest cluster in the new cluster list (newClustersUVect, newClustersVVect or newClustersWVect?
       std::map<int,const Cluster*> clustersForLeftoverHitsMap;
       if(hitType==TPC_VIEW_U) clustersForLeftoverHitsMap = newClustersUMap;
       else if(hitType==TPC_VIEW_V) clustersForLeftoverHitsMap = newClustersVMap;
       else if(hitType==TPC_VIEW_W) clustersForLeftoverHitsMap = newClustersWMap;
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
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, initialListName));

   }

   const PfoList *pNewPfoList(nullptr);
   std::string newPfoListName = "changedShowerParticles3D";
   PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewPfoList, newPfoListName));

   std::string originalClusterListName="InitialCluster";
   int iCluster(0);

   for(const Cluster *const pNewThreeDCluster : newThreeDClusters)
   {
           PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
           const bool isAvailableU((newClustersUMap.count(iCluster)) && newClustersUMap.at(iCluster)->IsAvailable());
           const bool isAvailableV((newClustersVMap.count(iCluster)) && newClustersVMap.at(iCluster)->IsAvailable());
           const bool isAvailableW((newClustersWMap.count(iCluster)) && newClustersWMap.at(iCluster)->IsAvailable());
           CaloHitList clusterUHits, clusterVHits, clusterWHits;
           if(isAvailableU)newClustersUMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterUHits);
           if(isAvailableV)newClustersVMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterVHits);
           if(isAvailableW)newClustersWMap.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterWHits);
           if(isAvailableU) pfoParameters.m_clusterList.push_back(newClustersUMap.at(iCluster));
           if(isAvailableV) pfoParameters.m_clusterList.push_back(newClustersVMap.at(iCluster));
           if(isAvailableW) pfoParameters.m_clusterList.push_back(newClustersWMap.at(iCluster));
           pfoParameters.m_clusterList.push_back(pNewThreeDCluster);

           pfoParameters.m_particleId = pPfoToRebuild->GetParticleId(); // SHOWER, placeholder for now... Are the new clusters all showers???
           pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
           pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
           pfoParameters.m_energy = 0.f;
           pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);

           const ParticleFlowObject *pNewPfo(nullptr);
           PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pNewPfo));

           ClusterList newClusterList2D;
           LArPfoHelper::GetTwoDClusterList(pNewPfo, newClusterList2D);

           iCluster++;
           
   }
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, newPfoListName, m_newPfosListNameAllAfterReclustering));
    return STATUS_CODE_SUCCESS;
}

StatusCode ThreeDReclusteringAlgorithm::FindBestThreeDClusters(ClusterList clusterList3D, std::string &reclusterListName, ClusterList &minimumFigureOfMeritClusterList)
{
    CaloHitList caloHitList3D;
    clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);
    float initialFigureOfMerit=this->GetFigureOfMerit(caloHitList3D);
    if(initialFigureOfMerit<0) {std::cout << "Could not calculate initial FOM!" << std::endl; return STATUS_CODE_FAILURE;} 

    //Call the reclustering algos that produce new cluster candidates
    float minimumFigureOfMerit(initialFigureOfMerit);
    minimumFigureOfMeritClusterList = clusterList3D;
    std::vector<float> mainClusterFractionVector,  initialFigureOfMeritVector, newFigureOfMeritVector, nHitsInitialFomVector, nFinalClustersVector;
    const ClusterList *pReclusterList = NULL;
    for (StringVector::const_iterator clusteringIter = m_clusteringAlgorithms.begin(), clusteringIterEnd = m_clusteringAlgorithms.end();
        clusteringIter != clusteringIterEnd; ++clusteringIter)
    {
        // Produce new cluster candidates
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, *clusteringIter, pReclusterList, reclusterListName));

        if (pReclusterList->empty())
            continue;
        
        //Loop over clusters and calculate new FOM including them if they have more than 10 hits
        std::vector<CaloHitList> newClustersCaloHitLists3D;
        ClusterList clustersFromReclustering;
        for(const Cluster *const pNewCluster : *pReclusterList)
        {
          CaloHitList newClusterCaloHitList3D;
          pNewCluster->GetOrderedCaloHitList().FillCaloHitList(newClusterCaloHitList3D);

          if((int)newClusterCaloHitList3D.size()<m_hitThresholdForNewPfo) continue; //This should remove empty clusters as well (isolated hits problem)

          newClustersCaloHitLists3D.push_back(newClusterCaloHitList3D);
          clustersFromReclustering.push_back(pNewCluster);
        }


        std::sort (newClustersCaloHitLists3D.begin(), newClustersCaloHitLists3D.end(), sortByCaloHits);
        if(!newClustersCaloHitLists3D.size()) continue; 
        float mainClusterFraction = (float)newClustersCaloHitLists3D.front().size()/caloHitList3D.size();
        //float newFigureOfMerit = this->GetFigureOfMerit(caloHitList3D, newClustersCaloHitLists3D); //Cheated FOM for group of new clusters is the minimum cheated FOM among those clusters
        float newFigureOfMerit = this->GetFigureOfMerit(newClustersCaloHitLists3D); //Cheated FOM for group of new clusters is the minimum cheated FOM among those clusters
        mainClusterFractionVector.push_back(mainClusterFraction);  ///watch out, these are in the loop over many algorithms! if I add more clustering algos I will need to differentiate the entries in these

        //Will print these for study/debug purposes
        initialFigureOfMeritVector.push_back(initialFigureOfMerit);
        newFigureOfMeritVector.push_back(newFigureOfMerit);
        nHitsInitialFomVector.push_back(caloHitList3D.size());
        nFinalClustersVector.push_back(newClustersCaloHitLists3D.size());

        if(newFigureOfMerit<minimumFigureOfMerit)
        {
            minimumFigureOfMerit=newFigureOfMerit;
            minimumFigureOfMeritClusterList=clustersFromReclustering; 
        }
    }
    return STATUS_CODE_SUCCESS;
}

float ThreeDReclusteringAlgorithm::GetFigureOfMerit(std::string figureOfMeritName, CaloHitList mergedClusterCaloHitList3D)
{
    float figureOfMerit(-999);
    if(figureOfMeritName=="cheated") figureOfMerit=this->GetCheatedFigureOfMerit(mergedClusterCaloHitList3D);
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


//float ThreeDReclusteringAlgorithm::GetFigureOfMerit(std::string figureOfMeritName, CaloHitList mergedClusterCaloHitList3D, std::vector<CaloHitList> newClustersCaloHitLists3D)
float ThreeDReclusteringAlgorithm::GetFigureOfMerit(std::string figureOfMeritName, std::vector<CaloHitList> newClustersCaloHitLists3D)
{
////    if(figureOfMeritName=="cheated")figureOfMerit=this->GetCheatedFigureOfMerit(mergedClusterCaloHitList3D, newClustersCaloHitLists3D);
////    return figureOfMerit;
      std::vector<float> newClustersFigureOfMeritVector;
      for(auto clusterCaloHitLists3D: newClustersCaloHitLists3D)
      {
        if(figureOfMeritName=="cheated")newClustersFigureOfMeritVector.push_back(this->GetCheatedFigureOfMerit(clusterCaloHitLists3D));
//        newClustersFigureOfMeritVector.push_back(this->GetFigureOfMerit(*iter,mergedClusterCaloHitList3D));
      }
      float figureOfMerit=*(std::min_element(newClustersFigureOfMeritVector.begin(), newClustersFigureOfMeritVector.end()));
      return figureOfMerit;
}

//float ThreeDReclusteringAlgorithm::GetFigureOfMerit(CaloHitList mergedClusterCaloHitList3D, std::vector<CaloHitList> newClustersCaloHitLists3D)
float ThreeDReclusteringAlgorithm::GetFigureOfMerit(std::vector<CaloHitList> newClustersCaloHitLists3D)
{
    std::vector<float> figureOfMeritVector;
    for (StringVector::const_iterator iter = m_figureOfMeritNames.begin(), iterEnd = m_figureOfMeritNames.end(); iter != iterEnd; ++iter)
    {
        //figureOfMeritVector.push_back(this->GetFigureOfMerit(*iter,mergedClusterCaloHitList3D,newClustersCaloHitLists3D)); 
        figureOfMeritVector.push_back(this->GetFigureOfMerit(*iter,newClustersCaloHitLists3D)); 
    }
    
    float figureOfMerit=*(std::min_element(figureOfMeritVector.begin(), figureOfMeritVector.end()));
    return figureOfMerit;
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

int ThreeDReclusteringAlgorithm::GetMCParticleHierarchyTier(const pandora::MCParticle *const pMCParticle)
{
    MCParticleVector mcParticleHierarchyVector;
    int mcParticleHierarchyTier(0);
    if(pMCParticle->GetParentList().size())
    {
        mcParticleHierarchyTier++;
        const pandora::MCParticle *pParentMCParticle = nullptr;
        pParentMCParticle=LArMCParticleHelper::GetParentMCParticle(pMCParticle);
        while(pParentMCParticle->GetParentList().size()){
            pParentMCParticle=LArMCParticleHelper::GetParentMCParticle(pParentMCParticle);
            mcParticleHierarchyTier++;
        }
    }
    return mcParticleHierarchyTier;
}



//float ThreeDReclusteringAlgorithm::GetCheatedFigureOfMerit(std::vector<CaloHitList> newClustersCaloHitLists3D)
//{
//    float minimumFigureOfMerit(999999);
//
//    for(CaloHitList newCaloHitList3D: newClustersCaloHitLists3D)
//    {
//        float currentFigureOfMerit = this->GetCheatedFigureOfMerit(newCaloHitList3D);
//        if(currentFigureOfMerit<minimumFigureOfMerit)minimumFigureOfMerit=currentFigureOfMerit;
//    }
//
//    return minimumFigureOfMerit;
//}

//At the moment I only want to select showers with reasonably high impurity; this is to produce a training sample and to test on the most useful scenario
//(In reality training stage I am even more restrictive: I ask that the 2nd most contributing particle contributes at least 30%. So the total impurity across all further contributions can be higher than this 30%. This is ok.
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

StatusCode ThreeDReclusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "FigureOfMeritNames", m_figureOfMeritNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VisualDisplaysOn", m_visualDisplaysOn));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "ClusteringAlgorithms", m_clusteringAlgorithms));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "HitThresholdForNewPfo", m_hitThresholdForNewPfo));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
