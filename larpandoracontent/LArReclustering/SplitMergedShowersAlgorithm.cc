/**
 *  @file   larpandoracontent/LArReclustering/SplitMergedShowersAlgorithm.cc
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

#include "larpandoracontent/LArReclustering/SplitMergedShowersAlgorithm.h"

const float convertADCToMeV = 0.0075;
const int atomicNumberArgon = 18;
const float atomicMassArgon = 39.948;
const float criticalEnergyArgon = 32.84;
const float moliereRadiusCmArgon = 9.043; //cm

bool sortByCaloHits (pandora::CaloHitList a, pandora::CaloHitList b) { return (a.size()>b.size()); }

using namespace pandora;

namespace lar_content
{

SplitMergedShowersAlgorithm::SplitMergedShowersAlgorithm():
    m_drawProfiles(false),
    m_hitThresholdForNewPfo(0)
{
}

SplitMergedShowersAlgorithm::~SplitMergedShowersAlgorithm()
{
}

StatusCode SplitMergedShowersAlgorithm::Run()
{
    std::cout << "--------------------------------------------- NEW EVENT ------------------------------------------------" << std::endl;
    if(m_drawProfiles)PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));

    // Get shower pfos and then find the 3D cluster in the shower pfo.
    const PfoList *pShowerPfoList(nullptr);
    std::string initialPfosListName;//, initialVertexListName;
    //std::string daughterVertexListName="DaughterVertices3D";
    std::string newPfosListNameAllAfterReclustering = "newShowerParticles3D";
    //std::string newVertexListNameAllAfterReclustering = "newVertices";
    std::string newPfoListName = "changedShowerParticles3D";
    //std::string newVertexListName = "changedVertices";
    std::string newPfosListNameUnchanged = "unchangedShowerParticles3D";
    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, "ShowerParticles3D", pShowerPfoList)) && pShowerPfoList)
    {
        std::cout << "In this event there are " << pShowerPfoList->size() << " shower pfos." << std::endl;
        if (pShowerPfoList->size()==0) return STATUS_CODE_NOT_FOUND;
        PfoList unchangedPfoList;   
        //VertexList unchangedVertexList;   
        std::vector<float> mainClusterFractionVector, initialFigureOfMeritVector, newFigureOfMeritVector, nHitsInitialFomVector, nFinalClustersVector;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Pfo>(*this, initialPfosListName));
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Vertex>(*this, initialVertexListName));

        int iPfo(0); //for debug
        for (const Pfo *const pShowerPfo : *pShowerPfoList)
        {
            std::cout << "iPfo = " << iPfo << std::endl;

            const PfoList *pNewPfoList(nullptr);
            //const VertexList *pNewVertexList(nullptr);
            ClusterList clusterList3D;
            LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);

            if(!this->PassesCutsForReclustering(pShowerPfo)) continue; // this function just checks it's a shower at the moment

            // Get the longitudinal and transverse shower profiles
            if (clusterList3D.empty())
                continue;

            CaloHitList caloHitList3D;
            clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);
            //float cheatedFOM = (float)caloHitList3D.front()->size()/caloHitList3D.size();

            //if(m_drawProfiles)
            //{
                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(),&caloHitList3D, "ShowerPfoCaloHitsBeforeReclustering", RED);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            //}


            //Quality cuts
            if (caloHitList3D.size() < 2)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            std::cout << "caloHitList3D size = " << caloHitList3D.size() << std::endl;

            //Some pfos are shower-like and yet include track-like 3D clusters. For the moment I don't want to deal with these.
            const ClusterList *pShowerClusters(nullptr);
            PandoraContentApi::GetList(*this, "ShowerClusters3D", pShowerClusters);
            if(pShowerClusters->end() == std::find(pShowerClusters->begin(), pShowerClusters->end(), clusterList3D.front())) continue;

            float initialFigureOfMerit=(1-this->GetCheatedFigureOfMerit(caloHitList3D));
            //float initialFigureOfMerit=this->GetTransverseProfileFigureOfMerit(caloHitList3D);
            std::cout << "InitialFigureOfMerit = " << initialFigureOfMerit << std::endl;

            //Free the hits in this cluster, so that they are available for reclustering!
            //Ask to remove the 3D cluster from the parent pfo, so that it's not owned any more
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pShowerPfo, clusterList3D.front()));
            std::cout << "deb0" << std::endl;
            //Pop this cluster in a local clusterlist
            const ClusterList reclusterClusterList(1, clusterList3D.front());
            const TrackList reclusterTrackList; //dummy track list

            //Also remove vertex from shower pfo
            //const Vertex *const pInitialVertex(LArPfoHelper::GetVertex(pShowerPfo));
            //std::cout << "deb1" << std::endl;
            //PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pShowerPfo, pInitialVertex));
            //std::cout << "deb2" << std::endl;
                        
            // Initialize reclustering with these local lists
            std::string currentClustersListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentClustersListName));
            //std::cout << "deb3" << std::endl;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, "ShowerClusters3D"));
            //std::cout << "deb4" << std::endl;


            // Specify clusters and tracks to be used in reclustering
            std::string originalClustersListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, reclusterTrackList, reclusterClusterList, originalClustersListName));
            //std::cout << "deb5" << std::endl;

            if(m_drawProfiles)
            {
                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&reclusterClusterList, "ClustersToBeReclustered", RED);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }
            //std::cout << "deb6" << std::endl;

            //Call the reclustering algos that produce new cluster candidates
            float minimumFigureOfMerit(initialFigureOfMerit);
            ClusterList minimumFigureOfMeritClusterList = clusterList3D;
            const ClusterList *pReclusterList = NULL;
            std::string reclusterListName;
            for (StringVector::const_iterator clusteringIter = m_clusteringAlgorithms.begin(), clusteringIterEnd = m_clusteringAlgorithms.end();
                clusteringIter != clusteringIterEnd; ++clusteringIter)
            {
                // Produce new cluster candidates
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, *clusteringIter, 
                    pReclusterList, reclusterListName));
   
                
                if (pReclusterList->empty())
                    continue;
                 
                if(m_drawProfiles)
                {
                    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),pReclusterList, "ReclusteredClusters", BLUE);
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                }

                
                //Loop over clusters and calculate new FOM including them if they have more than 10 hits
                std::vector<CaloHitList> newClustersCaloHitLists3D;
                ClusterList clustersFromReclustering;
                for(const Cluster *const pNewCluster : *pReclusterList)
                {
                  CaloHitList newClusterCaloHitList3D;
                  pNewCluster->GetOrderedCaloHitList().FillCaloHitList(newClusterCaloHitList3D);

                  if((int)newClusterCaloHitList3D.size()<m_hitThresholdForNewPfo) continue;

                  newClustersCaloHitLists3D.push_back(newClusterCaloHitList3D);
                  clustersFromReclustering.push_back(pNewCluster);
                 
                  if(m_drawProfiles)
                  {
                      PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &newClusterCaloHitList3D, "newClusterCaloHitList3D", AUTOITER));                 
                      PandoraMonitoringApi::ViewEvent(this->GetPandora());
                  }
                }

                //Draw calo hits from new clusters, for visual debug
                for(CaloHitList newClustersCaloHitList3D: newClustersCaloHitLists3D)
                {
                    //if(m_drawProfiles)
                    //{
                        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(),&newClustersCaloHitList3D, "ShowerPfoCaloHits", AUTOID);
                        PandoraMonitoringApi::ViewEvent(this->GetPandora());
                    //}
                }


                std::sort (newClustersCaloHitLists3D.begin(), newClustersCaloHitLists3D.end(), sortByCaloHits);
                std::cout << "newClustersCaloHitLists3D size = " << newClustersCaloHitLists3D.size() << std::endl;
                if(!newClustersCaloHitLists3D.size()) continue; 
                float mainClusterFraction = (float)newClustersCaloHitLists3D.front().size()/caloHitList3D.size();
                //float newFigureOfMerit = this->GetTransverseProfileFigureOfMerit(caloHitList3D, newClustersCaloHitLists3D); //I am also passing original caloHitList3D as the FOM uses the original PCA calculations
                float newFigureOfMerit = (1-this->GetCheatedFigureOfMerit(newClustersCaloHitLists3D)); //Cheated FOM for group of new clusters is the minimum cheated FOM among those clusters
                std::cout << "newFigureOfMerit = " << newFigureOfMerit << std::endl;
                mainClusterFractionVector.push_back(mainClusterFraction);  ///watch out, these are in the loop over many algorithms! if I add more clustering algos I will need to differentiate the entries in these

                //Will print these for study/debug purposes
                initialFigureOfMeritVector.push_back(initialFigureOfMerit);
                newFigureOfMeritVector.push_back(newFigureOfMerit);
                nHitsInitialFomVector.push_back(caloHitList3D.size());
                nFinalClustersVector.push_back(newClustersCaloHitLists3D.size());

                if(newFigureOfMerit<minimumFigureOfMerit)
                {
                    minimumFigureOfMerit=newFigureOfMerit;
                    //minimumFigureOfMeritClusterList=*pReclusterList; 
                    minimumFigureOfMeritClusterList=clustersFromReclustering; 
                }
            }

           //std::cout << "minimumFigureOfMerit = " << minimumFigureOfMerit << std::endl;
           if(minimumFigureOfMeritClusterList==clusterList3D)
           {
               std::cout << "NO CHANGE! Keep the original pfo" << std::endl;
               PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pShowerPfo, clusterList3D.front()));
               //PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pShowerPfo, pInitialVertex));
               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, originalClustersListName));
               unchangedPfoList.push_back(pShowerPfo);
               //unchangedVertexList.push_back(LArPfoHelper::GetVertex(pShowerPfo));
               
           }
           else
           {
               std::cout << "THE CLUSTER LIST CHANGED! I need to make " << minimumFigureOfMeritClusterList.size() << " new pfos." << std::endl;
               PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, reclusterListName));


               ClusterList clusterList2D;
               LArPfoHelper::GetTwoDClusterList(pShowerPfo, clusterList2D);
               
               std::vector<const Cluster*> newClustersUVect, newClustersVVect,newClustersWVect;
 
               for(const Cluster *const pTwoDCluster : clusterList2D)
               {

                   PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pShowerPfo, pTwoDCluster));

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
                   int iCluster(0);
                   for(const Cluster *const pNewCluster : minimumFigureOfMeritClusterList)
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
                                   {parameters.m_caloHitList.push_back(static_cast<const CaloHit *>(pCaloHit));}
                               }
                           }
                       }
                       const Cluster *pNewTwoDCluster(nullptr);
                       if (!parameters.m_caloHitList.empty())
                       {
                           PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewTwoDCluster));
                           CaloHitList clusterHits;
                           pNewTwoDCluster->GetOrderedCaloHitList().FillCaloHitList(clusterHits);

                       }
                       if(hitType==TPC_VIEW_U) newClustersUVect.push_back(pNewTwoDCluster);
                       else if(hitType==TPC_VIEW_V) newClustersVVect.push_back(pNewTwoDCluster);
                       else if(hitType==TPC_VIEW_W) newClustersWVect.push_back(pNewTwoDCluster);
                       iCluster++;
                   }
                   
                   PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
                   PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, initialListName));

               }
               PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewPfoList, newPfoListName));
               //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pNewVertexList, newVertexListName));

               std::string originalClusterListName="InitialCluster";
               int iCluster(0);
               for(const Cluster *const pNewThreeDCluster : minimumFigureOfMeritClusterList)
               {

                       PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;

                       const bool isAvailableU((NULL != newClustersUVect.at(iCluster)) && newClustersUVect.at(iCluster)->IsAvailable());
                       const bool isAvailableV((NULL != newClustersVVect.at(iCluster)) && newClustersVVect.at(iCluster)->IsAvailable());
                       const bool isAvailableW((NULL != newClustersWVect.at(iCluster)) && newClustersWVect.at(iCluster)->IsAvailable());
                       CaloHitList clusterUHits, clusterVHits, clusterWHits;
                       if(isAvailableU)newClustersUVect.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterUHits);
                       if(isAvailableV)newClustersVVect.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterVHits);
                       if(isAvailableW)newClustersWVect.at(iCluster)->GetOrderedCaloHitList().FillCaloHitList(clusterWHits);

                       if(isAvailableU) pfoParameters.m_clusterList.push_back(newClustersUVect.at(iCluster));
                       if(isAvailableV) pfoParameters.m_clusterList.push_back(newClustersVVect.at(iCluster));
                       if(isAvailableW) pfoParameters.m_clusterList.push_back(newClustersWVect.at(iCluster));
                       pfoParameters.m_clusterList.push_back(pNewThreeDCluster);

                       /*PandoraContentApi::Vertex::Parameters vertexParameters;
                       vertexParameters.m_position = pInitialVertex->GetPosition();
                       vertexParameters.m_vertexLabel = VERTEX_START;
                       vertexParameters.m_vertexType = VERTEX_3D;
                       const Vertex *pNewVertex(nullptr);
                       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Vertex::Create(*this, vertexParameters, pNewVertex));
                       pfoParameters.m_vertexList.push_back(pNewVertex); //For now, using the input shower pfo vertex list (but will need modifying)
                       */


                       pfoParameters.m_particleId = pShowerPfo->GetParticleId(); // SHOWER, placeholder for now... Are the new clusters all showers???
                       pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
                       pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
                       pfoParameters.m_energy = 0.f;
                       pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);

                       const ParticleFlowObject *pNewPfo(nullptr);
                       PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pNewPfo));
                       std::cout << "done! pfo n. = " << iCluster << std::endl;

                       ClusterList newClusterList2D;
                       LArPfoHelper::GetTwoDClusterList(pNewPfo, newClusterList2D);
                       //if(m_drawProfiles)
                       //{
                          PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                          PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&newClusterList2D, "Reclustered2DClusters", RED);
                          PandoraMonitoringApi::ViewEvent(this->GetPandora());
                       //}

                       iCluster++;
                       
               }
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, newPfoListName, newPfosListNameAllAfterReclustering));
                //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, newVertexListName, newVertexListNameAllAfterReclustering));

                const PfoList *pDebugNewPfos(nullptr);
                PandoraContentApi::GetList(*this, newPfosListNameAllAfterReclustering, pDebugNewPfos);

                //if(m_drawProfiles)
                //{
                   PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                   PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(),pDebugNewPfos, "newPfosListNameAllAfterReclustering", RED);
                   PandoraMonitoringApi::ViewEvent(this->GetPandora());
                //}

            }
            iPfo++;

            //if(m_drawProfiles && pNewPfoList->size()>0)
            //{
            //    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
            //    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(),pNewPfoList, "PfosToVisualizeAfter", RED);
            //    PandoraMonitoringApi::ViewEvent(this->GetPandora());
            //}

        }

        if(unchangedPfoList.size()>0) PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<PfoList>(*this, "ShowerParticles3D", newPfosListNameAllAfterReclustering,  unchangedPfoList));


        //DEBUGGING ONLY
        //const VertexList *pInitialVertexList;
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList<VertexList>(*this, daughterVertexListName, pInitialVertexList));
        //std::cout << "initialVertexListName = " << initialVertexListName << "pInitialVertexList size = " << pInitialVertexList->size() << std::endl;

        //if(unchangedVertexList.size()>0) PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<VertexList>(*this, daughterVertexListName, newVertexListNameAllAfterReclustering,  unchangedVertexList));
        const PfoList *pNewPfosListAllAfterReclustering;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList<PfoList>(*this, newPfosListNameAllAfterReclustering, pNewPfosListAllAfterReclustering));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, newPfosListNameAllAfterReclustering, "ShowerParticles3D"));
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Vertex>(*this, newVertexListNameAllAfterReclustering, daughterVertexListName));
        const PfoList *pShowerParticles3DDebug;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList<PfoList>(*this, "ShowerParticles3D", pShowerParticles3DDebug));

        //if(m_drawProfiles)
        //{
           PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
           PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(),pShowerParticles3DDebug, "ShowerParticles3DDebug", RED);
           PandoraMonitoringApi::ViewEvent(this->GetPandora());
         //}

/*        for(long unsigned int iFom=0; iFom<newFigureOfMeritVector.size(); iFom++)
        {
            std::cout << iFom << " " << mainClusterFractionVector.at(iFom) << " " << initialFigureOfMeritVector.at(iFom) << " " << newFigureOfMeritVector.at(iFom) << " " << nHitsInitialFomVector.at(iFom) << " " << nFinalClustersVector.at(iFom) << std::endl;
        }*/

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, initialPfosListName));
        //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Vertex>(*this, initialVertexListName));

    }

    return STATUS_CODE_SUCCESS;
}

//Lateral profile at the shower maximum
float SplitMergedShowersAlgorithm::GetLateralProfileAtShowerMaximum(float clusterEnergyInMeV, float radiusInCm){
    const float tau=1; //shower maximum
    float energy=clusterEnergyInMeV/criticalEnergyArgon;
    float radius=radiusInCm/moliereRadiusCmArgon;
    float z1 = 0.0251 + 0.00319*std::log(energy);
    float z2 = 0.1162 - 0.000381*atomicNumberArgon;
    float k1 = 0.659 - 0.00309*atomicNumberArgon;
    float k2 = 0.645;
    float k3 = -2.59;
    float k4 = 0.3585 + 0.0421*std::log(energy);
    float p1 = 2.632 - 0.00094*atomicNumberArgon;
    float p2 = 0.401 + 0.00187*atomicNumberArgon;
    float p3 = 1.313 - 0.0686*std::log(energy);
    float Rc = z1+ z2*tau;
    float Rt = k1*(std::exp(k3*(tau-k2))+std::exp(k4*(tau-k2)));
    float prob = p1*std::exp((p2-tau)/p3 - std::exp((p2-tau)/p3));
    float profile = 2*radius*(prob*Rc*Rc/std::pow((radius*radius+Rc*Rc),2)+(1-prob)*Rt*Rt/std::pow((radius*radius+Rt*Rt),2));
    return profile;

}


float SplitMergedShowersAlgorithm::GetCheatedFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
//    float mainClusterFraction = (float)mergedClusterCaloHitList3D.front().size()/mergedClusterCaloHitList3D.size();
    std::map<int,int> mainMcParticleMap;
    for (const CaloHit *const pCaloHit : mergedClusterCaloHitList3D)
    {
        const CaloHit *const pParentCaloHit = static_cast<const CaloHit *>(pCaloHit->GetParentAddress());
        //if (TPC_VIEW_W != pParentCaloHit->GetHitType()) 
        //    continue; 
        //Find main contributing MC particle Id 
        int mainMcParticleIndex = this->GetMainMcParticleIndex(pParentCaloHit);
        //std::cout << "This hit has mainMcParticleIndex = " << mainMcParticleIndex << std::endl;
        std::map<int, int>::iterator it = mainMcParticleMap.find(mainMcParticleIndex);
 
        if (it != mainMcParticleMap.end()) {
            it->second++;
        }
        else {
            mainMcParticleMap.insert(std::make_pair(mainMcParticleIndex, 1));
        }
    }
    //std::cout << "debug mainMcParticleMap size = " << mainMcParticleMap.size() << std::endl;
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
    MCParticleVector mcParticleVector(pMCParticleList->begin(),pMCParticleList->end());
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    std::cout << "Printing map for debug:" << std::endl;
    for (auto i = mainMcParticleMap.begin(); i!= mainMcParticleMap.end(); i++) {
      std::cout << i->first << " : "<< i->second << " pdg code: " << mcParticleVector.at(i->first)->GetParticleId() << " n. parents = " << mcParticleVector.at(i->first)->GetParentList().size() << "  n. daughters = "  << mcParticleVector.at(i->first)->GetDaughterList().size() << " hierarchy tier = " << this->GetMCParticleHierarchyTier(mcParticleVector.at(i->first)) << " momentum = " << mcParticleVector.at(i->first)->GetMomentum().GetMagnitude() << std::endl;
    }
    auto maxSharedHits = std::max_element(mainMcParticleMap.begin(), mainMcParticleMap.end(), [](const auto &x, const auto &y) {
                    return x.second < y.second;
                });
    float mainMcParticleFraction = (float)maxSharedHits->second/mergedClusterCaloHitList3D.size();
    //float mainMcParticleFraction = 1;
    //std::cout << "mainMcParticleFraction = " << mainMcParticleFraction << std::endl; 
    return mainMcParticleFraction;
}



int SplitMergedShowersAlgorithm::GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit)
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

int SplitMergedShowersAlgorithm::GetMCParticleHierarchyTier(const pandora::MCParticle *const pMCParticle)
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



float SplitMergedShowersAlgorithm::GetCheatedFigureOfMerit(std::vector<CaloHitList> newClustersCaloHitLists3D)
{
    float minimumFigureOfMerit(999999);

    for(CaloHitList newCaloHitList3D: newClustersCaloHitLists3D)
    {
        float currentFigureOfMerit = this->GetCheatedFigureOfMerit(newCaloHitList3D);
        if(currentFigureOfMerit<minimumFigureOfMerit)minimumFigureOfMerit=currentFigureOfMerit;
    }

    return minimumFigureOfMerit;
}


float SplitMergedShowersAlgorithm::GetTransverseProfileFigureOfMerit(CaloHitList mergedClusterCaloHitList3D, std::vector<CaloHitList> newClustersCaloHitLists3D)
{
    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));
    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());
    const CartesianVector orthoDirection2(axisDirection.GetCrossProduct(orthoDirection1).GetUnitVector());

    //Estimate total cluster energy
    float clusterEnergyInMeV(0);
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D) clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;

    //Observed transverse profile
    //int transverseProfileNBins = 1001;
    int transverseProfileNBins = 50;
    //float transverseProfileLow = -150.15;
    float transverseProfileLow = -50;
    float transverseProfileHigh = 50;
    float transverseProfileBinSize = (transverseProfileHigh-transverseProfileLow)/transverseProfileNBins;
    TwoDHistogram observedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
        const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
        const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
        observedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
    }

    //Scale observed profile to total cluster energy
    observedTransverseProfile.Scale(clusterEnergyInMeV/observedTransverseProfile.GetCumulativeSum());



    TwoDHistogram expectedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
    if(newClustersCaloHitLists3D.size()==1 || newClustersCaloHitLists3D.at(0)==mergedClusterCaloHitList3D)
    {
        //Expected tranvserse profile (Grindhammer parametrisation)
        for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
        {
            float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
            for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
            {
               float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
               float profileRadius=std::sqrt(profileX*profileX+profileY*profileY);
               float profileValue=GetLateralProfileAtShowerMaximum(clusterEnergyInMeV,profileRadius);
               expectedTransverseProfile.SetBinContent(iBinX, iBinY, profileValue);
            }
        }
        expectedTransverseProfile.Scale(clusterEnergyInMeV/expectedTransverseProfile.GetCumulativeSum());
    }

    //Now, if I passed a new set of cluster hit lists in the second argument, I should separately calculate all of their transverse profiles projecting them onto the same plane as merged shower.
    else if(newClustersCaloHitLists3D.size()!=1 and newClustersCaloHitLists3D.at(0)!=mergedClusterCaloHitList3D)
    {
        std::vector<TwoDHistogram> newObservedTransverseProfiles;
        std::vector<double> newClusterEnergies;
        std::vector<float> newClustersCenterPositionsX, newClustersCenterPositionsY;
        for(CaloHitList newCaloHitList3D: newClustersCaloHitLists3D)
        {
            TwoDHistogram newObservedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);;
             
            double clusterEnergy(0);
            for (const CaloHit *const pCaloHit3D : newCaloHitList3D)
            {
                const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
                const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
                const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
                newObservedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
                clusterEnergy+=pCaloHit3D->GetInputEnergy();
            }
            if(m_drawProfiles)
            { 
                PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), newObservedTransverseProfile, "COLZ");
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }
            newClusterEnergies.push_back(clusterEnergy);
            newClustersCenterPositionsX.push_back(newObservedTransverseProfile.GetMeanX());
            newClustersCenterPositionsY.push_back(newObservedTransverseProfile.GetMeanY());
        }

        //Expected tranvserse profile (Grindhammer parametrisation as a combination of N shower profiles)
        for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
        {
            float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
            for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
            {
               float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
               //float profileRadius=std::sqrt(profileX*profileX+profileY*profileY);
               float profileValue(0), shiftedRadius(0), shiftedX(0), shiftedY(0);
               for(std::vector<double>::size_type iCluster=0; iCluster<newClusterEnergies.size(); iCluster++){
                   shiftedX=profileX-newClustersCenterPositionsX.at(iCluster);
                   shiftedY=profileY-newClustersCenterPositionsY.at(iCluster);
                   shiftedRadius=std::sqrt(shiftedX*shiftedX+shiftedY*shiftedY);
                   profileValue+=GetLateralProfileAtShowerMaximum(newClusterEnergies.at(iCluster)*convertADCToMeV,shiftedRadius);
               }
               expectedTransverseProfile.SetBinContent(iBinX, iBinY, profileValue);
            }
        }
        expectedTransverseProfile.Scale(clusterEnergyInMeV/expectedTransverseProfile.GetCumulativeSum());
    }


    //Calculate figure of merit for this cluster
    float squaredDiffSum(0);
    for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
    {
        for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
        {
          float diff = expectedTransverseProfile.GetBinContent(iBinX, iBinY)-observedTransverseProfile.GetBinContent(iBinX,iBinY);
          float squaredDiff = diff*diff;     
          squaredDiffSum+=squaredDiff;
        }
    }
    float figureOfMerit = squaredDiffSum/clusterEnergyInMeV;

    if(m_drawProfiles) {
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedTransverseProfile, "COLZ");
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedTransverseProfile, "COLZ");
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
    
    return figureOfMerit;
}

float SplitMergedShowersAlgorithm::GetLongitudinalProfileFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

   // Observed longitudinal profile
    Histogram observedLongitudinalProfile(140, 0., 140.);

    const float convertADCToMeV(0.0075f); // (c) Maria
    const float convertGeVToMeV(1000.f);
    const float convertCmToX0(1.f / 14.f);
    float clusterEnergyInMeV(0.f);

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CaloHit *const pParentCaloHit(static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress()));

        if (TPC_VIEW_W != pParentCaloHit->GetHitType())
            continue;

        clusterEnergyInMeV += convertADCToMeV * pParentCaloHit->GetInputEnergy(); // Used later on

        const float longitudinalCoordInCm((pCaloHit3D->GetPositionVector() - axisIntercept).GetDotProduct(axisDirection));
        observedLongitudinalProfile.Fill(longitudinalCoordInCm * convertCmToX0, convertADCToMeV * pParentCaloHit->GetInputEnergy());
    }

    std::cout << "Observed longitudinal energy profile " << std::endl;
    if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedLongitudinalProfile, "");

    // Expected longitudinal profile
    Histogram expectedLongitudinalProfile(140, 0., 140.);

    const float clusterEnergyInGeV(clusterEnergyInMeV / convertGeVToMeV);
    const float longProfileCriticalEnergy(0.08f);
    const float longProfileParameter0(1.25f);
    const float longProfileParameter1(0.5f);

    const double a(longProfileParameter0 + longProfileParameter1 * std::log(clusterEnergyInGeV / longProfileCriticalEnergy));
    const double gammaA(std::exp(lgamma(a)));

    float t(0.f);
    for (int iBin = 0; iBin < expectedLongitudinalProfile.GetNBinsX(); ++iBin)
    {
        t += expectedLongitudinalProfile.GetXBinWidth();
        expectedLongitudinalProfile.Fill(t, convertGeVToMeV * clusterEnergyInGeV / 2. * std::pow(t / 2.f, static_cast<float>(a - 1.)) *
            std::exp(-t / 2.) * expectedLongitudinalProfile.GetXBinWidth() / gammaA);
    }

    std::cout << "Expected longitudinal energy profile " << std::endl;
    if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedLongitudinalProfile, "");

    return 1; //placeholder for fom

}


float SplitMergedShowersAlgorithm::GetTransverseProfileFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));
    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());
    const CartesianVector orthoDirection2(axisDirection.GetCrossProduct(orthoDirection1).GetUnitVector());

    //Estimate total cluster energy
    float clusterEnergyInMeV(0);
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D) clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;

    //Observed transverse profile
    //int transverseProfileNBins = 1001;
    int transverseProfileNBins = 50;
    //float transverseProfileLow = -150.15;
    float transverseProfileLow = -50;
    float transverseProfileHigh = 50;
    float transverseProfileBinSize = (transverseProfileHigh-transverseProfileLow)/transverseProfileNBins;
    TwoDHistogram observedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
        const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
        const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
        observedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
    }

    //Scale observed profile to total cluster energy
    observedTransverseProfile.Scale(clusterEnergyInMeV/observedTransverseProfile.GetCumulativeSum());



    TwoDHistogram expectedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
    //Expected tranvserse profile (Grindhammer parametrisation)
    for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
    {
        float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
        for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
        {
           float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
           float profileRadius=std::sqrt(profileX*profileX+profileY*profileY);
           float profileValue=GetLateralProfileAtShowerMaximum(clusterEnergyInMeV,profileRadius);
           expectedTransverseProfile.SetBinContent(iBinX, iBinY, profileValue);
        }
    }
    expectedTransverseProfile.Scale(clusterEnergyInMeV/expectedTransverseProfile.GetCumulativeSum());

    //Calculate figure of merit for this cluster
    float squaredDiffSum(0);
    for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
    {
        for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
        {
          float diff = expectedTransverseProfile.GetBinContent(iBinX, iBinY)-observedTransverseProfile.GetBinContent(iBinX,iBinY);
          float squaredDiff = diff*diff;     
          squaredDiffSum+=squaredDiff;
        }
    }
    float figureOfMerit = squaredDiffSum/clusterEnergyInMeV;

    if(m_drawProfiles) {
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedTransverseProfile, "COLZ");
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedTransverseProfile, "COLZ");
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
    
    return figureOfMerit;
}






bool SplitMergedShowersAlgorithm::PassesCutsForReclustering(const pandora::ParticleFlowObject *const pShowerPfo)
{
    if (LArPfoHelper::IsShower(pShowerPfo)) return true;
    return false;
}

StatusCode SplitMergedShowersAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DrawProfiles", m_drawProfiles));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "clusteringAlgorithms", m_clusteringAlgorithms));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "HitThresholdForNewPfo", m_hitThresholdForNewPfo));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
