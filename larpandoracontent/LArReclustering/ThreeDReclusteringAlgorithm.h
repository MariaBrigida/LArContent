/**
 *  @file   larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h
 *
 *  @brief  Header file for the 3D reclustering algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THREE_D_RECLUSTERING_ALGORITHM_H
#define LAR_THREE_D_RECLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

 class ClusteringTool; 
 
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief ThreeDReclusteringAlgorithm class
 */
class ThreeDReclusteringAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ThreeDReclusteringAlgorithm();

    /**
     *  @brief  Destructor
     */
    ~ThreeDReclusteringAlgorithm();

private:
    /**
     *  @brief  Run the algorithm
     */
    pandora::StatusCode Run();


    /**
     *  @brief Find the calo hit configuration that yields the best 3D clusters according to the chosen metric 
     *
     *  @param clusterList3D 
     *  @param minimumFigureOfMeritClusterListName
     *
     *  @return List of 
     */
    //pandora::StatusCode FindBestThreeDClusters(pandora::ClusterList clusterList3D, std::string minimumFigureOfMeritClusterListName);

    float GetFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);
    
    float GetFigureOfMerit(std::string figureOfMeritName, pandora::CaloHitList mergedClusterCaloHitList3D);

    //float GetFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D, std::vector<pandora::CaloHitList> newClustersCaloHitList3D);
    float GetFigureOfMerit(std::vector<pandora::CaloHitList*> newClustersCaloHitList3D);

    //float GetFigureOfMerit(std::string figureOfMeritName, pandora::CaloHitList mergedClusterCaloHitList3D, std::vector<pandora::CaloHitList> newClustersCaloHitLists3D);
    float GetFigureOfMerit(std::string figureOfMeritName, std::vector<pandora::CaloHitList*> newClustersCaloHitLists3D);

    /** 
     *  @brief Use this to find cheated FOM by finding the shower purity as the fraction of hits that the main Mc particle is contributing to
     * @param Merged cluster hit list
     * @return The figure of merit
     */
    float GetCheatedFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);
    
    //float GetCheatedFigureOfMerit(std::vector<pandora::CaloHitList> newClustersCaloHitList3D);



    /** 
     *  @brief Use this to find FOM by comparing observed transverse profiles to prediction for one shower
     * @param Merged cluster hit list
     * @return The figure of merit
     */
    float GetTransverseProfileFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);
    float GetLongitudinalProfileFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D);

    /** 
     *  @brief Use this to find FOM by comparing observed transverse profile (merged) to prediction produced as combination of a list of underlying clusters
     * @param Merged cluster hit list
     * @param Reclustered clusters hit lists vector

     * @return The figure of merit
     */
    float GetTransverseProfileFigureOfMerit(pandora::CaloHitList mergedClusterCaloHitList3D, std::vector<pandora::CaloHitList> newClustersCaloHitList3D); //Eventually, perhaps move each of these in its own class, and let choose which ones to calculate via XML, like the clustering algos themselves


    //Decide whether to try reclustering for this pfo
    bool PassesCutsForReclustering(const pandora::ParticleFlowObject *const pPfo);

    int GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit);
    int GetMCParticleHierarchyTier(const pandora::MCParticle *const pMCParticle);

    //This can then go in the SDK with the name/format PandoraContentApi::RebuildLArTPCPfo(Algorithm &, Pfo  *pPfoToRebuild, Cluster *pTemplate3DCluster, MapOfListNameInfo &)
     //pandora::StatusCode RebuildPfo(const pandora::Pfo *pPfoToRebuild, pandora::ClusterList *newThreeDClusters); 
     pandora::StatusCode RebuildPfo(const pandora::Pfo *pPfoToRebuild); 
     //pandora::StatusCode BuildNewTwoDClusters(const pandora::Pfo *pPfoToRebuild, pandora::ClusterList *newThreeDClusters); 
     pandora::StatusCode BuildNewTwoDClusters(const pandora::Pfo *pPfoToRebuild); 
     //pandora::StatusCode BuildNewPfos(const pandora::Pfo *pPfoToRebuild, pandora::ClusterList *newThreeDClusters); 
     pandora::StatusCode BuildNewPfos(const pandora::Pfo *pPfoToRebuild); 

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    typedef std::vector<ClusteringTool *> ClusteringToolVector;
    ClusteringToolVector m_algorithmToolVector; ///< The algorithm tool vector

    std::string m_pfoListName; ///< The input pfo list name (e.g. list of neutrino or testbeam pfos)
    pandora::StringVector m_figureOfMeritNames; ///what figure(s) of merit to use
    std::string m_newPfosListNameAllAfterReclustering;
    bool m_visualDisplaysOn; //Boolean to enable and disable displaying transverse profiles
    int m_hitThresholdForNewPfo; //Minimum nr. of hits to form new 3Dcluster and pfo
    pandora::StringVector   m_clusteringAlgorithms; ///< The ordered list of clustering algorithms to be used
    std::string m_mcParticleListName; ///< The mc particle list name 
    
    std::map<int, const pandora::Cluster*> m_newClustersUMap;
    std::map<int, const pandora::Cluster*> m_newClustersVMap;
    std::map<int, const pandora::Cluster*> m_newClustersWMap;
    pandora::ClusterList m_newClustersList;
};


//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  ClusteringTool class
 */



class ClusteringTool: public pandora::AlgorithmTool {

public:
    ClusteringTool() = default;
    virtual ~ClusteringTool() = default;
    virtual bool Run(const pandora::Algorithm *const pAlgorithm, std::vector<pandora::CaloHitList*> &newCaloHitListsVector) = 0;
};

} // namespace lar_content

#endif // #ifndef LAR_THREE_D_RECLUSTERING_ALGORITHM_H
