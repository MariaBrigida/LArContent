/**
 *  @file   larpandoracontent/LArThreeDReco/ThirdViewRecoveryAlgorithm.h
 *
 *  @brief  Header file for the clustering parent algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_THIRD_VIEW_RECOVERY_ALGORITHM_H
#define LAR_THIRD_VIEW_RECOVERY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
//#include "Objects/Vertex.h"

//#include "larpandoracontent/LArObjects/LArPfoObjects.h"


namespace lar_content
{

/**
 *  @brief  ThirdViewRecoveryAlgorithm class
 */
class ThirdViewRecoveryAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ThirdViewRecoveryAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::Cluster *, const pandora::ParticleFlowObject *> ClusterToPfoMap;
    /**
     *  @brief  Get all 3d clusters contained in the input pfo lists and a mapping from clusters to pfos
     *
     *  @param  clusters3D to receive the sorted list of 3d clusters
     *  @param  clusterToPfoMap to receive the mapping from 3d cluster to pfo
     */
    //void GetThreeDClusters(pandora::ClusterVector &clusters3D, ClusterToPfoMap &clusterToPfoMap) const;
    //void GetTwoViewPfos(pandora::PfoList &twoViewPfoList) const;
    void GetTwoViewPfos() const;

    /**
     *  @brief  
     *
     *  @param  
     *  @param  
     */
    //void GetThreeDClusters(pandora::ClusterVector &clusters3D, ClusterToPfoMap &clusterToPfoMap) const;
    //void AddThirdViewClusters(pandora::PfoList &twoViewPfoList, ClusterToPfoMap &clusterToPfoMap) const;


    pandora::StringVector m_inputPfoListNames; ///< The input pfo list names
    pandora::StringVector m_inputClusterListNames; ///< The input cluster list names
    int m_slidingLinearFitWindow;
    unsigned int m_minimumHitsToProject;
//    std::string m_clusteringAlgorithmName;  ///< The name of the clustering algorithm to run
//    std::string m_associationAlgorithmName; ///< The name of the topological association algorithm to run

//    std::string m_inputCaloHitListName; ///< The name of the input calo hit list, containing the hits to be clustered
//    bool m_replaceCurrentCaloHitList;   ///< Whether to permanently replace the original calo hit list as the "current" list upon completion

//    std::string m_clusterListName;    ///< The name under which to save the new cluster list
//    bool m_replaceCurrentClusterList; ///< Whether to subsequently use the new cluster list as the "current" list
};

} // namespace lar_content

#endif // #ifndef LAR_THIRD_VIEW_RECOVERY_ALGORITHM_H
