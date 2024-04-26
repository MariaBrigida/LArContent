/**
 *  @file   larpandoradlcontent/LArReclustering/DlGraphClusteringAlgorithm.h
 *
 *  @brief  Header file for the deep learning vertexing algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_DL_GRAPH_CLUSTERING_ALGORITHM_H
#define LAR_DL_GRAPH_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief   class
 */
class DlGraphClusteringAlgorithm : public pandora::Algorithm
{
public:
    typedef std::map<std::pair<int, int>, std::vector<const pandora::CaloHit *>> PixelToCaloHitsMap;

    /**
     *  @brief Default constructor
     */
    DlGraphClusteringAlgorithm();

    virtual ~DlGraphClusteringAlgorithm();

private:

    //const pandora::CartesianVector &GetPosition() const;
    //std::string ToString() const;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();
    pandora::StatusCode Test();
    int GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit);



    bool m_trainingMode;                      ///< Training mode
    std::string m_trainingOutputFile;         ///< Output file name for training examples
    int m_event;                              ///< The current event number
    bool m_writeTree;                         ///< Whether or not to write validation details to a ROOT tree
    std::string m_rootTreeName;               ///< The ROOT tree name
    std::string m_rootFileName;               ///< The ROOT file name
    //pandora::StringVector m_caloHitListNames; ///< Names of input calo hit lists
    std::mt19937 m_rng;                       ///< The random number generator
    std::string m_mcParticleListName;         ///< The mc particle list name 
    LArDLHelper::TorchModel m_model;



    /*
     *  @brief  Retrieve the map from MC to calo hits for reconstructable particles
     *
     *  @param  mcToHitsMap The map to populate
     *
     *  @return The StatusCode resulting from the function
     **/
    //pandora::StatusCode GetMCToHitsMap(LArMCParticleHelper::MCContributionMap &mcToHitsMap) const;

    /*
     *  @brief  Construct a list of the MC particles from the MC to calo hits map, completing the interaction hierarchy with the invisible
     *          upstream particles.
     *
     *  @param  mcToHitsMap The map of reconstructible MC particles to calo hits
     *  @param  mcHierarchy The output list of MC particles representing the interaction
     *
     *  @return The StatusCode resulting from the function
     **/
    //pandora::StatusCode CompleteMCHierarchy(const LArMCParticleHelper::MCContributionMap &mcToHitsMap, pandora::MCParticleList &mcHierarchy) const;

};

} // namespace lar_dl_content

#endif // LAR_DL_GRAPH_CLUSTERING_ALGORITHM_H
