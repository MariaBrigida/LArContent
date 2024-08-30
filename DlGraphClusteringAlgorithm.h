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
 *  @brief  DlGraphClusteringAlgorithm class
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

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    pandora::StatusCode Infer();

    /**
     *  @brief Get the MC index for the true particle that contributes most energy to the calo hit 
     *
     *  @param pCaloHit address of the calo hit
     *
     *  @return MCParticle Index
     */
    int GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit);

    /**
     *  @brief Get the pdg code for the true particle that contributes most energy to the calo hit 
     *
     *  @param pCaloHit address of the calo hit
     *
     *  @return pdg code
     */
    int GetMainMcParticlePdgCode(const pandora::CaloHit *const pCaloHit);



    bool m_trainingMode;                      ///< Training mode
    std::string m_trainingOutputFile;         ///< Output file name for training examples
    int m_event;                              ///< The current event number
    bool m_writeTree;                         ///< Whether or not to write validation details to a ROOT tree
    std::string m_rootTreeName;               ///< The ROOT tree name
    std::string m_rootFileName;               ///< The ROOT file name
    std::mt19937 m_rng;                       ///< The random number generator
    std::string m_mcParticleListName;         ///< The mc particle list name 
    LArDLHelper::TorchModel m_model;          ///< The model

};

} // namespace lar_dl_content

#endif // LAR_DL_GRAPH_CLUSTERING_ALGORITHM_H
