/**
 *  @file   larpandoradlcontent/LArTrackShowerId/DlTrackShowerAlgorithm.h
 *
 *  @brief  Header file for the graph neural network track shower Id.
 *
 *  $Log: $
 */
#ifndef LAR_DL_GRAPH_TRACK_SHOWER_ID_ALGORITHM_H
#define LAR_DL_GRAPH_TRACK_SHOWER_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

#include <random>

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  DeepLearningGraphTrackShowerIdAlgorithm class
 */
class DlGraphTrackShowerIdAlgorithm : public pandora::Algorithm
{
public:
    typedef std::map<std::pair<int, int>, std::vector<const pandora::CaloHit *>> PixelToCaloHitsMap;

    /**
     *  @brief Default constructor
     */
    DlGraphTrackShowerIdAlgorithm();

    virtual ~DlGraphTrackShowerIdAlgorithm();

private:

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    pandora::StatusCode PrepareTrainingSample();
    //pandora::StatusCode Infer();
    int GetMainMcParticleIndex(const pandora::CaloHit *const pCaloHit);
    int GetMainMcParticlePdgCode(const pandora::CaloHit *const pCaloHit);



    bool m_trainingMode;                      ///< Training mode
    std::string m_trainingOutputFile;         ///< Output file name for training examples
    int m_event;                              ///< The current event number
    std::mt19937 m_rng;                       ///< The random number generator
    std::string m_mcParticleListName;         ///< The mc particle list name 
    LArDLHelper::TorchModel m_model;


};

} // namespace lar_dl_content

#endif // LAR_DL_GRAPH_TRACK_SHOWER_ID_ALGORITHM_H
                                      
