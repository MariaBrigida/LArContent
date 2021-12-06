/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/ParticleConsistencyAlgorithm.h
 *
 *  @brief  Header file for the track recovery algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PARTICLE_CONSISTENCY_ALGORITHM_H
#define LAR_PARTICLE_CONSISTENCY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ParticleConsistencyAlgorithm class
 */
class ParticleConsistencyAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ParticleConsistencyAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_inputPfoListName; ///< The PFO list name
};

} // namespace lar_content

#endif // #ifndef LAR_PARTICLE_CONSISTENCY_ALGORITHM_H

