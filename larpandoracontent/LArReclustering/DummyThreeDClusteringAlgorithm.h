/**
 *  @file   larpandoracontent/LArReclustering/DummyThreeDClusteringAlgorithm.h
 *
 *  @brief  Header file for the reclustering algorithm that uses transverse calorimetric profiles.
 *
 *  $Log: $
 */

#ifndef LAR_DUMMY_THREE_D_CLUSTERING_ALGORITHM_H
#define LAR_DUMMY_THREE_D_CLUSTERING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{
/**
 *  @brief  RecursivePfoMopUpAlgorithm class
 */
class DummyThreeDClusteringAlgorithm : public pandora::Algorithm
{
public:

    DummyThreeDClusteringAlgorithm();

private:

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);
    bool m_drawProfiles; //Boolean to enable and disable displaying transverse profiles
};

} // namespace lar_content

#endif // #endif LAR_DUMMY_THREE_D_CLUSTERING_ALGORITHM_H
