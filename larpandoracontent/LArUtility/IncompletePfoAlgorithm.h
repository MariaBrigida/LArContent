/**
 *  @file   larpandoracontent/LArUtility/IncompletePfoAlgorithm.h
 *
 *  @brief  Header file for the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_INCOMPLETE_PFO_ALGORITHM_H
#define LAR_INCOMPLETE_PFO_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  IncompletePfoAlgorithm class
 */
class IncompletePfoAlgorithm : public pandora::Algorithm
{
public:
    IncompletePfoAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;       ///< The input cluster list names
    pandora::StringVector m_clusterListNames;   ///< The input cluster list names
    pandora::StringVector m_caloHitListNames;   ///< The input cluster list names
};

} // namespace lar_content

#endif // #ifndef LAR_INCOMPLETE_PFO_ALGORITHM_H
