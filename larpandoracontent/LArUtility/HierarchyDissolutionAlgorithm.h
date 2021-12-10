/**
 *  @file   larpandoracontent/LArUtility/HierarchyDissolutionAlgorithm.h
 *
 *  @brief  Header file for the pfo hit cleaning algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_HIERARCHY_DISSOLUTION_ALGORITHM_H
#define LAR_HIERARCHY_DISSOLUTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  HierarchyDissolutionAlgorithm class
 */
class HierarchyDissolutionAlgorithm : public pandora::Algorithm
{
private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StringVector m_pfoListNames;   ///< The input PFO list names
};

} // namespace lar_content

#endif // #ifndef LAR_HIERARCHY_DISSOLUTION_ALGORITHM_H
