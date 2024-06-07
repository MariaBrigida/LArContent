/*   larpandoradlcontent/LArControlFlow/DlThreeDReclusteringAlgorithm.h
 *
 *  @brief  Header file for the master algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_THREE_D_RECLUSTERING_ALGORITHM_H
#define LAR_DL_THREE_D_RECLUSTERING_ALGORITHM_H 1

#include "larpandoracontent/LArReclustering/ThreeDReclusteringAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  MasterAlgorithm class
 */
class DlThreeDReclusteringAlgorithm : public lar_content::ThreeDReclusteringAlgorithm
{

public:
    /**
     *  @brief  Default constructor
     */
    DlThreeDReclusteringAlgorithm() = default;

};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_THREE_D_RECLUSTERING_ALGORITHM_H
