/**
 *  @file   larpandoracontent/LArReclustering/SplitMergedShowersAlgorithm.h
 *
 *  @brief  Header file for the reclustering algorithm that runs other algs.
 *
 *  $Log: $
 */

#ifndef LAR_SPLIT_MERGED_SHOWERS_ALGORITHM_H
#define LAR_SPLIT_MERGED_SHOWERS_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{
/**
 *  @brief  RecursivePfoMopUpAlgorithm class
 */
class SplitMergedShowersAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    SplitMergedShowersAlgorithm();

   /**
    *  @brief  Destructor
    */
    ~SplitMergedShowersAlgorithm();

private:

    pandora::StatusCode Run();
    /**
     *  @brief Get the theoretical lateral profile value at the shower maximum for a given energy and radius (Grindhammer)
     *
     *  @param energy the cluster energy in MeV
     *  @param radiusInCm the radius at which to calculate the profile value
     *
     *  @return List of PfoMergeStats for each Pfo
     */
    float GetLateralProfileAtShowerMaximum(float clusterEnergyInMeV, float radiusInCm);

    float GetFigureOfMerit(pandora::CaloHitList caloHitList3D);

    //Decide whether to try reclustering for this pfo
    bool PassesCutsForReclustering(const pandora::ParticleFlowObject *const pPfo);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_pfoListName; ///< The input pfo list name (e.g. list of neutrino or testbeam pfos)
    bool m_drawProfiles; //Boolean to enable and disable displaying transverse profiles
};

} // namespace lar_content

#endif // #ifndef LAR_SPLIT_MERGED_SHOWERS_ALGORITHM_H
