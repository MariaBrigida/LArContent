/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the longitudinal association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/LongitudinalAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

LongitudinalAssociationAlgorithm::LongitudinalAssociationAlgorithm() :
    m_minClusterLayers(4),
    m_maxGapLayers(7),
    m_fitLayers(30),
    m_maxGapDistanceSquared(10.f),
    m_minCosRelativeAngle(0.985f),
    m_maxTransverseDisplacement(2.f),
    m_maxLongitudinalDisplacement(2.f),
    m_hitSizeZ(0.3f),
    m_hitSizeX(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        const Cluster *const pCluster = *iter;

        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < m_minClusterLayers)
            continue;

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerLayer);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // ATTN This method assumes that clusters have been sorted by layer
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterIEnd = clusterVector.end(); iterI != iterIEnd; ++iterI)
    {
        const Cluster *const pInnerCluster = *iterI;

        for (ClusterVector::const_iterator iterJ = iterI, iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
        {
            const Cluster *const pOuterCluster = *iterJ;

            if (pInnerCluster == pOuterCluster)
                continue;

            if (!this->AreClustersAssociated(pInnerCluster, pOuterCluster))
                continue;

            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
    const unsigned int currentLayer(isForward ? pCurrentCluster->GetOuterPseudoLayer() : pCurrentCluster->GetInnerPseudoLayer());
    const unsigned int testLayer(isForward ? pTestCluster->GetOuterPseudoLayer() : pTestCluster->GetInnerPseudoLayer());

    if (isForward && ((testLayer > currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    if (!isForward && ((testLayer < currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalAssociationAlgorithm::AreClustersAssociated(const Cluster *const pInnerCluster, const Cluster *const pOuterCluster) const
{
    if (pOuterCluster->GetInnerPseudoLayer() < pInnerCluster->GetInnerPseudoLayer())
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);

    // TODO Remove hardcoded numbers
    if ((pOuterCluster->GetInnerPseudoLayer() < pInnerCluster->GetInnerPseudoLayer() + 3) ||
        (pInnerCluster->GetOuterPseudoLayer() + 3 > pOuterCluster->GetOuterPseudoLayer()))
    {
        return false;
    }

    if ((pInnerCluster->GetOuterPseudoLayer() > pOuterCluster->GetInnerPseudoLayer() + 1) ||
        (pOuterCluster->GetInnerPseudoLayer() > pInnerCluster->GetOuterPseudoLayer() + m_maxGapLayers))
    {
        return false;
    }

    if ((2 * pInnerCluster->GetOuterPseudoLayer() < pOuterCluster->GetInnerPseudoLayer() + pInnerCluster->GetInnerPseudoLayer()) ||
        (pInnerCluster->GetOuterPseudoLayer() + pOuterCluster->GetOuterPseudoLayer() < 2 * pOuterCluster->GetInnerPseudoLayer()))
    {
        return false;
    }

    const CartesianVector innerEndCentroid(pInnerCluster->GetCentroid(pInnerCluster->GetOuterPseudoLayer()));
    const CartesianVector outerStartCentroid(pOuterCluster->GetCentroid(pOuterCluster->GetInnerPseudoLayer()));
    const CartesianVector outerEndCentroid(pOuterCluster->GetCentroid(pOuterCluster->GetOuterPseudoLayer()));

    if ((innerEndCentroid - outerStartCentroid).GetMagnitudeSquared() > m_maxGapDistanceSquared)
        return false;

    ClusterFitResult innerEndFit, outerStartFit;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterFitHelper::FitEnd(pInnerCluster, m_fitLayers, innerEndFit));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterFitHelper::FitStart(pOuterCluster, m_fitLayers, outerStartFit));

    ClusterList innerClusters, outerClusters;
    innerClusters.push_back(pInnerCluster);
    outerClusters.push_back(pOuterCluster);
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&innerClusters, "InnerClusters", RED);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&outerClusters, "OuterClusters", BLUE);

    std::cout << "Calling overloaded function..." <<std::endl;
    if (this->AreClustersAssociated(innerEndCentroid, outerStartCentroid, outerEndCentroid, innerEndFit, outerStartFit))
        {
	std::cout <<"CLUSTERS ARE ASSOCIATED!!!" << std::endl;
	PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return true;}

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalAssociationAlgorithm::AreClustersAssociated(const CartesianVector &innerClusterEnd, const CartesianVector &outerClusterStart,
    const CartesianVector &outerClusterEnd, const ClusterFitResult &innerFit, const ClusterFitResult &outerFit) const
{
    if (!innerFit.IsFitSuccessful() || !outerFit.IsFitSuccessful())
        return false;

    std::cout << "debug: opening angle = " << innerFit.GetDirection().GetCosOpeningAngle(outerFit.GetDirection()) << " min angle = " << m_minCosRelativeAngle << std::endl;
    if (innerFit.GetDirection().GetCosOpeningAngle(outerFit.GetDirection()) < m_minCosRelativeAngle)
        return false;

    const CartesianVector innerEndFit1(innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(innerClusterEnd - innerFit.GetIntercept())));
    const CartesianVector innerEndFit2(outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(innerClusterEnd - outerFit.GetIntercept())));

    const CartesianVector outerStartFit1(outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(outerClusterStart - outerFit.GetIntercept())));
    const CartesianVector outerStartFit2(innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(outerClusterStart - innerFit.GetIntercept())));
    const CartesianVector outerEndFit1(outerFit.GetIntercept() + outerFit.GetDirection() * (outerFit.GetDirection().GetDotProduct(outerClusterEnd - outerFit.GetIntercept())));
    const CartesianVector outerEndFit2(innerFit.GetIntercept() + innerFit.GetDirection() * (innerFit.GetDirection().GetDotProduct(outerClusterEnd - innerFit.GetIntercept())));

    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&innerEndFit1,"innerEndFit1",BLACK,2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&innerEndFit2,"innerEndFit2",GREEN,2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&outerStartFit1,"outerStartFit1",BLACK,2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&outerStartFit2,"outerStartFit2",GREEN,2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&outerEndFit1,"outerEndFit1",BLACK,3);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&outerEndFit2,"outerEndFit2",GREEN,3);

    ///New vars debug
    const CartesianVector fittedOuterEndSeparation(outerEndFit2 - outerEndFit1);
    std::cout << "debug: fittedOuterEndSeparation.GetX() " << std::fabs(fittedOuterEndSeparation.GetX()) <<  std::endl;
    std::cout << "debug: fittedOuterEndSeparation.GetZ() " << std::fabs(fittedOuterEndSeparation.GetZ()) <<  std::endl;
    //////

    const CartesianVector clusterSeparation(outerClusterStart - innerClusterEnd);

    std::cout << "debug: std::fabs(clusterSeparation.GetX()) " << std::fabs(clusterSeparation.GetX()) << " m_hitSizeX * m_maxTransverseDisplacement = " << m_hitSizeX * m_maxTransverseDisplacement << std::endl;
    std::cout << "debug: std::fabs(clusterSeparation.GetZ()) " << std::fabs(clusterSeparation.GetZ()) << " m_hitSizeZ * m_maxLongitudinalDisplacement = " << m_hitSizeZ * m_maxLongitudinalDisplacement << std::endl;
    std::cout << "condition = " << ((std::fabs(clusterSeparation.GetX()) < m_hitSizeX * m_maxTransverseDisplacement) && (std::fabs(clusterSeparation.GetZ()) < m_hitSizeZ * m_maxLongitudinalDisplacement)) << std::endl;
    if ((std::fabs(clusterSeparation.GetX()) < m_hitSizeX * m_maxTransverseDisplacement) &&
        (std::fabs(clusterSeparation.GetZ()) < m_hitSizeZ * m_maxLongitudinalDisplacement))
        return true;

    const CartesianVector fittedSeparation(outerStartFit1 - innerEndFit1);

    std::cout << "debug: fittedSeparation.GetX() " << std::fabs(fittedSeparation.GetX()) << " m_hitSizeX * m_maxTransverseDisplacement " << m_hitSizeX * m_maxTransverseDisplacement << std::endl;
    std::cout << "debug: fittedSeparation.GetZ() " << std::fabs(fittedSeparation.GetZ()) << " m_hitSizeZ * m_maxLongitudinalDisplacement " << m_hitSizeZ * m_maxLongitudinalDisplacement << std::endl;
    if ((std::fabs(fittedSeparation.GetX()) < m_hitSizeX * m_maxTransverseDisplacement) &&
        (std::fabs(fittedSeparation.GetZ()) < m_hitSizeZ * m_maxLongitudinalDisplacement))
        return true;

    const CartesianVector fittedInnerSeparation(innerEndFit2 - innerEndFit1);

    std::cout << "debug: fittedInnerSeparation.GetX() " << std::fabs(fittedInnerSeparation.GetX()) << " m_hitSizeX * m_maxTransverseDisplacement " << m_hitSizeX * m_maxTransverseDisplacement << std::endl;
    std::cout << "debug: fittedInnerSeparation.GetZ() " << std::fabs(fittedInnerSeparation.GetZ()) << " m_hitSizeZ * m_maxLongitudinalDisplacement " << m_hitSizeZ * m_maxLongitudinalDisplacement << std::endl;

    if ((std::fabs(fittedInnerSeparation.GetX()) < m_hitSizeX * m_maxTransverseDisplacement) &&
        (std::fabs(fittedInnerSeparation.GetZ()) < m_hitSizeZ * m_maxLongitudinalDisplacement))
        return true;

    const CartesianVector fittedOuterSeparation(outerStartFit2 - outerStartFit1);

    if ((std::fabs(fittedOuterSeparation.GetX()) < m_hitSizeX * m_maxTransverseDisplacement) &&
        (std::fabs(fittedOuterSeparation.GetZ()) < m_hitSizeZ * m_maxLongitudinalDisplacement))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxGapLayers", m_maxGapLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FitLayers", m_fitLayers));

    float maxGapDistance = std::sqrt(m_maxGapDistanceSquared);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxGapDistance", maxGapDistance));
    m_maxGapDistanceSquared = maxGapDistance * maxGapDistance;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosRelativeAngle", m_minCosRelativeAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTransverseDisplacement", m_maxTransverseDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxLongitudinalDisplacement", m_maxLongitudinalDisplacement));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HitSizeZ", m_hitSizeZ));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "HitSizeX", m_hitSizeX));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
