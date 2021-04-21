/**
 *  @file
 * larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewOvershootTracksTool.cc
 *
 *  @brief  Implementation of the two view overshoot tracks tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

#include "larpandoracontent/LArThreeDReco/LArTwoViewMatching/TwoViewOvershootTracksTool.h"

using namespace pandora;

namespace lar_content
{

TwoViewOvershootTracksTool::TwoViewOvershootTracksTool() :
    TwoViewThreeDKinkBaseTool(1),
    m_splitMode(true),
    m_cosThetaCutForKinkSearch(0.94f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewOvershootTracksTool::GetIteratorListModifications(
    TwoViewTransverseTracksAlgorithm *const pAlgorithm, const IteratorList &iteratorList, ModificationList &modificationList) const
{
    for (IteratorList::const_iterator iIter1 = iteratorList.begin(), iIter1End = iteratorList.end(); iIter1 != iIter1End; ++iIter1)
    {
        for (IteratorList::const_iterator iIter2 = iIter1; iIter2 != iIter1End; ++iIter2)
        {
            if (iIter1 == iIter2)
                continue;

            try
            {
                const unsigned int nMatchedReUpsampledSamplingPoints1((*iIter1)->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints());
                const unsigned int nMatchedReUpsampledSamplingPoints2((*iIter2)->GetOverlapResult().GetNMatchedReUpsampledSamplingPoints());
                IteratorList::const_iterator iIterA((nMatchedReUpsampledSamplingPoints1 >= nMatchedReUpsampledSamplingPoints2) ? iIter1 : iIter2);
                IteratorList::const_iterator iIterB((nMatchedReUpsampledSamplingPoints1 >= nMatchedReUpsampledSamplingPoints2) ? iIter2 : iIter1);

                Particle particle(*(*iIterA), *(*iIterB));
                const LArPointingCluster pointingClusterA(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA));
                const LArPointingCluster pointingClusterB(pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB));

                LArPointingCluster::Vertex vertexA, vertexB;
                LArPointingClusterHelper::GetClosestVerticesInX(pointingClusterA, pointingClusterB, vertexA, vertexB);

                if (!this->PassesVertexCuts(vertexA, vertexB))
                    continue;

                this->SetSplitPosition(vertexA, vertexB, particle);

                const bool isALowestInX(this->IsALowestInX(pointingClusterA, pointingClusterB));
                const bool isThreeDKink(m_majorityRulesMode ? true : this->IsThreeDKink(pAlgorithm, particle, isALowestInX));

                if (isThreeDKink != m_splitMode)
                    continue;

                // Construct the modification object
                Modification modification;

                if (m_splitMode)
                {
                    modification.m_splitPositionMap[particle.m_pCommonCluster].push_back(particle.m_splitPositionCommon);
                }
                else
                {
                    const bool vertexAIsLowX(vertexA.GetPosition().GetX() < vertexB.GetPosition().GetX());
                    const Cluster *const pLowXCluster(vertexAIsLowX ? particle.m_pClusterA : particle.m_pClusterB);
                    const Cluster *const pHighXCluster(vertexAIsLowX ? particle.m_pClusterB : particle.m_pClusterA);
                    modification.m_clusterMergeMap[pLowXCluster].push_back(pHighXCluster);
                }

                modification.m_affectedClusters.push_back(particle.m_pCommonCluster);
                modification.m_affectedClusters.push_back(particle.m_pClusterA);
                modification.m_affectedClusters.push_back(particle.m_pClusterB);

                modificationList.push_back(modification);
            }
            catch (StatusCodeException &statusCodeException)
            {
                if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
                    throw statusCodeException;

                continue;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewOvershootTracksTool::PassesVertexCuts(const LArPointingCluster::Vertex &vertexA, const LArPointingCluster::Vertex &vertexB) const
{
    float longitudinalAB(-std::numeric_limits<float>::max()), transverseAB(std::numeric_limits<float>::max());
    LArPointingClusterHelper::GetImpactParameters(vertexA, vertexB, longitudinalAB, transverseAB);

    float longitudinalBA(-std::numeric_limits<float>::max()), transverseBA(std::numeric_limits<float>::max());
    LArPointingClusterHelper::GetImpactParameters(vertexB, vertexA, longitudinalBA, transverseBA);

    if (std::min(longitudinalAB, longitudinalBA) < m_minLongitudinalImpactParameter)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TwoViewOvershootTracksTool::SetSplitPosition(
    const LArPointingCluster::Vertex &vertexA, const LArPointingCluster::Vertex &vertexB, Particle &particle) const
{
    particle.m_splitPositionAB = (vertexA.GetPosition() + vertexB.GetPosition()) * 0.5;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TwoViewOvershootTracksTool::IsThreeDKink(TwoViewTransverseTracksAlgorithm *const pAlgorithm, const Particle &particle, const bool isALowestInX) const
{
    try
    {
        const TwoDSlidingFitResult &lowXFitResult(isALowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA)
                                                               : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB));
        const TwoDSlidingFitResult &highXFitResult(isALowestInX ? pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterB)
                                                                : pAlgorithm->GetCachedSlidingFitResult(particle.m_pClusterA));
        const TwoDSlidingFitResult &fitResultCommon(pAlgorithm->GetCachedSlidingFitResult(particle.m_pCommonCluster));

        const float minusX(this->GetXSamplingPoint(particle.m_splitPositionAB, false, fitResultCommon, lowXFitResult));
        const float plusX(this->GetXSamplingPoint(particle.m_splitPositionAB, true, fitResultCommon, highXFitResult));

        CartesianVector minusAB(0.f, 0.f, 0.f), splitAB(particle.m_splitPositionAB), plusAB(0.f, 0.f, 0.f);
        CartesianVector minusCommon(0.f, 0.f, 0.f), splitCommon(particle.m_splitPositionCommon), plusCommon(0.f, 0.f, 0.f);

        if ((STATUS_CODE_SUCCESS != lowXFitResult.GetGlobalFitPositionAtX(minusX, minusAB)) ||
            (STATUS_CODE_SUCCESS != highXFitResult.GetGlobalFitPositionAtX(plusX, plusAB)) ||
            (STATUS_CODE_SUCCESS != fitResultCommon.GetGlobalFitPositionAtX(minusX, minusCommon)) ||
            (STATUS_CODE_SUCCESS != fitResultCommon.GetGlobalFitPositionAtX(plusX, plusCommon)))
        {
            return true; // majority rules, by default
        }

        // Extract results
        const HitType hitTypeAB(LArClusterHelper::GetClusterHitType(particle.m_pClusterA));
        const HitType hitTypeCommon(LArClusterHelper::GetClusterHitType(particle.m_pCommonCluster));

        CartesianVector minus(0.f, 0.f, 0.f), split(0.f, 0.f, 0.f), plus(0.f, 0.f, 0.f);
        float chi2Minus(std::numeric_limits<float>::max()), chi2Split(std::numeric_limits<float>::max()),
            chi2Plus(std::numeric_limits<float>::max());
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitTypeAB, hitTypeCommon, minusAB, minusCommon, minus, chi2Minus);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitTypeAB, hitTypeCommon, splitAB, splitCommon, split, chi2Split);
        LArGeometryHelper::MergeTwoPositions3D(this->GetPandora(), hitTypeAB, hitTypeCommon, plusAB, plusCommon, plus, chi2Plus);

        // Apply final cuts
        const CartesianVector minusToSplit((split - minus).GetUnitVector());
        const CartesianVector splitToPlus((plus - split).GetUnitVector());
        const float dotProduct(minusToSplit.GetDotProduct(splitToPlus));

        if (dotProduct > m_cosThetaCutForKinkSearch)
            return false;
    }
    catch (StatusCodeException &)
    {
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoViewOvershootTracksTool::Particle::Particle(const MatrixType::Element &elementA, const MatrixType::Element &elementB) :
    m_splitPositionAB(0.f, 0.f, 0.f),
    m_splitPositionCommon(0.f, 0.f, 0.f)
{
    m_pCommonCluster = (elementA.GetCluster1() == elementB.GetCluster1()) ? elementA.GetCluster1() : elementA.GetCluster2();
    m_pClusterA = (elementA.GetCluster1() == elementB.GetCluster1()) ? elementA.GetCluster2() : elementA.GetCluster1();
    m_pClusterB = (elementA.GetCluster1() == elementB.GetCluster1()) ? elementB.GetCluster2() : elementB.GetCluster1();

    if (m_pClusterA == m_pClusterB)
        throw StatusCodeException(STATUS_CODE_FAILURE);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TwoViewOvershootTracksTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SplitMode", m_splitMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "CosThetaCutForKinkSearch", m_cosThetaCutForKinkSearch));

    return TwoViewThreeDKinkBaseTool::ReadSettings(xmlHandle);
}

} // namespace lar_content
