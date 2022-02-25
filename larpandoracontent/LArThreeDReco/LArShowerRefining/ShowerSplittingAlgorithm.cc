/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerRefining/ShowerSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the shower hierarchy mop up algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArThreeDReco/LArShowerRefining/ShowerSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

StatusCode ShowerSplittingAlgorithm::Run()
{

    // Hacky location for new shower profile examination code!
    const PfoList *pPfoList(nullptr);
    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)) && pPfoList)
    {
        for (const Pfo *const pShowerPfo : *pPfoList)
        {
            ClusterList clusterList3D;
            LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);

            // Get the longitudinal and transverse shower profiles
            if (clusterList3D.empty())
                continue;

            CaloHitList caloHitList3D;
            clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

            if (caloHitList3D.size() < 2)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            // Begin with a PCA
            CartesianVector centroid(0.f, 0.f, 0.f);
            LArPcaHelper::EigenVectors eigenVecs;
            LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
            LArPcaHelper::RunPca(caloHitList3D, centroid, eigenValues, eigenVecs);

            // By convention, the primary axis has a positive z-component.
            const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

            // Place intercept at hit with minimum projection
            float minProjection(std::numeric_limits<float>::max());
            for (const CaloHit *const pCaloHit3D : caloHitList3D)
                minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

            const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

            // Now define ortho directions
            const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
                (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));
            const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());
            const CartesianVector orthoDirection2(axisDirection.GetCrossProduct(orthoDirection1).GetUnitVector());

            // Visualise axes
            std::cout << "axisDirection " << axisDirection << std::endl << "orthoDirection1 " << orthoDirection1 << std::endl << "orthoDirection2 " << orthoDirection2 << std::endl;
            PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);
            PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &clusterList3D, "ClusterList3D", GREEN);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &axisIntercept, "axisIntercept", RED, 1);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &centroid, "centroid", RED, 1);
            const CartesianVector centroidPlusAxisDir(centroid + axisDirection * 10.f);
            const CartesianVector centroidPlusOrthoDir1(centroid + orthoDirection1 * 10.f);
            const CartesianVector centroidPlusOrthoDir2(centroid + orthoDirection2 * 10.f);
            PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &centroid, &centroidPlusAxisDir, "axisDirection", RED, 1, 1);
            PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &centroid, &centroidPlusOrthoDir1, "orthoDirection1", RED, 1, 1);
            PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &centroid, &centroidPlusOrthoDir2, "orthoDirection2", RED, 1, 1);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            // Transverse profile
            TwoDHistogram transverseProfile(1001, -150.15, 150.15, 1001, -150.15, 150.15);

            for (const CaloHit *const pCaloHit3D : caloHitList3D)
            {
                const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
                const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
                const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
                transverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
            }

            std::cout << "Observed transverse energy profile " << std::endl;
            PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), transverseProfile, "COLZ");

            // Observed longitudinal profile
            Histogram observedLongitudinalProfile(140, 0., 140.);

            const float convertADCToMeV(0.0075f); // (c) Maria
            const float convertGeVToMeV(1000.f);
            const float convertCmToX0(1.f / 14.f);
            float clusterEnergyInMeV(0.f);

            for (const CaloHit *const pCaloHit3D : caloHitList3D)
            {
                const CaloHit *const pParentCaloHit(static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress()));

                if (TPC_VIEW_W != pParentCaloHit->GetHitType())
                    continue;

                clusterEnergyInMeV += convertADCToMeV * pParentCaloHit->GetInputEnergy(); // Used later on

                const float longitudinalCoordInCm((pCaloHit3D->GetPositionVector() - axisIntercept).GetDotProduct(axisDirection));
                observedLongitudinalProfile.Fill(longitudinalCoordInCm * convertCmToX0, convertADCToMeV * pParentCaloHit->GetInputEnergy());
            }

            std::cout << "Observed longitudinal energy profile " << std::endl;
            PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedLongitudinalProfile, "");

            // Expected longitudinal profile
            Histogram expectedLongitudinalProfile(140, 0., 140.);

            const float clusterEnergyInGeV(clusterEnergyInMeV / convertGeVToMeV);
            const float longProfileCriticalEnergy(0.08f);
            const float longProfileParameter0(1.25f);
            const float longProfileParameter1(0.5f);

            const double a(longProfileParameter0 + longProfileParameter1 * std::log(clusterEnergyInGeV / longProfileCriticalEnergy));
            const double gammaA(std::exp(lgamma(a)));

            float t(0.f);
            for (int iBin = 0; iBin < expectedLongitudinalProfile.GetNBinsX(); ++iBin)
            {
                t += expectedLongitudinalProfile.GetXBinWidth();
                expectedLongitudinalProfile.Fill(t, convertGeVToMeV * clusterEnergyInGeV / 2. * std::pow(t / 2.f, static_cast<float>(a - 1.)) *
                    std::exp(-t / 2.) * expectedLongitudinalProfile.GetXBinWidth() / gammaA);
            }

            std::cout << "Expected longitudinal energy profile " << std::endl;
            PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedLongitudinalProfile, "");
        }
    }

/*    const PfoList *pLeadingPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_leadingPfoListName, pLeadingPfoList));

    if (!pLeadingPfoList || pLeadingPfoList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerSplittingAlgorithm: unable to find pfos in provided list, " << m_leadingPfoListName << std::endl;

        return STATUS_CODE_SUCCESS;
    }

    PfoList parentShowerPfos;
    this->FindParentShowerPfos(pLeadingPfoList, parentShowerPfos);
    this->PerformPfoMerges(parentShowerPfos);
*/
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

/*void ShowerSplittingAlgorithm::FindParentShowerPfos(const PfoList *const pLeadingPfoList, PfoList &parentShowerPfos) const
{
    for (const Pfo *const pLeadingPfo : *pLeadingPfoList)
    {
        this->FindParentShowerPfos(pLeadingPfo, parentShowerPfos);
    }
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

/*void ShowerSplittingAlgorithm::FindParentShowerPfos(const Pfo *const pPfo, PfoList &parentShowerPfos) const
{
    if (LArPfoHelper::IsShower(pPfo))
    {
        if (pPfo->GetDaughterPfoList().empty())
            return;

        if (parentShowerPfos.end() != std::find(parentShowerPfos.begin(), parentShowerPfos.end(), pPfo))
            throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);

        parentShowerPfos.push_back(pPfo);
    }
    else
    {
        for (const Pfo *const pDaughterPfo : pPfo->GetDaughterPfoList())
            this->FindParentShowerPfos(pDaughterPfo, parentShowerPfos);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ShowerSplittingAlgorithm::PerformPfoMerges(const PfoList &parentShowerPfos) const
{
    for (const Pfo *const pParentShowerPfo : parentShowerPfos)
    {
        PfoList downstreamPfos;
        LArPfoHelper::GetAllDownstreamPfos(pParentShowerPfo, downstreamPfos);

        for (const Pfo *const pDownstreamPfo : downstreamPfos)
        {
            if (pDownstreamPfo != pParentShowerPfo)
                this->MergeAndDeletePfos(pParentShowerPfo, pDownstreamPfo);
        }
    }
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ShowerSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
