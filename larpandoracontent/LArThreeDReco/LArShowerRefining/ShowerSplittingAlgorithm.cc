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

ShowerSplittingAlgorithm::ShowerSplittingAlgorithm():
    m_drawProfiles(false),
    m_writeToTree(false),
    m_fileName("OutputFile.root"),
    m_treeName("OutputTree")
{
}


ShowerSplittingAlgorithm::~ShowerSplittingAlgorithm()
{
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "ShowerSplittingAlgorithm: Unable to write tree " << m_treeName << " to file " << m_fileName << std::endl;
        }
    }
}


StatusCode ShowerSplittingAlgorithm::Run()
{

    // Hacky location for new shower profile examination code!
    const PfoList *pPfoList(nullptr);
    int pfoId(0);
    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)) && pPfoList)
    {
        FloatVector position1Vect, position2Vect, energyVect;
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
            TwoDHistogram transverseProfile(501, -150.3, 150.3, 501, -150.3, 150.3);
            const float convertADCToMeV(0.0075f); // (c) Maria
            for (const CaloHit *const pCaloHit3D : caloHitList3D)
            {
                const CaloHit *const pParentCaloHit(static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress()));
                if (TPC_VIEW_W != pParentCaloHit->GetHitType())
                    continue;
                const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
                const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
                const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
                position1Vect.push_back(position1);
                position2Vect.push_back(position2);
                energyVect.push_back(pCaloHit3D->GetInputEnergy());
                transverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
            }
            std::cout << "pfoId = " << pfoId << " Write to tree = " << m_writeToTree << " drawProfiles = " << m_drawProfiles << std::endl;
            if(m_writeToTree){
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoId", pfoId);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "position1", &position1Vect);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "position2", &position2Vect);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "energy", &energyVect);
                PandoraMonitoringApi::FillTree(this->GetPandora(), m_treeName.c_str());

                position1Vect.clear();
                position2Vect.clear();
                energyVect.clear();
            }

            std::cout << "Observed transverse energy profile " << std::endl;
            if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), transverseProfile, "COLZ");
            // Observed longitudinal profile
            Histogram observedLongitudinalProfile(101, -0.2, 40.2);
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
            if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedLongitudinalProfile, "");

            // Expected longitudinal profile
            Histogram expectedLongitudinalProfile(101, -0.2, 40.2);

            const float clusterEnergyInGeV(clusterEnergyInMeV / convertGeVToMeV);
            const float longProfileCriticalEnergy(0.08f);
            const float longProfileParameter0(1.25f);
            const float longProfileParameter1(0.5f);
            const float longProfileMaxDifference(0.1f);

            const double a(longProfileParameter0 + longProfileParameter1 * std::log(clusterEnergyInGeV / longProfileCriticalEnergy));
            const double gammaA(std::exp(lgamma(a)));
            float t(0.f);
            for (int iBin = 0; iBin < expectedLongitudinalProfile.GetNBinsX(); ++iBin)
            {
                t += expectedLongitudinalProfile.GetXBinWidth();
                expectedLongitudinalProfile.Fill(t, convertGeVToMeV * clusterEnergyInGeV / 2. * std::pow(t / 2.f, static_cast<float>(a - 1.)) *
                    std::exp(-t / 2.) * expectedLongitudinalProfile.GetXBinWidth() / gammaA);
            }

            // Compare the observed and expected longitudinal profiles
            int binOffsetAtMinDifference(0);
            float minProfileDifference(std::numeric_limits<float>::max());

            for (int iBinOffset = 0; iBinOffset < expectedLongitudinalProfile.GetNBinsX(); ++iBinOffset)
            {
                float profileDifference(0.);

                for (int iBin = 0; iBin < observedLongitudinalProfile.GetNBinsX(); ++iBin)
                {
                    if (iBin < iBinOffset)
                    {
                        profileDifference += observedLongitudinalProfile.GetBinContent(iBin);
                    }
                    else
                    {
                        profileDifference += std::fabs(expectedLongitudinalProfile.GetBinContent(iBin - iBinOffset) - observedLongitudinalProfile.GetBinContent(iBin));
                    }
                }

                if (profileDifference < minProfileDifference)
                {
                    minProfileDifference = profileDifference;
                    binOffsetAtMinDifference = iBinOffset;
                }

                if (profileDifference - minProfileDifference > longProfileMaxDifference)
                    break;
            }

            const float profileStart(binOffsetAtMinDifference * expectedLongitudinalProfile.GetXBinWidth());
            const float profileDiscrepancy((clusterEnergyInMeV > 0.f) ? minProfileDifference / clusterEnergyInMeV : -1.f);

            std::cout << "ClusterEnergyInMeV " << clusterEnergyInMeV << ", profileStart " << profileStart << ", profileDiscrepancy " << profileDiscrepancy << std::endl;

            std::cout << "Observed longitudinal energy profile " << std::endl;
            if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedLongitudinalProfile, "");

            std::cout << "Expected longitudinal energy profile " << std::endl;
            if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedLongitudinalProfile, "");
            
            pfoId++;
        }
    }

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
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DrawProfiles", m_drawProfiles));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteToTree", m_writeToTree));
    if (m_writeToTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    }

    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
