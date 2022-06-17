/**
 *  @file   larpandoracontent/LArReclustering/SplitMergedShowersAlgorithm.cc
 *
 *  @brief  Implementation file for the reclustering algorithm that runs other algs.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArReclustering/SplitMergedShowersAlgorithm.h"

const float convertADCToMeV = 0.0075;
const int atomicNumberArgon = 18;
const float atomicMassArgon = 39.948;
const float criticalEnergyArgon = 32.84;
const float moliereRadiusCmArgon = 9.043; //cm

using namespace pandora;

namespace lar_content
{

SplitMergedShowersAlgorithm::SplitMergedShowersAlgorithm():
    m_drawProfiles(false)
{
}

SplitMergedShowersAlgorithm::~SplitMergedShowersAlgorithm()
{
}

StatusCode SplitMergedShowersAlgorithm::Run()
{
    // Get shower pfos and then find the 3D cluster in the shower pfo.
    const PfoList *pShowerPfoList(nullptr);
    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, "ShowerParticles3D", pShowerPfoList)) && pShowerPfoList)
    {
        for (const Pfo *const pShowerPfo : *pShowerPfoList)
        {
            ClusterList clusterList3D;
            LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);

            // Get the longitudinal and transverse shower profiles
            if (clusterList3D.empty())
                continue;

            CaloHitList caloHitList3D;
            clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

            //Quality cuts
            if (caloHitList3D.size() < 2)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            std::cout << "Figure of merit = " << this->GetFigureOfMerit(caloHitList3D) << std::endl;

            //Now let's free the hits in this cluster, so that they are available for reclustering!
            //ask to remove the 3D cluster from the parent pfo, so that it's not owned any more
            //ClusterToPfoMap &clusterToPfoMap;
            //const Pfo *const pPfo(clusterToPfoMap.at(pCluster3D));
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pShowerPfo, clusterList3D.front()));
            //pop this cluster in a local clusterlist
            ClusterList tmpClusterList;
            tmpClusterList.push_back(clusterList3D.front());
            //probably need a local dummy tracklist too

            //then you'd call InitializeReclustering
        }
    }


    return STATUS_CODE_SUCCESS;
}

//Lateral profile at the shower maximum
float SplitMergedShowersAlgorithm::GetLateralProfileAtShowerMaximum(float clusterEnergyInMeV, float radiusInCm){
    const float tau=1; //shower maximum
    float energy=clusterEnergyInMeV/criticalEnergyArgon;
    float radius=radiusInCm/moliereRadiusCmArgon;
    float z1 = 0.0251 + 0.00319*std::log(energy);
    float z2 = 0.1162 - 0.000381*atomicNumberArgon;
    float k1 = 0.659 - 0.00309*atomicNumberArgon;
    float k2 = 0.645;
    float k3 = -2.59;
    float k4 = 0.3585 + 0.0421*std::log(energy);
    float p1 = 2.632 - 0.00094*atomicNumberArgon;
    float p2 = 0.401 + 0.00187*atomicNumberArgon;
    float p3 = 1.313 - 0.0686*std::log(energy);
    float Rc = z1+ z2*tau;
    float Rt = k1*(std::exp(k3*(tau-k2))+std::exp(k4*(tau-k2)));
    float prob = p1*std::exp((p2-tau)/p3 - std::exp((p2-tau)/p3));
    float profile = 2*radius*(prob*Rc*Rc/std::pow((radius*radius+Rc*Rc),2)+(1-prob)*Rt*Rt/std::pow((radius*radius+Rt*Rt),2));
    return profile;

}

float SplitMergedShowersAlgorithm::GetFigureOfMerit(CaloHitList caloHitList3D)
{
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

    //Estimate total cluster energy
    float clusterEnergyInMeV(0);
    for (const CaloHit *const pCaloHit3D : caloHitList3D) clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;

    //Observed transverse profile
    //int transverseProfileNBins = 1001;
    int transverseProfileNBins = 50;
    //float transverseProfileLow = -150.15;
    float transverseProfileLow = -50;
    float transverseProfileHigh = 50;
    float transverseProfileBinSize = (transverseProfileHigh-transverseProfileLow)/transverseProfileNBins;
    TwoDHistogram observedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);

    for (const CaloHit *const pCaloHit3D : caloHitList3D)
    {
        const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
        const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
        const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
        observedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
    }

    //Scale observed profile to total cluster energy
    observedTransverseProfile.Scale(clusterEnergyInMeV/observedTransverseProfile.GetCumulativeSum());

    //Expected tranvserse profile (Grindhammer parametrisation)
    TwoDHistogram expectedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
    for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
    {
        float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
        for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
        {
           float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
           float profileRadius=std::sqrt(profileX*profileX+profileY*profileY);
           float profileValue=GetLateralProfileAtShowerMaximum(clusterEnergyInMeV,profileRadius);
           expectedTransverseProfile.SetBinContent(iBinX, iBinY, profileValue);
        }
    }

    //Calculate figure of merit for this cluster
    float squaredDiffSum(0);
    for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
    {
        for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
        {
          float diff = expectedTransverseProfile.GetBinContent(iBinX, iBinY)-observedTransverseProfile.GetBinContent(iBinX,iBinY);
          float squaredDiff = diff*diff;     
          squaredDiffSum+=squaredDiff;
        }
    }
    float figureOfMerit = squaredDiffSum/clusterEnergyInMeV;

    if(m_drawProfiles) {
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedTransverseProfile, "COLZ");
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedTransverseProfile, "COLZ");
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
    
    return figureOfMerit;
}

StatusCode SplitMergedShowersAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DrawProfiles", m_drawProfiles));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
