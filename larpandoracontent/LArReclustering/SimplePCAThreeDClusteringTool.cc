/**
 *  @file   larpandoracontent/LArReclustering/SimplePCAThreeDClusteringTool.cc
 *
 *  @brief  Implementation file for the reclustering algorithm that uses transverse calorimetric profiles.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArReclustering/SimplePCAThreeDClusteringTool.h"
#include "larpandoracontent/LArHelpers/LArPointingClusterHelper.h"

using namespace pandora;

namespace lar_content
{

SimplePCAThreeDClusteringTool::SimplePCAThreeDClusteringTool()
{
}

bool SimplePCAThreeDClusteringTool::Run(const Algorithm *const /*pAlgorithm*/, std::vector<pandora::CaloHitList*> &newCaloHitListsVector)
{
    if (newCaloHitListsVector.empty() || newCaloHitListsVector.size() != 1)
        return false;

    CaloHitList *initialCaloHitList = newCaloHitListsVector.at(0);
    newCaloHitListsVector.clear();

    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(*initialCaloHitList, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : *initialCaloHitList)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));
    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());

    CaloHitList *posCaloHitList = new CaloHitList;
    CaloHitList *negCaloHitList = new CaloHitList;

    for (const CaloHit *const pCaloHit : *initialCaloHitList)
    {
        const CartesianVector pCaloHitPosition = pCaloHit->GetPositionVector();

        if((pCaloHitPosition-centroid).GetDotProduct(centroid+orthoDirection1)<0)
		{
            negCaloHitList->push_back(pCaloHit);
        }
		else
		{ 
            posCaloHitList->push_back(pCaloHit);
		}
    }

    newCaloHitListsVector.push_back(posCaloHitList);
    newCaloHitListsVector.push_back(negCaloHitList);

    //Get pos list PCA
    CartesianVector centroidPos(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecsPos;
    LArPcaHelper::EigenValues eigenValuesPos(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(*posCaloHitList, centroidPos, eigenValuesPos, eigenVecsPos);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirectionPos(eigenVecsPos.at(0).GetZ() > 0.f ? eigenVecsPos.at(0) : eigenVecsPos.at(0) * -1.f);


    //Get neg list PCA
    CartesianVector centroidNeg(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecsNeg;
    LArPcaHelper::EigenValues eigenValuesNeg(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(*negCaloHitList, centroidNeg, eigenValuesNeg, eigenVecsNeg);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirectionNeg(eigenVecsNeg.at(0).GetZ() > 0.f ? eigenVecsNeg.at(0) : eigenVecsNeg.at(0) * -1.f);

    //Get intersection point of the two principal axes
    CartesianVector intersectionPoint(0.f, 0.f, 0.f);
    float displacementPos(0.f),displacementNeg(0.f);
    LArPointingClusterHelper::GetIntersection(centroidPos,axisDirectionPos,centroidNeg,axisDirectionNeg,intersectionPoint,displacementPos,displacementNeg);

    //Clear pos and neg lists
    posCaloHitList->clear();
    negCaloHitList->clear();

    //Loop over original hit list, check whether it is within smaller cone of pos or neg axis, attach to relevant list
    float cosConeAxisPos(0), cosConeAxisNeg(0);
    for (const CaloHit *const pCaloHit3D : *initialCaloHitList)
    {
        cosConeAxisPos = axisDirectionPos.GetCosOpeningAngle(pCaloHit3D->GetPositionVector()-intersectionPoint);
        cosConeAxisNeg = axisDirectionNeg.GetCosOpeningAngle(pCaloHit3D->GetPositionVector()-intersectionPoint);
        if(cosConeAxisPos>cosConeAxisNeg) posCaloHitList->push_back(pCaloHit3D);
        else negCaloHitList->push_back(pCaloHit3D);
    }

    return true;
}


StatusCode SimplePCAThreeDClusteringTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
