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

    //const Cluster *pClusterPos = nullptr;
    //const Cluster *pClusterNeg = nullptr;
    CaloHitList *posCaloHitList = new CaloHitList;
    CaloHitList *negCaloHitList = new CaloHitList;


    //std::cout << "debug before calo hit list loop" << std::endl;
    for (const CaloHit *const pCaloHit : *initialCaloHitList)
    {
        const CartesianVector pCaloHitPosition = pCaloHit->GetPositionVector();
        //std::cout << "projection = " << (pCaloHitPosition-centroid).GetDotProduct(centroid+orthoDirection1) << std::endl;

        if((pCaloHitPosition-centroid).GetDotProduct(centroid+orthoDirection1)<0)
		{
            negCaloHitList->push_back(pCaloHit);
            //std::cout << "this was neg, done" << std::endl;
        }
		else
		{ 
            posCaloHitList->push_back(pCaloHit);
            //std::cout << "this was pos, done" << std::endl;
		}

    }

    newCaloHitListsVector.push_back(posCaloHitList);
    newCaloHitListsVector.push_back(negCaloHitList);

    return true;
}


StatusCode SimplePCAThreeDClusteringTool::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
