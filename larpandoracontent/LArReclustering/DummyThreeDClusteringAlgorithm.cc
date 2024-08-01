/**
 *  @file   larpandoracontent/LArReclustering/DummyThreeDClusteringAlgorithm.cc
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

#include "larpandoracontent/LArReclustering/DummyThreeDClusteringAlgorithm.h"


using namespace pandora;

namespace lar_content
{

DummyThreeDClusteringAlgorithm::DummyThreeDClusteringAlgorithm()
{
}

StatusCode DummyThreeDClusteringAlgorithm::Run()
{

    //Access the clusters
    //Split the clusters using truth info
    //Return the new cluster list

    //InitializeReclustering in has set the hits to be reclustered as current
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));

    if (pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(*pCaloHitList, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : *pCaloHitList)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));
    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());

    const Cluster *pClusterPos = nullptr;
    const Cluster *pClusterNeg = nullptr;


    ClusterVector clusterVector;
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const CartesianVector pCaloHitPosition = pCaloHit->GetPositionVector();
        std::cout << "projection = " << (pCaloHitPosition-centroid).GetDotProduct(centroid+orthoDirection1) << std::endl;
        if((pCaloHitPosition-centroid).GetDotProduct(centroid+orthoDirection1)<0)
		{
            if(!pClusterNeg)
		    {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit); 
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterNeg));
                clusterVector.push_back(pClusterNeg);
		    }
		    else
			{ 
				PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterNeg, pCaloHit));
			}
        }
        else
        {
            if(!pClusterPos)
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit); 
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterPos));
                clusterVector.push_back(pClusterPos);
            }
            else
            {

				PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterPos, pCaloHit));
            }

        }
    }
    for(const Cluster* pNewCluster: clusterVector)
    { 
      ClusterList newClusters;
      newClusters.push_back(pNewCluster);  
      //if(m_drawProfiles)PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &newClusters, "newClusters", AUTOITER));
    }
//    if(m_drawProfiles)PandoraMonitoringApi::ViewEvent(this->GetPandora());

    return STATUS_CODE_SUCCESS;
}


StatusCode DummyThreeDClusteringAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DrawProfiles", m_drawProfiles));
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
