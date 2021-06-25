/**
 *  @file   larpandoracontent/LArThreeDReco/ThirdViewRecoveryAlgorithm.cc
 *
 *  @brief  Implementation of the clustering parent algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArThreeDReco/ThirdViewRecoveryAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

ThirdViewRecoveryAlgorithm::ThirdViewRecoveryAlgorithm():
    m_slidingLinearFitWindow(10),
    m_minimumHitsToProject(10)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThirdViewRecoveryAlgorithm::Run()
{

        //ClusterVector clusters3D;
        //this->GetThreeDClusters(clusters3D, clusterToPfoMap);
        //PfoList pTwoViewPfos;
        //ClusterToPfoMap clusterToPfoMap;
        //this->GetTwoViewPfos(pTwoViewPfos);
        this->GetTwoViewPfos();

        //this->AddThirdViewClusters(pTwoViewPfos, clusterToPfoMap);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------


//void ThirdViewRecoveryAlgorithm::GetThreeDClusters(ClusterVector &clusters3D, ClusterToPfoMap &clusterToPfoMap) const
void ThirdViewRecoveryAlgorithm::GetTwoViewPfos() const
{
    for (const std::string &pfoListName : m_inputPfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, pfoListName, pPfoList))
            continue;

        for (const Pfo *const pPfo : *pPfoList)
        {
//            ClusterList pfoClusters3D;
            ClusterList pfoClustersU;
            ClusterList pfoClustersV;
            ClusterList pfoClustersW;

            //LArPfoHelper::GetThreeDClusterList(pPfo, pfoClusters3D);
            LArPfoHelper::GetClusters(pPfo,pandora::TPC_VIEW_U,pfoClustersU);
            LArPfoHelper::GetClusters(pPfo,pandora::TPC_VIEW_V,pfoClustersV);
            LArPfoHelper::GetClusters(pPfo,pandora::TPC_VIEW_W,pfoClustersW);

            int ugood = (int)(pfoClustersU.size()>0);
            int vgood = (int)(pfoClustersV.size()>0);
            int wgood = (int)(pfoClustersW.size()>0);
            //int nGoodViews = (int)(pfoClustersU.size()>0) + (int)pfoClustersV.size()>0 + (int)(pfoClustersW.size()>0);
            int nGoodViews = ugood + vgood + wgood; 
            //std::cout << "test N GOOD VIEWS = " << nGoodViews << " U = " << ugood << " V = " << vgood << " W = " << wgood << std::endl;
            std::cout << "test N GOOD VIEWS = " << nGoodViews << std::endl;
            if (nGoodViews!=2) continue;

            ClusterList goodView1Clusters;
            ClusterList goodView2Clusters;
            std::string missingView="";

            //twoViewPfoList.push_back(pPfo); 

            if (pfoClustersU.size())
            {
                goodView1Clusters = pfoClustersU;
                if(pfoClustersV.size()) {goodView2Clusters = pfoClustersV; missingView="W";}
                else {goodView2Clusters = pfoClustersW; missingView="V";}
            }
            else
            {
                missingView="U";
                goodView1Clusters = pfoClustersV; 
                goodView2Clusters = pfoClustersW; 
            } 

            std::cout << "goodView1Clusters size = " << goodView1Clusters.size() << std::endl; 
            std::cout << "goodView2Clusters size = " << goodView2Clusters.size() << std::endl; 
            for (const Cluster *const pCluster1 : goodView1Clusters)
            {
                for(const Cluster *const pCluster2 : goodView2Clusters)
                {
                //try
                //{
                    std::cout << "---------------------------------------------------------------------" << std::endl;
                    std::cout << "pCluster1->GetNCaloHits() = " << pCluster1->GetNCaloHits() << " pCluster2->GetNCaloHits() = " << pCluster2->GetNCaloHits() << std::endl;
                    if(pCluster1->GetNCaloHits()<m_minimumHitsToProject || pCluster2->GetNCaloHits()<m_minimumHitsToProject) continue;

                   //Just for debugging purposes
//                    std::cout << "LArPfoHelper::IsTrack(pPfo) = " << LArPfoHelper::IsTrack(pPfo) << " LArPfoHelper::IsShower(pPfo) = " << LArPfoHelper::IsShower(pPfo) << " GetNumberOfTwoDHits = " << LArPfoHelper::GetNumberOfTwoDHits(pPfo) << std::endl;
//                    std::cout << "pfoClustersU.size() = " << pfoClustersU.size();
//                    std::cout << "pfoClustersV.size() = " << pfoClustersV.size();
//                    std::cout << "pfoClustersW.size() = " << pfoClustersW.size();
                    ClusterList firstCluster, secondCluster;
                    firstCluster.push_back(pCluster1);
                    secondCluster.push_back(pCluster2);
                    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&firstCluster, "cluster1", RED);
                    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&secondCluster, "cluster2", BLUE);
//                  PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&pfoClustersW, "ClustersW", GREEN);

                    const TwoDSlidingFitResult slidingFitResult1(pCluster1, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
                    const TwoDSlidingFitResult slidingFitResult2(pCluster2, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
                    float minX1(0), maxX1(0), minX2(0), maxX2(0), minZ1(0), maxZ1(0), minZ2(0), maxZ2(0);
                    //Call function to project sliding fit results in third view and return a list of clusters in third view that could be added in the pfo
                    slidingFitResult1.GetMinAndMaxX(minX1,maxX1);
                    slidingFitResult1.GetMinAndMaxZ(minZ1,maxZ1);
                    slidingFitResult2.GetMinAndMaxX(minX2,maxX2);
                    slidingFitResult2.GetMinAndMaxZ(minZ2,maxZ2);

                    //there should be some X overlap!
                    if(maxX1 < minX2 || maxX2 < minX1) continue;
                    std::cout << "debug0" << std::endl;
                    double commonXMin =  (minX1 < minX2) ? minX2 : minX1; 
                    double commonXMax =  (maxX1 > maxX2) ? maxX2 : maxX1;
                    double m_tolerance = 1;
                    commonXMin = commonXMin + m_tolerance;
                    commonXMax = commonXMax - m_tolerance;
                    std::cout << "minX1 = " << minX1 << " minX2 = " << minX2 << " maxX1 = " << maxX1 << " maxX2 = " << maxX2 << std::endl;
                    std::cout << "commonXMin = " << commonXMin << " commonXMax = " << commonXMax << std::endl;
                    //Project extremes in third view
                    //double minZ3_1(0), maxZ3_1(0), minZ3_2(0), maxZ3_2(0);
                    
                    pandora::CartesianVector view1AtCommonXMin(0,0,0), view1AtCommonXMax(0,0,0),  view2AtCommonXMin(0,0,0), view2AtCommonXMax(0,0,0), view3AtCommonXMin(0,0,0), view3AtCommonXMax(0,0,0);
                    if (STATUS_CODE_SUCCESS != slidingFitResult1.GetGlobalFitPositionAtX(commonXMin, view1AtCommonXMin)) continue;
                    std::cout << "debug1A: view1AtCommonXMin = " << view1AtCommonXMin << std::endl;
                    if (STATUS_CODE_SUCCESS != slidingFitResult1.GetGlobalFitPositionAtX(commonXMax, view1AtCommonXMax)) continue;
                    std::cout << "debug1B: view1AtCommonXMax = " << view1AtCommonXMax << std::endl;
                    if (STATUS_CODE_SUCCESS != slidingFitResult2.GetGlobalFitPositionAtX(commonXMin, view2AtCommonXMin)) continue;
                    std::cout << "debug1C: view2AtCommonXMin = " << view2AtCommonXMin << std::endl;
                    if (STATUS_CODE_SUCCESS != slidingFitResult2.GetGlobalFitPositionAtX(commonXMax, view2AtCommonXMax)) continue;
                    std::cout << "debug1D: view2AtCommonXMax = " << view2AtCommonXMax << std::endl;
                    float chi2min(0), chi2max(0);
                    if (missingView=="U") {
                        //double vAtCommonXMin(), vAtCommonXMax(), wAtCommonXMin(), wAtCommonXMax();
                        LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, view1AtCommonXMin, view2AtCommonXMin, view3AtCommonXMin, chi2min);
                        LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, view1AtCommonXMax, view2AtCommonXMax, view3AtCommonXMax, chi2max);
                        //thirdViewAtCommonXMin = uAtCommonXMin;
                        //thirdViewAtCommonXMax = uAtCommonXMax;
                        //minZ3=uAtCommonXMin.GetZ();
                        //maxZ3=uAtCommonXMax.GetZ();
                        //minZ3_1=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(commonXMin,minZ1);
                        //maxZ3_1=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(commonXMax,maxZ1);
                        //minZ3_2=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(commonXMin,minZ2);
                        //maxZ3_2=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->VWtoU(commonXMax,maxZ2);
                    }
                    else if (missingView=="V"){
                        LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, view1AtCommonXMin, view2AtCommonXMin, view3AtCommonXMin, chi2min);
                        LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_W, view1AtCommonXMax, view2AtCommonXMax, view3AtCommonXMax, chi2max);
                        //thirdViewAtCommonXMin = vAtCommonXMin;
                        //thirdViewAtCommonXMax = vAtCommonXMax;

                        //minZ3_1=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(commonXMin,minZ1);
                        //maxZ3_1=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(commonXMax,maxZ1);
                        //minZ3_2=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(commonXMin,minZ2);
                        //maxZ3_2=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->WUtoV(commonXMax,maxZ2);
                   }
                   else if (missingView=="W"){
                        LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, view1AtCommonXMin, view2AtCommonXMin, view3AtCommonXMin, chi2min);
                        LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, view1AtCommonXMax, view2AtCommonXMax, view3AtCommonXMax, chi2max);
                        //thirdViewAtCommonXMin = wAtCommonXMin;
                        //thirdViewAtCommonXMax = wAtCommonXMax;
                        //minZ3_1=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(commonXMin,minZ1);
                        //maxZ3_1=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(commonXMax,maxZ1);
                        //minZ3_2=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(commonXMin,minZ2);
                        //maxZ3_2=PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->UVtoW(commonXMax,maxZ2);
                   }

                   //const CartesianVector cluster3start(commonXMin,0,minZ3);
                   //const CartesianVector cluster3end(commonXMax,0,maxZ3);

                   PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&view3AtCommonXMin,"view3AtCommonXMin",BLACK,1);                   
                   PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&view3AtCommonXMax,"view3AtCommonXMax",BLACK,1);                   
                   //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&cluster3start_2,"cluster3start_2",BLUE,1);                   
                   //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(),&cluster3end_2,"cluster3end_2",BLUE,1);                   
                   PandoraMonitoringApi::ViewEvent(this->GetPandora());


                   //Loop over clusters in view 3
                   for (const std::string &clusterListName : m_inputClusterListNames)
                   {
                          const ClusterList *pClusterList(nullptr);

                          if (STATUS_CODE_SUCCESS != PandoraContentApi::GetList(*this, clusterListName, pClusterList))
                          continue;

                       for (const Cluster *const pCluster3 : *pClusterList)
                       {
                          if ((LArClusterHelper::GetClusterHitType(pCluster3)==LArClusterHelper::GetClusterHitType(pCluster1)) || (LArClusterHelper::GetClusterHitType(pCluster3)==LArClusterHelper::GetClusterHitType(pCluster2))) continue;
                          std::cout << "pCluster3->GetNCaloHits() = " << pCluster3->GetNCaloHits() << std::endl;
                          if(pCluster3->GetNCaloHits()<m_minimumHitsToProject) continue;
                          std::cout << "pCluster3->IsAvailable() = " << pCluster3->IsAvailable() << std::endl;
                          if(!pCluster3->IsAvailable()) continue;
                          std::cout << "LArClusterHelper::GetClusterHitType(pCluster3) = " << LArClusterHelper::GetClusterHitType(pCluster3) << std::endl;
                          ClusterList thirdCluster;
                          thirdCluster.push_back(pCluster3);
                          PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                          PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&thirdCluster, "cluster3", GREEN);
                          PandoraMonitoringApi::ViewEvent(this->GetPandora());

                          const TwoDSlidingFitResult slidingFitResult3(pCluster3, m_slidingLinearFitWindow, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
                          //Ask that there is an x overlap
                          //Define a direction for both this cluster and projected points
                          //Check if direction is within angular tolerance
                       } 
                   }




                //}
                //catch (StatusCodeException &)
                //{
                //}
            }
           
            


                //Just for debugging purposes
                /*std::cout << "---------------------------------------------------------------------" << std::endl;
                std::cout << "LArPfoHelper::IsTrack(pPfo) = " << LArPfoHelper::IsTrack(pPfo) << " LArPfoHelper::IsShower(pPfo) = " << LArPfoHelper::IsShower(pPfo) << " GetNumberOfTwoDHits = " << LArPfoHelper::GetNumberOfTwoDHits(pPfo) << std::endl;
                std::cout << "pfoClustersU.size() = " << pfoClustersU.size();
                std::cout << "pfoClustersV.size() = " << pfoClustersV.size();
                std::cout << "pfoClustersW.size() = " << pfoClustersW.size();
                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&pfoClustersU, "ClustersU", RED);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&pfoClustersV, "ClustersV", BLUE);
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&pfoClustersW, "ClustersW", GREEN);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());*/
                
                //Take every pair of clusters with a large time overlap?
                
//                for (const Cluster *const pCluster : pfoClustersU)
                //TwoDSlidingFitResultMap slidingFitResultMap;

//        try
//        {
//            (void)slidingFitResultMap.insert(
//                TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));
//        }
//        catch (StatusCodeException &)
//        {
//        }


            } 

//        pointVector.push_back(protoHit.GetPosition3D());

//    const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
//    const unsigned int layerWindow(m_slidingFitHalfWindow);
//        const ThreeDSlidingFitResult originalSlidingFitResult(&currentPoints3D, layerWindow, layerPitch);


//        const double rL(slidingFitResult.GetLongitudinalDisplacement(protoHit.GetPosition3D()));
//
//        if (STATUS_CODE_SUCCESS != slidingFitResult.GetGlobalFitPosition(rL, pointOnFit))
//            continue;

//        const double uFit(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoU(pointOnFit.GetY(), pointOnFit.GetZ()));
//        const double vFit(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoV(pointOnFit.GetY(), pointOnFit.GetZ()));
//        const double wFit(PandoraContentApi::GetPlugins(*this)->GetLArTransformationPlugin()->YZtoW(pointOnFit.GetY(), pointOnFit.GetZ()));



                //sort updated cluster list?
                //std::sort(clusters3D.begin(), clusters3D.end(), LArClusterHelper::SortByNHits);
                //pfo.m_clusterList = updatedClusterList


    



//            for (const Cluster *const pCluster3D : pfoClusters3D)
//            {
//                if (LArPfoHelper::IsTrack(pPfo) && (pCluster3D->GetNCaloHits() > m_maxHitsToConsider3DTrack))
//                    continue;
//
//                if (!clusterToPfoMap.insert(ClusterToPfoMap::value_type(pCluster3D, pPfo)).second)
//                    throw StatusCodeException(STATUS_CODE_ALREADY_PRESENT);
//                 
//                if()
//                clusters3D.push_back(pCluster3D);
                
//            }
        }
    }

//    std::sort(clusters3D.begin(), clusters3D.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//void ThirdViewRecoveryAlgorithm::AddThirdViewClusters(PfoList &twoViewPfoList, ClusterToPfoMap &clusterToPfoMap) const  //pass by reference? I don't have to modify the list... Or should I make this a private variable?
//{

//}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThirdViewRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputPfoListNames", m_inputPfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "InputClusterListNames", m_inputClusterListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinimumHitsToProject", m_minimumHitsToProject));


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
