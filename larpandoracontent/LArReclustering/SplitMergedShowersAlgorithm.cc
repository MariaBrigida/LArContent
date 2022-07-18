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

bool sortByCaloHits (pandora::CaloHitList a, pandora::CaloHitList b) { return (a.size()>b.size()); }

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
    if(m_drawProfiles)PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));

    // Get shower pfos and then find the 3D cluster in the shower pfo.
    const PfoList *pShowerPfoList(nullptr);
    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, "ShowerParticles3D", pShowerPfoList)) && pShowerPfoList)
    {   
          std::vector<float> mainClusterFractionVector, initialFigureOfMeritVector, newFigureOfMeritVector, nHitsInitialFomVector, nFinalClustersVector;

        for (const Pfo *const pShowerPfo : *pShowerPfoList)
        {
            ClusterList clusterList3D;
            LArPfoHelper::GetThreeDClusterList(pShowerPfo, clusterList3D);

            if(!this->PassesCutsForReclustering(pShowerPfo)) continue; // this just checks it's a shower at the moment

            // Get the longitudinal and transverse shower profiles
            if (clusterList3D.empty())
                continue;

            CaloHitList caloHitList3D;
            clusterList3D.front()->GetOrderedCaloHitList().FillCaloHitList(caloHitList3D);

            //Quality cuts
            if (caloHitList3D.size() < 2)
                throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

            float initialFigureOfMerit=this->GetTransverseProfileFigureOfMerit(caloHitList3D);
            //std::cout << "Figure of merit = " << initialFigureOfMerit << std::endl;

            //Now let's free the hits in this cluster, so that they are available for reclustering!
            //ask to remove the 3D cluster from the parent pfo, so that it's not owned any more
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pShowerPfo, clusterList3D.front()));
            //pop this cluster in a local clusterlist
            const ClusterList reclusterClusterList(1, clusterList3D.front());
            const TrackList reclusterTrackList; //dummy track list

                        
            // Initialize reclustering with these local lists
            std::string currentClustersListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentClustersListName));
            //std::cout << "debug: current list name before reclustering = " << currentClustersListName << std::endl;
            
            //Some pfos are shower-like and yet include track-like 3D clusters. For the moment I don't want to deal with these.
            const ClusterList *pShowerClusters(nullptr);
            PandoraContentApi::GetList(*this, "ShowerClusters3D", pShowerClusters);
            if(pShowerClusters->end() == std::find(pShowerClusters->begin(), pShowerClusters->end(), reclusterClusterList.front())) continue;

            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, "ShowerClusters3D"));


            // Specify clusters and tracks to be used in reclustering
            std::string originalClustersListName;
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeReclustering(*this, reclusterTrackList, reclusterClusterList, originalClustersListName));

            if(m_drawProfiles)
            {
                PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),&reclusterClusterList, "ClustersToBeReclustered", RED);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }
            //Call the reclustering algos that produce new cluster candidates
            for (StringVector::const_iterator clusteringIter = m_clusteringAlgorithms.begin(), clusteringIterEnd = m_clusteringAlgorithms.end();
                clusteringIter != clusteringIterEnd; ++clusteringIter)
            {
                // Produce new cluster candidates
                std::string reclusterListName;
                const ClusterList *pReclusterList = NULL;
                PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunClusteringAlgorithm(*this, *clusteringIter, 
                    pReclusterList, reclusterListName));
   
                //std::cout << "Did I manage to fill the recluster list? size = " << pReclusterList->size() << std::endl; //(Who is actually filling this? the reclustering algo itself doesn't set the current list to have the new clusters in!)
                
                if (pReclusterList->empty())
                    continue;
                 
                //std::cout << "debug x0" << std::endl; 
                //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_associationAlgorithmName));
                //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RunDaughterAlgorithm(*this, m_trackClusterAssociationAlgName));
                if(m_drawProfiles)
                {
                    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
                    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(),pReclusterList, "ReclusteredClusters", BLUE);
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                }
                //std::cout << "debug x1" << std::endl; 
                
                //Loop over clusters and calculate new FOM including them if they have more than 10 hits
                
                std::vector<CaloHitList> newClustersCaloHitLists3D;  //I'm doing this because the function calls for two arguments, but I am essentially passing it the same input in both
                for(const Cluster *const pNewCluster : *pReclusterList)
                {
                  CaloHitList newClusterCaloHitList3D;
                  pNewCluster->GetOrderedCaloHitList().FillCaloHitList(newClusterCaloHitList3D);

                  if(newClusterCaloHitList3D.size()<30) continue;

                  newClustersCaloHitLists3D.push_back(newClusterCaloHitList3D);
                 
                  if(m_drawProfiles)
                  {
                      PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &newClusterCaloHitList3D, "newClusterCaloHitList3D", AUTOITER));                 
                      PandoraMonitoringApi::ViewEvent(this->GetPandora());
                  }
                }
                //std::cout << "debug x2" << std::endl; 

                std::sort (newClustersCaloHitLists3D.begin(), newClustersCaloHitLists3D.end(), sortByCaloHits);
                //std::cout << "debug x3" << std::endl; 
                //std::cout << "newClustersCaloHitLists3D size = " << newClustersCaloHitLists3D.size() << std::endl;
                if(!newClustersCaloHitLists3D.size()) continue; 
                /*for(CaloHitList testList : newClustersCaloHitLists3D)
                {
                  std::cout << "test list size = " << testList.size() << std::endl;
                }*/
                //std::cout << "newClustersCaloHitLists3D.front().size() = " << newClustersCaloHitLists3D.front().size() << " caloHitList3D.size() = " << caloHitList3D.size() << std::endl; 
                float mainClusterFraction = (float)newClustersCaloHitLists3D.front().size()/caloHitList3D.size();
                //std::cout << "mainClusterFraction = " << mainClusterFraction << std::endl;
                //sortedNewClustersCaloHitLists3D 
                float newFigureOfMerit = this->GetTransverseProfileFigureOfMerit(caloHitList3D, newClustersCaloHitLists3D);
                //std::cout << "Initial figure of merit = " << initialFigureOfMerit << " new figure of merit = " << newFigureOfMerit << std::endl;                     
                mainClusterFractionVector.push_back(mainClusterFraction);  ///watch out, these are in the loop over many algorithms! if I add more clustering algos I will need to differentiate the entries in these
                initialFigureOfMeritVector.push_back(initialFigureOfMerit);
                newFigureOfMeritVector.push_back(newFigureOfMerit);
                nHitsInitialFomVector.push_back(caloHitList3D.size());
                nFinalClustersVector.push_back(newClustersCaloHitLists3D.size());
            }
            //Return to the original cluster list
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentClustersListName));

            // Choose the best recluster candidates, which may still be the originalsi. Watch out though, because if I chose the reclustered list, this won't contain the track-like clusters that I got rid of by only selecting shower clusters!!!!! Ah no, I don't think it matters, as I've specified the clusters to be re clustered in reclusterClusterList.
            //I suppose I need to calculate FOM in the loop, then make decision outside, here. How to access the different lists from here?
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndReclustering(*this, originalClustersListName));

        }
          
        //Histogram highPurityFOMratio(20, 0, 2);
        //Histogram lowPurityFOMratio(20, 0, 2);
        for(long unsigned int iFom=0; iFom<newFigureOfMeritVector.size(); iFom++)
        {
            //std::cout << "iFom " << iFom << " purity " << mainClusterFractionVector.at(iFom) << " initialFom = " << initialFigureOfMeritVector.at(iFom) << " finalFom = " << newFigureOfMeritVector.at(iFom) << " nHitsInitialFom " << nHitsInitialFomVector.at(iFom) << " nFinalClusters " << nFinalClustersVector.at(iFom) << std::endl;
            std::cout << iFom << " " << mainClusterFractionVector.at(iFom) << " " << initialFigureOfMeritVector.at(iFom) << " " << newFigureOfMeritVector.at(iFom) << " " << nHitsInitialFomVector.at(iFom) << " " << nFinalClustersVector.at(iFom) << std::endl;
        //    if(mainClusterFractionVector.at(iFom)>0.9)highPurityFOMratio.Fill(newFigureOfMeritVector.at(iFom)/initialFigureOfMeritVector.at(iFom));
        //    else lowPurityFOMratio.Fill(newFigureOfMeritVector.at(iFom)/initialFigureOfMeritVector.at(iFom));
        }
        //PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), highPurityFOMratio);
        //PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), lowPurityFOMratio);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());

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




float SplitMergedShowersAlgorithm::GetTransverseProfileFigureOfMerit(CaloHitList mergedClusterCaloHitList3D, std::vector<CaloHitList> newClustersCaloHitLists3D)
{
    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));
    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());
    const CartesianVector orthoDirection2(axisDirection.GetCrossProduct(orthoDirection1).GetUnitVector());

    //Estimate total cluster energy
    float clusterEnergyInMeV(0);
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D) clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;

    //Observed transverse profile
    //int transverseProfileNBins = 1001;
    int transverseProfileNBins = 50;
    //float transverseProfileLow = -150.15;
    float transverseProfileLow = -50;
    float transverseProfileHigh = 50;
    float transverseProfileBinSize = (transverseProfileHigh-transverseProfileLow)/transverseProfileNBins;
    TwoDHistogram observedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
        const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
        const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
        observedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
    }

    //Scale observed profile to total cluster energy
    observedTransverseProfile.Scale(clusterEnergyInMeV/observedTransverseProfile.GetCumulativeSum());



    TwoDHistogram expectedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
    if(newClustersCaloHitLists3D.size()==1 || newClustersCaloHitLists3D.at(0)==mergedClusterCaloHitList3D)
    {
        //Expected tranvserse profile (Grindhammer parametrisation)
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
        expectedTransverseProfile.Scale(clusterEnergyInMeV/expectedTransverseProfile.GetCumulativeSum());
    }

    //Now, if I passed a new set of cluster hit lists in the second argument, I should separately calculate all of their transverse profiles projecting them onto the same plane as merged shower.
    else if(newClustersCaloHitLists3D.size()!=1 and newClustersCaloHitLists3D.at(0)!=mergedClusterCaloHitList3D)
    {
        std::vector<TwoDHistogram> newObservedTransverseProfiles;
        std::vector<double> newClusterEnergies;
        std::vector<float> newClustersCenterPositionsX, newClustersCenterPositionsY;
        for(CaloHitList newCaloHitList3D: newClustersCaloHitLists3D)
        {
            TwoDHistogram newObservedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);;
             
            double clusterEnergy(0);
            for (const CaloHit *const pCaloHit3D : newCaloHitList3D)
            {
                const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
                const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
                const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
                newObservedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
                clusterEnergy+=pCaloHit3D->GetInputEnergy();
            }
            //newObservedTransverseProfiles.push_back(newObservedTransverseProfile);
            if(m_drawProfiles)
            { 
                PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), newObservedTransverseProfile, "COLZ");
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }
            newClusterEnergies.push_back(clusterEnergy);
            newClustersCenterPositionsX.push_back(newObservedTransverseProfile.GetMeanX());
            newClustersCenterPositionsY.push_back(newObservedTransverseProfile.GetMeanY());
        }

        //Expected tranvserse profile (Grindhammer parametrisation as a combination of N shower profiles)
        for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
        {
            float profileX=transverseProfileLow+iBinX*transverseProfileBinSize;
            for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
            {
               float profileY=transverseProfileLow+iBinY*transverseProfileBinSize;
               //float profileRadius=std::sqrt(profileX*profileX+profileY*profileY);
               float profileValue(0), shiftedRadius(0), shiftedX(0), shiftedY(0);
               for(std::vector<double>::size_type iCluster=0; iCluster<newClusterEnergies.size(); iCluster++){
                   shiftedX=profileX-newClustersCenterPositionsX.at(iCluster);
                   shiftedY=profileY-newClustersCenterPositionsY.at(iCluster);
                   shiftedRadius=std::sqrt(shiftedX*shiftedX+shiftedY*shiftedY);
                   profileValue+=GetLateralProfileAtShowerMaximum(newClusterEnergies.at(iCluster)*convertADCToMeV,shiftedRadius);
               }
               expectedTransverseProfile.SetBinContent(iBinX, iBinY, profileValue);
            }
        }
        expectedTransverseProfile.Scale(clusterEnergyInMeV/expectedTransverseProfile.GetCumulativeSum());
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

float SplitMergedShowersAlgorithm::GetLongitudinalProfileFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

   // Observed longitudinal profile
    Histogram observedLongitudinalProfile(140, 0., 140.);

    const float convertADCToMeV(0.0075f); // (c) Maria
    const float convertGeVToMeV(1000.f);
    const float convertCmToX0(1.f / 14.f);
    float clusterEnergyInMeV(0.f);

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
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
    if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedLongitudinalProfile, "");

    return 1; //placeholder for fom

}


float SplitMergedShowersAlgorithm::GetTransverseProfileFigureOfMerit(CaloHitList mergedClusterCaloHitList3D)
{
    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(mergedClusterCaloHitList3D, centroid, eigenValues, eigenVecs);

    // By convention, the primary axis has a positive z-component.
    const CartesianVector axisDirection(eigenVecs.at(0).GetZ() > 0.f ? eigenVecs.at(0) : eigenVecs.at(0) * -1.f);

    // Place intercept at hit with minimum projection
    float minProjection(std::numeric_limits<float>::max());
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
        minProjection = std::min(minProjection, axisDirection.GetDotProduct(pCaloHit3D->GetPositionVector() - centroid));

    const CartesianVector axisIntercept(centroid + (axisDirection * minProjection));

    // Now define ortho directions
    const CartesianVector seedDirection((axisDirection.GetX() < std::min(axisDirection.GetY(), axisDirection.GetZ())) ? CartesianVector(1.f, 0.f, 0.f) :
        (axisDirection.GetY() < std::min(axisDirection.GetX(), axisDirection.GetZ())) ? CartesianVector(0.f, 1.f, 0.f) : CartesianVector(0.f, 0.f, 1.f));
    const CartesianVector orthoDirection1(seedDirection.GetCrossProduct(axisDirection).GetUnitVector());
    const CartesianVector orthoDirection2(axisDirection.GetCrossProduct(orthoDirection1).GetUnitVector());

    //Estimate total cluster energy
    float clusterEnergyInMeV(0);
    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D) clusterEnergyInMeV += pCaloHit3D->GetInputEnergy()*convertADCToMeV;

    //Observed transverse profile
    //int transverseProfileNBins = 1001;
    int transverseProfileNBins = 50;
    //float transverseProfileLow = -150.15;
    float transverseProfileLow = -50;
    float transverseProfileHigh = 50;
    float transverseProfileBinSize = (transverseProfileHigh-transverseProfileLow)/transverseProfileNBins;
    TwoDHistogram observedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);

    for (const CaloHit *const pCaloHit3D : mergedClusterCaloHitList3D)
    {
        const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
        const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
        const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
        observedTransverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
    }

    //Scale observed profile to total cluster energy
    //std::cout << "clusterEnergyInMeV = " << clusterEnergyInMeV << " observedTransverseProfile.GetCumulativeSum() = " << observedTransverseProfile.GetCumulativeSum() << std::endl;
    observedTransverseProfile.Scale(clusterEnergyInMeV/observedTransverseProfile.GetCumulativeSum());



    TwoDHistogram expectedTransverseProfile(transverseProfileNBins, transverseProfileLow, transverseProfileHigh, transverseProfileNBins, transverseProfileLow, transverseProfileHigh);
    //Expected tranvserse profile (Grindhammer parametrisation)
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
    expectedTransverseProfile.Scale(clusterEnergyInMeV/expectedTransverseProfile.GetCumulativeSum());

    //Calculate figure of merit for this cluster
    float squaredDiffSum(0);
    for (int iBinX=0; iBinX<transverseProfileNBins; iBinX++) 
    {
        for (int iBinY=0; iBinY<transverseProfileNBins; iBinY++)
        {
          float diff = expectedTransverseProfile.GetBinContent(iBinX, iBinY)-observedTransverseProfile.GetBinContent(iBinX,iBinY);
          //std::cout << "diff = " << diff << " obs = " << observedTransverseProfile.GetBinContent(iBinX,iBinY) << " exp = " << expectedTransverseProfile.GetBinContent(iBinX,iBinY) << std::endl;
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






bool SplitMergedShowersAlgorithm::PassesCutsForReclustering(const pandora::ParticleFlowObject *const pShowerPfo)
{
    if (LArPfoHelper::IsShower(pShowerPfo)) return true;
    return false;
}

StatusCode SplitMergedShowersAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "DrawProfiles", m_drawProfiles));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmList(*this, xmlHandle, "clusteringAlgorithms", m_clusteringAlgorithms));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
