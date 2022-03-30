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
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArThreeDReco/LArShowerRefining/ShowerSplittingAlgorithm.h"

using namespace pandora;

const float convertGeVToMeV(1000.f);
const float convertCmToX0(1.f / 14.f);
const float criticalEnergyMeV(30.5f);
const float atomicNumberArgon(18);

const float longProfileCriticalEnergy(0.08f);
const float longProfileParameter0(1.25f);
const float longProfileParameter1(0.5f);
const float longProfileMaxDifference(0.1f);


//bool frequencysort(int i, int j) { return i > j; }
//bool frequencysort (int i,int j) { return (std::count_if(sortedMainMcPartVect.begin(), sortedMainMcPartVect.end(),[&](int const &i) {return (i == val);}))<std::count_if(sortedMainMcPartVect.begin(), sortedMainMcPartVect.end(),[&](int const &j) {return (j == val);}); }

bool mapValueComparison(const std::pair<int, int>& a,
         const std::pair<int, int>& b)
{
    return a.second < b.second;
}


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

/////////This comes from MyTrackShowerId algorithm///////
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));
   /* 
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // Get reconstructable MC
    LArMCParticleHelper::MCContributionMap primaryMCParticleToHitsMap;
    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_foldBackHierarchy = true;

    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters,
        LArMCParticleHelper::IsBeamNeutrinoFinalState, primaryMCParticleToHitsMap);
    MCParticleList primaryMCList;
    for (auto [ pMC, hits ] : primaryMCParticleToHitsMap)
    {   (void)hits;
        primaryMCList.emplace_back(pMC);
    }
*/
///////////////////////////////////////////////////////

    if ((STATUS_CODE_SUCCESS == PandoraContentApi::GetList(*this, m_pfoListName, pPfoList)) && pPfoList)
    {
        FloatVector position1Vect, position2Vect, energyVect, expectedTransverseProfileVect, expectedTransverseProfileRadiusVect;
        IntVector mainMcPartVect; int mainMcParticle; double mainMcParticleEnergy; 
        for (const Pfo *const pShowerPfo : *pPfoList)
        {       
            position1Vect.clear();
            position2Vect.clear();
            energyVect.clear();
            mainMcPartVect.clear();
            expectedTransverseProfileVect.clear();
            expectedTransverseProfileRadiusVect.clear();

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
            const int nTransverseProfileBins = 61;
            const float transverseProfileMin = -150.3;
            const float transverseProfileMax = 150.3;
            TwoDHistogram transverseProfile(nTransverseProfileBins, transverseProfileMin, transverseProfileMax, nTransverseProfileBins, transverseProfileMin, transverseProfileMax);
            TwoDHistogram expectedTransverseProfile(nTransverseProfileBins, transverseProfileMin, transverseProfileMax, nTransverseProfileBins, transverseProfileMin, transverseProfileMax);
            const float convertADCToMeV(0.0075f); // (c) Maria
            for (const CaloHit *const pCaloHit3D : caloHitList3D)
            {
                const CaloHit *const pParentCaloHit(static_cast<const CaloHit *>(pCaloHit3D->GetParentAddress()));
                if (TPC_VIEW_W != pParentCaloHit->GetHitType())
                    continue;
                //Find main contributing MC particle Id
                //MCParticleVector mcParticleVector;
                int mcParticleIndex(-999),iMcPart(0);
                for (const auto &weightMapEntry : pParentCaloHit->GetMCParticleWeightMap())
                {
                  if(weightMapEntry.second>0.5)
                  {
                    //std::cout << weightMapEntry.second << std::endl;
                    iMcPart=0; 
                    for(const MCParticle *const pMCParticle: *pMCParticleList)
                    {
                      if(pMCParticle==weightMapEntry.first) mcParticleIndex=iMcPart;
                      iMcPart++;
                    }
                  }
                  //mcParticleVector.push_back(weightMapEntry.first);
                //  std::cout << " weightMapEntry.first = " << weightMapEntry.first << " second = " << weightMapEntry.second << std::endl;
                }
                //std::cout << "largest contributing mc particle id = " << mcParticleIndex << std::endl;
                /*int mcParticleIndex(-999),iMcPart(0);
                for(const MCParticle *const pMCParticle: *pMCParticleList)
                {
                      //std::cout << *pMCParticle << std::endl;
                      if(pParentCaloHit->GetMCParticleWeightMap().at(*pMCParticle))
                      {
                        std::cout << "found" << std::endl;
                      }
                    const float weight(pParentCaloHit->GetMCParticleWeightMap().at(pMCParticle));
                    if(weight>0.5)mcParticleIndex=iMcPart;
                    iMcPart++;
                }
                */
                const CartesianVector hitCoordinate(pCaloHit3D->GetPositionVector() - centroid);
                const float position1(hitCoordinate.GetDotProduct(orthoDirection1));
                const float position2(hitCoordinate.GetDotProduct(orthoDirection2));
                position1Vect.push_back(position1);
                position2Vect.push_back(position2);
                energyVect.push_back(pCaloHit3D->GetInputEnergy());
                mainMcPartVect.push_back(mcParticleIndex);
                transverseProfile.Fill(position1, position2, pCaloHit3D->GetInputEnergy()); // Units: ADCs; note counting U, V and W parent hits here!
            }
            //Find main mc particle contributor to pfo
            IntVector sortedMainMcPartVect = mainMcPartVect;

            std::map<int, int> elementToFrequencyMap; 
            for (unsigned int i = 0; i < sortedMainMcPartVect.size(); i++) elementToFrequencyMap[sortedMainMcPartVect[i]]++;
            //int size = elementToFrequencyMap.size(); 
            mainMcParticle = std::max_element(std::begin(elementToFrequencyMap), std::end(elementToFrequencyMap), mapValueComparison)->first;
            int mainMcParticleFrequency = std::max_element(std::begin(elementToFrequencyMap), std::end(elementToFrequencyMap), mapValueComparison)->second;

//            std::sort(sortedMainMcPartVect.begin(), sortedMainMcPartVect.end(), frequencysort);
 //           mainMcParticle=sortedMainMcPartVect.at(0);
            
            std::cout << "MainMcParticle = " << mainMcParticle << " frequency = " << mainMcParticleFrequency << " total: " << sortedMainMcPartVect.size() << std::endl;
            //std::cout << "pfoId = " << pfoId << " Write to tree = " << m_writeToTree << " drawProfiles = " << m_drawProfiles << std::endl;
            MCParticleVector mcParticleVector;
            for (const MCParticle *pMCParticle : *pMCParticleList) //maybe move this outside pfo loop?
                 mcParticleVector.push_back(pMCParticle);

            mainMcParticleEnergy= mcParticleVector.at(mainMcParticle)->GetEnergy();  //how do we do this in Pandora without using iterators like this?

            std::cout << "Observed transverse energy profile " << std::endl;
            if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), transverseProfile, "COLZ");


            // Observed longitudinal profile
            Histogram observedLongitudinalProfile(101, -0.2, 40.2);
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

            std::cout << "mainMcParticleEnergy = " << mainMcParticleEnergy*1000 << " clusterEnergyInMeV = " << clusterEnergyInMeV << std::endl;
            std::cout << "Observed longitudinal energy profile " << std::endl;
            if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), observedLongitudinalProfile, "");

            //Expected transverse profile
            //const float showerDepthAtCentroid = (LArPfoHelper::GetVertex(pShowerPfo)->GetPosition() - centroid).GetMagnitude();//Distance between centroid and vertex
            const float showerDepthAtCentroid = (centroid-axisIntercept).GetMagnitude()*convertCmToX0;//Distance between centroid and vertex
            //const float showerDepthAtMaximum = std::log(clusterEnergyInMeV/criticalEnergyMeV)/std::log(2);
            const float showerDepthAtMaximum = std::log(mainMcParticleEnergy*1000/criticalEnergyMeV)/std::log(2);
            //const float showerDepthAtMaximum = std::log(mainMcParticleEnergy*1000/criticalEnergyMeV)-0.858;
            //const float showerDepthAtMaximum = std::log(mainMcParticleEnergy*1000/criticalEnergyMeV)-0.5;
            const float tau = showerDepthAtCentroid / showerDepthAtMaximum;
            std::cout << "showerDepthAtCentroid = " << showerDepthAtCentroid << " showerDepthAtMaximum = " << showerDepthAtMaximum << " tau = " << tau << std::endl;
            const float z1 = 0.0251 + 0.00319*std::log(clusterEnergyInMeV);
            const float z2 = 0.1162 - 0.000381*atomicNumberArgon;
            const float k1 = 0.659 - 0.00309*atomicNumberArgon;
            const float k2 = 0.645;
            const float k3 = -2.59;
            const float k4 = 0.3585 + 0.0421*std::log(clusterEnergyInMeV);
            const float p1 = 2.632 - 0.00094*atomicNumberArgon;
            const float p2 = 0.401 + 0.00187*atomicNumberArgon;
            const float p3 = 1.313 - 0.0686*std::log(clusterEnergyInMeV);
            const float Rc = z1+ z2*tau;
            const float Rt = k1*(std::exp(k3*(tau-k2))+std::exp(k4*(tau-k2)));
            const float prob = p1*std::exp((p2-tau)/p3 - std::exp((p2-tau)/p3));
        
            float radius(0.f), xDist(0.f), yDist(0.f), profile(0.f);
            //std::cout << "expectedTransverseProfile.GetNBinsX() = " << expectedTransverseProfile.GetNBinsX() << " nTransverseProfileBins = " << nTransverseProfileBins << std::endl;
            for(int iBin = 0; iBin < expectedTransverseProfile.GetNBinsX(); ++iBin)
            {
                for(int jBin = 0; jBin < expectedTransverseProfile.GetNBinsX(); ++jBin)
                {   
                    //std::cout << "deb1" << std::endl;
                    xDist = transverseProfileMin+iBin*expectedTransverseProfile.GetXBinWidth(); 
                    //std::cout << "deb2" << std::endl;
                    yDist = transverseProfileMin+jBin*expectedTransverseProfile.GetXBinWidth();   
                    //std::cout << "deb3" << std::endl;
                    radius=std::sqrt(xDist*xDist+yDist*yDist);
                    //std::cout << "deb4" << std::endl;
                    profile = 2*radius*(prob*Rc*Rc/std::pow((radius*radius+Rc*Rc),2)+(1-prob)*Rt*Rt/std::pow((radius*radius+Rt*Rt),2)); 
                    //std::cout << "deb5" << std::endl;
                    expectedTransverseProfile.SetBinContent(iBin,jBin,profile);
                    //std::cout << "deb6" << std::endl;
                    //std::cout << "iBin = " << iBin << " jBin = " << jBin << " xDist = " << xDist << " yDist = " << yDist << " radius = " << radius << " profile = " << profile << std::endl;
                }
            }
            
            //Fill vectors with expected tranvserse profile in bins of radius
            for(int iBin = 0; iBin < expectedTransverseProfile.GetNBinsX(); ++iBin){
                    radius = std::sqrt(2)*(transverseProfileMin+iBin*expectedTransverseProfile.GetXBinWidth());   
                    profile = 2*radius*(prob*Rc*Rc/std::pow((radius*radius+Rc*Rc),2)+(1-prob)*Rt*Rt/std::pow((radius*radius+Rt*Rt),2)); 
                    expectedTransverseProfileVect.push_back(profile); 
                    expectedTransverseProfileRadiusVect.push_back(radius);
            }

            if(m_writeToTree){
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "pfoId", pfoId);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "position1", &position1Vect);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "position2", &position2Vect);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "energy", &energyVect);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mainMcParticle", &mainMcPartVect);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "expectedTransverseProfileRadius", &expectedTransverseProfileRadiusVect);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "expectedTransverseProfile", &expectedTransverseProfileVect);
                PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "expectedTransverseProfileRadius", &expectedTransverseProfileRadiusVect);
                PandoraMonitoringApi::FillTree(this->GetPandora(), m_treeName.c_str());
            }


            float maximumExpectedTransverseProfile(0.f);
            int maximumExpectedTransverseProfileBinX(0);
            int maximumExpectedTransverseProfileBinY(0);
            expectedTransverseProfile.GetMaximum( maximumExpectedTransverseProfile,maximumExpectedTransverseProfileBinX,maximumExpectedTransverseProfileBinY);
            expectedTransverseProfile.Scale(1/maximumExpectedTransverseProfile);
            if(m_drawProfiles) PandoraMonitoringApi::DrawPandoraHistogram(this->GetPandora(), expectedTransverseProfile, "COLZ");

            // Expected longitudinal profile
            Histogram expectedLongitudinalProfile(101, -0.2, 40.2);


            const float clusterEnergyInGeV(clusterEnergyInMeV / convertGeVToMeV);
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
////////////////////////////
            //1)Get calo hit list for this pfo
            CaloHitList pShowerCaloHits;
            LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, pShowerCaloHits);
            //LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_W, pShowerCaloHits);
            //2)Fill McToCaloHitList map for these calo hits
/*            LArMCParticleHelper::MCRelationMap mcToTargetMCMap;
            LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToTargetMCMap);
            LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
            LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
            LArMCParticleHelper::GetMCParticleToCaloHitMatches(&pShowerCaloHits,mcToTargetMCMap,hitToMCMap,mcToTrueHitListMap);
            //3)Count how many MCparticles in this map contribute to at least 5 hits per view
            std::cout << "Number of true particles contributing to this pfo = " << mcToTrueHitListMap.size() << " mcToTargetMCMap size = " << mcToTargetMCMap.size() << " hitToMCMap size = " << hitToMCMap.size() << " pShowerCaloHits size = " << pShowerCaloHits.size() <<  std::endl;
            //4)Save this number to output tree
*/
    //CaloHitToCaloHitMap ThreeDHitToTwoDHitMap;

/*    for (const pandora::CaloHit* const pCaloHit : pShowerCaloHits) {
      if (pCaloHit->GetHitType() != pandora::TPC_3D)
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER); //found a non-3D hit in the input list*/
    CaloHitList pTwoDHitList;
    const pandora::CaloHit* parentAddress;
    for (const pandora::CaloHit* const pCaloHit : pShowerCaloHits)
    {
        parentAddress=static_cast<const pandora::CaloHit*>(pCaloHit->GetParentAddress());
        if(static_cast<const pandora::CaloHit*>(parentAddress)->GetHitType()==pandora::TPC_VIEW_W)pTwoDHitList.insert(pTwoDHitList.end(),static_cast<const pandora::CaloHit*>(parentAddress));
       
    }


/*      // ATTN get the 2D calo hit from the 3D calo hit
      if (!ThreeDHitToTwoDHitMap
             .insert(CaloHitToCaloHitMap::value_type(
               pCaloHit, static_cast<const pandora::CaloHit*>(pCaloHit->GetParentAddress())))
             .second)
        {std::cout << "Found repeated input hit" << std::endl; throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);} //found repeated input hits
    }

    //Loop over map to check what we've filled it with
    for (auto const& x : ThreeDHitToTwoDHitMap)
    {
        std::cout << "key hit type = " << x.first->GetHitType()  // string (key)
              << " value hit type = " 
              << x.second->GetHitType() // string's value 
              << std::endl;
    }*/
/*
            //2)Fill McToCaloHitList map for these calo hits
    LArMCParticleHelper::MCRelationMap mcToTargetMCMap;
    //LArMCParticleHelper::GetMCPrimaryMap(pMCParticleList, mcToTargetMCMap); //folded hierarchy
    LArMCParticleHelper::GetMCToSelfMap(pMCParticleList, mcToTargetMCMap);  //Unfolded hierarchy
    LArMCParticleHelper::CaloHitToMCMap hitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&pTwoDHitList,mcToTargetMCMap,hitToMCMap,mcToTrueHitListMap);

   MCParticleList pPfoContributingMCParticleList;
   std::vector<int> pPfoContributingMCParticleNHits;    
   for ( const auto &pair : mcToTrueHitListMap ) {
        std::cout << "mc particle energy = " << pair.first->GetEnergy() << " number of hits  = " << pair.second.size() << std::endl;
        pPfoContributingMCParticleList.insert(pPfoContributingMCParticleList.end(),pair.first);
        pPfoContributingMCParticleNHits.push_back(pair.second.size());
    }
         PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), false, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
         PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(),&pTwoDHitList, "2DCaloHitList", BLUE);
         PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(),&pShowerCaloHits, "3DCaloHitList", GREEN);
         PandoraMonitoringApi::VisualizeMCParticles(this->GetPandora(),&pPfoContributingMCParticleList, "ContributingMCParticles", RED);
         PandoraMonitoringApi::ViewEvent(this->GetPandora());

            std::cout << "Number of true particles contributing to this pfo = " << mcToTrueHitListMap.size() << " mcToTargetMCMap size = " << mcToTargetMCMap.size() << " hitToMCMap size = " << hitToMCMap.size() << " pShowerCaloHits size = " << pShowerCaloHits.size() << " pTwoDHitList size = " << pTwoDHitList.size() <<  std::endl;

/////////////////////////
         if(m_writeToTree)
          { 
            PandoraMonitoringApi::SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PfoCpntributingMCParticleHits", &pPfoContributingMCParticleNHits);
            PandoraMonitoringApi::FillTree(this->GetPandora(), m_treeName.c_str());
          }*/
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
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
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
