/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoRecovery/ParticleConsistencyAlgorithm.cc
 *
 *  @brief  Implementation of the particle recovery algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPcaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArThreeDReco/LArPfoRecovery/ParticleConsistencyAlgorithm.h"

#include <algorithm>
#include <numeric>

using namespace pandora;

namespace lar_content
{

// To be moved to CaloUtility file later

// p is mom in MeV, M is particle mass in MeV
float bethe_bloch(const float p, const float M)
{
    const float K{0.307075f};       // MeV g^-1 cm^2
    const float z{1.f};             // incident charge
    const float zOverA{0.450586f};  // Argon z = 18, A = 39.948
    const float Me{0.511f};        // electron mass
    const float rho{1.3954f};       // liquid argon density (1.3954 g cm^-3 at boiling point)
    const float pOverM{p / M};
    const float gamma{std::sqrt(1.f + pOverM * pOverM)};
    const float beta{std::sqrt(1.f - 1.f / (gamma * gamma))};
    const float betaGamma{pOverM};
    const float massRatio{Me / M};
    const float Tmax{2.f * Me * betaGamma * betaGamma / (1 + 2 * gamma * massRatio + massRatio * massRatio)};
    const float hOmegaP = std::sqrt(rho * zOverA) * 28.816 * 1e-6;  // Units are eV, so convert to MeV - see PDG review
    
    const float coeff{K * z * z * zOverA / (beta * beta)};
    const float logParam{2.f * Me * Tmax / (hOmegaP * hOmegaP)};

    const float dedx{coeff * (0.5f * std::log(logParam) - (beta * beta) + 0.5f)};

    return dedx;
}

ParticleConsistencyAlgorithm::ParticleConsistencyAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleConsistencyAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_inputPfoListName, pPfoList));

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        CaloHitList caloHitList;
        LArPfoHelper::GetCaloHits(pPfo, HitType::TPC_3D, caloHitList);
        if (caloHitList.empty())
            continue;

        const CartesianVector &vertex{LArPfoHelper::GetVertex(pPfo)->GetPosition()};
        // Get the PCA axis in the 3D, project the hits onto it and produce longitudinal and transverse energy profiles
        // Compare to theoretical profiles for tracks and showers. Also look for rphi style peaks to find evidence for
        // merged particles in one pfo
        CartesianVector centroid(0.f, 0.f, 0.f);
        LArPcaHelper::EigenValues eigenVals(0.f, 0.f, 0.f);
        LArPcaHelper::EigenVectors eigenVecs;
        LArPcaHelper::RunPca(caloHitList, centroid, eigenVals, eigenVecs);

        // Having more than one non-negligible eigen value could indicate either shower or mis-merged track
        const CartesianVector lAxis(eigenVecs[0]);
        const CartesianVector tAxis(-lAxis.GetZ(), lAxis.GetY(), lAxis.GetX());
        const float lVertex{vertex.GetDotProduct(lAxis)};
        const float tVertex{vertex.GetDotProduct(tAxis)};
        // Aside - tracks will likely have very asymmetric projected transverse hit positions
        std::vector<std::pair<const float, const float>> lProjEnergyVec;
        std::vector<std::pair<const float, const float>> tProjEnergyVec;
        float lMin{std::numeric_limits<float>::max()}, lMax{-std::numeric_limits<float>::max()};
        float tMin{std::numeric_limits<float>::max()}, tMax{-std::numeric_limits<float>::max()};
        // Project hits onto the longitudinal and transverse axes
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            const float lProj{pos.GetDotProduct(lAxis)};
            const float tProj{pos.GetDotProduct(tAxis) - tVertex};
            if (lProj < lMin)
                lMin = lProj;
            if (lProj > lMax)
                lMax = lProj;
            if (tProj < tMin)
                tMin = tProj;
            if (tProj > tMax)
                tMax = tProj;
            tProjEnergyVec.emplace_back(std::make_pair(tProj, pCaloHit->GetInputEnergy()));
        }
        // Determine which end of the longitudinal projection is closest to the vertex
        float start{0.f};
        if (std::abs(lMin - lVertex) <= std::abs(lMax - lVertex))
            start = lMin;
        else
            start = lMax;
        // Adjust projection so that the extremal hit closest to the vertex is longitudinally projected onto zero
        float totalEnergy{0.f};
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            const float lProj{std::abs(pos.GetDotProduct(lAxis) - start)};
            const float hitEnergy{pCaloHit->GetElectromagneticEnergy()};
            lProjEnergyVec.emplace_back(std::make_pair(lProj, hitEnergy));
            totalEnergy += hitEnergy;
        }
        // GeV to MeV
        totalEnergy *= 1e3;
        const float lBinMax{std::max(std::abs(lMax - start), std::abs(lMin - start)) + std::numeric_limits<float>::epsilon()};
        const float radiationLength{14.f};
        const float lBinSize{0.5f * radiationLength};
        const int nBinsL{static_cast<int>(std::ceil(lBinMax / lBinSize))};
        FloatVector binnedLongEnergy(nBinsL);
        std::cout << "N: " << nBinsL << std::endl;
        for (const auto point : lProjEnergyVec)
        {
            const int bin{static_cast<int>(point.first / lBinSize)};
            binnedLongEnergy[bin] += point.second;
        }

        FloatVector binnedMuonEnergy(nBinsL);
        FloatVector binnedPionEnergy(nBinsL);
        FloatVector binnedProtonEnergy(nBinsL);
        float residualEnergy{totalEnergy};
        for (int i = 0; i < nBinsL; ++i)
        {
            std::cout << "Residual: " << residualEnergy << std::endl;
            if (residualEnergy > 0.f)
            {
                binnedMuonEnergy[i] = bethe_bloch(residualEnergy, 105.658f) * lBinSize;
                binnedPionEnergy[i] = bethe_bloch(residualEnergy, 139.570f) * lBinSize;
                binnedProtonEnergy[i] = bethe_bloch(residualEnergy, 938.272f) * lBinSize;
            }
            else
            {
                binnedMuonEnergy[i] = 0.f;
                binnedPionEnergy[i] = 0.f;
                binnedProtonEnergy[i] = 0.f;
            }

            residualEnergy -= (binnedLongEnergy[i] * 1e3);
            std::cout << i << ": " << (binnedLongEnergy[i] * 1e3) << " dedx(mu): " << binnedMuonEnergy[i] <<
                " dedx(pi): " << binnedPionEnergy[i] << " dedx(p): " << binnedProtonEnergy[i] << std::endl;
        }
        std::cout << std::endl;

        // Sum up energy, get theory profiles for e and gamma (and tracks too)
        // Compare, also generally assess the shape etc to look for problems
        // Do an r-phi based on the vertex and see how many peaks there are
        // Can child particles be folded into parent?

        PANDORA_MONITORING_API(VisualizeCaloHits(this->GetPandora(), &caloHitList, "PFO", AUTOITER));
        PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ParticleConsistencyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "InputPfoListName", m_inputPfoListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

