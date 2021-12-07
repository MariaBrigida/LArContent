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

void get_longitudinal_energy_profile(const float e0, const FloatVector &positions, const float cj, const float binSize,
    FloatVector &fractionalEnergies)
{
    const unsigned int N{static_cast<unsigned int>(positions.size())};
    const float eCrit{32.84f};
    const float b{0.5f};
    float a{0.f};
    if (e0 > eCrit)
        a = (1 + b * cj) + b * std::log(e0 / eCrit);
    else
        a = 1 + b * cj;
    float gammaA{std::tgamma(a)};
    float cumulativeEnergy{0.f};
    for (unsigned int i = 0; i < N; ++i)
    {
        const float t{positions[i]};
        fractionalEnergies[i] = b * std::pow(b * t * binSize, a - 1) * std::exp(-b * t) / gammaA;
        cumulativeEnergy += fractionalEnergies[i];
    }
    if (cumulativeEnergy > std::numeric_limits<float>::epsilon())
    {
        for (unsigned int i = 0; i < N; ++i)
            fractionalEnergies[i] /= cumulativeEnergy;
    }
}

void get_photon_longitudinal_energy_profile(const float e0, const FloatVector &positions, const float binSize,
    FloatVector &fractionalEnergies)
{
    get_longitudinal_energy_profile(e0, positions, +0.5f, binSize, fractionalEnergies);
}

void get_electron_longitudinal_energy_profile(const float e0, const FloatVector &positions, const float binSize,
    FloatVector &fractionalEnergies)
{
    get_longitudinal_energy_profile(e0, positions, -0.5f, binSize, fractionalEnergies);
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
            tProjEnergyVec.emplace_back(std::make_pair(tProj, 1000.f * pCaloHit->GetElectromagneticEnergy()));
        }
        // Determine which end of the longitudinal projection is closest to the vertex
        float start{0.f};
        if (std::abs(lMin - lVertex) <= std::abs(lMax - lVertex))
            start = lMin;
        else
            start = lMax;
        // Adjust projection so that the extremal hit closest to the vertex is longitudinally projected onto zero
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const CartesianVector &pos{pCaloHit->GetPositionVector()};
            const float lProj{std::abs(pos.GetDotProduct(lAxis) - start)};
            const float hitEnergy{pCaloHit->GetElectromagneticEnergy()};
            lProjEnergyVec.emplace_back(std::make_pair(lProj, 1000.f * hitEnergy));
        }
        const float lBinMax{std::max(std::abs(lMax - start), std::abs(lMin - start)) + std::numeric_limits<float>::epsilon()};
        const float radiationLength{14.f};
        const float lBinSize{0.5f * radiationLength};
        const int nBinsL{static_cast<int>(std::ceil(lBinMax / lBinSize))};

        if (nBinsL <= 1)
            continue;

        FloatVector binnedLongEnergy(nBinsL);
        std::cout << "N: " << nBinsL << std::endl;
        for (const auto point : lProjEnergyVec)
        {
            const int bin{static_cast<int>(point.first / lBinSize)};
            binnedLongEnergy[bin] += point.second;
        }

        // Get theoretical shower profiles
        const float binSize{0.5f};
        FloatVector binCentres(nBinsL);
        for (int i = 0; i < nBinsL; ++i)
            binCentres[i] = 0.5f * (2 * i + 1) * binSize;
        const float e0{std::accumulate(binnedLongEnergy.begin(), binnedLongEnergy.end(), 0.f)};

//        for (int i = 0; i < nBinsL; ++i)
//            binnedLongEnergy[i] /= e0;

        FloatVector electronProfile(nBinsL), photonProfile(nBinsL);
        get_electron_longitudinal_energy_profile(e0, binCentres, binSize, electronProfile);
        get_photon_longitudinal_energy_profile(e0, binCentres, binSize, photonProfile);

        for (int i = 0; i < nBinsL; ++i)
        {
            electronProfile[i] *= e0;
            photonProfile[i] *= e0;
        }

        // Assess profiles
        float photonChi2{0.f}, electronChi2{0.f};
        // Only look at the first 6 bins
        const int nShowerBins{std::min(nBinsL, 6)};
        for (int i = 0; i < nShowerBins; ++i)
        {
            float dE{binnedLongEnergy[i] - electronProfile[i]};
            electronChi2 += dE * dE / electronProfile[i];
            dE = binnedLongEnergy[i] - photonProfile[i];
            photonChi2 += dE * dE / photonProfile[i];

            std::cout << i << " - Actual: " << binnedLongEnergy[i] << " " << electronProfile[i] << " " << photonProfile[i] << std::endl;
        }
        electronChi2 /= nShowerBins - 1;
        photonChi2 /= nShowerBins - 1;

        std::cout << "Chi2(e): " << electronChi2 << " Chi2(gamma): " << photonChi2 << std::endl;

        // Look for approx constant energy deposition prior to a potential Bragg peak (essentially the last bin)
        const int nConstantBins{nBinsL - 1};
        FloatVector gradient(nConstantBins - 1);
        int drop{0};
        for (int i = 0; i < nConstantBins - 1; ++i)
        {
            const float scale{std::max(std::abs(binnedLongEnergy[i + 1]), std::abs(binnedLongEnergy[i]))};
            if (scale < std::numeric_limits<float>::epsilon())
            {
                // Both bins are zero, drop the bin from the calculation
                ++drop;
                continue;
            }
            gradient[i] = (binnedLongEnergy[i + 1] - binnedLongEnergy[i]) / (scale * binSize);
            gradient[i] *= gradient[i];
            std::cout << i << " - Grad: " << gradient[i] << std::endl;
        }
        const float flatness{drop < (nConstantBins - 1) ? std::sqrt(std::accumulate(gradient.begin(), gradient.end(), 0.f) / (nConstantBins - (drop + 1))) : 0.f};
        std::cout << "Flatness: " << flatness << std::endl;

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

