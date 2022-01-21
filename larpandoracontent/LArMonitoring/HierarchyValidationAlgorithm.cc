/**
 *  @file   larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.cc
 *
 *  @brief  Implementation of the particle visualisation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArMonitoring/HierarchyValidationAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"

using namespace pandora;

namespace lar_content
{

HierarchyValidationAlgorithm::HierarchyValidationAlgorithm() :
    m_event{-1},
    m_writeTree{false},
    m_foldToPrimaries{false},
    m_foldDynamic{false},
    m_foldToLeadingShowers{false},
    m_validateEvent{false},
    m_validateMC{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

HierarchyValidationAlgorithm::~HierarchyValidationAlgorithm()
{
    if (m_writeTree)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treename.c_str(), m_filename.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyValidationAlgorithm::Run()
{
    ++m_event;
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    LArHierarchyHelper::FoldingParameters foldParameters;
    if (m_foldToPrimaries)
        foldParameters.m_foldToTier = true;
    else if (m_foldDynamic)
        foldParameters.m_foldDynamic = true;
    else if (m_foldToLeadingShowers)
        foldParameters.m_foldToLeadingShowers = true;
    LArHierarchyHelper::MCHierarchy mcHierarchy;
    LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, foldParameters, mcHierarchy);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FillRecoHierarchy(*pPfoList, foldParameters, recoHierarchy);
    LArHierarchyHelper::MatchInfo matchInfo;
    LArHierarchyHelper::MatchHierarchies(mcHierarchy, recoHierarchy, matchInfo);
    //matchInfo.Print(mcHierarchy);

    if (m_validateEvent)
        this->EventValidation(matchInfo);
    else if (m_validateMC)
        this->MCValidation(matchInfo);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::EventValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeTree)
    {
        const LArHierarchyHelper::MCMatchesVector &matches{matchInfo.GetMatches()};
        MCParticleSet primaryMCSet;
        std::set<const LArHierarchyHelper::MCHierarchy::Node *> trackNodeSet, showerNodeSet;
        int nGoodTrackMatches{0}, nGoodShowerMatches{0};
        int nGoodMatches{0}, nPoorMatches{0}, nUnmatched{0};
        int nGoodTier1Matches{0}, nTier1Nodes{0};
        int nGoodTier1TrackMatches{0}, nTier1TrackNodes{0};
        int nGoodTier1ShowerMatches{0}, nTier1ShowerNodes{0};
        int hasLeadingMuon{0}, hasLeadingElectron{0}, isLeadingLeptonCorrect{0};
        // ATTN: Probably want quality cuts here for "good match" definition
        for (const LArHierarchyHelper::MCMatches &mcMatch : matches)
        {
            const LArHierarchyHelper::MCHierarchy::Node *pNode{mcMatch.GetMC()};
            const MCParticle *const pMC{pNode->GetLeadingMCParticle()};
            primaryMCSet.insert(LArMCParticleHelper::GetPrimaryMCParticle(pMC));
            const int nReco{static_cast<int>(mcMatch.GetRecoMatches().size())};
            const bool isQuality{mcMatch.IsQuality(matchInfo.GetQualityCuts())};
            if (isQuality)
                ++nGoodMatches;
            else if (nReco >= 1)
                ++nPoorMatches;
            else
                ++nUnmatched;
            if (pNode->GetHierarchyTier() == 1)
            {
                ++nTier1Nodes;
                if (isQuality)
                    ++nGoodTier1Matches;
            }

            const int pdg{std::abs(pNode->GetParticleId())};
            if (pNode->IsLeadingLepton())
            {
                if (pdg == MU_MINUS)
                    hasLeadingMuon = 1;
                else if (pdg == E_MINUS)
                    hasLeadingElectron = 1;
                isLeadingLeptonCorrect = nReco == 1 ? 1 : 0;
            }

            if (pdg == PHOTON || pdg == E_MINUS)
            {
                showerNodeSet.insert(pNode);
                if (isQuality)
                {
                    ++nGoodShowerMatches;
                    if (pNode->GetHierarchyTier() == 1)
                        ++nGoodTier1ShowerMatches;
                }
                if (pNode->GetHierarchyTier() == 1)
                    ++nTier1ShowerNodes;
            }
            else
            {
                trackNodeSet.insert(pNode);
                if (isQuality)
                {
                    ++nGoodTrackMatches;
                    if (pNode->GetHierarchyTier() == 1)
                        ++nGoodTier1TrackMatches;
                }
                if (pNode->GetHierarchyTier() == 1)
                    ++nTier1TrackNodes;
            }
        }

        MCParticleList primaryMCList;
        for (const MCParticle *const pMC : primaryMCSet)
            primaryMCList.emplace_back(pMC);
        const int interactionType{static_cast<int>(LArInteractionTypeHelper::GetInteractionType(primaryMCList))};
        const int nNodes{static_cast<int>(matchInfo.GetNMCNodes())};

        const int nTrackNodes{static_cast<int>(trackNodeSet.size())}, nShowerNodes{static_cast<int>(showerNodeSet.size())};
        const CartesianVector &trueVertex{matchInfo.GetMCNeutrino()->GetVertex()};
        const CartesianVector &recoVertex{LArPfoHelper::GetVertex(matchInfo.GetRecoNeutrino())->GetPosition()};
        const float vtxDx{recoVertex.GetX() - trueVertex.GetX()};
        const float vtxDy{recoVertex.GetY() - trueVertex.GetY()};
        const float vtxDz{recoVertex.GetZ() - trueVertex.GetZ()};
        const float vtxDr{std::sqrt(vtxDx * vtxDx + vtxDy * vtxDy + vtxDz * vtxDz)};

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "interactionType", interactionType));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodMatches", nGoodMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nPoorMatches", nPoorMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nUnmatched", nUnmatched));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nNodes", nNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTier1Matches", nGoodTier1Matches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTier1Nodes", nTier1Nodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTrackMatches", nGoodTrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodShowerMatches", nGoodShowerMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTrackNodes", nTrackNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nShowerNodes", nShowerNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTier1TrackMatches", nGoodTier1TrackMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTier1TrackNodes", nTier1TrackNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nGoodTier1ShowerMatches", nGoodTier1ShowerMatches));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nTier1ShowerNodes", nTier1ShowerNodes));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hasLeadingMuon", hasLeadingMuon));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "hasLeadingElectron", hasLeadingElectron));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLeptonCorrect", isLeadingLeptonCorrect));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
        PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::MCValidation(const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    if (m_writeTree)
    {
        for (const LArHierarchyHelper::MCMatches &match : matchInfo.GetMatches())
            this->Fill(match, matchInfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HierarchyValidationAlgorithm::Fill(const LArHierarchyHelper::MCMatches &matches, const LArHierarchyHelper::MatchInfo &matchInfo) const
{
    const LArHierarchyHelper::MCHierarchy::Node *pMCNode{matches.GetMC()};
    const int isTestBeam{pMCNode->IsTestBeamParticle() ? 1 : 0};
    const int isCosmicRay{!isTestBeam && pMCNode->IsCosmicRay() ? 1 : 0};
    const int isNeutrinoInt{!(isTestBeam || isCosmicRay) ? 1 : 0};
    const int mcId{pMCNode->GetId()};
    const int pdg{pMCNode->GetParticleId()};
    const int tier{pMCNode->GetHierarchyTier()};
    const int mcHits{static_cast<int>(pMCNode->GetCaloHits().size())};
    const int isLeadingLepton{pMCNode->IsLeadingLepton() ? 1 : 0};
    const int isQuality{matches.IsQuality(matchInfo.GetQualityCuts()) ? 1: 0};

    const MCParticleList &parentList{pMCNode->GetLeadingMCParticle()->GetParentList()};
    const int isElectron{std::abs(pMCNode->GetLeadingMCParticle()->GetParticleId()) == E_MINUS ? 1 : 0};
    const int hasMuonParent{parentList.size() == 1 && std::abs(parentList.front()->GetParticleId()) == MU_MINUS ? 1 : 0};
    const int isMichel{isElectron && hasMuonParent && LArMCParticleHelper::IsDecay(pMCNode->GetLeadingMCParticle()) ? 1 : 0};

    const LArHierarchyHelper::RecoHierarchy::NodeVector &nodeVector{matches.GetRecoMatches()};
    const int nMatches{static_cast<int>(nodeVector.size())};
    IntVector recoIdVector, nRecoHitsVector, nSharedHitsVector;
    float purity{0.f}, completeness{0.f};
    float purityAdc{0.f}, completenessAdc{0.f};
    float purityU{0.f}, purityV{0.f}, purityW{0.f}, completenessU{0.f}, completenessV{0.f}, completenessW{0.f};
    float purityAdcU{0.f}, purityAdcV{0.f}, purityAdcW{0.f}, completenessAdcU{0.f}, completenessAdcV{0.f}, completenessAdcW{0.f};
    const CartesianVector &trueVertex{pMCNode->GetLeadingMCParticle()->GetVertex()};
    float vtxDx{0.f}, vtxDy{0.f}, vtxDz{0.f}, vtxDr{0.f};
    // Identify the best matching reco node (most complete) and use that to determine per particle metrics
    const LArHierarchyHelper::RecoHierarchy::Node *pBestRecoNode{nullptr};
    float maxCompleteness{0.f};
    for (const LArHierarchyHelper::RecoHierarchy::Node *pRecoNode : nodeVector)
    {
        recoIdVector.emplace_back(pRecoNode->GetParticleId());
        nRecoHitsVector.emplace_back(static_cast<int>(pRecoNode->GetCaloHits().size()));
        nSharedHitsVector.emplace_back(static_cast<int>(matches.GetSharedHits(pRecoNode)));

        float recoCompleteness{matches.GetCompleteness(pRecoNode)};
        if (recoCompleteness > maxCompleteness)
        {
            pBestRecoNode = pRecoNode;
            maxCompleteness = recoCompleteness;
        }
    }
    if (pBestRecoNode)
    {
        purity = matches.GetPurity(pBestRecoNode);
        completeness = matches.GetCompleteness(pBestRecoNode);
        purityAdc = matches.GetPurity(pBestRecoNode, true);
        completenessAdc = matches.GetCompleteness(pBestRecoNode, true);
        purityU = matches.GetPurity(pBestRecoNode, TPC_VIEW_U);
        purityV = matches.GetPurity(pBestRecoNode, TPC_VIEW_V);
        purityW = matches.GetPurity(pBestRecoNode, TPC_VIEW_W);
        completenessU = matches.GetCompleteness(pBestRecoNode, TPC_VIEW_U);
        completenessV = matches.GetCompleteness(pBestRecoNode, TPC_VIEW_V);
        completenessW = matches.GetCompleteness(pBestRecoNode, TPC_VIEW_W);
        purityAdcU = matches.GetPurity(pBestRecoNode, TPC_VIEW_U, true);
        purityAdcV = matches.GetPurity(pBestRecoNode, TPC_VIEW_V, true);
        purityAdcW = matches.GetPurity(pBestRecoNode, TPC_VIEW_W, true);
        completenessAdcU = matches.GetCompleteness(pBestRecoNode, TPC_VIEW_U, true);
        completenessAdcV = matches.GetCompleteness(pBestRecoNode, TPC_VIEW_V, true);
        completenessAdcW = matches.GetCompleteness(pBestRecoNode, TPC_VIEW_W, true);

        // Temporary - need to identify the PFO vertex for the node - actually slightly tricky
        const CartesianVector recoVertex(0.f, 0.f, 0.f);
        //const CartesianVector &recoVertex{LArPfoHelper::GetVertex(pBestRecoNode)->GetPosition()};
        vtxDx = recoVertex.GetX() - trueVertex.GetX();
        vtxDy = recoVertex.GetY() - trueVertex.GetY();
        vtxDz = recoVertex.GetZ() - trueVertex.GetZ();
        vtxDr = std::sqrt(vtxDx * vtxDx + vtxDy * vtxDy + vtxDz * vtxDz);
    }

    // Would like to add information on hierarchy matching. Needs some thought, it's extremely complicated

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "event", m_event));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcId", mcId));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcPDG", pdg));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcTier", tier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "mcNHits", mcHits));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isNuInteration", isNeutrinoInt));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isCosmicRay", isCosmicRay));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isTestBeam", isTestBeam));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isLeadingLepton", isLeadingLepton));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isMichel", isMichel));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nMatches", nMatches));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "isQuality", isQuality));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "recoIdVector", &recoIdVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nRecoHitsVector", &nRecoHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "nSharedHitsVector", &nSharedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purity", purity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completeness", completeness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdc", purityAdc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdc", completenessAdc));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityU", purityU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityV", purityV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityW", purityW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessU", completenessU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessV", completenessV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessW", completenessW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcU", purityAdcU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcV", purityAdcV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "purityAdcW", purityAdcW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcU", completenessAdcU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcV", completenessAdcV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "completenessAdcW", completenessAdcW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDx", vtxDx));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDy", vtxDy));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDz", vtxDz));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treename.c_str(), "vtxDr", vtxDr));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treename.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HierarchyValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    if (m_pfoListName.empty())
        m_pfoListName = "RecreatedPfos";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidateEvent", m_validateEvent));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidateMC", m_validateMC));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteTree", m_writeTree));
    if (m_writeTree)
    {
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_filename));
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treename));
        if (!(m_validateEvent || m_validateMC))
        {
            std::cout << "Error: WriteTree requested but no tree names found" << std::endl;
            return STATUS_CODE_NOT_FOUND;
        }
        else if (m_validateEvent && m_validateMC)
        {
            std::cout << "Error: Both event-level and MC-level validation requested simulataneously" << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToPrimaries", m_foldToPrimaries));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldDynamic", m_foldDynamic));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FoldToLeadingShowers", m_foldToLeadingShowers));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
