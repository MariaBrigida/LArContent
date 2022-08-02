/**
 *  @file   larpandoracontent/LArCheating/MyCheatingVertexSelectionAlgorithm.h
 *
 *  @brief  Header file for the cheating vertex selection algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MY_CHEATING_VERTEX_SELECTION_ALGORITHM_H
#define LAR_MY_CHEATING_VERTEX_SELECTION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  MyCheatingVertexSelectionAlgorithm class
 */
class MyCheatingVertexSelectionAlgorithm : public pandora::Algorithm
{
public:

    MyCheatingVertexSelectionAlgorithm();

    //virtual ~MyCheatingVertexSelectionAlgorithm();

private:
//    void GetVertexScoreList(const pandora::VertexVector &vertexVector, const BeamConstants &beamConstants, HitKDTree2D &kdTreeU,
//        HitKDTree2D &kdTreeV, HitKDTree2D &kdTreeW, VertexScoreList &vertexScoreList) const;

    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_candidateVertexListName;
    std::string m_outputVertexListName;
    bool m_replaceCurrentVertexList;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_VERTEX_SELECTION_ALGORITHM_H
