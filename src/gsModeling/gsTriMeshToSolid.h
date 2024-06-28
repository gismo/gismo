/** @file gsTriMeshToSolid.h

    @brief Provides declaration of gsTriMeshToSolid: a triangle mesh to
    gsSolid convertor class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mayer, M. Pauley
*/

#pragma once

#include <gsUtils/gsMesh/gsMesh.h>

namespace gismo
{

/**
    @brief Class gsTriMeshToSolid object.

    Construct an instance of this class from a gsMesh, and then use it
    to perform "CAD model reconstruction" to produce a gsSolid.

    \ingroup Modeling
*/
template<class T>
class gsTriMeshToSolid
{
public:
    // typedefs
    typedef gsMeshElement<T> MeshElement;
    typedef typename MeshElement::scalar_t scalar_t;
    typedef typename MeshElement::gsVertexHandle VertexHandle;
    typedef typename MeshElement::gsFaceHandle FaceHandle;
    typedef typename MeshElement::gsEdgeHandle EdgeHandle;
    typedef gsEdge<T>                         Edge;
    typedef gsVertex<T>                       Vertex;

    // constructors
    gsTriMeshToSolid(gsMesh<T> *sourceMesh):
        mesh(sourceMesh),
        vertex(sourceMesh->m_vertex),
        face(sourceMesh->m_face),
        edge(sourceMesh->m_edge)
    {}

    // public functions

    /** \brief Computes data describing the patch structure of the mesh. Combines getFeatures and getFaces.
     *
     * \param angle - defines the angle between 2 triangles before the
     *                edge between them counts as sharp. By increasing
     *                this parameter, the function will identify a
     *                smaller number of patches of the mesh.
     *
     * \param innerAngle - after the initial patch generation there is
     *                     the possibility to subdivide the larger
     *                     patches by a smaller angle.  Larger is
     *                     defined by patchAreaWeight. \e innerAngle
     *                     should be smaller than angle, else it will
     *                     not have any effect.
     *
     * \param patchAreaWeight - patchAreaWeight determines how large a
     * patch has to be such that it is subdivided by using the
     * parameter \e innerAngle.
     *
     * \param mergeSmallPatches - by using a low angle a lot of small
     *                            patches might be produced. \e
     *                            mergeSmallPatches gives the
     *                            opportunity to merge the small
     *                            patches. \e mergeSmallPatches
     *                            decides how small the patches are
     *                            allowed to be before they are
     *                            merged.  A value of 1 will merge
     *                            patches if the area of them is
     *                            smaller than the average area.
     *
     * \param[out] iPoints - the data about the inner points of a
     *                       patch is stored in an unsorted
     *                       vector. For each patch such a vector is
     *                       generated.  The parameter \e iPoints is a
     *                       vector of these generated vectors.
     *
     * \param[out] oPoints - the data about the boundary points of a
     *                       patch is stored in a counter clockwise
     *                       sorted vector. For each patch such a
     *                       vector is generated.  The parameter \e
     *                       oPoints is a vector of these generated
     *                       vectors.
     *
     * \param[out] innerBdrys - the data about each hole of a patch is
     *                          stored in a clockwise sorted
     *                          vector. For each hole of a patch such
     *                          a vector is generated.  For each patch
     *                          a vector of its holes is
     *                          generated. The parameter \e innerBdrys
     *                          is a vector of these generated
     *                          vectors.
     *
     * \param[out] innerBdrysMassP - the data about the mass points of
     *                               each hole is stored here.
     *
     * \param[out] oPointsConvexFlag - each boundary point can be
     * convex of not, depending on the turning angles of the adjacent
     * edges on the boundary.
     *
     * \param filenameFeatures - the path to the text file where the
     *                           manually provided features (sharp edges) can be found.
     *
     * \param useFeatures - if set to 1 the text file is used to add
     *                      the features from the text file. If set to 2 these features are
     *                      used exclusively.
     *
     * \ingroup Modeling
     */
    void getPatchData(T angle, T innerAngle, T patchAreaWeight, T mergeSmallPatches,
                      std::vector<std::vector<VertexHandle> > & iPoints,
                      std::vector<std::vector<VertexHandle> > & oPoints,
                      std::vector< std::vector<std::vector<VertexHandle> > > & innerBdrys,
                      std::vector< std::vector<Vertex>  > & innerBdrysMassP,
                      std::vector<std::vector<bool> > & oPointsConvexFlag,
                      std::string filenameFeatures,
                      int useFeatures);


    /** \brief generates edges for a mesh consisting of vertices and faces. Determines if these edges are sharp,
     *         depending on the angle between the two adjacent faces of the edges.
     *  \param angleGrad - if the angle (degree) between the two adjacent faces of an edge is bigger than angleGrad the edge is marked as sharp.
     *  \param[out] bWarnNonManifold - set to one if an edge is detected with more than 2 adjacent faces.
     *  \param[out] bWarnBorders - set to one if an edge is detected with 1 adjacent face.
     */
    void getFeatures(T angleGrad,bool& bWarnNonManifold,bool& bWarnBorders);

    /** \brief Each face obtains a patch number, faces of the same number belong to the same patch.
     */
    void calcPatchNumbers();

    /** \brief Sets sharp edges according to the value of useFeatures
    * \param featEdges - a vector of edges which can be used to manually add features.
    * \param useFeatures - if set to 1 the text file is used to add the features from the text file. If set to 2 these features are used exclusively.
    */
    void setSharpEdges(std::vector< gsEdge<T> > & featEdges, int useFeatures);


    /** \brief Store each edge's neighboring faces **/
    void storeNeighboringFaces();

    /** \brief Improve surface segmentation by using more complex rules
     * \param innerAngle - after the initial patch generation there is the possibility to subdivide the larger patches by a smaller angle.
     *                     Larger is defined by patchAreaWeight. \e innerAngle should be smaller than angle, else it will not have
     *                     any effect.
     * \param patchAreaWeight - patchAreaWeight determines how large a patch has to be such that it is subdivided by using the parameter \e innerAngle.
     * \param mergeSmallPatches - by using a low angle a lot of small patches might be produced. \e mergeSmallPatches gives the opportunity
     *                            to merge the small patches. \e mergeSmallPatches decides how \"small\" the patches are allowed to be before they are merged.
     *                            A value of 1 will merge patches if the area of them is smaller than the average area. A value of 0 turns this option off.
     */
    void divideAndMergePatches(T innerAngle, T patchAreaWeight, T mergeSmallPatches);

    /** \brief Computes data describing the patch structure of the mesh.
     * \param[out] iPoints - the data about the inner points of a patch is stored in an unsorted vector. For each patch such a vector is generated.
     *                       The parameter \e iPoints is a vector of these generated vectors.
     * \param[out] oPoints - the data about the boundary points of a patch is stored in a counter clockwise sorted vector. For each patch such a vector is generated.
     *                       The parameter \e oPoints is a vector of these generated vectors.
     * \param[out] innerBdrys - the data about each hole of a patch is stored in a clockwise sorted vector. For each hole of a patch such a vector is generated.
     *                          For each patch a vector of its holes is generated. The parameter \e innerBdrys is a vector of these generated vectors.
     * \param[out] innerBdrysMassP - the data about the mass points of each hole is stored here.
     * \param[out] oPointsConvexFlag - each boundary point can be convex of not, depending on the turning angles of the adjacent edges on the boundary.
    */
    void getFaces(std::vector<std::vector<VertexHandle> > & iPoints, std::vector<std::vector<VertexHandle> > & oPoints,
                  std::vector< std::vector<std::vector<VertexHandle> > > & innerBdrys, std::vector< std::vector<Vertex>  > & innerBdrysMassP,
                  std::vector<std::vector<bool> > & oPointsConvexFlag);


    /** \brief Parametrized a number of patches given by \e iPoints, \e oPoints, \e innerBdrys, \e innerBdrysMassP and \e oPointsConvexFlag,
     *         fits B-spline surfaces to them and trimmes the resulting surfaces. The trimmed surfaces are then added to the empty solid \e sl
     *         and a half-edge structure for \e sl is generated. If desired, in paraview plotable meshes to visualize the solid are generated.
     * \param [out] sl - The empty solid \e sl is designed by adding trimmed surface patches and incorporating a half-edge structure.
     * \param iPoints - the data about the inner points of a patch is stored in an unsorted vector. For each patch such a vector is generated.
     *                       The parameter \e iPoints is a vector of these generated vectors.
     * \param oPoints - the data about the boundary points of a patch is stored in a counter clockwise sorted vector. For each patch such a vector is generated.
     *                       The parameter \e oPoints is a vector of these generated vectors.
     * \param innerBdrys - the data about each hole of a patch is stored in a clockwise sorted vector. For each hole of a patch such a vector is generated.
     *                          For each patch a vector of its holes is generated. The parameter \e innerBdrys is a vector of these generated vectors.
     * \param innerBdrysMassP - the data about the mass points of each hole is stored here.
     * \param oPointsConvexFlag - each boundary point can be convex of not, depending on the turning angles of the adjacent edges on the boundary.
     * \param [out] paraMeshes - a vector of meshes visualizing the parametrizations of the patches.
     * \param [out] fitMeshes - a vector of meshes visualizing the trimmed surfaces of the patches.
     * \param [out] patchMeshes - a vector of meshes visualizing the patches.
     * \param kvOuterPoints - determines the multiplicity at the two end knots of the knot vector needed to approximate the surface
     *                        of each face.
     * \param kvAdditionalInnerPoints - determines the number of additionl inner Points in the knot vector (equally distributed). The
     *                                  sqrt of the number of corners of a face are taken as a base value.
     * \param plot - if set to 1, it will generate meshes of each of the patches, its parametrization and also of the trimmed surfaces.
     * \param meshPoints - the accuracy of the visualization for the _fit plots.
     * \param moreInner - if additional interior points shall be added between large Edges. The number of points added per Edge are
     *                   calculated by length(edge)/h. h is the length of one of the smaller Edges. (The Edge such that 10% of the
     *                   other Edges are smaller.)
     * \param wE - weights for edge points during fitting the B-spline surfaces.
     * \param wI - weights for interior points during fitting the B-spline surfaces.
     * \param closeBoundary - if additional points on the edges shall be added, to reduce the possible gaps between two patches.
     * \param noSmooth - if the trimming curves shall not be smoothed out.
     *
     * \ingroup Modeling
     */
    void toSolid(gsSolid<T> & sl, std::vector<std::vector<VertexHandle> > & iPoints,
                 std::vector<std::vector<VertexHandle> > & oPoints,
                 std::vector< std::vector<std::vector<VertexHandle> > > & innerBdrys,
                 std::vector< std::vector<Vertex>  > & innerBdrysMassP,
                 std::vector<std::vector<bool> > & oPointsConvexFlag, std::vector<gsMesh<T> *> &paraMeshes,
                 std::vector<gsMesh<T> *> &fitMeshes, std::vector<gsMesh<T> *> &patchMeshes, int kvOuterPoints,
                 int kvAdditionalInnerPoints,
                 bool plot, int meshPoints, bool moreInner=true,
                 T wE=5, T wI=1, int closeBoundary=0, bool noSmooth=false);

    /** \brief calculates a curve loop consisting of B-spline curves of degree one and two from a vector of vertices and the information about the convexity of these vertices.
     *  \param outerPoints - the vertices representing the boundary.
     *  \param isCorner - information about the convexity of the vertices.
     *  \param noSmooth - if set to 1, only curves of degree one will be used for the curve loop.
     *  \return the resuting curve loop.
     */
    gsCurveLoop<T> * calculateLoop(std::vector<Vertex> outerPoints, std::vector<bool > const & isCorner, bool noSmooth = false);


    /// Check if a given normal vector is consistent with a collection of
    /// triangles. Return +1 if the inner product of globalNormal with each
    /// triangle's normal is positive. Return -1 if all the inner products are
    /// negative. Return 0 if mixed.
    static int normalMult(gsVector3d<T> globalNormal,
                          std::vector<FaceHandle> & face, int bigFaceIdx);

    /** \brief using Heron's Formula for the area of a triangle.
     *  \param f1 - the face of  which the area shall be calculated.
     *  \return the area of the input face \e f1.
     */
    static T calcArea(FaceHandle f1);

    /** \brief calculates a weight between 2 vertices used in Floater's algorithm.
     *  \param v1 - first vertex.
     *  \param v2 - second vertex.
     *  \param vertexFaceSet - used to check if both vertices are in the same patch.
     *  \return the weight.
     */
    static T calcWeight(VertexHandle v1,VertexHandle v2,
                        std::set<VertexHandle>
                        const & vertexFaceSet);

    /** \brief checks if two edges are very close to each other.
     *  \param e1 - first edge.
     *  \param e2 - second edge.
     *
     *  \return true if the points of the edges are closer to each
     *  other than one percent of the length of the first edge.
     */
    static bool approxEqual(const gsEdge<T> & e1,const gsEdge<T> & e2);


    /** \brief calculates the conditioned angle between 2 edges.
     *  \param e1 - first edge.
     *  \param e2 - second edge.
     *  \param faceNum - number of the patch, used to compute the normal.
     *  \return the conditioned angle.
     */
    static T calcAngle(EdgeHandle e1,EdgeHandle e2, int faceNum);


    /** \brief Adds up the lengths between 2 neighboring vertices of a vector of vertices.
     *  \param vec - a vector of vertices describing the boundary.
     *  \return the length of the input vector \e vec.
     */
    static T calcBdryLength(std::vector<VertexHandle > vec);


    /** \brief calculates the mass point of a vector of vertices.
     *  \param vec - a vector of vertices from which the mass points is calculated.
     *  \return the mass point
     */
    static gsVertex<T> getMassP(std::vector<VertexHandle > vec);


    /** \brief calculates the distance between 2 vertices.
     *  \param v1 - first vertex.
     *  \param v2 - second vertex.
     *  \return the distance.
     */
    static T calcDist(VertexHandle v1,VertexHandle v2);


    // Documentation guessed by S.K.
    /** \brief checks whether vertices are corners.
     *  \param vertexVec3d vector of vertices.
     *  \return vector indicating whether the corresponding vertex is a corner or not.
     *
     *  \todo Verify that this description is correct.
     */
    static std::vector<bool> isCorner(std::vector<VertexHandle > const & vertexVec3d);


    /** \brief calculates a B-spline of degree one from 2 vertices.
     *  \param v1 - first vertex.
     *  \param v2 - second vertex.
     *  \return the B-spline.
     */
    gsBSpline<T> * calcTCurve(Vertex v1,Vertex v2);


    /** \brief calculates a B-spline of degree two from 3 vertices.
     *  \param v1 - first vertex.
     *  \param v2 - second vertex.
     *  \param v3 - third vertex.
     *  \return the B-spline.
     */
    gsBSpline<T> * calcTCurve(Vertex v1,Vertex v2,Vertex v3);


    /** \brief calculates the midpoint of two vertices.
     *  \param v1 - first vertex.
     *  \param v2 - second vertex.
     *  \return the Vertex in the middle of \e v1 and \e v2.
     */
    static Vertex giveMidpoint(Vertex v1,Vertex v2);


    /** \brief reads a text file consisting of lines of 6 values, each line representing an edge.
     *  \param fn - path to text file.
     *  \param[out] edges - vector of the resulting edges.
     */
    void readEdges( std::string const & fn, std::vector<gsEdge <T> > & edges );


    // public members
    gsMesh<> *mesh;

    int numEdges;
    int numBigFaces;
    std::vector<VertexHandle > &vertex;
    std::vector<FaceHandle > &face;
    gsSortedVector<Edge> &edge;
};

}


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTriMeshToSolid.hpp)
#endif
