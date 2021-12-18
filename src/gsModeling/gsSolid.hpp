/** @file gsSolid.hpp

    @brief Provides implementation of gsSolid class

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, D.-M. Nguyen, M. Pauley
*/

#pragma once

#include <queue>
#include <set>

#include <gsCore/gsMultiPatch.h>
#include <gsModeling/gsCurveLoop.h>
#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsPlanarDomain.h>

#include <gsNurbs/gsKnotVector.h>
//#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsTensorBSpline.h>

#include <gsNurbs/gsBSpline.h>



namespace gismo {  
  
template <class T>
std::ostream &gsSolid<T>::print(std::ostream &os) const 
{
  os<<"\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";       
  os<<"A free-form solid modeled as a gsSolid:\n"; 
  os<<"Number of volumes   : "<< numVolumes <<"\n";
  os<<"Number of half faces: "<< numHalfFaces <<"\n"; 
  os<<"Number of half edges: "<< numHalfEdges <<"\n"; 
  os<<"Number of vertices  : "<< numVertices <<"\n"; 
  // print a face
  os<< *face[0];
  os<<"vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n";       
  return os;       
}


template <class T>
void gsSolid<T>::addHeVertex(scalar_t const& x, scalar_t const& y, scalar_t const& z)
{
    gsSolidHeVertexHandle v = this->makeHeVertex(x,y,z);
    vertex.push_back(v );
    v->setId(numVertices++);   
    
    if(!initialized) {
	bb.pMax.x() = bb.pMin.x() = x;
	bb.pMax.y() = bb.pMin.y() = y;
	bb.pMax.z() = bb.pMin.z() = z;
	initialized = true;
    }	
    else {
	if (bb.pMax.x() < x)
	    bb.pMax.x() = x;
	if (bb.pMin.x() > x)
	    bb.pMin.x() = x;
	
	if (bb.pMax.y() < y)
	    bb.pMax.y() = y;
	if (bb.pMin.y() > y)
	    bb.pMin.y() = y;
	
	if (bb.pMax.z() < z)
	    bb.pMax.z() = z;
	if (bb.pMin.z() > z)
	    bb.pMin.z() = z;
    };  
}

template <class T>
typename gsSolid<T>::gsSolidHalfFaceHandle gsSolid<T>::addFace(std::vector<gsSolidHeVertexHandle> V)
{
    if (V.size() != 4) //{ gsDebug<<"\n The following initialization of trimmed surf does not work for less or more than 4 corners yet\n";};
    {
        return addFace_PlanarPolygon(V);
    }

    gsMatrix<T> corner = gsMatrix<T>(4,3);
    for (unsigned i=0; i<=3; i++)
    { 
      for (unsigned xi=0;xi<=2;xi++)
      {
        corner(i,xi) = (V[i]->coords)[xi];
      }
    }
    //gsTrimSurface<T>* tsurf = new  gsTrimSurface<T>(corner, 3, 3, 3);
    gsTrimSurface<T>* tsurf = new  gsTrimSurface<T>(corner, 1, 1, 1);
    return addFace(V,tsurf);
}

template<class T>
gsSolid<T>::~gsSolid()
{
    //apaga os vertices
    freeAll( vertex );

    //apaga as arestas
    freeAll( edge );
    
    //apaga as faces
    freeAll( face );
    
    //apaga as faces
    freeAll( volume );
}


template<class T>     
typename gsSolid<T>::gsSolidHalfFaceHandle gsSolid<T>::addFace(std::vector<gsSolidHeVertexHandle> V, gsTrimSurface<T> * tsurf)  
{
    gsSolidHalfFaceHandle f = new gsSolidHalfFace<T>();
    f->surf = tsurf;
    std::vector<gsSolidHalfEdgeHandle> H;
    for ( typename std::vector<gsSolidHeVertexHandle>::iterator 
	      it = V.begin(); it!= V.end(); ++it)
    {
	// initialize the half-edge associated with the vertex *it and the face f
	gsSolidHalfEdgeHandle tempHe = this->makeSolidHalfEdge(*it,f);
	tempHe->setId(numHalfEdges++);
	H.push_back( tempHe );
	edge.push_back( tempHe );
	// assign the incident halfedge associated with the vertex *it to *it
	// fist check if hed is already assigned to avoid doubling assignments
	if (!(*it)->hed)
	(*it)->hed = tempHe; 
	
    }          
    f->setBoundary( H );
    face.push_back(f);
    f->setId(numHalfFaces++);
    return f;
}

template<class T>
typename gsSolid<T>::gsSolidHalfFaceHandle gsSolid<T>::addFace(std::vector< std::vector<gsSolidHeVertexHandle> > loopV, gsTrimSurface<T> * tsurf)
{
    gsSolidHalfFaceHandle f = new gsSolidHalfFace<T>();
    f->surf = tsurf;
    std::vector<gsSolidHalfEdgeHandle> H;

    std::vector<gsSolidHeVertexHandle>  V;
    for ( size_t il=0; il<loopV.size(); il++)
    {
        V=loopV[il];
        H.clear();
        for ( typename std::vector<gsSolidHeVertexHandle>::iterator
              it = V.begin(); it!= V.end(); ++it)
        {
            // initialize the half-edge associated with the vertex *it and the face f
            gsSolidHalfEdgeHandle tempHe = this->makeSolidHalfEdge(*it,f);
            tempHe->setId(numHalfEdges++);
            H.push_back( tempHe );
            edge.push_back( tempHe );
            // assign the incident halfedge associated with the vertex *it to *it
            // fist check if hed is already assigned to avoid doubling assignments
            if (!(*it)->hed)
                (*it)->hed = tempHe;
        }
        f->setBoundary( H );
    }

    face.push_back(f);
    f->setId(numHalfFaces++);
    return f;
}

/// Add a planar face, specified by a list of vertices that should all be in the same plane.
template<class T>
typename gsSolid<T>::gsSolidHalfFaceHandle gsSolid<T>::addFace_PlanarPolygon(std::vector<gsSolidHeVertexHandle> V)
{
  // make sure there are at least 3 vertices.
  int faceVerts = V.size();
  assert(faceVerts >= 3);
  
  // compute normal to plane and make sure all vertices are in the plane.
  // (find a vertex with minimal x-coordinate, and minimal y in the case of a tie.
  // this ensures that the vertex is convex.)
  int minVert = faceVerts - 1;
  for(int i = 0; i < faceVerts - 1; i++)
  {
    if(V[i]->coords(0) < V[minVert]->coords(0) ||
      (V[i]->coords(0) <= V[minVert]->coords(0) && V[i]->coords(1) < V[minVert]->coords(1)))
      {
        //gsDebug << "\n" << V[i]->coords(0) << "," << V[i]->coords(1) << "," << V[i]->coords(2);
        minVert = i;
      }
  }
  gsVector3d<T> prevVector = V[(minVert + faceVerts - 1) % faceVerts]->coords - V[minVert]->coords;
  
  gsVector3d<T> normal(0, 0, 0);
  int nextVertex = minVert;
  while(normal.squaredNorm() == 0)
  {
    nextVertex = (nextVertex + 1) % faceVerts;
    assert(nextVertex != minVert); // failed to find a vertex that is not collinear with minVert and the one before it.
    normal = (V[nextVertex]->coords - V[minVert]->coords).cross(prevVector);
  }
  normal.normalize();
  for(int i = 3; i < faceVerts; i++)
  {
    assert((V[i]->coords - V[0]->coords).dot(normal) == 0);
  }
  T normalCoord = normal.dot(V[0]->coords);
  //gsDebug << "\nNormal coefs: " << normal << std::endl << normalCoord;
  
  // compute a good transformation that we can use to generate
  // corners in the 2D domain. We make an
  // arbitrary vector b1 orthogonal to normal, then cross product it with normal
  // to produce a second vector b2. {b1, b2} form a basis for the plane.
  // use the dot product of a vertex's coordinates with b1 and b2 to produce
  // 2D coordinates.
  
  gsVector3d<T> a1(1, 0, 0), a2(0, 1, 0);
  a1 = a1.cross(normal);
  a2 = a2.cross(normal);
  gsVector3d<T> b1, b2;
  if(a1.squaredNorm() > a2.squaredNorm())
  {
    b1 = a1;
  }
  else
  {
    b1 = a2;
  }
  b1.normalize();
  b2 = normal.cross(b1);
  
  //gsDebug << "\nbasis for plane. b1=" << b1 << ", b2="<< b2;
  
  // make a matrix containing a list of (transposed) corners in the domain.
  // compute a good rectangle for the domain of the base surface. do this by
  // finding the bounding box and enlarging it by a fixed amount (margin) on
  // each side.
  T margin(0.1);
  T x, y;
  gsMatrix<T> DomCor(faceVerts + 1, 2);
  for(int i = 0; i < faceVerts; i++)
  {
    // compute x and y coords and insert them into the matrix DomCor
    x = b1.dot(V[i]->coords);
    y = b2.dot(V[i]->coords);
    //gsDebug << "\nvertex " << i << " coords in domain: " << x << "," << y;
    DomCor(i, 0) = x;
    DomCor(i, 1) = y;
  }
  DomCor(faceVerts, 0) = DomCor(0, 0);
  DomCor(faceVerts, 1) = DomCor(0, 1);
  
  // find min and max x and y values (before resizing the polygon)
  T minx = DomCor.col(0).minCoeff();
  T maxx = DomCor.col(0).maxCoeff();
  T miny = DomCor.col(1).minCoeff();
  T maxy = DomCor.col(1).maxCoeff();
  T rangex = maxx - minx, rangey = maxy - miny;
  T midx = (minx + maxx) / 2, midy = (miny + maxy) / 2;
  T halfMaxRange = std::max(rangex, rangey) * (margin * 2 + 1) / 2;
  minx -= rangex * margin;
  maxx += rangex * margin;
  miny -= rangey * margin;
  maxy += rangey * margin;
  
  // resize and shift to be inside [0,1]x[0,1]
  gsCurveLoop<T>::adjustPolygonToUnitSquare(DomCor, margin);
  
  //gsDebug << "\nDomain rectangular bounds: [" << minx << "," << maxx << "] x [" << miny << "," << maxy << "]";
  
  // create a trimming loop and construct a planar domain from it.
  gsCurveLoop<T> * tloop = new gsCurveLoop<T>();
  for(int i = 0; i < faceVerts; i++)
  {
    gsMatrix<T> tcp(2, 2);
    tcp << DomCor(i, 0), DomCor(i, 1), DomCor(i + 1, 0), DomCor(i + 1, 1);
    gsBSpline<T> * tcurve = new gsBSpline<T>( 0, 1, 0, 1, give(tcp) );
    tloop->insertCurve(tcurve);
  }
  gsPlanarDomain<T> * domain1= new gsPlanarDomain<T>(tloop);
  
  // create the basis
  
  gsKnotVector<T> KV1 = gsKnotVector<T>(0, 1, 0, 2);//multiplicity at ends==1+1
  gsKnotVector<T> KV2 = gsKnotVector<T>(0, 1, 0, 2);//multiplicity at ends==1+1
  
  // unused variables:
  //gsBSplineBasis<T> * Bu= new gsBSplineBasis<T>(KV1);
  //gsBSplineBasis<T> * Bv= new gsBSplineBasis<T>(KV2);
  
  // produce a 4x3 matrix with the control points
  // to map the 2D domain linearly into 3D.
  gsMatrix<T> pcp(4, 3);
  for(int xi = 0; xi < 3; xi++) // loop over coordinates
  {
    pcp(0, xi) = normal[xi] * normalCoord + (midx - halfMaxRange) * b1[xi] + (midy - halfMaxRange) * b2[xi];
    pcp(1, xi) = normal[xi] * normalCoord + (midx + halfMaxRange) * b1[xi] + (midy - halfMaxRange) * b2[xi];
    pcp(2, xi) = normal[xi] * normalCoord + (midx - halfMaxRange) * b1[xi] + (midy + halfMaxRange) * b2[xi];
    pcp(3, xi) = normal[xi] * normalCoord + (midx + halfMaxRange) * b1[xi] + (midy + halfMaxRange) * b2[xi];
  }
  //gsDebug << "\nControl point matrix: " << pcp << "\n";
  //gsDebug << "\nPolygon points: " << DomCor << "\n";
  typename gsTensorBSpline<2,T>::Ptr tp1(new gsTensorBSpline<2,T>(KV1, KV2, give(pcp)));
  
  // instantiate and add the half-face object.
  
  gsTrimSurface<T> * ts = new gsTrimSurface<T>(tp1, domain1);
  return addFace(V, ts);
  
}

template<class T>
void gsSolid<T>::setHeMate()
{      
  // Method in the meantime: consider each pair of halfedges and see if they are mates
  // Todo: consider each pair of faces 
  unsigned int noMate(0); // number of mates 
  for (typename std::vector<gsSolidHalfEdgeHandle>::iterator it1=edge.begin(); it1!=edge.end()-1; ++it1)
  {
    for (typename std::vector<gsSolidHalfEdgeHandle>::iterator it2=it1+1; it2!=edge.end(); ++it2)
    {
      // a pair of half edges are mates iff the source of one of them is the target of the orther
      gsSolidHeVertexHandle source1, source2, target1, target2;
      source1 = (*it1)->source;
      source2 = (*it2)->source;
      target1 = (*it1)->next->source;
      target2 = (*it2)->next->source;
      if (source1 == target2 && source2 == target1)
      {
	noMate++;
	(*it1)->mate = *it2;
	(*it2)->mate = *it1;
      }
    }
  }	  
  // check if the number of mates is the same as the number of assignments
  if (2*noMate!=edge.size())
  {
//      gsDebug <<"\n"<<"The number of assignments of HE mates (="<< noMate <<") is NOT equal to number of edges (not halfedges) (="<< edge.size()/2 <<"), this is most likely because of the wrong order of the vertices of a face "<<"\n";

      gsWarn << "The number of assignments of HE mates (="<< noMate <<") is NOT equal to number of edges (not halfedges) (="<< edge.size()/2 <<"), this is most likely because of the wrong order of the vertices of a face, or the model is not a manifold\n";
  }
}

template<class T>
std::vector< typename gsSolid<T>::gsSolidHalfEdgeHandle > gsSolid<T>::detectNonConvexEdges(std::vector<int> const & ncEdgeV1, std::vector<int> const & ncEdgeV2)
{
  // Todo: detect nonconvex edges automatically
  size_t nEdge = ncEdgeV1.size();
  gsSolidHalfEdgeHandle he;
  assert(nEdge==ncEdgeV2.size());
  std::vector<gsSolidHalfEdgeHandle> nce;
  for (size_t i=0;i<nEdge;i++)
  {
      he = this->getVertexFromID(ncEdgeV1[i])->getHalfEdge( this->getVertexFromID(ncEdgeV2[i]) );
      he->is_convex = false;
      he->mate->is_convex = false;
      nce.push_back(he);
  }
  return nce;
}




template <class T>
typename gsSolid<T>::gsSolidHalfFaceHandle gsSolid<T>::splitFace(
  gsSolidHalfFaceHandle f,
  gsSolidHeVertexHandle startVertex, gsSolidHeVertexHandle endVertex, gsBSpline<T> *domainSpline)
{
  checkStructure();

  int numEdges = f->surf->domain().outer().size();
  // search for a half-edge on this face whose source is endVertex
  gsSolidHalfEdgeHandle nextEdge = endVertex->getHalfEdgeOnFace(f, false);
  gsSolidHalfEdgeHandle matesPrev = nextEdge->prev;
  int nextEdgeIdx = f->indexOfEdge(nextEdge);
  // search for a half-edge whose end vertex (next->source) is startVertex
  gsSolidHalfEdgeHandle prevEdge = startVertex->getHalfEdgeOnFace(f, true);
  gsSolidHalfEdgeHandle matesNext = prevEdge->next;
  int prevEdgeIdx = f->indexOfEdge(prevEdge);
  // create a new half-edge
  gsSolidHalfEdgeHandle newHE = this->makeSolidHalfEdge(startVertex, f);
  newHE->face = f;
  newHE->setId(numHalfEdges++);
  
  // close off the first part of this face off along the new half-edge
  newHE->prev = prevEdge;
  newHE->next = nextEdge;
  prevEdge->next = newHE;
  nextEdge->prev = newHE;
  // create the half-edge's mate
  gsSolidHalfEdgeHandle mate = this->makeSolidHalfEdge(endVertex, f);
  newHE->mate = mate;
  mate->mate = newHE;
  mate->setId(numHalfEdges++);
  
  // close off the second part
  mate->prev = matesPrev;
  mate->next = matesNext;
  matesPrev->next = mate;
  matesNext->prev = mate;
  // build a collection of edges for the new face
  std::vector<gsSolidHalfEdgeHandle> newFaceBoundary;
  newFaceBoundary.push_back(mate);
  gsSolidHalfEdgeHandle tempEdge;
  for(tempEdge = matesNext; tempEdge != mate; tempEdge = tempEdge->next)
  {
    newFaceBoundary.push_back(tempEdge);
  }

  // create a reverse spline of domainSpline
  std::vector<T> origKnots = domainSpline->knots();
  std::vector<T> mateKnots;
  mateKnots.reserve( origKnots.size() );
  T flipPoint = (*origKnots.begin()) + (*origKnots.rbegin());
  for(int i = origKnots.size() - 1; i >= 0; i--)
  {
    mateKnots.push_back(flipPoint - origKnots[i]);
  }
  gsKnotVector<T> mateKV(domainSpline->knots().degree(), mateKnots.begin(), mateKnots.end());
  gsMatrix<T> mateCoefs = domainSpline->coefs().colwise().reverse();
  gsBSpline<T> * reverseSpline = new gsBSpline<T>(mateKV, give(mateCoefs));
  
  // split the domain for the new face off from the old one
  typename gsPlanarDomain<T>::uPtr newDomain =
      f->surf->domain().split((prevEdgeIdx + 1) % numEdges, nextEdgeIdx, domainSpline, reverseSpline);
  if(prevEdgeIdx > nextEdgeIdx)
  {
    // if startVertex is the source of the first edge, or is later in the cycle than endVertex,
    // then the new face will contain loop[0]. The easiest thing to do is to reset loop[0] to
    // the new start of f's cycle.
    f->loop[0] = nextEdge;
  }
  typename gsSurface<T>::Ptr newBaseSurface = f->surf->getTP();
  
  // create the new face
  gsTrimSurface<T> * newTS = new gsTrimSurface<T>(newBaseSurface , newDomain.release());
  
  // add the new face using existing half-edges. (calling addFace would re-add the edges)
  // gsSolidHalfFace<T> * newFace = addFace(newFaceBoundary, newTS);
  edge.push_back( newHE );
  edge.push_back( mate );
  
  gsSolidHalfFaceHandle newFace = new gsSolidHalfFace<T>();
  newFace->surf = newTS;
  newFace->setBoundary( newFaceBoundary );
  face.push_back(newFace);
  newFace->setId(numHalfFaces++);
  newFace->vol = f->vol;
  newFace->vol->face.push_back(newFace);
  
  // update all the half-edges on the new face
  tempEdge = mate;
  do
  {
    tempEdge->face = newFace;
    tempEdge = tempEdge->next;
  } while(tempEdge != mate);

  checkStructure();
  
  return newFace;
}


template <class T>
gsMultiPatch<T> gsSolid<T>::plotEdgeGraph()
{
    gsMultiPatch<T> mp;
    std::vector<gsCurve<T>*>  loopv;	
    for (unsigned i=0;i!=numHalfFaces;i++)
    {
        loopv = face[i]->surf->domain().outer().curves();
        for (unsigned j=0;j!=loopv.size();j++)
        {
            mp.addPatch(*loopv[j]);
        }
    }
    return mp;
}
  
template <class T>
gsVolumeBlock<T> *gsSolid<T>::newVolume(gsSolidHalfFaceHandle startingFace)
{
  // set up the new volume
  gsVolumeHandle newVol = new gsVolumeBlock<T>;
  volume.push_back(newVol);
  newVol->setId(numVolumes++);
  // TODO: set the new volume's faces, and update the old volume's faces.
  //newVol->face.push(startingFace);
  // start the queue off with a single face
  std::queue<gsSolidHalfFaceHandle> unprocessed;
  unprocessed.push(startingFace);
  std::set<gsSolidHalfFaceHandle> seen;
  seen.insert(startingFace);
  // process the queue
  while(!unprocessed.empty())
  {
    // grab first face off the queue
    gsSolidHalfFaceHandle nextFace = unprocessed.front();
    unprocessed.pop();
    // set the face's volume
    nextFace->vol = newVol;
    newVol->face.push_back(nextFace);
    size_t numLoops = nextFace->loop.size();
    // add all neighbours to the queue, if they haven't already been seen
    for(size_t i = 0; i < numLoops; i++)
    {
      gsSolidHalfEdgeHandle firstEdge = nextFace->loop[i];
      gsSolidHalfEdgeHandle nextEdge = firstEdge;
      do
      {
        gsSolidHalfFaceHandle neighbFace = nextEdge->mate->face;
        if(seen.find(neighbFace) == seen.end())
        {
          unprocessed.push(neighbFace);
          seen.insert(neighbFace);
        }
        nextEdge = nextEdge->next;
      } while(nextEdge != firstEdge);
    }
  }

  // process every volume's list of faces, keeping only faces with matching volume.
  unsigned foundFaces = 0;
  GISMO_ASSERT(numVolumes == volume.size(), "Number of volumes does not match");
  for(size_t idxV = 0; idxV < numVolumes; idxV++)
  {
    size_t outF = 0;
    size_t numF = volume[idxV]->face.size();
    for(size_t inF = 0; inF < numF; inF++)
    {
      if(volume[idxV]->face[inF]->vol == volume[idxV])
      {
        volume[idxV]->face[outF] = volume[idxV]->face[inF];
        outF++;
        foundFaces++;
      }
    }
    volume[idxV]->face.resize(outF);
  }
  GISMO_ASSERT(foundFaces == face.size() && foundFaces == numHalfFaces,
               "Failed to assign a volume for all faces.");
  return newVol;
}

template <class T>
void gsSolid<T>::checkStructure(bool checkVerts) const
{
    GISMO_UNUSED(checkVerts);
    size_t numEdges = edge.size();
    for(size_t idxE = 0; idxE < numEdges; idxE++)
    {
        gsSolidHalfEdgeHandle e = edge[idxE];
        if (e->mate->mate != e) gsWarn<< "Inconsistent solid graph.\n";
        GISMO_ASSERT(e->prev->target() == e->source, "Inconsistent solid graph");
        GISMO_ASSERT(e->mate->face->vol == e->face->vol, "Inconsistent solid graph");
        GISMO_ASSERT(e->next->face == e->face, "Inconsistent solid graph");
        GISMO_ASSERT(e->next->prev == e, "Inconsistent solid graph");
        GISMO_ASSERT(e->prev->next == e, "Inconsistent solid graph");
    }

    // check that the trimmed surface data matches the info we have about the face
    for(size_t idxF = 0; idxF < numHalfFaces; idxF++)
    {
        gsSolidHalfFace<T> *f = face[idxF];
        gsTrimSurface<T> *surf = f->surf;
        size_t numLoops = surf->domain().numLoops();
        GISMO_ASSERT(numLoops == f->loop.size(), "Number of holes in face does not match corresponding trimmed surface");
        for(size_t loopNum = 0; loopNum < numLoops; loopNum++)
        {
            gsSolidHalfEdge<T> *e = f->loop[loopNum];
            const gsCurveLoop<T> & curveLoop = surf->domain().loop(loopNum);
            size_t numCurves = curveLoop.size();
            for(size_t curveNum = 0; curveNum < numCurves; curveNum++)
            {
                gsMatrix<T> surfCoord = surf->vertexCoord(loopNum, curveNum);
                GISMO_ASSERT(!checkVerts ||
                             (surfCoord - e->source->coords).norm() < 0.0001, "Vertices in surface do not match vertex coordinates");

                const gsMatrix<T>& thisCurveCoefs = curveLoop.curve(curveNum).coefs();
                gsMatrix<T> thisCurveEnd = thisCurveCoefs.row(thisCurveCoefs.rows() - 1);
                const gsMatrix<T>& nextCurveCoefs = curveLoop.curve((curveNum + 1) % numCurves).coefs();
                gsMatrix<T> nextCurveStart = nextCurveCoefs.row(0);
                GISMO_ASSERT((thisCurveEnd - nextCurveStart).norm() < 0.0001, "gsCurveLoop did not close up");
                e = e->next;
            }
            GISMO_ASSERT(e == f->loop[loopNum], "Face's loop did not close up properly");
        }
    }
}

template <class T>
gsSolidHalfFace<T> *gsSolid<T>::addFaceWithMate(const std::vector<gsSolidHeVertexHandle> &verts, gsTrimSurface<T> *surf)
{
  GISMO_ASSERT(verts.size() >= 2, "Unexpected number of vertices in new face");
  checkStructure();

  std::vector<gsSolidHeVertexHandle> reorderedVerts;
  reorderedVerts.push_back(verts[verts.size() - 1]);
  for(size_t iV = 0; iV < verts.size() - 1; iV++)
  {
    reorderedVerts.push_back(verts[iV]);
  }

  surf->cleanEndpoints(0.0000001);

  // special check to make sure the corners of the new trimmed surface actually
  // match the coordinates of the vertex.
  size_t numLoops = surf->domain().numLoops();
  for(size_t loopNum = 0; loopNum < numLoops; loopNum++)
  {
    size_t numCurves = surf->domain().loop(loopNum).size();
    for(size_t curveNum = 0; curveNum < numCurves; curveNum++)
    {
      gsMatrix<T> surfCoord = surf->vertexCoord(loopNum, curveNum);
      GISMO_ASSERT((surfCoord - reorderedVerts[curveNum]->coords).norm() < 0.0001, "Vertices in surface do not match vertex coordinates");
    }
  }

  typedef typename gsTensorBSplineBasis<2,T>::GeometryType MasterSurface;
  // Create a new master surface by taking this surface's one and flipping one coordinate.
  gsTrimSurface<T> *surfReverse = new gsTrimSurface<T>(surf->getTP()->clone().release(), surf->domain().clone().release());
  gsMatrix<T> &surfRevCoefs = surfReverse->getTP()->coefs();
  MasterSurface *genGeom = dynamic_cast<MasterSurface *>(surfReverse->getTP().get());
  GISMO_ASSERT(genGeom != NULL, "This procedure requires a gsGenericGeometry");
  gsMatrix<T> revSupp = surfReverse->getTP()->support();
  GISMO_ASSERT(revSupp.cols() == 2, "Unexpected support");

  gsMatrix<T> testPt(2,4);
  testPt << 0.75 * revSupp(0, 0) + 0.25 * revSupp(0, 1), 0.75 * revSupp(0, 0) + 0.25 * revSupp(0, 1), 0.25 * revSupp(0, 0) + 0.75 * revSupp(0, 1), 0.25 * revSupp(0, 0) + 0.75 * revSupp(0, 1),
            0.75 * revSupp(1, 0) + 0.25 * revSupp(1, 1), 0.25 * revSupp(1, 0) + 0.75 * revSupp(1, 1), 0.75 * revSupp(1, 0) + 0.25 * revSupp(1, 1), 0.25 * revSupp(1, 0) + 0.75 * revSupp(1, 1);
  gsMatrix<T> testVal1;
  genGeom->eval_into(testPt, testVal1);

  unsigned totalCPs = surfRevCoefs.rows();
  unsigned dim = surfRevCoefs.cols();
  unsigned blockSize = genGeom->basis().component(0).size();
  for(unsigned blockStart = 0; blockStart < totalCPs; blockStart += blockSize) // note step size
  {
    surfRevCoefs.block(blockStart, 0, blockSize, dim) =
      surfRevCoefs.block(blockStart, 0, blockSize, dim).colwise().reverse().eval();
  }

  // test that the surface has been correctly reversed
  gsMatrix<T> testVal2;
  testPt << 0.25 * revSupp(0, 0) + 0.75 * revSupp(0, 1), 0.25 * revSupp(0, 0) + 0.75 * revSupp(0, 1), 0.75 * revSupp(0, 0) + 0.25 * revSupp(0, 1), 0.75 * revSupp(0, 0) + 0.25 * revSupp(0, 1),
            0.75 * revSupp(1, 0) + 0.25 * revSupp(1, 1), 0.25 * revSupp(1, 0) + 0.75 * revSupp(1, 1), 0.75 * revSupp(1, 0) + 0.25 * revSupp(1, 1), 0.25 * revSupp(1, 0) + 0.75 * revSupp(1, 1);

  genGeom->eval_into(testPt, testVal2);
  gsMatrix<T> origCoefs = surf->getTP()->coefs();
  gsMatrix<T> testRevCoefs = surfReverse->getTP()->coefs();
  GISMO_ASSERT((testVal2 - testVal1).norm() < 0.0001, "Surface was not correctly reversed");

  // Create new trimming loops by starting with this surfaces's one and
  //  (i) reversing the whole loop
  //  (ii) flipping one coordinate
  numLoops = surfReverse->domain().numLoops();
  for(size_t loopNum = 0; loopNum < numLoops; loopNum++)
  {
    gsCurveLoop<T> & revCurveLoop = surfReverse->domain().loop(loopNum);
    size_t numCurves = revCurveLoop.size();
    for(size_t curveNum = 0; curveNum < numCurves; curveNum++)
    {
      gsGeometry<T> & thisCurve = revCurveLoop.curve(curveNum);
      gsMatrix<T> &curveCoefs = thisCurve.coefs();
      curveCoefs = curveCoefs.colwise().reverse().eval();
      size_t rows = curveCoefs.rows();
      for(size_t cpIdx = 0; cpIdx < rows; cpIdx++)
      {
        curveCoefs(cpIdx, 0) = revSupp(0, 0) + revSupp(0, 1) - curveCoefs(cpIdx, 0);
      }
    }

    // reverse order of curves in the curve loop
    std::reverse( revCurveLoop.curves().begin(), revCurveLoop.curves().end() );
  }

  // the following loop is only for debugging purposes: check that curves are joined
  for(size_t loopNum = 0; loopNum < numLoops; loopNum++)
  {
    const gsCurveLoop<T> & revCurveLoop = surfReverse->domain().loop(loopNum);
    size_t numCurves = revCurveLoop.size();
    for(size_t curveNum = 0; curveNum < numCurves; curveNum++)
    {
      // make sure each curve starts where the previous one ends
      const gsGeometry<T> & thisCurve = revCurveLoop.curve(curveNum);
      const gsGeometry<T> & nextCurve = revCurveLoop.curve((curveNum + 1) % numCurves);
      gsMatrix<T> thisSupp = thisCurve.support();
      gsMatrix<T> nextSupp = nextCurve.support();
      gsMatrix<T> thisEndPt = nextCurve.eval(nextSupp.col(0));
      gsMatrix<T> nextStartPt = thisCurve.eval(thisSupp.col(1));
      GISMO_ASSERT((thisEndPt - nextStartPt).norm() < 0.0001, "Invalid curve loop");
      // make sure corresponding points on the opposite faces have the same location in space
      gsMatrix<T> thisCoord = surfReverse->vertexCoord(loopNum, curveNum);
      gsMatrix<T> otherCoord = surf->vertexCoord(loopNum, (numCurves - curveNum) % numCurves);
      GISMO_ASSERT((otherCoord - thisCoord).norm() < 0.0001, "Mate faces do not have matching vertices");
      GISMO_ASSERT((otherCoord - reorderedVerts[(numCurves - curveNum) % numCurves]->coords).norm() < 0.0001, "Vertices in surface do not match vertex coordinates");
    }
  }

  // record the edge where we will start the updating process
  gsSolidHalfEdgeHandle startEdge = reorderedVerts[0]->getHalfEdge(reorderedVerts[1]);
  // duplicate all of the relevant vertices
  std::vector<gsSolidHeVertexHandle> newVerts, newVertsReverse;
  for(typename std::vector<gsSolidHeVertexHandle>::const_iterator iter = reorderedVerts.begin(); iter != reorderedVerts.end(); iter++)
  {
    addHeVertex((*iter)->coords(0), (*iter)->coords(1), (*iter)->coords(2));
    newVerts.push_back(vertex[numVertices - 1]);
  }
  // reverse vertices start with vertex 0 and go back the other way
  newVertsReverse.push_back(newVerts[0]);
  for(int i = newVerts.size() - 1; i > 0; i--) newVertsReverse.push_back(newVerts[i]);
  // create the first new face using the already existing vertices
  gsSolidHalfFaceHandle frontFace = addFace(reorderedVerts, surf);
  // create its mate using the duplicate vertices and the reverse trim surface
  gsSolidHalfFaceHandle backFace = addFace(newVertsReverse, surfReverse);
  backFace->vol = NULL;
  // set faces as each others' mates
  frontFace->mate = backFace;
  backFace->mate = frontFace;
  // make sure none of the vertices along the loop have their "hed" pointing
  // into the new volume
  size_t vertNum, numVerts = reorderedVerts.size();
  for(vertNum = 0; vertNum < numVerts; vertNum++)
  {
    gsSolidHeVertex<T> *thisVert = reorderedVerts[vertNum];
    gsSolidHeVertex<T> *nextVert = reorderedVerts[(vertNum + 1) % numVerts];
    while(thisVert->hed->next->source != nextVert)
    {
      thisVert->hed = thisVert->hed->mate->next;
      GISMO_ASSERT(thisVert->hed->source == thisVert, "Inconsistent solid graph");
    }
    thisVert->hed = thisVert->hed->mate->next;
    GISMO_ASSERT(thisVert->hed->source == thisVert, "Inconsistent solid graph");
  }
  // loop over half-edges "he" in the loop
  gsSolidHalfEdge<T> *he = startEdge;
  vertNum = 0;
  gsSolidHalfEdgeHandle newFrontFaceEdge = frontFace->loop[0];
  while(newFrontFaceEdge->source != reorderedVerts[0]) newFrontFaceEdge = newFrontFaceEdge->next;
  gsSolidHalfEdgeHandle newBackFaceEdge = backFace->loop[0];
  while(newBackFaceEdge->source != newVerts[1]) newBackFaceEdge = newBackFaceEdge->next;
  do
  {
    // set mates for the new front face
    he->mate->mate = newFrontFaceEdge;
    newFrontFaceEdge->mate = he->mate;
    // set mates for the new back face
    he->mate = newBackFaceEdge;
    newBackFaceEdge->mate = he;
    // set volume for the new back face (new front face's volume will
    // get set later when we chase round the volume)
    GISMO_ASSERT(backFace->vol == NULL || backFace->vol == he->face->vol, "Inconsistent solid graph");
    backFace->vol = he->face->vol;
    // set source
    GISMO_ASSERT(newFrontFaceEdge->source == reorderedVerts[vertNum], "Inconsistent solid graph");
    GISMO_ASSERT(newBackFaceEdge->source == newVerts[(vertNum + 1) % numVerts], "Inconsistent solid graph");
    //if(he->source->hed == he) he->source->hed = newFrontFaceEdge;
    //he->source = newVerts[vertNum];
    // set up for the next edge
    newFrontFaceEdge = newFrontFaceEdge->next;
    newBackFaceEdge = newBackFaceEdge->prev;
    vertNum++;
    gsSolidHeVertexHandle nextVert = reorderedVerts[(vertNum + 1) % numVerts];
    // loop round the edges coming out of this vertex, updating their source to the
    // new vertex and looking for the next edge in the loop.
    while(true) {
      he = he->next;
      GISMO_ASSERT(he->source->hed->source == he->source, "Inconsistent solid graph");
      if(he->source->hed == he) he->source->hed = newFrontFaceEdge;
      he->source = newVerts[vertNum % numVerts];
      if(he->next->source == nextVert || he == startEdge) break;
      GISMO_ASSERT(he->source->hed->source == he->source, "Inconsistent solid graph");
      GISMO_ASSERT(he->prev->target() == he->source, "Inconsistent solid graph");
      he = he->mate;
    }
  } while(he != startEdge);
  // add the newly created back face to its volume's list (the front face will
  // be assigned to a new volume when we call newVolume)
  GISMO_ASSERT(backFace->vol != NULL, "Could not find the volume for mate face");
  backFace->vol->face.push_back(backFace);

  // chase all the faces in the new volume, setting their volume member correctly
  newVolume(frontFace);

  checkStructure();


  // TODO: the following not always true, only do this temporarily
  gsSolidHalfEdge<T>* he0;
  he0 = frontFace->loop[0];
  he = he0;
  for (size_t i=0;i<this->nHalfEdges();i++)
  {
      he->is_convex = true;
      he->mate->is_convex = true;
      if (he->target()==he0->source) break;
      he = he->next;
  }

  he0 = frontFace->mate->loop[0];
  he = he0;
  for (size_t i=0;i<this->nHalfEdges();i++)
  {
      he->is_convex = true;
      he->mate->is_convex = true;
      if (he->target()==he0->source) break;
      he = he->next;
  }

  return frontFace;
}

template <class T>
std::vector<typename gsSolid<T>::gsSolidHalfEdgeHandle > gsSolid<T>::impedingEdges(gsSolidHalfEdgeHandle he) const
{
    //int vol = he->face->vol->getId();
    // collect vertices of the two faces incident to the HE
    std::vector<gsSolidHeVertexHandle> v0,v1;
    gsSolidHalfFaceHandle face0 = he->face;
    gsSolidHalfFaceHandle face1 = he->mate->face;    
    v0 = face0->getVertices();
    v1 = face1->getVertices();
    // remove two vertices of the HE from the two vertex sets
    typename std::vector<gsSolidHeVertexHandle>::iterator it,it0,it1;
    for (it = v0.begin();it!=v0.end();++it)
    {
        if ( ((*it)==he->source) || ((*it)==he->target()) )
        {it0 = it;break;}
    }
    v0.erase(it0, it0 + 2);
    for (it = v1.begin();it!=v1.end();++it)
    {
        if ( ((*it)==he->source) || ((*it)==he->target()) )
        {it0 = it;break;}
    }
    v1.erase(it0, it0 + 2);
    // check if each vertex in v0 is connected by a HE to another in v1
    std::vector<gsSolidHalfEdgeHandle> iHE,vHE;
    typename std::vector<gsSolidHalfEdgeHandle>::iterator ite;
    for (it0 = v0.begin();it0!=v0.end();++it0)
    {
        vHE = (*it0)->halfEdges();
        for (ite = vHE.begin();ite!=vHE.end();++ite)
        {
            for (it1 = v1.begin();it1!=v1.end();++it1)
            {
                if ( (*ite)->target()== *it1) iHE.push_back(*ite);
            }

        }
    }
    return iHE;
}

template <class T>
void gsSolid<T>::insertNewVertex(gsSolidHalfEdgeHandle he)
{
    checkStructure();
    gsSolidHalfEdgeHandle hem = he->mate;
    gsSolidHalfFaceHandle face0,face1;
    face0 = he->face;
    face1 = he->mate->face;
    size_t nFace = this->nHalfFaces();
    bool convex = he->is_convex;
    int heLoopN = he->loopN();
    int hemLoopN = hem->loopN();
    // adapting trimming structure of the two incident faces: first, face0:
    size_t he_tCurveID = face0->indexOfEdge(he);
    size_t hem_tCurveID = face1->indexOfEdge(hem);
    gsMatrix<T> mid0 = face0->surf->splitCurve(he->loopN(), he_tCurveID);
    gsMatrix<T> space0;
    face0->surf->getTP()->eval_into(mid0.transpose(), space0);
    // TODO: deal with possibility that the new vertex is on an inner loop
    T nearestParam = face1->surf->nearestPoint(0, hem_tCurveID, 10, 10, space0);
    gsMatrix<T> mid1 = face1->surf->splitCurve(hem->loopN(), hem_tCurveID, nearestParam);
    gsMatrix<T> space1;
    face1->surf->getTP()->eval_into(mid1.transpose(), space1);
    //GISMO_ASSERT((space0 - space1).norm() < 0.0001, "Midpoints of curves do not match");
    //denote:
    // he->source----------newV-----------he->target()
    //          s-----------n-------------t
    gsSolidHeVertexHandle s,n,t;
    // create new HeVertex located at the midpoint between the two new trimmed surface corners
    n = new gsSolidHeVertex<T>(
                0.5 * (space0(0, 0) + space1(0, 0)),
                0.5 * (space0(1, 0) + space1(1, 0)),
                0.5 * (space0(2, 0) + space1(2, 0)),
                this->nVertices());
    // add n to the vertex list
    (this->vertex).push_back(n);

    // add four new HEs: sn,ns,nt,tn
    s = he->source;
    t = he->target();
    gsSolidHalfEdgeHandle sn = new gsSolidHalfEdge<T>(s,face0,nFace,convex,heLoopN);
    gsSolidHalfEdgeHandle nt = new gsSolidHalfEdge<T>(n,face0,nFace+1,convex,heLoopN);
    gsSolidHalfEdgeHandle tn = new gsSolidHalfEdge<T>(t,face1,nFace+2,convex,hemLoopN);
    gsSolidHalfEdgeHandle ns = new gsSolidHalfEdge<T>(n,face1,nFace+3,convex,hemLoopN);
    //
    he->prev->next = sn;
    sn->next = nt;
    nt->next = he->next;
    he->next->prev = nt;
    nt->prev = sn;
    sn->prev = he->prev;
    //
    hem->prev->next = tn;
    tn->next = ns;
    ns->next = hem->next;
    hem->next->prev = ns;
    ns->prev = tn;
    tn->prev = hem->prev;
    //
    sn->mate = ns;
    ns->mate = sn;
    nt->mate = tn;
    tn->mate = nt;
    //
    n->hed = nt;
    this->edge.push_back(sn);
    this->edge.push_back(nt);
    this->edge.push_back(tn);
    this->edge.push_back(ns);

    // reassign pointers to each loop of a face
    if ((face0->loop).at(he->loopN())==he) (face0->loop).at(he->loopN()) = sn;
    if ((face1->loop).at(he->mate->loopN())==he->mate) (face1->loop).at(he->mate->loopN()) = tn;

    // reassign HE pointers of s and t
    if (s->hed == he) s->hed = he->prev->mate;
    if (t->hed == hem) t->hed = hem->prev->mate;

    // remove the HE with larger ID
    int idmax = he->getId();
    if (he->mate->getId()>idmax) idmax=he->mate->getId();
    int idmin = he->getId();
    if (he->mate->getId()<idmin) idmin=he->mate->getId();
    // recalculate ID for HEs
    for (size_t i = idmax+1; i < this->nHalfEdges(); i++)
    {
        this->edge[i]->setId(i-1);
    }
    this->edge.erase(this->edge.begin() + idmax);
    // remove the HE with smaller ID
    for (size_t i = idmin+1; i < this->nHalfEdges(); i++)
    {
        this->edge[i]->setId(i-1);
    }
    this->edge.erase(this->edge.begin() + idmin);
    //
    this->numVertices = this->nVertices();
    this->numHalfEdges = this->nHalfEdges();
    this->numHalfFaces = this->nHalfFaces();
    this->numVolumes = this->nVolumes();
    gsDebug<<"A new vertex is added to the halfedge with source: "<<*he->source<<
               " and target: "<<*he->target()<<std::endl;
    delete he;
    delete hem;
    checkStructure();
}

template <class T>
void gsSolid<T>::handleImpedingEdges(gsSolidHalfEdgeHandle he)
{
    typename std::vector<gsSolidHalfEdgeHandle> iHe;
    iHe=this->impedingEdges(he);
    typename std::vector<gsSolidHalfEdgeHandle>::const_iterator it;
    if (iHe.empty()==false)
    {
        for (it=iHe.begin();it!=iHe.end();++it)
        {
            this->insertNewVertex(*it);
        }
    }
    checkStructure();
}

} // namespace
