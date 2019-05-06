/** @file gsTriMeshToSolid.hpp

    @brief Provides implementation gsTriMeshToSolid class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mayer, M. Pauley
*/

#include <gsModeling/gsModelingUtils.hpp>
#include <gsModeling/gsCurveLoop.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsModeling/gsTrimSurface.h>
 #include <gsModeling/gsSolid.h>
 #include <gsNurbs/gsBSpline.h>

#include <fstream>

namespace gismo
{

template <class T>
void gsTriMeshToSolid<T>::calcPatchNumbers()
{
    // determine which faces form a big face
    std::vector<int> added;
    for (std::vector<int>::size_type i=0;i!=face.size();i++)
        added.push_back(0);
    numBigFaces=0;
    std::stack<FaceHandle> faceCollect;
    FaceHandle tempFace;
    for( typename std::vector<FaceHandle >::iterator it(face.begin());it!=face.end();++it)
    {
        if (added[(**it).getId()]==0&&
                *((**it).vertices[0])!=*((**it).vertices[1])&&
                *((**it).vertices[2])!=*((**it).vertices[1])&&
                *((**it).vertices[0])!=*((**it).vertices[2]))
        {
            numBigFaces++;
            faceCollect.push(*it);
            while(!faceCollect.empty())
            {
                (*(faceCollect.top())).faceIdentity=numBigFaces;
                added[(*(faceCollect.top())).getId()]=1;
                tempFace=faceCollect.top();
                faceCollect.pop();

                for (size_t i =0;i<(*tempFace).nFaces.size();i++)
                {
                    if(added[tempFace->nFaces[i]->getId()]==0)
                    {
                        int k=0;
                        while((*(tempFace->nFaces[i])!=*(edge[tempFace->nEdges[k]->getId()].nFaces[0])&&
                               ((edge[tempFace->nEdges[k]->getId()].nFaces.size()==1)||
                              *(tempFace->nFaces[i])!=*(edge[tempFace->nEdges[k]->getId()].nFaces[1]))))
                            k++;
                        if (edge[tempFace->nEdges[k]->getId()].sharp==0)
                            faceCollect.push((*tempFace).nFaces[i]);
                    }
                }
            }
        }
    }
}


template <class T>
void gsTriMeshToSolid<T>::getFeatures(T angleGrad,bool& bWarnNonManifold,bool& bWarnBorders)
{

    gsDebug<<"Getting the features..."<<"\n";
    bWarnNonManifold=false;
    bWarnBorders=false;
    //create 3 edges for each face
    int fsize=face.size();
    for( int it=0;it<fsize;++it)
    {
        for(int i=0;i<3;i++)
        {
            VertexHandle p0 = face[it]->vertices[i];
            VertexHandle p1 = face[it]->vertices[(i+1)%3];
            if (*p0!=*p1)
            {
                if ( Xless<T>(p0,p1) ) std::swap(p0,p1);
                //Edge tempEdge;
                //tempEdge=Edge(p0,p1);
                //insertObject<Edge >(tempEdge,edge);
                edge.push_sorted_unique( Edge(p0,p1) );
            }
            else
                gsWarn<<"face "<<it<<" has 2 common vertices"<<"\n"<<*p0<<*p1<<"\n";
        }
    }

    numEdges=edge.size();
    //number the edges
    int iterId=0;
    for(typename std::vector<Edge>::iterator iter(edge.begin());iter!=edge.end();++iter)
    {
        iter->setId(iterId);
        iterId++;
    }
    //determine neighboring faces of each edge
    for(size_t it=0;it<face.size();++it)
    {
        if(*(face[it]->vertices[0])!=*(face[it]->vertices[1])&&
           *(face[it]->vertices[2])!=*(face[it]->vertices[1])&&
           *(face[it]->vertices[0])!=*(face[it]->vertices[2]))
        {
            for(int i=0;i<3;i++)
            {
                VertexHandle p0 = face[it]->vertices[i];
                VertexHandle p1 = face[it]->vertices[(i+1)%3];
                if ( Xless<T>(p0,p1) ) std::swap(p0,p1);
                //Edge tempEdge;
                //tempEdge=Edge(p0,p1);
                //edge[getIndex<Edge>(tempEdge,edge)].nFaces.push_back(face[it]);
                Edge * eit = edge.find_ptr_or_fail( Edge(p0,p1) );

                eit->nFaces.push_back(face[it]);

                //face[it]->nEdges.push_back(&edge[getIndex<Edge>(tempEdge,edge)]);
                face[it]->nEdges.push_back(eit);
            }
        }
    }

    // Extract those edges whose adjacent triangles form a large angle
    for(typename std::vector<Edge>::iterator iter(edge.begin());iter!=edge.end();++iter)
    {
        std::vector<FaceHandle > vT;
        for (size_t i=0;i<iter->nFaces.size();i++)
        {
            vT.push_back(iter->nFaces[i]);
        }
        if(vT.size()==1)
        {
            bWarnBorders=true;
            iter->sharp=1;
            continue;
        }
        if(vT.size()>2)
        {
            bWarnNonManifold=true;
            gsWarn<<"non manifold edge"<<"\n";
            iter->sharp=1;
            continue;
        }
        GISMO_ASSERT(vT.size()==2, "Edge must belong to two triangles, got "<<vT.size() );
        gsVector3d<T> nv0((*vT[0]).orthogonalVector());
        nv0 = nv0/(math::sqrt(nv0.squaredNorm()));
        gsVector3d<T> nv1((*vT[1]).orthogonalVector());
        nv1 = nv1/(math::sqrt(nv1.squaredNorm()));
        T cosPhi( nv0.dot(nv1) );
        // Numerical robustness
        if(cosPhi>1.0) cosPhi=1.0;
        else if(cosPhi<-1.0) cosPhi=-1.0;

        const T PI_(3.14159);
        T phiGrad(math::acos(cosPhi)/PI_*180);
        if(phiGrad>=angleGrad)
            iter->sharp=1;
        else
            iter->sharp=0;

        (*face[(*vT[0]).getId()]).nFaces.push_back(face[(*vT[1]).getId()]);
        (*face[(*vT[1]).getId()]).nFaces.push_back(face[(*vT[0]).getId()]);
    }

    if(bWarnNonManifold)
    {
        gsDebug<<"Attention: There were non-manifold edges. These are always features"<<"\n";
    }
    if(bWarnBorders)
    {
        gsDebug<<"Attention: There were border edges. These are always features"<<"\n";
    }
}

template <class T>
void gsTriMeshToSolid<T>::setSharpEdges(std::vector< gsEdge<T> > & featEdges, int useFeatures)
{
    if(useFeatures==1||useFeatures==2)
    {
        if(useFeatures==2)
        {
            for(size_t i=0;i<edge.size();i++)
            {
                edge[i].sharp=0;
            }
        }
        for(size_t i=0;i<featEdges.size();i++)
        {
            bool foundEdge=0;
            for(size_t j=0;j<edge.size();j++)
            {
                if(approxEqual(featEdges[i],edge[j])==1)
                {
                    edge[j].sharp=1;
                    foundEdge=1;
                }
            }
            if (foundEdge==0)
                gsWarn<<"could not find the feature number"<<i<<",review the input data"<<"\n";
        }
    }
    int i1=0;
    for(typename std::vector<Edge>::iterator iter(edge.begin());iter!=edge.end();++iter)
    {
        if(iter->sharp==1)
            i1++;
    }
    gsDebug<<"Feature lines: "<<i1<<"\n";
}


template <class T>
void gsTriMeshToSolid<T>::storeNeighboringFaces()
{
    //store information about the neighboring faces of each edge in the edges
    for (typename std::vector<Edge >::iterator it(edge.begin());it!=edge.end();++it)
    {
        if ((*it).nFaces[0]->faceIdentity==(*it).nFaces[1]->faceIdentity)
            (*it).numPatches.push_back((*it).nFaces[0]->faceIdentity);
        else
        {
            (*it).numPatches.push_back((*it).nFaces[0]->faceIdentity);
            (*it).numPatches.push_back((*it).nFaces[1]->faceIdentity);
        }
    }
}


template <class T>
void gsTriMeshToSolid<T>::divideAndMergePatches(T innerAngle, T patchAreaWeight, T mergeSmallPatches)
{
    //calculate areas of the patches
    std::vector<T > areas;
    for (int i=0;i<numBigFaces;i++)
        areas.push_back(0);
    for (size_t i=0;i<face.size();i++)
    {
        if(face[i]->faceIdentity>0&&face[i]->faceIdentity<=numBigFaces)
            areas[face[i]->faceIdentity-1]+=calcArea(face[i]);
    }
    T totalArea=0;
    for(size_t i=0;i<areas.size();i++)
    {
        //gsDebug<<"area of ith patch: "<<areas[i]<<"\n";
        totalArea+=areas[i];
    }
    T averageArea=totalArea/areas.size();
    for(size_t j=0;j<edge.size();j++)
    {
        if(edge[j].numPatches.size()==1)
        {

            if(areas[edge[j].numPatches[0]-1]>averageArea*patchAreaWeight)
            {

                gsVector3d<T> nv0(edge[j].nFaces[0]->orthogonalVector());
                nv0 = nv0/(math::sqrt(nv0.squaredNorm()));
                gsVector3d<T> nv1(edge[j].nFaces[1]->orthogonalVector());
                nv1 = nv1/(math::sqrt(nv1.squaredNorm()));
                T cosPhi( nv0.dot(nv1) );
                // Numerical robustness
                if(cosPhi>1.0) cosPhi=1.0;
                else if(cosPhi<-1.0) cosPhi=-1.0;

                const T PI_(3.14159);
                T phiGrad(math::acos(cosPhi)/PI_*180);
                if(phiGrad>=innerAngle)
                    edge[j].sharp=1;
                else
                    edge[j].sharp=0;
            }
        }
    }
    int i2=0;
    for(typename std::vector<Edge>::iterator iter(edge.begin());iter!=edge.end();++iter)
    {
        if(iter->sharp==1)
            i2++;
    }
    gsDebug<<"Feature lines after dividing patches: "<<i2<<"\n";
    if(mergeSmallPatches!=0)
    {
        for(size_t j=0;j<edge.size();j++)
        {
            if(edge[j].numPatches.size()==2)
            {
                if(areas[edge[j].numPatches[0]-1]<averageArea*mergeSmallPatches&&areas[edge[j].numPatches[1]-1]<averageArea*mergeSmallPatches)
                    edge[j].sharp=0;
            }
        }
    }
    int i3=0;
    for(typename std::vector<Edge>::iterator iter(edge.begin());iter!=edge.end();++iter)
    {
        if(iter->sharp==1)
            i3++;
    }
    gsDebug<<"Feature lines after merging patches: "<<i3<<"\n";
}

template <class T>
void gsTriMeshToSolid<T>::getFaces(std::vector<std::vector<VertexHandle> > & iPoints, std::vector<std::vector<VertexHandle> > & oPoints,
              std::vector< std::vector<std::vector<VertexHandle> > > & innerBdrys, std::vector< std::vector<Vertex>  > & innerBdrysMassP,
              std::vector<std::vector<bool> > & oPointsConvexFlag)
{
    gsDebug<<"Getting the faces..."<<"\n";
    this->calcPatchNumbers();

    //eliminate sharp edges in the interior of big faces
    for (typename std::vector<Edge >::iterator it(edge.begin());it!=edge.end();++it)
    {
        if ((*it).sharp==1&&(*it).nFaces[0]->faceIdentity==(*it).nFaces[1]->faceIdentity)
            (*it).sharp=0;
    }
    //determine sharp attribute for vertices
    for (size_t i=0;i<edge.size();i++)
    {
        edge[i].source->nVertices.push_back(edge[i].target);
        edge[i].target->nVertices.push_back(edge[i].source);
        if (edge[i].sharp==1)
        {
            edge[i].source->sharp=1;
            edge[i].target->sharp=1;
        }
        else
        {
            if(edge[i].source->sharp!=1)
                edge[i].source->sharp=0;
            if(edge[i].target->sharp!=1)
                edge[i].target->sharp=0;
        }
    }
    // put all sharp edges of a big face in a multimap
    std::multimap<int,EdgeHandle> mmIE;
    for (typename std::vector<Edge >::iterator it(edge.begin());it!=edge.end();++it)
    {
        if ((*it).sharp==1&&(*it).nFaces[0]->faceIdentity!=(*it).nFaces[1]->faceIdentity)
        {
            if ((*it).nFaces.size()!=2)
                gsDebug<<"edge "<<*it<<" has "<<(*it).nFaces.size()<<" neighboring faces"<<"\n";
            mmIE.insert(std::make_pair((*it).nFaces[0]->faceIdentity,&(*it)));
            mmIE.insert(std::make_pair((*it).nFaces[1]->faceIdentity,&(*it)));
        }
    }
    for (std::vector<int>::size_type i=0;i!=vertex.size();i++)
        vertex[i]->numEdges=0;
    for( typename std::vector<Edge >::iterator it(edge.begin());it!=edge.end();++it)
    {
        if((*it).nFaces.size()==2&&((*it).nFaces[0]->faceIdentity!=(*it).nFaces[1]->faceIdentity))
        {
            vertex[(*it).source->getId()]->numEdges++;
            vertex[(*it).target->getId()]->numEdges++;
        }
    }

    gsDebug<<"Getting interior points..."<<"\n";
    //calculating interior points of each big face
    std::vector< std::set< VertexHandle> > iPointsSet;
    std::set< VertexHandle > vertexSet;

    // first set up iPointsSet[i] for each i, which is the set of all
    // vertices that are on some triangle belonging to big face i
    for (int i=0;i<numBigFaces;i++)
    {
        iPointsSet.push_back(vertexSet);
    }
    for( typename std::vector<FaceHandle >::iterator it(face.begin());it!=face.end();++it)
    {
        if ( *((**it).vertices[0])!=*((**it).vertices[1])&&
             *((**it).vertices[2])!=*((**it).vertices[1])&&
             *((**it).vertices[0])!=*((**it).vertices[2]))
        {
            for (int i=0;i<3;i++)
            {
                iPointsSet[(**it).faceIdentity-1].insert((**it).vertices[i]);
            }
        }
    }

    // now, for each i, set up iPoints[i] which consists of those elements
    // of iPointsSet[i] that are sharp vertices.
    for(int i=0;i<numBigFaces;i++)
    {
        std::vector< VertexHandle > tempVec;
        for (typename  std::set<VertexHandle> ::iterator it(iPointsSet[i].begin());it!=iPointsSet[i].end();++it)
        {
            if((**it).sharp==0)
                tempVec.push_back(*it);
        }
        iPoints.push_back(tempVec);
    }


    gsDebug<<"Getting boundary points..."<<"\n";

    int sourcePos=0;
    int targetPos=0;

    for(int i=1;i<numBigFaces+1;i++)
    {
        //check if all boundaries of a face are used
        bool allEdgesCovered=false;
        std::vector<bool> edgeAdded;
        for(size_t j=0;j<mmIE.count(i);j++)
            edgeAdded.push_back(false);
        T maxLength=0;
        T bdryLength=0;
        std::vector< std::vector<VertexHandle> > innerBdryHelpVec;
        std::vector<Vertex> innerBdryMassPHelpVec;
        while (allEdgesCovered==false)
        {
            std::vector< bool> isConvex;//required for mapping to a u,v plane
            std::vector< T> angle; //to calculate isConvex
            std::vector< VertexHandle> vertexVec; //required for mapping to a u,v plane

            //take the first edge
            EdgeHandle firstEdge=mmIE.find(i)->second;
            int EdgeCount=0;//determine where the first not already used edge is.
            while (edgeAdded[EdgeCount]==true)
                EdgeCount++;
            edgeAdded[EdgeCount]=true;
            int helpCount=0;
            for(typename std::multimap<int,EdgeHandle>::iterator edgeIter = mmIE.find(i);helpCount<=EdgeCount;edgeIter++)
            {
                helpCount++;
                firstEdge=(*edgeIter).second;
            }

            //use a neighboring face of the first edge to determine the direction of the boundary.
            FaceHandle firstFace=NULL;
            if (firstEdge->nFaces[0]->faceIdentity==i)
                firstFace=firstEdge->nFaces[0];
            else if(firstEdge->nFaces[1]->faceIdentity==i)
                firstFace=firstEdge->nFaces[1];
            else
                gsWarn<<"edge has wrong neighboring faces"<<"\n";
            for (int j=0;j<3;j++)
            {
                if (*(firstFace->vertices[j])==*(firstEdge->source))
                    sourcePos=j;
                if (*(firstFace->vertices[j])==*(firstEdge->target))
                    targetPos=j;
            }
            if ((sourcePos-targetPos)==1||sourcePos-targetPos==-2)
            {
                vertexVec.push_back((firstEdge->target));
                vertexVec.push_back((firstEdge->source));
            }
            else if((sourcePos-targetPos)==2||sourcePos-targetPos==-1)
            {
                vertexVec.push_back((firstEdge->source));
                vertexVec.push_back((firstEdge->target));
            }
            else
                gsWarn<<"Edge not found"<<"\n";
            //look for an edge with a common vertex and add it to the boundary.
            EdgeHandle currentEdge=firstEdge;
            int k=0;
            while (k==0||(vertexVec[0])!=(vertexVec[vertexVec.size()-1]))
            {
                k++;
                bool edgeNotFound=1;
                GISMO_UNUSED(edgeNotFound);
                int l=0;
                for (typename std::multimap<int,EdgeHandle>::iterator it(mmIE.find(i));(*it).first==i;++it)
                {
                    if (*(*it).second!=*currentEdge&&*(*it).second->target==*vertexVec[vertexVec.size()-1]&&edgeAdded[l]==false)
                    {
                        vertexVec.push_back(((*it).second->source));
                        //calculate Angle between currentEdge and *it.second
                        angle.push_back(calcAngle(currentEdge,(*it).second,i));
                        currentEdge=(*it).second;
                        edgeNotFound=0;
                        edgeAdded[l]=true;
                        break;
                    }
                    else if (*(*it).second!=*currentEdge&&*(*it).second->source==*vertexVec[vertexVec.size()-1]&&edgeAdded[l]==false)
                    {
                        vertexVec.push_back(((*it).second->target));
                        //calculate Angle between currentEdge and *it.second
                        angle.push_back(calcAngle(currentEdge,(*it).second,i));
                        currentEdge=(*it).second;
                        edgeNotFound=0;
                        edgeAdded[l]=true;
                        break;
                    }
                    l++;

                }
                //calculate angle between first and last edge.
                if (*(vertexVec[0])==*(vertexVec[vertexVec.size()-1]))
                {
                    typename std::vector<T>::iterator it=angle.begin();
                    angle.insert(it,calcAngle(currentEdge, firstEdge,i));
                }
                GISMO_ASSERT(edgeNotFound==0,"edge not found, could not create a closed boundary of sharp edges to identify a face");
            }

            for (size_t j=0;j<angle.size();j++)
            {
                isConvex.push_back(angle[j]<EIGEN_PI);
            }
            //first vertex is added 2 times -> delete last one
            vertexVec.pop_back();

            //store the data of the outer boundary.
            bdryLength=calcBdryLength(vertexVec);
            if (maxLength==0)
            {
                oPoints.push_back(vertexVec);
                maxLength=bdryLength;
                oPointsConvexFlag.push_back(isConvex);
            }
            else if (bdryLength>maxLength)
            {
                innerBdryHelpVec.push_back(oPoints.back());
                innerBdryMassPHelpVec.push_back(getMassP(oPoints.back()));
                oPoints.pop_back();
                oPoints.push_back(vertexVec);
                maxLength=bdryLength;
                oPointsConvexFlag.pop_back();
                oPointsConvexFlag.push_back(isConvex);
            }
            else if(bdryLength<=maxLength)
            {
                innerBdryHelpVec.push_back(vertexVec);
                innerBdryMassPHelpVec.push_back(getMassP(vertexVec));

            }
            //check if all Edges are used yet.
            allEdgesCovered=true;
            for (size_t j=0;j<mmIE.count(i);j++)
            {
                if (edgeAdded[j]==false)
                    allEdgesCovered=false;
            }
        }
        innerBdrys.push_back(innerBdryHelpVec);
        innerBdrysMassP.push_back(innerBdryMassPHelpVec);
    }
    //establish connections between boundary of a hole in a face and generated point in its interior
    for(size_t i=0;i<innerBdrys.size();i++)
    {
        for(size_t j=0;j<innerBdrys[i].size();j++)
        {
            for(size_t k=0;k<innerBdrys[i][j].size();k++)
            {
                innerBdrys[i][j][k]->nVertices.push_back(&innerBdrysMassP[i][j]);
                innerBdrysMassP[i][j].nVertices.push_back(innerBdrys[i][j][k]);
            }
        }
    }

}

template <class T>
void gsTriMeshToSolid<T>::toSolid(gsSolid<T> & sl, std::vector<std::vector<VertexHandle> > & iPoints,
                        std::vector<std::vector<VertexHandle> > & oPoints,
                        std::vector< std::vector<std::vector<VertexHandle> > > & innerBdrys,
                        std::vector< std::vector<Vertex>  > & innerBdrysMassP, // articifial points
                        std::vector<std::vector<bool> > & oPointsConvexFlag,
                        std::vector<gsMesh<T> *> & paraMeshes,
                        std::vector<gsMesh<T> *> & fitMeshes,
                        std::vector<gsMesh<T> *> & patchMeshes,
                        int kvOuterPoints, int kvAdditionalInnerPoints,
                        bool plot, int meshPoints,
                        bool moreInner, T wE, T wI, int closeBoundary, bool noSmooth)
{

    GISMO_ASSERT( (static_cast<int>(iPoints.size()           ) == numBigFaces) &&
                  (static_cast<int>(oPoints.size()           ) == numBigFaces) &&
                  (static_cast<int>(innerBdrys.size()        ) == numBigFaces) &&
                  (static_cast<int>(innerBdrysMassP.size()   ) == numBigFaces) &&
                  (static_cast<int>(oPointsConvexFlag.size() ) == numBigFaces),
                  "expecting the same number of big faces everywhere");

    gsDebug<<"mapping to plane..."<<"\n";
    std::vector<std::vector<Vertex> > iPoints2D;
    std::vector<std::vector<Vertex> > oPoints2D;
    std::vector< std::vector<std::vector<Vertex> > > innerBdrys2D;

    // Create pointers to gsVector3d out of oPoints so that the
    // gsCurveLoop constructor may be used
    for (size_t i=0;i<oPoints.size();i++) // for all patches (big faces, given as collection of vertices)
    {

        std::vector<gsVector3d<T>*> vertices;
        std::vector<VertexHandle> vertexVec=oPoints[i];

        std::vector<bool> isConvex=oPointsConvexFlag[i];
        for (size_t j=0;j<oPoints[i].size();j++)  // for all boundary points of i-th patch
        {
            vertices.push_back(&(*vertexVec[j]));
        }

        // ---- Construct the boundary curve loop out of the points
        // start construction
        gsVector3d<T> faceNormal;
        gsCurveLoop<T> loop(vertices,isConvex,0.01,&faceNormal);
        // check curve loop's normal is consistent with triangle normals,
        // try again with the alternate method if necessary.
        int nm = normalMult(faceNormal, face, (int)(i + 1));
        if(nm == 0)
        {
            bool success = loop.initFrom3DByAngles(vertices, isConvex, 0.01);
            // if THAT failed, last ditch effort which will always produce something
            // but it might not be very good.
            if(!success)
            {
                loop.initFromIsConvex(isConvex, 0.01);
            }
        }
        else if(nm == -1) loop.flip1(); // consistent once we flip
        // finish construction

        // Get the parameter pre-images ( 2D points - coefficients of the line segments of the parameter loops)
        std::vector<gsVertex<T> > vector2D;
        for(size_t j=0;j<loop.curves().size();j++)
        {
            vector2D.push_back(gsVertex<T>(loop.curve(j).coefs()(0,0),loop.curve(j).coefs()(0,1),0));
        }
        //calculate the inner 2D Points
        oPoints2D.push_back(vector2D);
    }
    gsDebug<<"generated loops"<<'\n';
    // calculate areas of the Patches - in case we need to add extra points
    std::vector<T > areas;
    for (int i=0;i<numBigFaces;i++)
        areas.push_back(0);
    for (size_t i=0;i<face.size();i++)
    {
        if(face[i]->faceIdentity>0&&face[i]->faceIdentity<=numBigFaces)
            areas[face[i]->faceIdentity-1]+=calcArea(face[i]);
    }

    T maxArea=0;
    for (size_t i=0;i<areas.size();i++)
        if(areas[i]>maxArea)
            maxArea=areas[i];
    // Store pointers to all edges of a face, in order to improve the fitting later
    std::vector< std::vector< EdgeHandle > > faceEdges;
    for (int i=0;i<numBigFaces;i++) // create container of edges for each patch
    {
        std::vector< EdgeHandle > helpVec;
        faceEdges.push_back(helpVec);
    }
    for(size_t i=0;i<edge.size();i++) // for all patches
    {
//        GISMO_ASSERT(edge[i].nFaces.size()==2,"each edge has to have 2 nFaces");
        if(edge[i].nFaces.size()==2 &&
                ((edge[i].nFaces[0]->faceIdentity)==(edge[i].nFaces[1]->faceIdentity))) // edge in the interior of the patch ?
        {
            faceEdges[(edge[i].nFaces[0]->faceIdentity)-1].push_back(&edge[i]);
        }
        else // edge on boundary ?
        {
            faceEdges[(edge[i].nFaces[0]->faceIdentity)-1].push_back(&edge[i]);
            if(edge[i].nFaces.size()==2)
            {
                faceEdges[(edge[i].nFaces[1]->faceIdentity)-1].push_back(&edge[i]);
            }
        }
    }
    std::vector<bool> isCylinder; // true if a patch is cylindric shaped
    std::vector<T> improveCylinderParametrization; // true if a patch is cylindric shaped
    for(size_t i=0;i<oPoints.size();i++) // For all patches
    {
        if(innerBdrys[i].size()!=1) // different than one hole
        {
            isCylinder.push_back(0);
            improveCylinderParametrization.push_back(0);
        }
        else // Has one hole ? (then it is cylindrical)
        {
            T diameterOut=0; // Measure outer diameter
            for(size_t j=0;j<oPoints[i].size()-1;j++)
            {
                for(size_t k=j+1;k<oPoints[i].size();k++)
                {
                    T tempDiameter=((oPoints[i][j]->x()-oPoints[i][k]->x())*(oPoints[i][j]->x()-oPoints[i][k]->x()))+
                            ((oPoints[i][j]->y()-oPoints[i][k]->y())*(oPoints[i][j]->y()-oPoints[i][k]->y()))+
                            ((oPoints[i][j]->z()-oPoints[i][k]->z())*(oPoints[i][j]->z()-oPoints[i][k]->z()));
                    if(tempDiameter>diameterOut)
                        diameterOut=tempDiameter;
                }
            }
            T diameterIn=0; // measure inner diameter
            for(size_t j=0;j<innerBdrys[i][0].size()-1;j++)
            {
                for(size_t k=j+1;k<innerBdrys[i][0].size();k++)
                {
                    T tempDiameter=((innerBdrys[i][0][j]->x()-innerBdrys[i][0][k]->x())*(innerBdrys[i][0][j]->x()-innerBdrys[i][0][k]->x()))+
                            ((innerBdrys[i][0][j]->y()-innerBdrys[i][0][k]->y())*(innerBdrys[i][0][j]->y()-innerBdrys[i][0][k]->y()))+
                            ((innerBdrys[i][0][j]->z()-innerBdrys[i][0][k]->z())*(innerBdrys[i][0][j]->z()-innerBdrys[i][0][k]->z()));
                    if(tempDiameter>diameterIn)
                        diameterIn=tempDiameter;
                }
            }
            T distance=0; // measure height of the cylinder
            for(size_t j=0;j<oPoints[i].size()-1;j++)
            {
                for(size_t k=0;k<innerBdrys[i][0].size();k++)
                {
                    T tempDist=((oPoints[i][j]->x()-innerBdrys[i][0][k]->x())*(oPoints[i][j]->x()-innerBdrys[i][0][k]->x()))+
                            ((oPoints[i][j]->y()-innerBdrys[i][0][k]->y())*(oPoints[i][j]->y()-innerBdrys[i][0][k]->y()))+
                            ((oPoints[i][j]->z()-innerBdrys[i][0][k]->z())*(oPoints[i][j]->z()-innerBdrys[i][0][k]->z()));
                    if(tempDist<distance||distance==0)
                        distance=tempDist;
                }
            }
            if(diameterOut/4<distance&&diameterIn/4<distance) // Is this patch a cylinder ??
            {
                isCylinder.push_back(1);
                improveCylinderParametrization.push_back(0.1);// scalar factor for Floater's weight for artificial vertex inside the hole
            }

            else
            {
                isCylinder.push_back(0);
                improveCylinderParametrization.push_back(0);
            }
        }
    }

    // At this point, we have a parametrized boundary point cloud
    // for each patch, and we go on to obtain a parameterization
    // of the inner points as well

    gsDebug<<"allocating surfaces"<<'\n';

    //allocating the trimmed surfaces
    std::vector<gsTrimSurface<T> *> tSurfVec;
    for (size_t i=0;i<iPoints.size();i++) // run through all the patches (big faces)
    {
        size_t iPsize=iPoints[i].size();
        // will contain: TOTAL number of points to map for this
        // patch ( inner boundary/hole points, artificial points (mass points)
        // and interior (non boundary) point cloud inside the face
        size_t n = iPsize+innerBdrysMassP[i].size();
        size_t innerBdrysSize=0;
        for (size_t m=0;m<innerBdrys[i].size();m++)
        {
            innerBdrysSize+=innerBdrys[i][m].size();
        }
        n+=innerBdrysSize;

        // Construct the 3 important point clouds
        std::set<VertexHandle> vertexFaceSet; // union of vertexFaceSetBdry and vertexFaceSetInner
        std::set<VertexHandle> vertexFaceSetBdry; // outer boundary points
        std::set<VertexHandle> vertexFaceSetInner; // all rest of points apart from outer boundary

        // --- start Fill up point sets
        for (size_t j=0;j<iPoints[i].size();j++)
        {
            vertexFaceSet.insert(iPoints[i][j]);
            vertexFaceSetInner.insert(iPoints[i][j]);
        }
        for (size_t j=0;j<innerBdrys[i].size();j++)
        {
            vertexFaceSet.insert(&innerBdrysMassP[i][j]);
            vertexFaceSetInner.insert(&innerBdrysMassP[i][j]);
            for (size_t k=0;k<innerBdrys[i][j].size();k++)
            {
                vertexFaceSet.insert(innerBdrys[i][j][k]);
                vertexFaceSetInner.insert(innerBdrys[i][j][k]);
            }
        }
        for (size_t j=0;j<oPoints[i].size();j++)
        {
            vertexFaceSet.insert(oPoints[i][j]);
            vertexFaceSetBdry.insert(oPoints[i][j]);
        }
        // --- end Fill up point sets

        // ----------------------- Start Floater's algorithm
        //set up a linear system of equations and solve it in
        //order to receive a parametrization of the inner points
        //of a face
        // u : u-parameter values
        // v : v-parameter values
        // b1, b2 : right hand sides of the linera system ( sums of boundary veterx weights or zeros)
        gsVector<T> u(n) ,v(n) ,b1(n) ,b2(n);
        // Linear system
        Eigen::SparseMatrix<T, ColMajor> A(n,n);
        gsSparseEntries<T> coefficients;
        typename gsSparseSolver<T>::LU solver;

        // idea: Pre-define a map (std::map) from the
        // vertexhandles to 0...n+k-1, also containing additional
        // info such as: for instance is it artificial THEN from
        // vertexhandle we can get the column index, and assemble the
        // matrix and rhs accordingly without searching every time

        for (size_t j=0;j<n;j++) //run through all inner points of a single face -- Rows of matrix A
        {
            // initialize rhs vector position
            b1(j)=0;
            b2(j)=0;
            coefficients.add(j,j,1.0); // A has ones on the diagonal

            // --- Local coefficients contributing to point j in the matrix and rhs
            std::vector<T> rhsCoefs; // weights of the boundary points which are connected to point j
            std::vector<VertexHandle> rhs; // vertices of right hand side (correst. rhsCoefs)
            std::vector<T> matCoefs; // weights of the Interiod points which are connected to point j
            std::vector<VertexHandle> mat; // vertices indexing columns of mat (correst. matCoefs)

            T check=0;//make sure the weights sum up to 1
            T normCoef=0; // normalization coefficient, ensures sum of weights at j is equal to 1
            T weight=0;// temporary
            if (j<iPsize) // is j an interior point ?
            {
                for (size_t k=0;k<iPoints[i][j]->nVertices.size();k++) // searching neighbors of point j
                {
                    // is it a boundary neighbor ?
                    if(vertexFaceSetBdry.find(iPoints[i][j]->nVertices[k])!=vertexFaceSetBdry.end())
                    {
                        weight=calcWeight(iPoints[i][j],iPoints[i][j]->nVertices[k],vertexFaceSet);
                        rhsCoefs.push_back(weight);
                        rhs.push_back(iPoints[i][j]->nVertices[k]);
                        normCoef+=weight;
                    }
                    else // is it an interior neighbor
                    {
                        weight=calcWeight(iPoints[i][j],iPoints[i][j]->nVertices[k],vertexFaceSet);
                        matCoefs.push_back(weight);
                        mat.push_back(iPoints[i][j]->nVertices[k]);
                        normCoef+=weight;
                    }
                }
            }
            else if(j<iPsize+innerBdrysMassP[i].size()) // is j an artificial point ?
            {
                for (size_t k=0;k<innerBdrys[i][j-iPsize].size();k++) // searching neighors of point j
                {
                    // get weight of artificial point and put it in the matrix
                    weight=calcWeight(&innerBdrysMassP[i][j-iPsize],innerBdrys[i][j-iPsize][k],vertexFaceSet);
                    matCoefs.push_back(weight);
                    mat.push_back(innerBdrys[i][j-iPsize][k]);
                    normCoef+=weight;
                }
            }
            else  // is it a point on a hole ? (inner boundary)
            {
                // allocate innerbdry points, first all neighbors, last inner point

                // j ----> k,l (: hole number, point number on k hole)

                //derive correct indices
                size_t innerIndex=j-iPsize-innerBdrysMassP[i].size();
                size_t c=0;
                size_t k=0;
                size_t l=0;
                while (innerIndex>=c)
                {
                    c+=innerBdrys[i][k].size();
                    k++;
                }
                k--;
                c-=innerBdrys[i][k].size();
                l=innerIndex-c;
                // searching neighors of point j
                for (size_t i2=0;i2<innerBdrys[i][k][l]->nVertices.size();i2++)
                {
                    // is it a boundary neighbor ?
                    if(vertexFaceSetBdry.find(innerBdrys[i][k][l]->nVertices[i2])!=vertexFaceSetBdry.end())
                    {
                        weight=calcWeight(innerBdrys[i][k][l],(innerBdrys[i][k][l]->nVertices[i2]),vertexFaceSet);
                        rhsCoefs.push_back(weight);
                        rhs.push_back(innerBdrys[i][k][l]->nVertices[i2]);
                        normCoef+=weight;
                    }
                    // is it an interior neighbor
                    else if(vertexFaceSetInner.find(innerBdrys[i][k][l]->nVertices[i2])!=vertexFaceSetInner.end())
                    {
                        weight=calcWeight(innerBdrys[i][k][l],(innerBdrys[i][k][l]->nVertices[i2]),vertexFaceSet);
                        // treating cylindrical patch
                        if(isCylinder[i]==1&&*(innerBdrys[i][k][l]->nVertices[i2])==innerBdrysMassP[i][0])
                        {
                            weight*=improveCylinderParametrization[i];
                        }
                        matCoefs.push_back(weight);
                        mat.push_back(innerBdrys[i][k][l]->nVertices[i2]);
                        normCoef+=weight;
                    }
                }
            }

            GISMO_ASSERT(normCoef > 0, "normCoef should be positive");

            for (size_t k=0;k<rhsCoefs.size();k++) // for all boundary connections of point j
            {
                rhsCoefs[k]=rhsCoefs[k]/normCoef; // normalize coefficient
                check+=rhsCoefs[k];
                int l=0;
                while (*oPoints[i][l]!=*rhs[k]) // locate outer boundary neighbors of point j
                {
                    l++;
                }
                // Fill in right hand side b
                b1(j)+=(oPoints2D[i][l].operator[](0))*rhsCoefs[k];
                b2(j)+=(oPoints2D[i][l].operator[](1))*rhsCoefs[k];

            }

            for (size_t k=0;k<matCoefs.size();k++) // for all inner/artificial/inner boundary connections of point j
            {
                matCoefs[k]=matCoefs[k]/normCoef;  // normalize coefficient
                check+=matCoefs[k];

                //-- start Locating NON-outer boundary neighbors of j-th point
                for (size_t l=0;l<iPsize;l++) // search in inner points
                {
                    if (*mat[k]==*iPoints[i][l])
                        coefficients.add(j,l,-matCoefs[k]);
                }

                for (size_t l=0;l<innerBdrysMassP[i].size();l++) // search in artifical points
                {
                    if(*mat[k]==innerBdrysMassP[i][l])
                        coefficients.add(j,l+iPsize,-matCoefs[k]);
                }

                int l=0;
                for (size_t i2=0;i2<innerBdrys[i].size();i2++) // search in inner boundaries (holes)
                {
                    for(size_t i3=0;i3<innerBdrys[i][i2].size();i3++)
                    {
                        if(*mat[k]==*innerBdrys[i][i2][i3])
                            coefficients.add(j,l+iPsize+innerBdrysMassP[i].size(),-matCoefs[k]);

                        l++;
                    }
                }
                //-- end Locating NON-outer boundary neighbors of j-th point
            }
            if(check<1-0.001||check>1+0.001)
                gsWarn<<"something might have gone wrong with norming: "<<check<<"!=1\n";
        }
        //build A from the entries
        A.setFromTriplets(coefficients.begin(), coefficients.end());
        //gsDebug<<A<<'\n';

        if(A.rows()!=0) // If there are interior points
        {
            // Compute the ordering permutation vector from the structural pattern of A
            solver.analyzePattern(A);
            // Compute the numerical factorization
            solver.factorize(A);
            //Use the factors to solve the linear system
            u = solver.solve(b1);
            v = solver.solve(b2);
        }
        // ----------------------- End Floater's algorithm

        // ----------------------- Start Fitting

        std::vector<gsVertex<T> > helpvec;
        for (size_t j=0;j<iPsize;j++)
        {
            helpvec.push_back(gsVertex<T>(u(j),v(j),0));
        }
        iPoints2D.push_back(helpvec);
        std::vector<std::vector<gsVertex<T> > > holeVecFace;
        int index=innerBdrys[i].size();//ignore the virtual points in the middle of the inner boundarys index=#holes
        for(size_t j=0;j<innerBdrys[i].size();j++)
        {
            std::vector<gsVertex<T> > holeVec;
            for(size_t k=0;k<innerBdrys[i][j].size();k++)
            {
                holeVec.push_back(gsVertex<T>(u(iPsize+index),v(iPsize+index),0));
                index++;
            }
            holeVecFace.push_back(holeVec);
        }
        innerBdrys2D.push_back(holeVecFace);

        //altering data to gsMatrix-form, in order to use a fitting function
        int nCorners=0;
        for (size_t j=0;j<oPoints[i].size();j++)
        {
            if(oPoints[i][j]->numEdges>2)
            {
                nCorners++;
            }
        }
        for (size_t j=0;j<innerBdrys[i].size();j++)
        {
            for(size_t k=0;k<innerBdrys[i][j].size();k++)
            {
                if(innerBdrys[i][j][k]->numEdges>2)
                {
                    nCorners++;
                }
            }
        }
        gsMatrix<T> Corners2d(2,nCorners);
        gsMatrix<T> Corners3d(3,nCorners);
        int nEdgePts=(oPoints[i].size()+innerBdrysSize)*(closeBoundary+1)-nCorners;

        gsMatrix<T> EdgePts2d(2,nEdgePts);
        gsMatrix<T> EdgePts3d(3,nEdgePts);
        int corNum=0;
        int edgNum=0;
        for (size_t j=0;j<oPoints[i].size();j++)
        {

            if(oPoints[i][j]->numEdges>2)
            {
                Corners2d(0,corNum)=oPoints2D[i][j].x();
                Corners2d(1,corNum)=oPoints2D[i][j].y();
                Corners3d(0,corNum)=oPoints[i][j]->x();
                Corners3d(1,corNum)=oPoints[i][j]->y();
                Corners3d(2,corNum)=oPoints[i][j]->z();
                corNum++;
            }
            else
            {
                EdgePts2d(0,edgNum)=oPoints2D[i][j].x();
                EdgePts2d(1,edgNum)=oPoints2D[i][j].y();
                EdgePts3d(0,edgNum)=oPoints[i][j]->x();
                EdgePts3d(1,edgNum)=oPoints[i][j]->y();
                EdgePts3d(2,edgNum)=oPoints[i][j]->z();
                edgNum++;
            }
        }
        for(size_t j=0;j<oPoints[i].size()-1;j++)
        {
            for(int k=0;k<closeBoundary;k++)
            {
                EdgePts2d(0,edgNum)=oPoints2D[i][j].x()*(k+1)/(closeBoundary+1)+oPoints2D[i][j+1].x()*(closeBoundary-k)/(closeBoundary+1);
                EdgePts2d(1,edgNum)=oPoints2D[i][j].y()*(k+1)/(closeBoundary+1)+oPoints2D[i][j+1].y()*(closeBoundary-k)/(closeBoundary+1);
                EdgePts3d(0,edgNum)=oPoints[i][j]->x()*(k+1)/(closeBoundary+1)+oPoints[i][j+1]->x()*(closeBoundary-k)/(closeBoundary+1);
                EdgePts3d(1,edgNum)=oPoints[i][j]->y()*(k+1)/(closeBoundary+1)+oPoints[i][j+1]->y()*(closeBoundary-k)/(closeBoundary+1);
                EdgePts3d(2,edgNum)=oPoints[i][j]->z()*(k+1)/(closeBoundary+1)+oPoints[i][j+1]->z()*(closeBoundary-k)/(closeBoundary+1);
                edgNum++;
            }
        }
        for(int k=0;k<closeBoundary;k++)
        {
            EdgePts2d(0,edgNum)=oPoints2D[i][0].x()*(k+1)/(closeBoundary+1)+oPoints2D[i].back().x()*(closeBoundary-k)/(closeBoundary+1);
            EdgePts2d(1,edgNum)=oPoints2D[i][0].y()*(k+1)/(closeBoundary+1)+oPoints2D[i].back().y()*(closeBoundary-k)/(closeBoundary+1);
            EdgePts3d(0,edgNum)=oPoints[i][0]->x()*(k+1)/(closeBoundary+1)+oPoints[i].back()->x()*(closeBoundary-k)/(closeBoundary+1);
            EdgePts3d(1,edgNum)=oPoints[i][0]->y()*(k+1)/(closeBoundary+1)+oPoints[i].back()->y()*(closeBoundary-k)/(closeBoundary+1);
            EdgePts3d(2,edgNum)=oPoints[i][0]->z()*(k+1)/(closeBoundary+1)+oPoints[i].back()->z()*(closeBoundary-k)/(closeBoundary+1);
            edgNum++;
        }

        int nInteriorPts=iPsize;
        //calculate number of additional inner Points
        std::multiset<T> edgeLengths;
        std::vector<int> nAddIntPointsPerEdge;
        int nAddIntPoints=0;
        if(moreInner)
        {
            T h=0;
            for(size_t j=0;j<faceEdges[i].size();j++)
            {
                edgeLengths.insert(calcDist(faceEdges[i][j]->source,faceEdges[i][j]->target));
            }
            //for(typename std::multiset<T >::iterator it=edgeLengths.begin();j!=edgeLengths.size()/10;it++)
            typename std::multiset<T >::iterator it=edgeLengths.begin();
            for(size_t j=0;j<=edgeLengths.size()/10;j++)
            {
                h=*it;
                it++;
            }
            for(size_t j=0;j<faceEdges[i].size();j++)
            {
                T l=calcDist(faceEdges[i][j]->source,faceEdges[i][j]->target);
                int add= math::max( cast<T,int>(l/h)-1, 0 );
                nAddIntPointsPerEdge.push_back(add);
                nAddIntPoints+=add;
            }
        }
        gsMatrix<T> interiorPts2d(2,nInteriorPts+nAddIntPoints);
        gsMatrix<T> interiorPts3d(3,nInteriorPts+nAddIntPoints);
        for (size_t j=0;j<iPoints2D[i].size();j++)
        {
            {
                interiorPts2d(0,j)=iPoints2D[i][j].x();
                interiorPts2d(1,j)=iPoints2D[i][j].y();
                interiorPts3d(0,j)=iPoints[i][j]->x();
                interiorPts3d(1,j)=iPoints[i][j]->y();
                interiorPts3d(2,j)=iPoints[i][j]->z();
            }
        }
        int intNum=iPsize;
        if(moreInner)
        {
            for(size_t j=0;j<faceEdges[i].size();j++)
            {
                if(nAddIntPointsPerEdge[j]>0)
                {
                    VertexHandle v1=NULL;
                    VertexHandle v2=NULL;
                    VertexHandle v1_2d=NULL;
                    VertexHandle v2_2d=NULL;

                    //find source and target of edge
                    for(size_t k=0;k<oPoints[i].size();k++)
                    {
                        if(oPoints[i][k]==faceEdges[i][j]->source)
                        {
                            v1=oPoints[i][k];
                            v1_2d=&oPoints2D[i][k];
                        }
                        else if(oPoints[i][k]==faceEdges[i][j]->target)
                        {
                            v2=oPoints[i][k];
                            v2_2d=&oPoints2D[i][k];
                        }
                    }
                    for(size_t k=0;k<iPoints[i].size();k++)
                    {
                        if(iPoints[i][k]==faceEdges[i][j]->source)
                        {
                            v1=iPoints[i][k];
                            v1_2d=&iPoints2D[i][k];
                        }
                        else if(iPoints[i][k]==faceEdges[i][j]->target)
                        {
                            v2=iPoints[i][k];
                            v2_2d=&iPoints2D[i][k];
                        }
                    }
                    for(size_t k=0;k<innerBdrys[i].size();k++)
                    {
                        for(size_t l=0;l<innerBdrys[i][k].size();l++)
                        {
                            if(innerBdrys[i][k][l]==faceEdges[i][j]->source)
                            {
                                v1=innerBdrys[i][k][l];
                                v1_2d=&innerBdrys2D[i][k][l];
                            }
                            else if(innerBdrys[i][k][l]==faceEdges[i][j]->target)
                            {
                                v2=innerBdrys[i][k][l];
                                v2_2d=&innerBdrys2D[i][k][l];
                            }
                        }
                    }
                    GISMO_ASSERT(v1!=NULL&&v2!=NULL,"could not find source or target of the edge in the face");
                    for (int k=0;k<nAddIntPointsPerEdge[j];k++)
                    {
                        interiorPts2d(0,intNum)=v1_2d->x()*(k+1)/(nAddIntPointsPerEdge[j]+1)+
                                v2_2d->x()*(nAddIntPointsPerEdge[j]-k)/(nAddIntPointsPerEdge[j]+1);
                        interiorPts2d(1,intNum)=v1_2d->y()*(k+1)/(nAddIntPointsPerEdge[j]+1)+
                                v2_2d->y()*(nAddIntPointsPerEdge[j]-k)/(nAddIntPointsPerEdge[j]+1);
                        interiorPts3d(0,intNum)=v1->x()*(k+1)/(nAddIntPointsPerEdge[j]+1)+
                                v2->x()*(nAddIntPointsPerEdge[j]-k)/(nAddIntPointsPerEdge[j]+1);
                        interiorPts3d(1,intNum)=v1->y()*(k+1)/(nAddIntPointsPerEdge[j]+1)+
                                v2->y()*(nAddIntPointsPerEdge[j]-k)/(nAddIntPointsPerEdge[j]+1);
                        interiorPts3d(2,intNum)=v1->z()*(k+1)/(nAddIntPointsPerEdge[j]+1)+
                                v2->z()*(nAddIntPointsPerEdge[j]-k)/(nAddIntPointsPerEdge[j]+1);
                        intNum++;
                    }

                }
            }
        }
        for(size_t j=0;j<innerBdrys[i].size();j++)
        {
            for(size_t k=0;k<innerBdrys[i][j].size();k++)
            {
                if(innerBdrys[i][j][k]->numEdges>2)
                {
                    Corners2d(0,corNum)=innerBdrys2D[i][j][k].x();
                    Corners2d(1,corNum)=innerBdrys2D[i][j][k].y();
                    Corners3d(0,corNum)=innerBdrys[i][j][k]->x();
                    Corners3d(1,corNum)=innerBdrys[i][j][k]->y();
                    Corners3d(2,corNum)=innerBdrys[i][j][k]->z();
                    corNum++;
                }
                else
                {
                    EdgePts2d(0,edgNum)=innerBdrys2D[i][j][k].x();
                    EdgePts2d(1,edgNum)=innerBdrys2D[i][j][k].y();
                    EdgePts3d(0,edgNum)=innerBdrys[i][j][k]->x();
                    EdgePts3d(1,edgNum)=innerBdrys[i][j][k]->y();
                    EdgePts3d(2,edgNum)=innerBdrys[i][j][k]->z();
                    edgNum++;
                }
            }
        }
        for(size_t it=0;it<innerBdrys[i].size();it++)
        {
            for(size_t j=0;j<innerBdrys[i][it].size()-1;j++)
            {
                for(int k=0;k<closeBoundary;k++)
                {
                    EdgePts2d(0,edgNum)=innerBdrys2D[i][it][j].x()*(k+1)/(closeBoundary+1)+innerBdrys2D[i][it][j+1].x()*(closeBoundary-k)/(closeBoundary+1);
                    EdgePts2d(1,edgNum)=innerBdrys2D[i][it][j].y()*(k+1)/(closeBoundary+1)+innerBdrys2D[i][it][j+1].y()*(closeBoundary-k)/(closeBoundary+1);
                    EdgePts3d(0,edgNum)=innerBdrys[i][it][j]->x()*(k+1)/(closeBoundary+1)+innerBdrys[i][it][j+1]->x()*(closeBoundary-k)/(closeBoundary+1);
                    EdgePts3d(1,edgNum)=innerBdrys[i][it][j]->y()*(k+1)/(closeBoundary+1)+innerBdrys[i][it][j+1]->y()*(closeBoundary-k)/(closeBoundary+1);
                    EdgePts3d(2,edgNum)=innerBdrys[i][it][j]->z()*(k+1)/(closeBoundary+1)+innerBdrys[i][it][j+1]->z()*(closeBoundary-k)/(closeBoundary+1);
                    edgNum++;
                }
            }
            for(int k=0;k<closeBoundary;k++)
            {
                EdgePts2d(0,edgNum)=innerBdrys2D[i][it][0].x()*(k+1)/(closeBoundary+1)+innerBdrys2D[i][it].back().x()*(closeBoundary-k)/(closeBoundary+1);
                EdgePts2d(1,edgNum)=innerBdrys2D[i][it][0].y()*(k+1)/(closeBoundary+1)+innerBdrys2D[i][it].back().y()*(closeBoundary-k)/(closeBoundary+1);
                EdgePts3d(0,edgNum)=innerBdrys[i][it][0]->x()*(k+1)/(closeBoundary+1)+innerBdrys[i][it].back()->x()*(closeBoundary-k)/(closeBoundary+1);
                EdgePts3d(1,edgNum)=innerBdrys[i][it][0]->y()*(k+1)/(closeBoundary+1)+innerBdrys[i][it].back()->y()*(closeBoundary-k)/(closeBoundary+1);
                EdgePts3d(2,edgNum)=innerBdrys[i][it][0]->z()*(k+1)/(closeBoundary+1)+innerBdrys[i][it].back()->z()*(closeBoundary-k)/(closeBoundary+1);
                edgNum++;
            }
        }
        gsMatrix<T> appxNormalPoints(2,0);
        gsMatrix<T> appxNormals(3,0);
        int innerPts=cast<T,int>(math::sqrt(cast<int,T>(nCorners+kvAdditionalInnerPoints)));
        if (innerPts<0)
            innerPts=0;

        gsKnotVector<T> kv1 (0, 1, innerPts, kvOuterPoints) ;

        typename gsTensorBSpline<2,T>::Ptr spline = gsInterpolateSurface(
                Corners2d, Corners3d,
                EdgePts2d, EdgePts3d,
                interiorPts2d, interiorPts3d,
                appxNormalPoints, appxNormals,
                T(wE),  T(wI),  T(0),  T(0.01),
                kv1,kv1,
                false
                );
        // ----------------------- End Fitting

        // ----------- start Calculate trimmed surfaces
        std::vector< gsCurveLoop<T> *> loops;
        gsCurveLoop<T> * loop = calculateLoop(oPoints2D[i], isCorner(oPoints[i]), noSmooth);

        loops.push_back(loop);
        for (size_t j=0;j<innerBdrys2D[i].size();j++)
        {
            gsCurveLoop<T> * innerloop = calculateLoop(innerBdrys2D[i][j], isCorner(innerBdrys[i][j]), noSmooth);
            loops.push_back(innerloop);
        }
        gsPlanarDomain<T> * domain=  new gsPlanarDomain<T>(loops);
        gsTrimSurface<T> * cface= new gsTrimSurface<T> (spline,domain);
        tSurfVec.push_back(cface);

        if(plot) // for debugging ( plot trimmed surfaces
        {
            typename gsMesh<T>::uPtr m;
            int nPoints=cast<T,int>(
                meshPoints*math::sqrt(math::sqrt(areas[i]/maxArea)));
            if (nPoints<10)
                nPoints=10;
            m = cface->toMesh(nPoints);
            fitMeshes.push_back(m.release());
        }

        // ----------- end Calculate trimmed surfaces

        //delete cface;

    } // End main loop over patches / big faces

    //allocate  gsSolid
    //std::vector<gsVertex<T> > vertVec; // ---sorted
    gsSortedVector<gsVertex<T> > vertVec;
    for (size_t i=0;i<oPoints.size();i++)
    {
        for (size_t j=0;j<oPoints[i].size();j++)
        {
            //insertObject(*oPoints[i][j],vertVec);
            vertVec.push_sorted_unique(*oPoints[i][j]);
        }
    }
    for (size_t i=0;i<innerBdrys.size();i++)
    {
        for (size_t j=0;j<innerBdrys[i].size();j++)
        {
            for(size_t k=0;k<innerBdrys[i][j].size();k++)
                //insertObject(*innerBdrys[i][j][k],vertVec);
                vertVec.push_sorted_unique(*innerBdrys[i][j][k]);
        }
    }
    for (size_t i=0;i<vertVec.size();i++)
        sl.addHeVertex(vertVec[i][0],vertVec[i][1],vertVec[i][2]);
    for (size_t i=0;i<oPoints.size();i++)
    {
        std::vector<std::vector<gsSolidHeVertex<T>* > > faceConstruct;
        std::vector<gsSolidHeVertex<T>* > faceHelp;
        for (size_t j=0;j<oPoints[i].size();j++)
        {
            //faceHelp.push_back(sl.vertex[getIndex(*oPoints[i][j],vertVec)]);
            faceHelp.push_back( sl.vertex[ vertVec.getIndex(*oPoints[i][j])] );
        }
        faceConstruct.push_back(faceHelp);
        for (size_t j=0;j<innerBdrys[i].size();j++)
        {
            faceHelp.clear();
            for (size_t k=0;k<innerBdrys[i][j].size();k++)
                //faceHelp.push_back(sl.vertex[getIndex(*innerBdrys[i][j][k],vertVec)]);
                faceHelp.push_back( sl.vertex[ vertVec.getIndex(*innerBdrys[i][j][k])] );
            faceConstruct.push_back(faceHelp);
        }
        sl.addFace( faceConstruct, tSurfVec[i]);

    }

    // Now: only members 'mate' of halfedge is missing, we will set them up now

    sl.setHeMate();
    // Add volume

    sl.addVolume(sl.face);
    // test what the faces look like
    if(plot) // for debugging, plot original meshes
    {
        for(int facenumber=1;facenumber<=numBigFaces;facenumber++) // for all patches
        {

            gsMesh<T> * mface = new gsMesh<T>();
            for (size_t it=0;it<this->face.size();it++)
            {
                if(this->face[it]->faceIdentity==facenumber)
                {

                    mface->addVertex(this->face[it]->vertices[0]->operator[](0),
                            this->face[it]->vertices[0]->operator[](1),
                            this->face[it]->vertices[0]->operator[](2));

                    mface->addVertex(this->face[it]->vertices[1]->operator[](0),
                            this->face[it]->vertices[1]->operator[](1),
                            this->face[it]->vertices[1]->operator[](2));

                    mface->addVertex(this->face[it]->vertices[2]->operator[](0),
                            this->face[it]->vertices[2]->operator[](1),
                            this->face[it]->vertices[2]->operator[](2));
                }
            }

            for (size_t i=0;i!=mface->numVertices();i=i+3)
            {
                (mface)->addFace(&mface->vertex(i), &mface->vertex(i+1), &mface->vertex(i+2));
            }
            patchMeshes.push_back(mface);


            gsMesh<T> * paraface=new gsMesh<T>();
            for(size_t it=0;it<this->face.size();it++)
            {
                if(this->face[it]->faceIdentity==facenumber)
                {
                    for(size_t vertIt=0;vertIt<3;vertIt++)
                    {
                        bool found=0;
                        for(size_t j=0;j<oPoints[facenumber-1].size();j++)
                        {
                            if(*this->face[it]->vertices[vertIt]==*oPoints[facenumber-1][j]&&found==0)
                            {

                                paraface->addVertex(oPoints2D[facenumber-1][j][0],oPoints2D[facenumber-1][j][1],0);
                                found=1;
                            }
                        }
                        for(size_t j=0;j<iPoints[facenumber-1].size();j++)
                        {
                            if(*this->face[it]->vertices[vertIt]==*iPoints[facenumber-1][j]&&found==0)
                            {

                                paraface->addVertex(iPoints2D[facenumber-1][j][0],iPoints2D[facenumber-1][j][1],0);
                                found=1;

                            }
                        }
                        for(size_t j=0;j<innerBdrys[facenumber-1].size();j++)
                        {
                            for(size_t k=0;k<innerBdrys[facenumber-1][j].size();k++)
                            {
                                if(*this->face[it]->vertices[vertIt]==*innerBdrys[facenumber-1][j][k]&&found==0)
                                {

                                    paraface->addVertex(innerBdrys2D[facenumber-1][j][k][0],innerBdrys2D[facenumber-1][j][k][1],0);
                                    found=1;
                                }
                            }
                        }
                    }
                }
            }

            for (size_t i=0; i!=paraface->numVertices(); i+=3)
            {
                (paraface)->addFace(&paraface->vertex(i), &paraface->vertex(i+1), &paraface->vertex(i+2));
            }

            paraMeshes.push_back(paraface);
        }
    }
}

template <class T>
void gsTriMeshToSolid<T>::getPatchData(T angle, T innerAngle,T patchAreaWeight,T mergeSmallPatches,
                        std::vector<std::vector<VertexHandle> > & iPoints,
                        std::vector<std::vector<VertexHandle> >  & oPoints,
                        std::vector< std::vector<std::vector<VertexHandle> > > & innerBdrys,
                        std::vector< std::vector<Vertex>  > & innerBdrysMassP,

                        std::vector<std::vector<bool> > & oPointsConvexFlag,
                        std::string filenameFeatures,
                        int useFeatures
                        )
{
    bool non_manifold, warning_borders;

    // compute the features
    this->getFeatures(angle, non_manifold, warning_borders);
    this->mesh->cleanMesh();

    // read features from a file, if required
    std::vector<gsEdge<T> > featEdges;
    if (useFeatures!=0)
        this->readEdges(filenameFeatures, featEdges);
    this->setSharpEdges(featEdges, useFeatures);

    // give every face a patch number
    this->calcPatchNumbers();

    // improve quality by further dividing and merging patches according to more complex rules
    storeNeighboringFaces();
    divideAndMergePatches(innerAngle, patchAreaWeight, mergeSmallPatches);

    // recompute patch numbers
    this->calcPatchNumbers();

    this->getFaces(iPoints,oPoints,innerBdrys,innerBdrysMassP,oPointsConvexFlag);

}


template <class T>
gsCurveLoop<T> * gsTriMeshToSolid<T>::calculateLoop(std::vector<Vertex> outerPoints, std::vector<bool > const & isCorner, bool noSmooth)
{
    GISMO_ASSERT(outerPoints.size()==isCorner.size(),"the vertices have to be defines as corners or edges by a vector of bools");
    gsCurveLoop<T> * loop = new gsCurveLoop<T>();
        //handle first point
        if (noSmooth || isCorner.at(0)==1)
        {
            if(noSmooth || isCorner.at(1)==1)
            {
                gsBSpline<T> * tcurve = calcTCurve(outerPoints[0],outerPoints[1]);
                loop->insertCurve(tcurve);
            }
        }
        else
        {
            Vertex v1(0,0,0);
            Vertex v2=outerPoints[0];
            Vertex v3(0,0,0);
            if(isCorner.back()==1)
                v1=outerPoints.back();
            else
                v1=giveMidpoint(outerPoints.back(),outerPoints[0]);
            if(isCorner.at(1)==1)
                v3=outerPoints[1];
            else
                v3=giveMidpoint(outerPoints[0],outerPoints[1]);
            gsBSpline<T> * tcurve = calcTCurve(v1,v2,v3);
            loop->insertCurve(tcurve);
        }
        //handle points in between
        for (size_t i=1;i<outerPoints.size()-1;i++)
        {
            if (noSmooth || isCorner.at(i)==1)
            {
                if(noSmooth || isCorner.at(i+1)==1)
                {
                    gsBSpline<T> * tcurve = calcTCurve(outerPoints[i],outerPoints[i+1]);
                    loop->insertCurve(tcurve);
                }
            }
            else
            {
                Vertex v1(0,0,0);
                Vertex v2=outerPoints[i];
                Vertex v3(0,0,0);
                if(isCorner.at(i-1)==1)
                    v1=outerPoints[i-1];
                else
                    v1=giveMidpoint(outerPoints[i-1],outerPoints[i]);
                if(isCorner.at(i+1)==1)
                    v3=outerPoints[i+1];
                else
                    v3=giveMidpoint(outerPoints[i],outerPoints[i+1]);
                gsBSpline<T> * tcurve = calcTCurve(v1,v2,v3);
                loop->insertCurve(tcurve);
            }
        }
        //handle last point
        if (noSmooth || isCorner.back()==1)
        {
            if(noSmooth || isCorner.at(0)==1)
            {
                gsBSpline<T> * tcurve = calcTCurve(outerPoints.back(),outerPoints[0]);
                loop->insertCurve(tcurve);
            }
        }
        else
        {
            Vertex v1(0,0,0);
            Vertex v2=outerPoints.back();
            Vertex v3(0,0,0);
            if(isCorner.at(isCorner.size()-2)==1)
                v1=outerPoints[isCorner.size()-2];
            else
                v1=giveMidpoint(outerPoints[isCorner.size()-2],outerPoints.back());
            if(isCorner.at(0)==1)
                v3=outerPoints[0];
            else
                v3=giveMidpoint(outerPoints.back(),outerPoints[0]);
            gsBSpline<T> * tcurve = calcTCurve(v1,v2,v3);
            loop->insertCurve(tcurve);
        }

    return loop;
}

template<class T>
T gsTriMeshToSolid<T>::calcWeight(VertexHandle v1,VertexHandle v2,
                    std::set<VertexHandle>
                    const & vertexFaceSet)
{
    T weight = 0;
    const gsVector3d<T> vec1 = *v2 - *v1;

    for(size_t i=0;i<v1->nVertices.size();i++)
    {
        for(size_t j=0;j<v2->nVertices.size();j++)
        {
            //check if neighboring vertices coincide and if point found is part of the face
            if ( *(v1->nVertices[i])==*(v2->nVertices[j]) &&
                 vertexFaceSet.find(v2->nVertices[j]) != vertexFaceSet.end()
                )
            {

                const gsVector3d<T> vec2 = *v2->nVertices[j] - *v1;
                weight+=math::tan(conditionedAngle( vec1,  vec2)/2);
            }
        }
    }
    weight /= calcDist(v1,v2);

    GISMO_ASSERT(weight > 0, "Weight should be positive");

    return weight;
}


template <class T>
int gsTriMeshToSolid<T>::normalMult(gsVector3d<T> globalNormal,
                        std::vector<FaceHandle> & face,
                        int bigFaceIdx)
{
    int result = 0;
    bool foundResult = false;
    // loop over faces
    size_t nf = face.size();
    // AM: A better way: use am "int count" variable to count the positive inner products
    // Then set the result if count==nf or count==0
    for(size_t i = 0; i < nf; i++)
    {
        // check it's a triangle
        GISMO_ASSERT(face[i]->vertices.size() == 3, "Expected triangle mesh");
        // only look at triangles inside the surface
        if(face[i]->faceIdentity != bigFaceIdx) continue;
        // compute the normal of this face
        const gsVector3d<T> v0 = *face[i]->vertices[0];
        const gsVector3d<T> off1 = *face[i]->vertices[1] - v0;
        const gsVector3d<T> off2 = *face[i]->vertices[2] - v0;
        const gsVector3d<T> thisNormal = off1.cross(off2);
        // compare against n and record
        T thisResult = thisNormal.dot(globalNormal);
        if((foundResult && thisResult * result <= 0) || thisResult == 0)
        {
            result = 0;
        }
        else result = (thisResult > 0)?1:-1;
        foundResult = true;
    }
    GISMO_ASSERT(foundResult, "Could not examine any triangles to get normal information");
    return result;
}

template<class T>
bool gsTriMeshToSolid<T>::approxEqual(const gsEdge<T> & e1,const gsEdge<T> & e2)
{
    const T epsilon= calcDist(e1.source, e1.target ) * 0.01 ;

    return ( (*e1.source - *e2.source).norm() < epsilon &&
             (*e1.target - *e2.target).norm() < epsilon );
    /*
      bool result=0;
      if((e1.source->x()>(e2.source->x()-epsilon))&&((e1.source->x()-epsilon)<e2.source->x())&&
      (e1.source->y()>(e2.source->y()-epsilon))&&((e1.source->y()-epsilon)<e2.source->y())&&
      (e1.source->z()>(e2.source->z()-epsilon))&&((e1.source->z()-epsilon)<e2.source->z())&&
      (e1.target->x()>(e2.target->x()-epsilon))&&((e1.target->x()-epsilon)<e2.target->x())&&
      (e1.target->y()>(e2.target->y()-epsilon))&&((e1.target->y()-epsilon)<e2.target->y())&&
      (e1.target->z()>(e2.target->z()-epsilon))&&((e1.target->z()-epsilon)<e2.target->z()))
      result=1;
      return result;
    */
}

template <class T>
T gsTriMeshToSolid<T>::calcAngle(EdgeHandle e1,EdgeHandle e2, int faceNum)
{
    VertexHandle anglePoint;
    if (e1->source==e2->source||e1->source==e2->target)
        anglePoint=e1->source;
    else if (e1->target==e2->target||e1->target==e2->source)
        anglePoint=e1->target;
    else
    {
        gsDebug<<"Edges are not neighbors"<<"\n";
        return 0;
    }
    gsVector3d<T> vec1(0,0,0);
    if (anglePoint==e1->source)
    {
        gsVector3d<T> helpvec(e1->target->operator[](0)-e1->source->operator[](0),e1->target->operator[](1)-e1->source->operator[](1),e1->target->operator[](2)-e1->source->operator[](2));
        vec1 = helpvec;
    }
    else
    {
        gsVector3d<T> helpvec(e1->source->operator[](0)-e1->target->operator[](0),e1->source->operator[](1)-e1->target->operator[](1),e1->source->operator[](2)-e1->target->operator[](2));
        vec1 = helpvec;
    }
    gsVector3d<T> vec2(0,0,0);
    if (anglePoint==e2->source)
    {
        gsVector3d<T> helpvec(e2->target->operator[](0)-e2->source->operator[](0),e2->target->operator[](1)-e2->source->operator[](1),e2->target->operator[](2)-e2->source->operator[](2));
        vec2 = helpvec;
    }
    else
    {
        gsVector3d<T> helpvec(e2->source->operator[](0)-e2->target->operator[](0),e2->source->operator[](1)-e2->target->operator[](1),e2->source->operator[](2)-e2->target->operator[](2));
        vec2=helpvec;
    }
    FaceHandle vec1Face=NULL;
    FaceHandle vec2Face=NULL;

    if (e1->nFaces[0]->faceIdentity==faceNum)
        vec1Face=e1->nFaces[0];
    else if (e1->nFaces[1]->faceIdentity==faceNum)
        vec1Face=e1->nFaces[1];
    else
        gsDebug<<"selected Edge has no valid neighboring face"<<"\n";
    if (e2->nFaces[0]->faceIdentity==faceNum)
        vec2Face=e2->nFaces[0];
    else if (e2->nFaces[1]->faceIdentity==faceNum)
        vec2Face=e2->nFaces[1];
    else
        gsDebug<<"selected Edge has no valid neighboring face"<<"\n";
    gsVector3d<T> normal1=vec1Face->orthogonalVector();
    gsVector3d<T> normal2=vec2Face->orthogonalVector();
    gsVector3d<T> normal=(normal1+normal2)/2;
    T angle=conditionedAngle(vec1,vec2,normal);

    return angle;

}


template<class T>
T gsTriMeshToSolid<T>::calcBdryLength(std::vector<VertexHandle > vec)
{
    T length=0;
    for (size_t i=0;i<vec.size()-1;i++)
        length+=calcDist(vec[i],vec[i+1]);
    length+=calcDist(vec[0],vec.back());
    return length;
}


template<class T>
gsVertex<T> gsTriMeshToSolid<T>::getMassP(std::vector<VertexHandle > vec)
{
    T x=0,y=0,z=0;
    for (size_t i=0;i<vec.size();i++)
    {
        x+=vec[i]->x();
        y+=vec[i]->y();
        z+=vec[i]->z();
    }
    x/=vec.size();
    y/=vec.size();
    z/=vec.size();
    return gsVertex<T>(x,y,z);
}


template<class T>
T gsTriMeshToSolid<T>::calcDist(VertexHandle v1,VertexHandle v2)
{
    return (*v1 - *v2).norm();
}


template<class T>
std::vector<bool> gsTriMeshToSolid<T>::isCorner(std::vector<VertexHandle > const & vertexVec3d)
{
    std::vector<bool> isCorner;
    for(size_t i=0;i<vertexVec3d.size();i++)
    {
        if(vertexVec3d[i]->numEdges>2)
            isCorner.push_back(1);
        else
            isCorner.push_back(0);
    }
    return isCorner;
}


template<class T>
gsBSpline<T> * gsTriMeshToSolid<T>::calcTCurve(Vertex v1,Vertex v2)
{
    gsMatrix<T> tcp(2, 2);
    tcp << v1.x(), v1.y(), v2.x(), v2.y();
    gsBSpline<T> * tcurve = new gsBSpline<T>( 0, 1, 0, 1, give(tcp) );
    return tcurve;
}


template<class T>
gsBSpline<T> * gsTriMeshToSolid<T>::calcTCurve(Vertex v1,Vertex v2,Vertex v3)
{
    gsMatrix<T> tcp(3, 2);
    tcp << v1.x(), v1.y(), v2.x(), v2.y(), v3.x(), v3.y();
    gsBSpline<T> * tcurve = new gsBSpline<T>( 0, 1, 0, 2, give(tcp) );
    return tcurve;
}


template<class T>
typename gsTriMeshToSolid<T>::Vertex gsTriMeshToSolid<T>::giveMidpoint(Vertex v1,Vertex v2)
{
    return gsVertex<T>((v1[0]+v2[0])/2,
                       (v1[1]+v2[1])/2,
                       (v1[2]+v2[2])/2);
}


template <class T>
void gsTriMeshToSolid<T>::readEdges( std::string const & fn, std::vector<gsEdge <T> > & edges )
{
    std::string line;
    std::ifstream myfile (fn.c_str());
    if (myfile.is_open())
    {
      int count=1;
      T x1=0,x2=0,y1=0,y2=0,z1=0,z2=0;
      while(myfile >> line)
      {

          if(count==1)
              x1 = atof(line.c_str());
          else if(count==2)
              y1 = atof(line.c_str());
          else if(count==3)
              z1 = atof(line.c_str());
          else if(count==4)
              x2 = atof(line.c_str());
          else if(count==5)
              y2 = atof(line.c_str());
          else if(count==6)
          {
              z2 = atof(line.c_str());
              gsVertex<T> * v1 = new gsVertex<T>(x1,y1,z1);
              gsVertex<T> * v2 = new gsVertex<T>(x2,y2,z2);

              if(Xless<T>(v1,v2))std::swap(v1,v2);
              gsEdge<T> featureEdge(v1,v2);
              edges.push_back(featureEdge);

              count-=6;
          }
          count++;
      }

      myfile.close();
    }
}


template <class T>
T gsTriMeshToSolid<T>::calcArea(FaceHandle f1)
{
    T d1=calcDist(f1->vertices[0],f1->vertices[1]);
    T d2=calcDist(f1->vertices[0],f1->vertices[2]);
    T d3=calcDist(f1->vertices[1],f1->vertices[2]);
    //p=perimeter/2
    T p=(d1+d2+d3)/2;
    T area=math::sqrt(p*(p-d1)*(p-d2)*(p-d3));
    return area;
}

}
