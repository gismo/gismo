/** @file gsPlanarDomain.hpp

    @brief Provides implementation of the gsPlanarDomain class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Falini, A. Mantzaflaris, D. Mayer, D.-M. Nguyen, M. Pauley
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsConstantFunction.h>
#include <gsCore/gsCurve.h>

#include <gsNurbs/gsNurbsCreator.h>

#include <gsUtils/gsMesh/gsMesh.h>

#include <gsModeling/gsTraceCurve.hpp>
#include <gsModeling/gsTemplate.h>

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{


template <class T>
gsPlanarDomain<T>::gsPlanarDomain( std::vector< gsCurveLoop<T> *> const & loops)
{
    if(!loops[0]->is_ccw())
        loops[0]->reverse();
    m_loops.push_back(loops[0]);
    for(size_t i=1; i<loops.size();i++ )
    {

        if( loops[i]->is_ccw() )

            loops[i]->reverse();
        m_loops.push_back(loops[i]);
    }
    updateBoundingBox();
}

/*
template<class T>
gsMatrix<T> gsPlanarDomain<T>::averageValue( std::vector<gsFunction<T>*> const &f,
                                             std::vector<T> const & breaks)
{
    gsMatrix<T> ngrid;
    gsVector<T> wgrid;
    gsMatrix<T> xev,B (f.size(),1);
    B.setZero();
    iteratedGaussRule(ngrid, wgrid, 2, breaks);
    for(size_t i=0; i<f.size();i++)
    {
        for(int k=0; k<ngrid.cols();k++) //for all quadrature points
        {
            const gsMatrix<T> & uc = ngrid.col(k);

            f[i]->eval_into(uc,xev);
            B(i,0)  += wgrid[k]*xev(0,0);
        }
    }
    gsDebugVar( B.transpose() );
    return B;
}
*/

template<class T>
bool gsPlanarDomain<T>::inDomain( gsMatrix<T> const & u, int direction)
{
    GISMO_ASSERT(u.cols() == 1, "Waiting for a single point");

    // compute intersections with line x=u(0,0);
    std::vector<T> tmp = m_loops[0]->lineIntersections(direction, u(direction,0) );

    //   gsDebug<< "intersections:\n";
    //   for( size_t i = 0; i!=tmp.size(); ++i)
    //     gsDebug<< " "<< tmp[i];
    //   gsDebug<<"\n";

    if ( tmp.empty() ) // point outside the outer loop
        return false;

    gsAsMatrix<T> xx(tmp);
    gsMatrix<T> e;
    m_loops[0]->curve(0).eval_into( xx, e );

    int count = 0; // count even means out of the domain

    //checking if abscissae of intersection are greater then point's abscissa
    for(index_t i=0; i!=e.cols(); i++)
    {
        if( (e( !direction, i )> u(!direction,0))  )
            count++;
    }
    if ( (count % 2) == 0 ) //this means point is outside the outer loop
        return false;

    for(index_t v = 1; v<this->numLoops(); v++)
    {

        count=0;
        tmp= m_loops[v]->lineIntersections( direction, u(direction,0) );

        if(tmp.size()!=0) // intersections detected with hole loop
        {
            gsAsMatrix<T> Ev(tmp);
            m_loops[v]->curve(0).eval_into( Ev, e );

            for(index_t i=0; i!=e.cols(); i++)
            {
                // gsDebug<<"e(1,i) "<< e(1,i)<<"\n";
                if( e(!direction,i) > u(!direction,0) )
                    count++;

            }
        }
        if( (count % 2) != 0 ) //this means point is inside hole v
            return false;
    }

    return  true; // if we get here then the point is in the domain
}


template<class T>
bool gsPlanarDomain<T>::onBoundary(gsMatrix<T> const & u)
{
    for(index_t v=0; v< this->numLoops();v++)
    {
        // true iff point u on boundary
        T parValue;
        if( m_loops[v]->isOn(u, parValue,1e-5 ) )
        {
            return true;
        }
    }
    return false;
}


template <class T>
std::ostream& gsPlanarDomain<T>::print(std::ostream &os) const
{
    os << "Outer boundary: " << *m_loops[0];
    if ( m_loops.size()>1 )
    {
        os << "Holes: ";
        for ( typename std::vector< gsCurveLoop<T> *>::const_iterator it =
              m_loops.begin()+1;  it != m_loops.end(); ++it)
            os << **it;
    }
    os << "Bounding box: ["<< m_bbox.col(0).transpose()<< ", "
       << m_bbox.col(1).transpose() << "]";
    return os;
}


/// linearly discriti
template <class T>
void gsPlanarDomain<T>::sampleLoop_into( int loopID, int npoints, int numEndPoints, gsMatrix<T> & u )
{
    assert( (loopID>=0) && (loopID < numLoops()) );

    int np; // new number of points
    switch (numEndPoints)
    {
    case (0):
        np = npoints-2;
        break;
    case (1):
        np = npoints-1;
        break;
    case (2):
        np = npoints;
        break;
    default:
        np = 0;
        break;
    }

    u.resize(2, (m_loops[loopID]->curves()).size() * np);
    int i=0;
    typename std::vector< gsCurve<T> *>::const_iterator it;
    gsMatrix<T> pts;
    gsMatrix<T> uCols;

    int firstInd=0;
    int secondInd=np-1;
    if (numEndPoints==0) {firstInd=1;secondInd=npoints-2;}

    for ( it= (m_loops[loopID]->curves()).begin(); it!= (m_loops[loopID]->curves()).end(); ++it )
    {
        //gsMatrix<T> * interval = (*it)->parameterRange();
        //gsMatrix<T> *  pts = gsPointGrid( interval->at(0), interval->at(1), npoints );
        pts.resize(1,np);
        for (int ii=firstInd;ii<=secondInd;ii++) pts(0,ii-firstInd)= T(ii)/(npoints-1);
        uCols.resize(2,np);
        (*it)->eval_into( pts, uCols );
        u.middleCols( i * np,np ) = uCols;

        (*it)->eval_into( pts, uCols );
        u.middleCols( i * np,np ) = uCols;
        ++i;
    }
}
template <class T>
T getDistance(gsVertex<T>* v1,gsVertex<T>* v2)  // todo: move as member of gsVertex 
{
    T x1=(v1->x());
    T x2=(v2->x());
    T y1=(v1->y());
    T y2=(v2->y());
    T z1=(v1->z());
    T z2=(v2->z());
    T dist=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
    return  dist;
}

template <class T>
void gsPlanarDomain<T>::sampleCurve_into( int loopID, int curveID, int npoints, gsMatrix<T> & u )
{
    assert( (loopID>=0) && (loopID < numLoops()) );
    assert( (curveID>=0) && (curveID < loop(loopID).size() ) );
    u.resize(2,npoints);
    typename std::vector< gsCurve<T> *>::const_iterator it;
    for ( it= (m_loops[loopID]->curves()).begin()+curveID; it!= (m_loops[loopID]->curves()).begin()+curveID+1 ; ++it )
    {
        //gsMatrix<T> * interval = (*it)->parameterRange();
        //gsMatrix<T> *  pts = gsPointGrid( interval->at(0), interval->at(1), npoints );
        gsMatrix<T> pts(1,npoints);
        for (int ii=0;ii<=npoints-1;ii++)
        { pts(0,ii)= T(ii)/(npoints-1); }; // todo: use gsPointGrid
        (*it)->eval_into( pts, u );
    }
}

/// Return a triangulation of the planar domain
template <class T>
memory::unique_ptr<gsMesh<T> > gsPlanarDomain<T>::toMesh(int npoints) const     // FOR NOW ONLY ONE LOOP
{
    typename gsMesh<T>::uPtr m(new gsMesh<T>());
    // Typedefs
    typedef typename gsVertex<T>::gsVertexHandle VertexHandle;
    //typedef std::deque<VertexHandle>            VertexList;
    //typedef typename std::deque<VertexHandle>::iterator  VertexListIter;

    GISMO_ASSERT(npoints > 0, "Number of points must be positive.");

#if FALSE
    T h_x = (m_bbox(0,1) - m_bbox(0,0) ) / (npoints-1) ,
            h_y = (m_bbox(1,1) - m_bbox(1,0) ) / (npoints-1) ;

    //========== 1. Compute y-parallel segments
    VertexList x_seg_begin; // Starting points of y-parallel segments
    VertexList x_seg_end;   // Endpoints  of y-parallel segments

    T xi = m_bbox(0,1);
    std::vector<T> tmp = m_loops[0]->lineIntersections(0, xi );// segments
    gsDebug<<"----- INTs:  "<<tmp.size()<<"( "<<xi<<") \n";
    gsAsMatrix<T> xxi(tmp);
    gsMatrix<T> xev;
    m_loops[0]->curve(0).eval_into(xxi, xev );// Push forward onto the curve: xev.row(1)==xi
    gsVector<T> x_seg = xev.row(1) ;// xev.row(0)==xi
    std::sort( x_seg.data(), x_seg.data()+x_seg.size() );

    for ( index_t j= 0; j!= x_seg.size(); j+=2 )
    {
        VertexHandle v = m->addVertex( xi, x_seg[j] );
        if ( x_seg[j] < x_seg[j+1] )
        {
            x_seg_begin.push_back( v );
            v = m->addVertex( xi, x_seg[j+1] );
            x_seg_end.push_back ( v );
        }
        else
        {
            x_seg_end.push_back( v );
            v = m->addVertex( xi, x_seg[j+1] );
            x_seg_begin.push_back ( v );
        }
    }

    xi = m_bbox(0,0);
    for ( int i = 0; i!= npoints-1; ++i )
    {
        tmp = m_loops[0]->lineIntersections(0, xi );
        gsDebug<<"----- INTs:  "<<tmp.size()<<"( "<<xi<<") \n";
        //gsDebug<<"----- INTs x:  "<<tmp.size()<<", i="<<i<<"\n";
        gsAsMatrix<T> xxi(tmp);
        m_loops[0]->curve(0).eval_into(xxi, xev );
        x_seg = xev.row(1) ;// xev.row(0)==xi
        std::sort( x_seg.data(), x_seg.data()+x_seg.size() );

        for ( index_t j= 0; j!= x_seg.size(); j+=2 )
        {
            VertexHandle v = m->addVertex( xi, x_seg[j] );
            if ( x_seg[j] < x_seg[j+1] )
            {
                x_seg_begin.push_back( v );
                v = m->addVertex( xi, x_seg[j+1] );
                x_seg_end.push_back ( v );
            }
            else
            {
                x_seg_end.push_back( v );
                v = m->addVertex( xi, x_seg[j+1] );
                x_seg_begin.push_back ( v );
            }
        }
        xi += h_x;
    }

    //========== 2. Sort the segments wrt both endpoints
    std::sort( x_seg_begin.begin(), x_seg_begin.end(), Yless<T> );
    std::sort( x_seg_end.begin()  , x_seg_end.end()  , Yless<T> );

    gsDebug<<"x-segments: "<<x_seg_begin.size()<<", "<<x_seg_end.size()<<"\n";

    //========== 3. March on x-parallel segments
    bool SegStart(false), SegEnd(false);
    VertexList line0; // x-sorted

    VertexListIter ss = x_seg_begin.begin(),
            es = x_seg_end.begin();
    T yi = m_bbox(1,0);
    VertexList Yprev;

    for ( int i = 0; i!= npoints; ++i )
    {
        tmp = m_loops[0]->lineIntersections(1, yi );
        //gsDebug<<"----- INTs y: "<<tmp.size()<<"\n";

        if ( ! tmp.empty() ) // Intersection event
        {
            gsDebug<<"---  Intersection event "<<i<<"( "<<yi<<" )\n";

            // Compute intersections with y=yi
            gsAsMatrix<T> yyi(tmp);
            gsMatrix<T> yev;
            m_loops[0]->curve(0).eval_into(yyi, yev );
            gsVector<T> y_seg = yev.row(0) ;// yev.row(1)==yi
            std::sort( y_seg.data(), y_seg.data()+y_seg.size() );
            VertexHandle v;
            VertexList line1;

            // Look for SegEnd events
            while( es != x_seg_end.end() &&  (*es)->y() < yi )
            {// Property: every line has unique points wrt x-coord
                gsDebug<<"SegEnd: "<< **es ;
                SegEnd=true;

                for ( VertexListIter it= line0.begin(); it!= line0.end(); ++it )
                    if (  (*it)->x() == (*es)->x() )
                    {
                        gsDebug<<"SegEnd: removed "<< **it ;
                        line0.erase( it );
                        break;
                    }
                es++;
            }

            // Mirror previous lines
            for ( VertexListIter it= line0.begin(); it!= line0.end(); ++it )
            {
                v = m->addVertex( (*it)->x(),  yi );
                line1.push_back(v);
            }

            // Look for SegStart events
            while( ss != x_seg_begin.end() &&  (*ss)->y() < yi )
            {
                SegStart=true;
                gsDebug<<"SegStart: added "<< (*ss)->x() <<".\n";
                line0.push_back(*(ss) );
                v = m->addVertex( (*ss)->x(),  yi );
                line1.push_back(v);
                ss++;
            }

            std::sort( line0.begin(), line0.end(), Xless<T> );
            std::sort( line1.begin(), line1.end(), Xless<T> );

            gsDebug<<"line0  has "<< line0.size()  <<" points\n";
            gsDebug<<"line1  has "<< line1.size()  <<" points\n";

            // add faces
            if ( ! line1.empty() && line1.size() == line0.size() )
            {
                gsDebug<<"Tiling row "<<i<<".\n";
                VertexListIter it0=line0.begin();
                for ( VertexListIter it1= line1.begin(); it1!= line1.end()-1; ++it1, ++it0 )
                    m->addFace( *it0, *it1,  *(it1+1), *(it0+1) );
            }
            else
            {
                gsDebug<<"Trouble..\n";
                VertexListIter it0=line0.begin();
                for ( VertexListIter it1= line1.begin(); it1!= line1.end(); ++it1, ++it0 )
                    gsDebug<< (*it0)->x() <<" - "<< (*it1)->x()  <<" \n";
                while ( it0 != line0.end() )
                    gsDebug<< (*it0++)->x() <<" -     * " <<" \n";

            }

            gsDebug<<"Connecting endpoints.\n";

            VertexList Yseg;
            // Make boundary vertices on y=yi
            Yseg.push_back( m->addVertex( y_seg[0]             , yi ) );
            Yseg.push_back( m->addVertex( y_seg[y_seg.size()-1], yi ) );

            if (SegStart )
            {
                // connect x-boundary point to level 0
                if ( Yprev.size() )
                {
                    m->addFace( Yprev[0],  line0.front(), line0[1] );
                    m->addFace( Yprev[1],  line0.back() , line0[line0.size()-2] );
                }
                // connect x-boundary point to level 1
                m->addFace( Yseg[0],  line1.front(), line0.front() );
                m->addFace( Yseg[1],  line1.back(), line0.back() );
            }
            else if ( Yprev.size() )
            {
                // Connect y-boundary points
                m->addFace( line1.front(),  line0.front(), Yprev[0], Yseg[0] );
                m->addFace( line1.back() ,  line0.back() , Yprev[1], Yseg[1] );
            }

            Yprev = Yseg;
            line0 = line1;
            SegStart=SegEnd=false;
        }
        else // No intersection event
        {
            Yprev.clear();
            line0.clear();
        }

        // next y-line
        yi += h_y;
    }

    // repeat for yi=m_bbox(1,1)

    return m;
#endif
    T bbox_length_y =m_bbox(1,1)-m_bbox(1,0);
    int yPoints = cast<T,int>( npoints*bbox_length_y) ;
    int lb_yPoints=25;
    if(yPoints<lb_yPoints)
    {
        if ( yPoints > 0 )
            npoints *= lb_yPoints / yPoints;
        else
            npoints *= lb_yPoints / 2;
        
        yPoints  = lb_yPoints;
    }

    // vector of y-coords of npoints lines parrallel to x-axis
    // IDEA: Also use these as x_sample_guides !!!
    gsMatrix<T> y_samples = gsPointGrid<T>( m_bbox(1,0), m_bbox(1,1), yPoints);

    // x-coords of sampling points on the line y=yi (if any)
    // int: yi index,
    //gsVector *: a pointer to a vector of npoints x-coords of sample points on the line y=yi
    //std::vector<index_t>: a vector of indices that indicate
    //the line segments of y=yi that intesect the domain

    std::vector<std::vector< VertexHandle > >samples;
    std::vector<std::vector< VertexHandle > >intersections;

    //    std::map<int,
    //            std::pair<
    //            std::vector< VertexHandle >,
    //            std::vector<T>
    //            >
    //            > x_samples;
    for ( int i = 0; i!= yPoints; ++i )
    {
        std::vector<T> x_all;
        //gsDebug<<" --- before intersections  " << i <<", y="<< y_samples(0,i)  <<"\n";
        for (size_t j=0;j<m_loops.size();j++)
        {
            std::vector<T> x = m_loops[j]->lineIntersections(1, y_samples(0,i) );
            if ( ! x.empty() )
            {
                typename gsCurve<T>::uPtr curve = m_loops[j]->singleCurve() ; // TO BE REMOVED later

                if ( x.size() == 1 )
                {
                    x.push_back(x[0]);
                    gsWarn<<"lineIntersection yielded a tangent"<<'\n';
                }
                else if (x.size()%2==1)
                {
                    x.push_back(x[0]);
                    gsWarn<<"lineIntersection yielded a tangent and more intersections, try another number of points"<<'\n';
                }

                gsAsMatrix<T> xx(x);
                gsMatrix<T> e = curve->eval( xx );
                for (size_t k=0;k<x.size();k++)
                {
                    x_all.push_back(e(0,k));
                }
            }
        }
        std::sort(x_all.begin(), x_all.end() );
        gsVector<int> numPoints( x_all.size() / 2 );
        int k(0);
        for ( typename std::vector<T>::const_iterator it = x_all.begin();
              it < x_all.end(); it += 2 )
            numPoints[k++] = cast<T,int>( (*(it+1) - *it)* npoints );

        std::vector<VertexHandle> x_line;
        std::vector<VertexHandle> intersection_vec;
        for(index_t j=0;j<numPoints.size();j++)
        {
            x_line.push_back(m->addVertex(x_all[j*2],y_samples(0,i)));
            intersection_vec.push_back(m->addVertex(x_all[j*2],y_samples(0,i)));

            for(k=0; k<numPoints[j]; k++)
            {
                x_line.push_back(m->addVertex(x_all[j*2]*(numPoints[j]-k)/(numPoints[j]+1)+x_all[j*2+1]*(k+1)/(numPoints[j]+1),y_samples(0,i)));
            }
            x_line.push_back(m->addVertex(x_all[j*2+1],y_samples(0,i)));
            intersection_vec.push_back(m->addVertex(x_all[j*2+1],y_samples(0,i)));
        }
        samples.push_back(x_line);
        intersections.push_back(intersection_vec);
    }
    T checkDist=(y_samples(0,1)-y_samples(0,0))*5;
    for(size_t i=0;i<intersections.size()-1;i++)
    {
        size_t currentTop=0;
        size_t currentBot=0;

        T topToBotDist=0;
        T botToTopDist=0;
        if(samples[i].size()>0&&samples[i+1].size()>0)
        {
            while (currentTop!=samples[i].size()-1||currentBot!=samples[i+1].size()-1)
            {
                if(currentTop==samples[i].size()-1)
                {
                    VertexHandle v1=samples[i][currentTop];
                    VertexHandle v2=samples[i+1][currentBot];
                    VertexHandle v3=samples[i+1][currentBot+1];
                    if(getDistance(v1,v2)<checkDist&&getDistance(v1,v3)<checkDist&&getDistance(v2,v3)<checkDist)
                        m->addFace(v1,v2,v3);
                    currentBot++;
                }
                else if(currentBot==samples[i+1].size()-1)
                {
                    VertexHandle v1=samples[i][currentTop];
                    VertexHandle v2=samples[i+1][currentBot];
                    VertexHandle v3=samples[i][currentTop+1];
                    if(getDistance(v1,v2)<checkDist&&getDistance(v1,v3)<checkDist&&getDistance(v2,v3)<checkDist)
                        m->addFace(v1,v2,v3);
                    currentTop++;
                }
                else
                {
                    topToBotDist=getDistance(samples[i][currentTop],samples[i+1][currentBot+1]);
                    botToTopDist=getDistance(samples[i][currentTop+1],samples[i+1][currentBot]);
                    if (topToBotDist<botToTopDist)
                    {
                        VertexHandle v1=samples[i][currentTop];
                        VertexHandle v2=samples[i+1][currentBot];
                        VertexHandle v3=samples[i+1][currentBot+1];
                        if(getDistance(v1,v2)<checkDist&&getDistance(v1,v3)<checkDist&&getDistance(v2,v3)<checkDist)
                            m->addFace(v1,v2,v3);
                        currentBot++;
                    }
                    else
                    {
                        VertexHandle v1=samples[i][currentTop];
                        VertexHandle v2=samples[i+1][currentBot];
                        VertexHandle v3=samples[i][currentTop+1];
                        if(getDistance(v1,v2)<checkDist&&getDistance(v1,v3)<checkDist&&getDistance(v2,v3)<checkDist)
                            m->addFace(v1,v2,v3);
                        currentTop++;
                    }
                }
            }
        }

    }


    //            // generate the Mesh vertices.
    //            k = npoints - x.size(); // points to distribute to the (interior of) the segments
    //            std::vector< VertexHandle > x_line;
    //            x_line.reserve(npoints);
    //            for (size_t v=0; v< x.size(); v+=2 )
    //            {
    //                x_line.push_back( m->addVertex(x[v], y_samples(0,i)) );
    //                // TO DO
    //                // ( ll[v/2] / ltotal )*k samples here
    //                for( int j = 1; j!= k+1; ++j)
    //                    x_line.push_back( m->addVertex(x[v] + T(j) * ll[v/2] / T(k+1), y_samples(0,i)) );

    //                x_line.push_back( m->addVertex(x[v+1], y_samples(0,i) ) );
    //            }
    //            //gsDebug<<"i: "<< i<< ", sz="<< x_line.size() <<", ints="<< x.size() <<"\n";
    //            //gsDebug<<"x= "<< x[0] <<", "<< x[1] <<"(y="<< y_samples(0,i) <<"\n";

    //            x_samples.insert( make_pair(i,make_pair( x_line, x ) ) ); // i-th sample line


    //    }

    //gsDebug<<" --- mesh DONE  "<< *m  <<"\n";
    // for ( int i = 0; i!= 25; ++i )
    //gsDebug<<"--"<<* m->vertex[i] ;

    // Zig-zag connection
    //    for ( int i = 0; i!= npoints-1; ++i )
    //    {
    //        typename std::map<int,std::pair<std::vector<VertexHandle>, std::vector<T> > >::const_iterator
    //            it0 = x_samples.find(i);

    //        if ( it0 != x_samples.end() ) // Check whether intersection with line 0 exists
    //        {
    //            typename std::map<int,std::pair<std::vector<VertexHandle>, std::vector<T> > >:: const_iterator
    //                it1 = x_samples.find(i+1);
    //            if ( it1 != x_samples.end() ) // Check whether intersection with line 1 exists
    //            {
    //                //typename std::vector<T>::const_iterator s0 = it0->second.second.begin();
    //                //typename std::vector<T>::const_iterator s1 = it1->second.second.begin();
    //                typename std::vector<VertexHandle>::const_iterator f1 = it1->second.first.begin();
    //                // Run through the sample points at line 0
    //                for ( typename std::vector<VertexHandle>::const_iterator f0 =
    //                        it0->second.first.begin();
    //                        f0 != it0->second.first.end()-1 ; ++f0 )
    //                {
    //                    // TO DO
    //                    //if ( *f0 == *s0 )// start segment level 0
    //                    // if ( *f1 == *s1 )// start segment level 1
    //                    m->addFace( *f0    , *f1    , *(f1+1), *(f0+1) ) ;
    //                    f1++;
    //                }
    //            }
    //        }
    //    }
    //gsDebug<<*m;
    return m;

}


}
