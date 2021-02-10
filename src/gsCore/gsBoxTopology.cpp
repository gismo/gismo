/** @file gsBoxTopology.cpp

    @brief Provides declaration of the BoxTopology class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): F. Buchegger, A. Mantzaflaris
*/

#include <gsCore/gsConfig.h>

#include <gsCore/gsBoxTopology.h>

namespace gismo
{

std::ostream & gsBoxTopology::print(std::ostream &os) const
{
    if ( nboxes > 0 )
    {
        os << "gsBoxTopology (" << nboxes << ").\n";
    }
    else
        os << "gsBoxTopology ( empty! ).\n";
    
    os << "Boundaries:";
    if ( m_boundary.size() )
    {
        for( std::vector< patchSide >::const_iterator bit =
                 m_boundary.begin(); bit != m_boundary.end(); ++bit )
        {
            os <<"\n"<< *bit <<" ";
        }
    }
    else
        os << " (none)";
    os << "\nInterfaces:";
    if ( m_interfaces.size() )
    {
        for( std::vector< boundaryInterface >::const_iterator iit =
                 m_interfaces.begin(); iit != m_interfaces.end(); ++iit )
        {
            os <<"\n"<< *iit <<" ";
        }
    }
    else
        os << " (none)";
    return os;
}


void gsBoxTopology::addAutoBoundaries()
{
    if (nboxes == 0)
    {
        return;
    }
    for (index_t b = 0; b < nboxes; ++b)
    {
        for (boxSide bs = boxSide::getFirst(m_dim); bs < boxSide::getEnd(m_dim); ++bs)
        {
            patchSide ps(b, bs);
            if (!isBoundary(ps) && !isInterface(ps))
            {
                addBoundary(ps);
            }
        }
    }
}

bool gsBoxTopology::isInterface( const patchSide& ps ) const
{
    for ( size_t i = 0; i < m_interfaces.size(); ++i )
        if ( m_interfaces[i].first() == ps || m_interfaces[i].second() == ps ) {
            return true;
        }
    return false;
}

void gsBoxTopology::checkConsistency() const
{
    GISMO_ASSERT(m_dim >= 0, "m_dim asserted to be positive");
    const size_t numSides = nboxes * 2 * m_dim;      // an n-D cube has 2*d sides
    if ( numSides != 2 * nInterfaces() + nBoundary() ) {
        std::cerr << "*** WARNING *** gsBoxTopology has inconsistent interfaces or boundaries!" <<
                     std::endl;
        std::cerr << "                " << nboxes << " patches with " << numSides << " sides" << std::endl;
        std::cerr << "                " << nInterfaces() << " declared interfaces" << std::endl;
        std::cerr << "                " << nBoundary() << " declared boundaries" << std::endl;
        std::cerr << "                this leaves " << numSides - 2 * nInterfaces() - nBoundary() <<
                     " sides unaccounted for" << std::endl;
    }
    for ( const_biterator i = bBegin(); i != bEnd(); ++i )
        if ( i->patch >= nboxes ) {
            std::cerr << "*** WARNING *** gsBoxTopology: box index " << i->patch << " in boundary out of range."
                      << std::endl;
        }
    for ( const_iiterator i = iBegin(); i != iEnd(); ++i )
        if ( i->first().patch >= nboxes || i->second().patch >= nboxes ) {
            std::cerr << "*** WARNING *** gsBoxTopology: box index " << i->first().patch << " or " << i->second().patch
                      << " in interface out of range." << std::endl;
        }
}

bool gsBoxTopology::getNeighbour(const patchSide& ps ,patchSide& result, int & ii) const
{
    for ( size_t i = 0; i < m_interfaces.size(); ++i )
    {
       if ( m_interfaces[i].first() == ps )
       {
           result = m_interfaces[i].second();
           ii     = i;
           return true;
       }
       else if ( m_interfaces[i].second() == ps )
       {
           result = m_interfaces[i].first();
           ii     = i;
           return true;
       }
    }
    return false;
}

bool gsBoxTopology::getNeighbour(const patchSide& ps ,patchSide& result) const
{
    int a;
    return getNeighbour(ps, result, a);
}

const boundaryInterface * gsBoxTopology::findInterface(const index_t b1, const index_t b2) const
{
    for ( size_t i = 0; i < m_interfaces.size(); ++i )
    {
        if ( (m_interfaces[i].first() .patch == b1  && 
              m_interfaces[i].second().patch == b2) ||
             (m_interfaces[i].first() .patch == b2  && 
              m_interfaces[i].second().patch == b1) )
            return & m_interfaces[i];
    }
    return NULL;
}

bool gsBoxTopology::getCornerList(const patchCorner& start,std::vector<patchCorner> & cornerList) const
{
    bool innerVertex=true;
    cornerList.clear();
    cornerList.push_back(start);
    std::vector<patchSide> visitedSides;
    std::vector<patchSide> psides;     // psides and vertices relate to each other
    std::vector<patchCorner> vertices;
    start.getContainingSides(m_dim,psides);
    for(size_t i = 0;i<psides.size();i++)
        vertices.push_back(start);
    patchSide ps,psNeighbour;
    patchCorner pc,pcNeighbour;
    boundaryInterface boundIf;
    do
    {
        ps = psides.back();
        pc = vertices.back();
        psides.pop_back();
        vertices.pop_back();
        if(std::find(visitedSides.begin(), visitedSides.end(), ps)==visitedSides.end())
            visitedSides.push_back(ps);
        if(!getNeighbour(ps,psNeighbour))
        {
            innerVertex=false;
            continue;
        }
        if(std::find(visitedSides.begin(), visitedSides.end(), psNeighbour)!=visitedSides.end())
            continue;
        visitedSides.push_back(psNeighbour);
        getInterface(ps,boundIf);
        pcNeighbour = boundIf.mapCorner(pc);
        if(pcNeighbour==pc)
            continue;
        std::vector<patchSide> neighbourSides;
        pcNeighbour.getContainingSides(m_dim,neighbourSides);
        for(size_t i=0;i<neighbourSides.size();++i)
        {
            psides.push_back(neighbourSides[i]);
            vertices.push_back(pcNeighbour);
        }
        if(std::find(cornerList.begin(), cornerList.end(), pcNeighbour)==cornerList.end())
            cornerList.push_back(pcNeighbour);
    }while(!psides.empty());
    return innerVertex;
}

int gsBoxTopology::getMaxValence() const
{
    patchCorner start;
    std::vector<patchCorner> cornerList;
    int valence,maxValence=-1;
    for(index_t i = 0;i<nboxes;++i)
    {
        for(int j = 1;j<=( 1 << m_dim );++j)
        {
            start=patchCorner(i,j);
            cornerList.clear();
            getCornerList(start,cornerList);
            valence=cornerList.size();
            if(valence>maxValence)
                maxValence=valence;
        }
    }
    return maxValence;
}

namespace {

// @brief Returns a canonic representation for the given patch corner
//
// This function determines all patchCorners that conicide with the given
// one and returns the one withe the smallest patch number.
patchCorner getCanonicCorner( const patchCorner& c, const gsBoxTopology& bt )
{
    std::vector< patchCorner > corners;
    bt.getCornerList(c,corners);
    return *std::min_element(corners.begin(), corners.end());
}

// @brief Returns the canonic corners for the given corners
//
// The result is sorted. Thus, it can be used to uniquely characterize all kinds
// of patchComponents (edges, faces, ...)
std::vector<patchCorner> getCanonicCorners( const std::vector<patchCorner>& c, const gsBoxTopology& bt )
{
    const index_t sz = c.size();
    std::vector< patchCorner > corners;
    corners.reserve(sz);
    for (index_t i=0; i<sz; ++i)
        corners.push_back( getCanonicCorner(c[i],bt) );
    std::sort(corners.begin(),corners.end());
    return corners;
}

// @brief Converts the given corners in unique corner ids
std::vector<index_t> getCornerIndices( const std::vector<patchCorner>& corner, index_t dim )
{
    const index_t sz = corner.size();
    std::vector<index_t> result(sz);
    for (index_t i=0; i<sz; ++i)
        result[i] = corner[i].patch*(1ull<<(dim)) + corner[i].m_index;
    return result;
}

} // end namespace

std::vector< std::vector<patchComponent> > gsBoxTopology::allComponents(bool combineCorners) const
{
    const size_t nPatches = nboxes;
    const short_t dim = m_dim;
    const short_t cnr = math::ipow(3, dim);

    typedef std::vector<patchComponent>                         component_coll_t;
    typedef std::map< std::vector<index_t>, component_coll_t>   map_t;

    std::vector<map_t> comps(dim+1);

    for (size_t i = 0; i<nPatches; ++i)
    {
        for (short_t j = 0; j<cnr; ++j)
        {
            patchComponent pc(i, j, dim);
            const index_t d = pc.dim();
            std::vector< patchCorner > crns = getCanonicCorners(pc.containedCorners(),*this);
            component_coll_t& g = comps[d][getCornerIndices(crns, dim)];
            g.push_back(pc);
        }
    }
    index_t sz = 0;
    for (short_t i=0; i<dim+1; ++i)
        sz += comps[i].size();

    std::vector<component_coll_t> result;
    result.reserve(sz);
    for (short_t i=dim; i>=0; --i)
    {
        if (!combineCorners || i>0)
            for( map_t::iterator it = comps[i].begin(); it != comps[i].end(); ++it )
                result.push_back( it->second );
        else
        {
            component_coll_t last;
            for( map_t::iterator it = comps[i].begin(); it != comps[i].end(); ++it )
            {
                const index_t nrcp = it->second.size();
                for (index_t j=0; j<nrcp; ++j)
                    last.push_back(it->second[j]);
            }
            result.push_back(last);
        }
    }
    return result;
}

void gsBoxTopology::getEVs(std::vector<std::vector<patchCorner> > & cornerLists) const
{
    GISMO_ASSERT(m_dim==2,"works only for 2D");
    cornerLists.clear();
    std::vector<patchCorner> cornerList;
    patchCorner c;
    for(index_t i = 0;i<nboxes;++i)
    {
        for(int j=1;j<=4;++j)
        {
            c=patchCorner(i,j);
            bool isCycle = getCornerList(c,cornerList);
            bool alreadyReached = false;
            for(size_t k = 0;k<cornerList.size();++k)
                if(cornerList[k].patch<i)
                    alreadyReached = true;
            if(isCycle&&cornerList.size()!=4&&!alreadyReached)
                cornerLists.push_back(cornerList);
        }
    }
}

void gsBoxTopology::getOVs(std::vector<std::vector<patchCorner> > & cornerLists) const
{
    GISMO_ASSERT(m_dim==2,"works only for 2D");
    cornerLists.clear();
    std::vector<patchCorner> cornerList;
    patchCorner c;
    for(index_t i = 0;i<nboxes;++i)
    {
        for(int j=1;j<=4;++j)
        {
            c=patchCorner(i,j);
            bool isCycle = getCornerList(c,cornerList);
            bool alreadyReached = false;
            for(size_t k = 0;k<cornerList.size();++k)
                if(cornerList[k].patch<i)
                    alreadyReached = true;
            if(isCycle&&cornerList.size()==4&&!alreadyReached)
                cornerLists.push_back(cornerList);
        }
    }
}
}
