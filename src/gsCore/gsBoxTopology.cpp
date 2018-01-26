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
    if ( nboxes == 0 ) {
        return;
    }
    for (index_t b=0; b<nboxes; ++b)
    {
        for (boxSide bs=boxSide::getFirst(m_dim); bs<boxSide::getEnd(m_dim); ++bs)
        {
            patchSide ps(b,bs);
            if ( !isBoundary( ps ) && !isInterface( ps ) ) {
                addBoundary( ps );
            }
        }
    }
}

bool gsBoxTopology::isInterface( const patchSide& ps ) const
{
    for ( unsigned i = 0; i < m_interfaces.size(); ++i )
        if ( m_interfaces[i].first() == ps || m_interfaces[i].second() == ps ) {
            return true;
        }
    return false;
}

void gsBoxTopology::checkConsistency() const
{
    const int numSides = nboxes * 2 * m_dim;      // an n-D cube has 2*d sides
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
    for ( unsigned i = 0; i < m_interfaces.size(); ++i ) 
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

const boundaryInterface * gsBoxTopology::findInterface(const int b1, const int b2) const
{
    for ( unsigned i = 0; i < m_interfaces.size(); ++i ) 
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
    for(unsigned i = 0;i<psides.size();i++)
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
        if(std::find(visitedSides.begin(), visitedSides.end(), psNeighbour)==visitedSides.end())
            visitedSides.push_back(psNeighbour);
        getInterface(ps,boundIf);
        pcNeighbour = boundIf.mapCorner(pc);
        if(pcNeighbour==pc)
            continue;
        std::vector<patchSide> neighbourSides;
        pcNeighbour.getContainingSides(m_dim,neighbourSides);
        for(unsigned i=0;i<neighbourSides.size();++i)
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
    for(int i = 0;i<nboxes;++i)
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

void gsBoxTopology::getEVs(std::vector<std::vector<patchCorner> > & cornerLists) const
{
    GISMO_ASSERT(m_dim==2,"works only for 2D");
    cornerLists.clear();
    std::vector<patchCorner> cornerList;
    patchCorner c;
    for(int i = 0;i<nboxes;++i)
    {
        for(int j=1;j<=4;++j)
        {
            c=patchCorner(i,j);
            bool isCycle = getCornerList(c,cornerList);
            bool alreadyReached = false;
            for(unsigned k = 0;k<cornerList.size();++k)
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
    for(int i = 0;i<nboxes;++i)
    {
        for(int j=1;j<=4;++j)
        {
            c=patchCorner(i,j);
            bool isCycle = getCornerList(c,cornerList);
            bool alreadyReached = false;
            for(unsigned k = 0;k<cornerList.size();++k)
                if(cornerList[k].patch<i)
                    alreadyReached = true;
            if(isCycle&&cornerList.size()==4&&!alreadyReached)
                cornerLists.push_back(cornerList);
        }
    }
}
}
