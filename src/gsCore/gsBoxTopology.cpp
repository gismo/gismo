
#include <gsCore/gsConfig.h>

#include <gsCore/gsBoxTopology.h>

namespace gismo
{
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
        if ( m_interfaces[i].ps1 == ps || m_interfaces[i].ps2 == ps ) {
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
        std::cerr << "                " << size() << " patches with " << numSides << " sides" << std::endl;
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
        if ( i->ps1.patch >= nboxes || i->ps2.patch >= nboxes ) {
            std::cerr << "*** WARNING *** gsBoxTopology: box index " << i->ps1.patch << " or " << i->ps2.patch
                      << " in interface out of range." << std::endl;
        }
}

bool gsBoxTopology::getNeighbour(const patchSide& ps ,patchSide& result, int & ii) const
{
    for ( unsigned i = 0; i < m_interfaces.size(); ++i ) 
    {
       if ( m_interfaces[i].ps1 == ps ) 
       {
           result = m_interfaces[i].ps2;
           ii     = i;
           return true;
       }
       else if ( m_interfaces[i].ps2 == ps ) 
       {
           result = m_interfaces[i].ps1;
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

bool gsBoxTopology::getCornerList(const patchCorner& start,std::vector<patchCorner> & cornerList) const
{
    bool innerVertex=true;
    cornerList.resize(0);
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
