
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

bool gsBoxTopology::getNeighbour(const patchSide& ps ,patchSide& result, boundaryInterface* iface) const
{
    for ( unsigned i = 0; i < m_interfaces.size(); ++i ) {
        if ( m_interfaces[i].ps1 == ps ) {
            result = m_interfaces[i].ps2;
            if(iface) *iface=m_interfaces[i];
            return true;
        }
        else if ( m_interfaces[i].ps2 == ps ) {
            result = m_interfaces[i].ps1;
            if(iface) *iface=m_interfaces[i];
            return true;
        }
    }
    return false;
}

bool gsBoxTopology::getCornerList(const patchCorner& start,std::vector<patchCorner> & cornerList) const
{
    GISMO_ASSERT(m_dim==2,"works only for 2D");
    cornerList.resize(0);
    std::vector<patchSide> psides;
    start.getContainingSides(m_dim,psides);
    GISMO_ASSERT(psides.size()==static_cast<size_t>(m_dim),"there should always be two patchsides on each patchCorner");
    patchSide curSide = psides[0];
    patchSide endSide = psides[1];
    patchSide neighbour;
    boundaryInterface interface;
    patchCorner curCorner = start;
    patchCorner newCorner;
    bool firstTurn = true;
    do
    {
        cornerList.push_back(curCorner);
        if(!getNeighbour(curSide,neighbour, &interface))
        {
            if(firstTurn)
            {
                patchSide tempSide=curSide;
                curSide=endSide;
                endSide=tempSide;
                curCorner = start;
                firstTurn = false;
                if(!getNeighbour(curSide,neighbour, &interface))
                    break;
            }
            else
                break;
        }
        newCorner = interface.mapCorner(curCorner);
        newCorner.getContainingSides(m_dim,psides);
        if(neighbour == psides[0])
            curSide = psides[1];
        else if(neighbour == psides[1])
            curSide = psides[0];
        else
            GISMO_ERROR("one of the two sides has to be the neighbour.");
        curCorner=newCorner;
    }while(!(curCorner==start));
    return firstTurn;
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
