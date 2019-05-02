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

#define watch(x) gsInfo << (#x) << " is " << (x) << "\n";

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

bool gsBoxTopology::getAllNeighbours(const patchSide& ps ,std::vector<patchSide> & result) const
{
    bool found = false;

    for ( unsigned i = 0; i < m_interfaces.size(); ++i )
    {
        if ( m_interfaces[i].first() == ps )
        {
            result.push_back(m_interfaces[i].second());
            found = true;
        }
        else if ( m_interfaces[i].second() == ps )
        {
            result.push_back(m_interfaces[i].first());
            found = true;
        }
    }
    return found;
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

bool gsBoxTopology::getNonMatchingCornerList(const patchCorner& start,std::vector<patchCorner> & cornerList) const
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
        gsInfo << "\n side patch " << ps.patch << " and " << ps.side() << "\n";
        gsInfo << "corner patch " << pc.patch << " and " << pc.m_index << "\n";
        GISMO_ERROR("stop");
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

namespace {

// @brief Returns a canonic representation for the given patch corner
//
// This function determines all patchCorners that conicide with the given
// one and returns the one withe the smallest patch number.
/*
 * // the original function
patchCorner getCanonicCorner( const patchCorner& c, const gsBoxTopology& bt )
{
    std::vector< patchCorner > corners;
    bt.getCornerList(c,corners);
    return *std::min_element(corners.begin(), corners.end());
}
*/

patchCorner getCanonicCorner( const patchCorner& c, const gsBoxTopology& bt )
{
    std::vector< patchCorner > corners;
    bt.getCornerList(c,corners);
    //gsInfo << "patch: " << c.patch << "\n";
    //for(int k = 0; k < corners.size(); k++)
    //{
    //    gsInfo << "corner: " << corners[k].m_index << " with patch " << corners[k].patch << "\n";
    //}
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
        result[i] = corner[i].patch*(1u<<(dim)) + corner[i].m_index;
    return result;
}

} // end namespace

std::vector< std::vector<patchComponent> > gsBoxTopology::allComponents(bool combineCorners) const
{
    const index_t nPatches = nboxes;
    const index_t dim = m_dim;

    index_t cnr = 1;
    for (index_t i=0; i<dim; ++i) cnr *= 3;

    typedef std::vector<patchComponent>                         component_coll_t;
    typedef std::map< std::vector<index_t>, component_coll_t>   map_t;

    std::vector<map_t> comps(dim+1);

    for (index_t i = 0; i<nPatches; ++i)
    {
        for (index_t j = 0; j<cnr; ++j)
        {
            patchComponent pc(i, j, dim);
            const index_t d = pc.dim();
            std::vector< patchCorner > crns = getCanonicCorners(pc.containedCorners(),*this);
            component_coll_t& g = comps[d][getCornerIndices(crns, dim)];
            g.push_back(pc);
        }
    }
    index_t sz = 0;
    for (index_t i=0; i<dim+1; ++i)
        sz += comps[i].size();

    std::vector<component_coll_t> result;
    result.reserve(sz);
    for (index_t i=dim; i>=0; --i)
    {
        if (!combineCorners || i>0)
            for( typename map_t::iterator it = comps[i].begin(); it != comps[i].end(); ++it )
                result.push_back( it->second );
        else
        {
            component_coll_t last;
            for( typename map_t::iterator it = comps[i].begin(); it != comps[i].end(); ++it )
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

std::vector< std::vector<patchComponent> > gsBoxTopology::allNonMatchingComponents(bool combineCorners) const
{
    const index_t nPatches = nboxes;
    const index_t dim = m_dim;

    index_t cnr = 1;
    for (index_t i=0; i<dim; ++i) cnr *= 3;

    typedef std::vector<patchComponent>                         component_coll_t;
    typedef std::map< std::vector<index_t>, component_coll_t>   map_t;

    std::vector<map_t> comps(dim+1);

    // first iterate over all interfaces and look for sides which have 2 or more neighbours
    // add all participating sides to one component
    // add all vertices between the neighbouring patches to the component
    std::vector<std::vector<std::vector<patchSide>> > neighbours;
    neighbours.resize(this->nBoxes());
    for(std::vector<std::vector<std::vector<patchSide>> >::iterator it = neighbours.begin(); it != neighbours.end(); ++it)
        it->resize(4);

    for(std::vector<boundaryInterface>::const_iterator it = m_interfaces.begin(); it != m_interfaces.end(); ++it)
    {
        boundaryInterface bi = *it;
        patchSide p1 = bi.first();
        patchSide p2 = bi.second();

        std::vector<patchSide> temp1, temp2;
        temp1.push_back(p2);
        temp2.push_back(p1);
        neighbours[p1.patch][p1.side().index()-1].push_back(p2);
        neighbours[p2.patch][p2.side().index()-1].push_back(p1);
    }

    component_coll_t extension;
    component_coll_t boundaryCorners;
    std::vector<component_coll_t> coll_boundaryCorners;
    std::vector<patchSide> processedSide;
    std::vector<patchCorner> processedCorner;

    for(size_t i = 0; i < neighbours.size(); ++i)
    {
        //gsInfo << "patch: " << i << "\n";
        for(size_t j = 0; j < neighbours[i].size(); ++j)
        {
            //gsInfo << " with side: " << j+1 << "\n";
            if(neighbours[i][j].size() >= 2)
            {
                // add the side of the patch with 2 or more neighbours
                patchSide ps1(i, j+1);
                //processedSide.push_back(ps1);

                if((std::find(processedSide.begin(), processedSide.end(), ps1)) == processedSide.end())
                {
                    if (ps1.index() == 4)
                        extension.push_back(patchComponent(ps1.patch, 6, m_dim));
                    else
                        extension.push_back(patchComponent(ps1.patch, ps1.index(), m_dim));
                }

                processedSide.push_back(ps1);
                size_t nAdj = neighbours[i][j].size();
                //gsInfo << "adjacent " << nAdj << "\n";
                std::vector<index_t> patches, patches2;
                for(size_t h = 0; h < neighbours[i][j].size(); ++h)
                {
                    // add the sides of the adjacent patches
                    //gsInfo << " this is patch side " << neighbours[i][j][h] << "\n";
                    patchSide ps = neighbours[i][j][h];
                    if(std::find(processedSide.begin(), processedSide.end(), ps)==processedSide.end())
                    {
                        if (ps.index() == 4)
                            extension.push_back(patchComponent(ps.patch, 6, m_dim));
                        else
                            extension.push_back(patchComponent(ps.patch, ps.index(), m_dim));
                    }

                    // collect the adjacent patch numbers for the knots
                    patches.push_back(ps.patch);
                    processedSide.push_back(ps);
                }

                patches2 = patches;

                // add the common corners of adjacent boxes
                for(size_t h = 0; h < nAdj-1; ++h)
                {
                    int p1 = patches.back();
                    patches.pop_back();
                    int p2 = patches.back();
                    gsInfo << "p1: " << p1 << ", p2: " << p2 << "\n";

                    std::vector<patchCorner> crn;
                    std::vector<patchSide> pss;
                    const boundaryInterface* bI_ptr = findInterface(p1, p2);
                    bI_ptr->first() < bI_ptr->second() ? bI_ptr->first().getContainedCorners(m_dim, crn) : bI_ptr->second().getContainedCorners(m_dim, crn);
                    processedSide.push_back(bI_ptr->first());
                    processedSide.push_back(bI_ptr->second());

                    for(std::vector<patchCorner>::iterator it = crn.begin(); it != crn.end(); ++it)
                    {
                        it->getContainingSides(m_dim, pss);
                        for(std::vector<patchSide>::iterator sit = pss.begin(); sit != pss.end(); ++sit)
                        {
                            std::vector<patchSide> tempNeighbour;
                            patchCorner other;
                            if(getAllNeighbours(*sit, tempNeighbour))
                            {
                                gsInfo << "p1: " << *sit << " with neighbour: " << tempNeighbour[0] << " and " << tempNeighbour[1] << "\n";
                                watch(i);
                                if((std::find(tempNeighbour.begin(), tempNeighbour.end(), ps1)) != tempNeighbour.end()) // found the correct vertex
                                {
                                    // add the vertex of patch 0
                                    if(it->m_index == 1)
                                        extension.push_back(patchComponent(it->patch, 3*(it->m_index)+1, m_dim));
                                    else if (it->m_index == 2 || it->m_index == 3)
                                        extension.push_back(patchComponent(it->patch, (3*(it->m_index) - (it->m_index-1)), m_dim));
                                    else
                                        extension.push_back(patchComponent(it->patch, 8, m_dim));

                                    // add the vertex of patch 1
                                    other = bI_ptr->mapCorner(*it);
                                    // !!! delete the following if then else statement if the corner is required only for the smallest participating patch
                                    if(other.m_index == 1)
                                        extension.push_back(patchComponent(other.patch, 3*(other.m_index)+1, m_dim));
                                    else if(other.m_index == 2 || other.m_index == 3)
                                        extension.push_back(patchComponent(other.patch, (3*(other.m_index) - (other.m_index-1)), m_dim));
                                    else
                                        extension.push_back(patchComponent(other.patch, 8, m_dim));

                                    // push both to the vector of processed corners
                                    processedCorner.push_back(*it);
                                    processedCorner.push_back(other);
                                }
                            }
                        }
                    }

                }

                // check if the vertices of the interface are on the boundary
                patchSide ps(i, j+1), nghbour;
                patchCorner pc;
                std::vector<patchCorner> crn;
                std::vector<patchSide> pss, pss2;
                int boundaries = 0, boundaries2 = 0;

                ps.getContainedCorners(m_dim, crn);

                for(std::vector<patchCorner>::iterator it = crn.begin(); it != crn.end(); ++it)
                {
                    //watch(*it);
                    boundaries = 0; boundaries2 = 0;
                    it->getContainingSides(m_dim, pss);
                    for (std::vector<patchSide>::iterator sit = pss.begin(); sit != pss.end(); ++sit)
                    {
                        if (isBoundary(*sit))
                            boundaries++;
                    }


                    if(boundaries == 1)
                    {
                        for (std::vector<int>::iterator sit = patches2.begin(); sit != patches2.end(); ++sit)
                        {
                            const boundaryInterface* bi = findInterface(i, *sit);
                            patchCorner temp_pc;
                            pc = bi->mapCorner(*it);
                            //watch(*sit);
                            //watch(pc);
                            if(std::find(processedCorner.begin(), processedCorner.end(), pc)==processedCorner.end())
                            {
                                pc.getContainingSides(m_dim, pss2);

                                for (std::vector<patchSide>::iterator sit2 = pss2.begin(); sit2 != pss2.end(); ++sit2)
                                {
                                    if (isBoundary(*sit2))
                                    {
                                        boundaries2++;
                                        //pc = temp_pc;
                                    }
                                }
                            }
                            if (boundaries == 1 && boundaries2 == 1)
                                break;
                        }
                    }

                    if (boundaries == 1 && boundaries2 == 1)
                    {
                        if (it->m_index == 1)
                            boundaryCorners
                                .push_back(patchComponent(it->patch, 4, m_dim));
                        else if (it->m_index == 2)
                            boundaryCorners
                                .push_back(patchComponent(it->patch, 5, m_dim));
                        else if (it->m_index == 3)
                            boundaryCorners
                                .push_back(patchComponent(it->patch, 7, m_dim));
                        else if (it->m_index == 4)
                            boundaryCorners
                                .push_back(patchComponent(it->patch, 8, m_dim));


                        if (pc.m_index == 1)
                            boundaryCorners.push_back(patchComponent(pc.patch, 3*(pc.m_index)+1, m_dim));
                        else if (pc.m_index == 2 || pc.m_index == 3)
                            boundaryCorners
                                .push_back(patchComponent(pc.patch, (3 * (pc.m_index) - (pc.m_index - 1)), m_dim));
                        else
                            boundaryCorners.push_back(patchComponent(pc.patch, 8, m_dim));

                        processedCorner.push_back(*it);
                        processedCorner.push_back(pc);
                    }

                    if(!boundaryCorners.empty())
                    {
                        coll_boundaryCorners.push_back(boundaryCorners);
                        boundaryCorners.clear();
                    }
                }
                //coll_boundaryCorners.push_back(boundaryCorners);
                //boundaryCorners.clear();
            }
        }
    }

    /*
    gsInfo << "number of elements: " << extension.size() << "\n";
    for(component_coll_t::iterator it = extension.begin(); it != extension.end(); ++it)
    {
        gsInfo << "Index: " << it->index() << "\n";
    }
     */

    for (index_t i = 0; i<nPatches; ++i)
    {
        for (index_t j = 0; j<cnr; ++j)
        {
            patchSide p;
            patchCorner c;
            p.patch = i;
            c.patch = i;

            if(j == 1 || j == 2 || j == 3 )
                p.side() = j;
            else
                if(j == 6)
                    p.side() = 4;

            if(j == 4 || j == 5)
                c.m_index = j - 3;
            else if(j == 7 || j == 8)
                c.m_index = j - 4;

            if((std::find(processedSide.begin(), processedSide.end(), p)==processedSide.end()) &&
                (std::find(processedCorner.begin(), processedCorner.end(), c)==processedCorner.end()))
            {
                patchComponent pc(i, j, dim);
                const index_t d = pc.dim();
                std::vector< patchCorner > crns = getCanonicCorners(pc.containedCorners(),*this);
                component_coll_t& g = comps[d][getCornerIndices(crns, dim)];
                g.push_back(pc);
            }
            else
                continue;
        }
    }

    index_t sz = 0;
    for (index_t i=0; i<dim+1; ++i)
        sz += comps[i].size();

    std::vector<component_coll_t> result;
    result.reserve(sz);
    for (index_t i=dim; i>=0; --i)
    {
        if (!combineCorners || i>0)
            for( typename map_t::iterator it = comps[i].begin(); it != comps[i].end(); ++it )
                result.push_back( it->second );
        else
        {
            component_coll_t last;
            for( typename map_t::iterator it = comps[i].begin(); it != comps[i].end(); ++it )
            {
                const index_t nrcp = it->second.size();
                for (index_t j=0; j<nrcp; ++j)
                    last.push_back(it->second[j]);
            }
            result.push_back(last);
        }
    }

    //then push the extensions...
    result.push_back(extension);

    //... and the boundary corners of non matching interfaces if existent
    if(coll_boundaryCorners.size() > 0)
        for(size_t i = 0; i < coll_boundaryCorners.size(); ++i)
            result.push_back(coll_boundaryCorners[i]);

    return result;
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
