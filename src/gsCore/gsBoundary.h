/** @file gsBoundary.h

    @brief Provides structs and classes related to interfaces and boundaries.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): F. Buchegger, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsAffineFunction.h>

namespace gismo
{
/** 
    Struct that defines the boundary sides and types. 
    The sides are numbered as follows:

    2D CASE                              |    3D CASE
    -------------------------------------|----------------------------------------
    Edge 1, {(u,v) : u = 0} :  \c west   |    Face 1, {(u,v,w) : u = 0}:  \c west
    Edge 2, {(u,v) : u = 1} :  \c east   |    Face 2, {(u,v,w) : u = 1}:  \c east
    Edge 3, {(u,v) : v = 0} :  \c south  |    Face 3, {(u,v,w) : v = 0}:  \c south
    Edge 4, {(u,v) : v = 1} :  \c north  |    Face 4, {(u,v,w) : v = 1}:  \c north
    &nbsp;                |    Face 5, {(u,v,w) : w = 0}:  \c front
    &nbsp;                |    Face 6, {(u,v,w) : w = 1}:  \c back

    \c none is a special compatibility value used to denote that this is not a boundary.

    The corners are numbered as follows:

    2D CASE                                            |    3D CASE
    ---------------------------------------------------|----------------------------------------
    Corner 1, {(u,v) : u = 0, v = 0} :  \c southwest   |    Corner 1, {(u,v,w) : u = 0, v = 0, w = 0}:  \c southwestfront
    Corner 2, {(u,v) : u = 1, v = 0} :  \c southeast   |    Corner 2, {(u,v,w) : u = 1, v = 0, w = 0}:  \c southeastfront
    Corner 3, {(u,v) : u = 0, v = 1} :  \c northwest   |    Corner 3, {(u,v,w) : u = 0, v = 1, w = 0}:  \c northwestfront
    Corner 4, {(u,v) : u = 1, v = 1} :  \c northeast   |    Corner 4, {(u,v,w) : u = 1, v = 1, w = 0}:  \c northeastfront
    &nbsp;                |    Corner 5, {(u,v,w) : u = 0, v = 0, w = 1}:  \c southwestback
    &nbsp;                |    Corner 6, {(u,v,w) : u = 1, v = 0, w = 1}:  \c southeastback
    &nbsp;                |    Corner 7, {(u,v,w) : u = 0, v = 1, w = 1}:  \c northwestback
    &nbsp;                |    Corner 8, {(u,v,w) : u = 1, v = 1, w = 1}:  \c northeastback
*/
struct boundary
{
    /// Identifiers for topological sides.
    enum side { west = 1, east = 2, south = 3, north= 4, front=5, back=6 , 
                left = 1, right= 2, down  = 3, up   = 4 , none = 0 };

    /// Types of boundary conditions.
    enum type {dirichlet=0, neumann=1, robin=2 };
    //mixed BD means: there are both dirichlet and neumann sides
    //robin: a linear combination of value and derivative
    //cauchy: there are 2 conditions (value+deriv) defined on the same side

    /// Identifiers for topological corners.
    enum corner { southwestfront = 1, southeastfront = 2, northwestfront = 3, northeastfront = 4,
                  southwestback  = 5, southeastback  = 6, northwestback  = 7, northeastback  = 8,
                  southwest      = 1, southeast      = 2, northwest      = 3, northeast      = 4 };
};

// functions for iterating over all sides for a given dimension
inline void firstSide(boundary::side& result);
inline bool nextSide(int dim, boundary::side& result);
  
/// Get a side from it's integer value
inline boundary::side getSide(int i)
{
    GISMO_ASSERT( i>0, "Requested side of none boundary.\n");
    return static_cast<boundary::side>(i);
}

// functions for iterating over all corners for a given dimension
inline void firstCorner(boundary::corner& result);
inline bool nextCorner(int dim, boundary::corner& result);

/// Get a corner from it's integer value
inline boundary::corner getCorner(int i)
{
    GISMO_ASSERT( i>0, "Requested corner of none boundary.\n");
    return static_cast<boundary::corner>(i);
}
    
/** 
    Struct patch_side represents a side of a patch
*/  
struct patch_side 
{
public:
    patch_side() { }

    patch_side(int p, boundary::side s)
    : patch(p), side(s) { }

    patch_side(int p, int s)
    : patch(p), side( getSide(s) ) { }

    int patch;              ///< The index of the patch.
    boundary::side side;    ///< The side of the patch.

    bool operator== (const patch_side& other) const
    { return patch == other.patch && side == other.side; }
};

/**
    Struct patch_corner represents a corner of a patch
*/
struct patch_corner
{
public:
    patch_corner() { }

    patch_corner(int p,boundary::corner c)
    : patch (p),corner(c) { }

    patch_corner(int p, int c)
    : patch(p), corner( getCorner(c) ) { }

    int patch;               ///< The index of the patch.
    boundary::corner corner; ///< The corner of the patch.

    bool operator== (const patch_corner& other) const
    { return patch == other.patch && corner == other.corner; }
};


/** 
    Struct boundaryInterface represents an interface between two patches
    by storing the two sides which are joined.
*/  
struct boundaryInterface
{
public:
    typedef std::pair<index_t,bool> dirInfo;

public:

    boundaryInterface() { }

    boundaryInterface(int p1, boundary::side s1,
                      int p2, boundary::side s2,
                      bool o1, bool o2)
    : ps1(p1, s1), ps2(p2, s2)
    {
        orient.resize(2);
        orient[0]= o1;
        orient[1]= o2;
    }

    boundaryInterface(int p1, boundary::side s1,
                      int p2, boundary::side s2,
                      bool o1)
    : ps1(p1, s1), ps2(p2, s2)
    {
        orient.resize(1);
        orient[0]= o1;
    }

    boundaryInterface(patch_side const & _ps1,
                      patch_side const & _ps2,
                      bool o1, bool o2)
    : ps1(_ps1), ps2(_ps2)
    {
        orient.resize(2);
        orient[0]= o1;
        orient[1]= o2;
    }

    boundaryInterface(patch_side const & _ps1,
                      patch_side const & _ps2,
                      bool o1)
    : ps1(_ps1), ps2(_ps2)
    {
        orient.resize(1);
        orient[0]= o1;
    }

    boundaryInterface(patch_side const & _ps1,
                      patch_side const & _ps2,
                      gsVector<bool> const & orient_flags)
    : ps1(_ps1), ps2(_ps2), orient(orient_flags)
    {  }

    boundaryInterface(int p1, boundary::side s1,
                      int p2, boundary::side s2,
                      gsVector<bool> const & orient_flags)
    : ps1(p1, s1), ps2(p2, s2), orient(orient_flags) { }

    boundaryInterface(gsVector<int>  const & psps,
                      gsVector<bool> const & orient_flags)
    : ps1(psps[0], psps[1]), ps2(psps[2], psps[3]), orient(orient_flags) { }
    
    bool operator== (const boundaryInterface& other) const
    {
        return ( ps1 == other.ps1 && ps2 == other.ps2 && orient == other.orient ) ||
               ( ps1 == other.ps2 && ps2 == other.ps1 && orient == other.orient );
    }
    
    patch_side  operator [] (size_t i) const
    { 
        if (i==0)
            return ps1;
        else if (i==1)
            return ps2;
        else
            GISMO_ERROR("Invalid index "<<i<<": Interface has 2 elements(sides).");
    }

    patch_side & operator [] (size_t i)
    { 
        if (i==0)
            return ps1;
        else if (i==1)
            return ps2;
        else
            GISMO_ERROR("Invalid index "<<i<<": Interface has 2 elements(sides).");
    }

public:

    patch_side ps1; ///< The first patch side.
    patch_side ps2; ///< The second patch side.

    /// dirInfo: m_iFaceInfo[i] is the pair (direction,bool) that
    /// denotes that direction \a i on \a ps1 is mapped to direction
    /// \a m_iFaceInfo[i].first with orientation m_iFaceInfo[i].second
    /// Orientation is a boolean that decides if the orientation is
    /// the same (1) or opposite (0).
    /// The entry \a m_iFaceInfo[ direction(ps1.side)] is included for
    /// compatibility.
    std::vector<dirInfo> m_iFaceInfo;

    // orient[i] is True iff ps1 and ps2 have the same orientation
    // with respect to parameter direction i
    // To be removed
    gsVector<bool> orient;

};

/// Returns the side corresponding to \a dir and \a par
inline boundary::side sideOf (int dir, int par)
{
    if ( par )
        return static_cast<boundary::side>(2*dir+2);
    else
        return static_cast<boundary::side>(2*dir+1);
}

/// Computes the orientation of the boundary side \a s with respect
/// to the interior of the parameter domain. The result is either +1
/// (counter-clockwise) or -1 (clock-wise).
///
/// \param s integer corresponding to a boundary side
inline int sideOrientation(int s)
{
    // The orientation sign is determined by the formula:
    // (direction + parameter) + 1 mod 2
    // where 0 means "-" (CW) and 1 means "+" (CCW)
    return ( ( s + (s+1) / 2 ) % 2 ? 1 : -1 );
}
  
/// \brief Returns the parametric direction that corresponds to side s.
/// \param[in] s integer corresponding to a boundary side
/// \return Integer which says which parameter has to be fixed
/// in order to get the boundary.
///
/// <b>Example:</b>\n
/// In 2D, let the parameter domain be defined by \f$(u,v)\in [0,1]^2\f$.
/// Since the side with index \em 3 corresponds to "south", i.e. to \f$ \{ (u,v):\ v = 0 \} \f$,
/// calling parameter(3) will return <em>1</em>, because \em v (i.e., parameter direction with index \em 1) is fixed/set to zero.\n
///
inline int direction (int s)
{
    GISMO_ASSERT( s>0, "Requested direction of none boundary.\n");
    return (s-1) / 2 ;
}  


/// \brief Returns the parameter value (false=0=start, true=1=end) that corresponds to side s
///
/// \param[in] s integer corresponding to a boundary side
/// \return \em false, if side \em s is defined by setting the
/// corresponding parameter to 0, and\n
/// \em true, if it is defined by setting the corresponding
/// parameter to 1.
///
/// <b>Example:</b>\n
/// In 2D, let the parameter domain be defined by \f$(u,v)\in [0,1]^2\f$.
/// Since the side with index \em 3 corresponds to "south", i.e. to \f$ \{ (u,v):\ v = 0 \} \f$,
/// calling parameter(3) will return <em>0=false</em>.
///
inline bool parameter (int s)
{
    GISMO_ASSERT( s>0, "Requested parameter of none boundary.\n");
    return ( (s+1) % 2 == 0 ? false : true ) ;
}

/**
 * takes a boundary side \a s and a dimension \a dim and fills the
 * vector corners with all the corners that belong to that side
 * @param[in] s : side, which is looked at
 * @param[in] dim : dimensions, 2 for 2D and 3 for 3D
 * @param[out] corners : vector of corners, will be filled with the corners of this side \a s
 */
inline void getCorners(const boundary::side& s,unsigned dim,std::vector<boundary::corner>& corners)
{
    GISMO_ASSERT(dim==2||dim==3,"only 2D and 3D allowed.");
    corners.resize(0);
    switch(s)
    {
    case boundary::south:
        corners.push_back(boundary::southwest);
        corners.push_back(boundary::southeast);
        break;
    case boundary::west:
        corners.push_back(boundary::southwest);
        corners.push_back(boundary::northwest);
        break;
    case boundary::north:
        corners.push_back(boundary::northwest);
        corners.push_back(boundary::northeast);
        break;
    case boundary::east:
        corners.push_back(boundary::southeast);
        corners.push_back(boundary::northeast);
        break;
    default: ;
    }
    if(dim==3)
    {
        switch(s)
        {
        case boundary::south:
            corners.push_back(boundary::southwestback);
            corners.push_back(boundary::southeastback);
            break;
        case boundary::west:
            corners.push_back(boundary::southwestback);
            corners.push_back(boundary::northwestback);
            break;
        case boundary::north:
            corners.push_back(boundary::southwestback);
            corners.push_back(boundary::northwestback);
            break;
        case boundary::east:
            corners.push_back(boundary::southeastback);
            corners.push_back(boundary::northeastback);
            break;
        case boundary::front:
            corners.push_back(boundary::southwestfront);
            corners.push_back(boundary::southeastfront);
            corners.push_back(boundary::northwestfront);
            corners.push_back(boundary::northeastfront);
            break;
        case boundary::back:
            corners.push_back(boundary::southwestback);
            corners.push_back(boundary::southeastback);
            corners.push_back(boundary::northwestback);
            corners.push_back(boundary::northeastback);
            break;
        default: ;
        }
    }
}

/**
 * takes a patch side \a ps and a dimension \a dim and fills the
 * vector pcorners with all the patch_corners that belong to that side
 * @param[in] ps : side, which is looked at
 * @param[in] dim : dimensions, 2 for 2D and 3 for 3D
 * @param[out] pcorners : vector of patch_corners, will be filled with the corners of this side \a ps
 */
inline void getPatchCorners(const patch_side& ps,unsigned dim,std::vector<patch_corner>& pcorners)
{
    std::vector<boundary::corner> corners;
    getCorners(ps.side,dim,corners);
    pcorners.resize(corners.size());
    for(unsigned i=0;i<corners.size();++i)
        pcorners[i]=patch_corner(ps.patch,corners[i]);
}

/// Takes a corner, a patch_side and the dimension and returns a gsVector defining
/// the parameters of that corners regarding the patch_side.
/// \param[in] c : patch_corner, which is looked at
/// \param[in] s : patch_side, where the patch_corner is on
/// \param[in] dim : dimensions, 2 for 2D and 3 for 3D
/// \param[out] pars : gsVector<bool> parameters, the (up to) two directions in ascending order
///                    (for north this would be (u) in 2D or (u,w) in 3D,
///                     for east this would be (v) in 2D or (v,w) in 3D)
inline void getParsOnSide(const patch_corner& c,const patch_side& s,unsigned dim,gsVector<bool>& pars)
{
    GISMO_ASSERT(dim==2||dim==3,"only 2D and 3D allowed.");
    std::vector<boundary::corner> corners;
    getCorners(s.side,dim,corners);
    GISMO_ASSERT(dim==2||(dim==3&&corners.size()==4),"3D case needs four corners for each patch_side");
    if(dim==2)
    {
        GISMO_ASSERT(corners.size()==2,"2D case needs two corners for each patch_side");
        pars.resize(1);
        if(c.corner==corners[0])
            pars(0)=0;
        else if(c.corner==corners[1])
            pars(0)=1;
        else
            GISMO_ERROR("one of the corners has to be the corner c");
    }
    else
    {
        GISMO_ASSERT(corners.size()==4,"3D case needs four corners for each patch_side");
        pars.resize(2);
        pars.setZero();
        if(c.corner==corners[0])
            return;
        else if(c.corner==corners[1])
            pars(0)=1;
        else if(c.corner==corners[2])
            pars(1)=1;
        else if(c.corner==corners[3])
        {
            pars(0)=1;
            pars(1)=1;
        }
        else
            GISMO_ERROR("one of the corners has to be the corner c");
    }
}

/// Returns the corner of a given \a s defined by
/// the gsVector<bool> \a pars, which gives the parameters of
/// the (up to) two directions in ascending order
/// (for north this would be (u) in 2D or (u,w) in 3D,
///  for east this would be (v) in 2D or (v,w) in 3D)
inline boundary::corner getCorner(const boundary::side s,const gsVector<bool>& pars)
{
    std::vector<bool> allPars(2);
    allPars[0] = pars(0);
    allPars[1] = pars.rows()==2 ? pars(1) : 0;
    std::vector<bool>::iterator it = allPars.begin();
    switch(s)
    {
    case boundary::west :
        allPars.insert(it,0);
        break;
    case boundary::east :
        allPars.insert(it,1);
        break;
    case boundary::south :
        allPars.insert(it+1,0);
        break;
    case boundary::north :
        allPars.insert(it+1,1);
        break;
    case boundary::front :
        allPars.insert(it+2,0);
        break;
    case boundary::back :
        allPars.insert(it+2,1);
        break;
    default: ;
    }
    return getCorner(allPars[0]+2*allPars[1]+4*allPars[2]+1);
}

/// Returns the patch_corner of a given \a ps defined by
/// the gsVector<bool> \a pars, which gives the parameters of
/// the (up to) two directions in ascending order
/// (for north this would be (u) in 2D or (u,w) in 3D,
///  for east this would be (v) in 2D or (v,w) in 3D)
inline patch_corner getPatchCorner(const patch_side& ps,const gsVector<bool>& pars)
{
    return patch_corner(ps.patch,getCorner(ps.side,pars));
}

/**
 * takes as input a corner \a c and a dimension \a dim and fills
 * the vector \a sides with all the sides connected to the corner \a c
 * @param[in] c : corner
 * @param[in] dim : dimensions, 2 for 2D and 3 for 3D
 * @param[out] sides : vector of sides, output variable
 */
inline void getSides(const boundary::corner& c,unsigned dim,std::vector<boundary::side>& sides)
{
    sides.resize(0);
    switch(c)
    {
    case boundary::southwest:
    case boundary::southeast:
    case boundary::southwestback:
    case boundary::southeastback:
        sides.push_back(boundary::south);
        break;
    case boundary::northwest:
    case boundary::northeast:
    case boundary::northwestback:
    case boundary::northeastback:
        sides.push_back(boundary::north);
        break;
    }
    switch(c)
    {
    case boundary::southwest:
    case boundary::northwest:
    case boundary::southwestback:
    case boundary::northwestback:
        sides.push_back(boundary::west);
        break;
    case boundary::southeast:
    case boundary::northeast:
    case boundary::southeastback:
    case boundary::northeastback:
        sides.push_back(boundary::east);
        break;
    }
    if(dim==3)
    {
        switch(c)
        {
        case boundary::southwestfront:
        case boundary::southeastfront:
        case boundary::northwestfront:
        case boundary::northeastfront:
            sides.push_back(boundary::front);
            break;
        case boundary::southwestback:
        case boundary::southeastback:
        case boundary::northwestback:
        case boundary::northeastback:
            sides.push_back(boundary::back);
            break;
        }
    }
}

/**
 * takes as input a patch_corner \a c and a dimension \a dim and fills
 * the vector \a psides with all the patch_sides connected to the corner \a c
 * @param[in] c : corner
 * @param[in] dim : dimensions, 2 for 2D and 3 for 3D
 * @param[out] psides : vector of sides, output variable
 */
inline void getPatchSides(const patch_corner& c,unsigned dim,std::vector<patch_side>& psides)
{
    std::vector<boundary::side> sides;
    getSides(c.corner,dim,sides);
    psides.resize(sides.size());
    for(unsigned i=0;i<sides.size();++i)
        psides[i]=patch_side(c.patch,sides[i]);
}

/// Returns in variable \a result the first side of a box
void firstSide(boundary::side& result)
{
    result = static_cast<boundary::side>(1);
}

/// Returns true if the next side after \a result of a box in
/// dimensio \a dim is existis, and updates \a result to that
/// side, false otherwise.  When returning false the \a result is
/// an invalid side
bool nextSide(int dim, boundary::side& result)
{
    result = static_cast<boundary::side>(result+1);
    return ( ( static_cast<int>(result) <= 2*dim)  && ( static_cast<int>(result) > 0) );
}

/// Returns in variable \a result the first corner of a box
void firstCorner(boundary::corner& result)
{
    result = static_cast<boundary::corner>(1);
}

/// Returns true if the next corner after \a result of a box in
/// dimensio \a dim is existis, and updates \a result to that
/// corner, false otherwise.  When returning false the \a result is
/// an invalid corner
bool nextCorner(int dim, boundary::corner& result)
{
    result = static_cast<boundary::corner>(result+1);
    return ( ( static_cast<int>(result) <= (1<<(dim-1)) )  && ( static_cast<int>(result) > 0) );
}

/// Print (as string) a boundary side
inline std::ostream &operator<<(std::ostream &os, const boundary::side& o)
{
    switch (o) 
    {
    case 0:
        os<<"none ";
        break;
    case 1:
        os<< "west ";
        break;
    case 2:
        os<< "east ";
        break;
    case 3:
        os<< "south";
        break;
    case 4:
        os<< "north";
        break;
    case 5:
        os<< "front";
        break;
    case 6:
        os<< "back ";
        break;
    default:
        os<< "side ";
        break;
    };           
    os<<"("<<static_cast<int>(o)<<")"; 
    return os;
};


/// Print (as string) an interface
inline std::ostream &operator<<(std::ostream &os, boundaryInterface const & i)
{
/*
  os<<"interface: first: "<< i.ps1.patch<<", "<< int(i.ps1.side)<<" ie. " <<i.ps1.side
  << " -  second: "     << i.ps2.patch<<", "<< int(i.ps2.side)<<" ie. " <<i.ps2.side
  <<", orientation = [";
  for ( index_t j = 0; j!=i.orient.size(); ++j)
  os << i.orient[0]<< ( j == i.orient.size() -1 ? "" : ", ") ;
  os<< "].\n";
*/
    os<< i.ps1.patch<<" "<< i.ps1.side<<" ~ "<<i.ps2.patch<<" "<<i.ps2.side<<" [";
    for ( index_t j = 0; j!=i.orient.size(); ++j) os << i.orient[0] ;
    os <<"]";
    return os; 
}

/// Print (as string) a patch side
inline std::ostream &operator<<(std::ostream &os, patch_side const & i)
{
    //os<<"Side: patch="<< i.patch<<", side="<< int(i.side)<<", ie. "<< i.side << ".\n";
    os<<i.patch<<" "<< int(i.side) ;
    return os; 
}

///Print (as string) a boundary type
inline std::ostream &operator<<(std::ostream &os, const boundary::type& o)
{
    switch (o) 
    {
    case 0:
        os<< "Dirichlet";
    case 1:
        os<< "Neumann";
    case 2:
        os<< "Mixed";
    default:
        gsInfo<<"boundary type not known.\n";
    };        
    return os;    
}

/// Returns the face of the box corresponding to the side
template <typename T>
gsMatrix<T> getFace (const boundary::side side, const gsMatrix<T> box)
{
    GISMO_ASSERT(direction(side)< box.rows(), "the specified side is not appropriate for the domain dimension");
    gsMatrix<T> temp=box;
    const index_t dir=direction(side);
    if (side%2) // bottom face
        temp(dir,1)=box(dir,0);
    else    // top face
        temp(dir,0)=box(dir,1);
    return temp;
}

/// construct an affine function mapping face s1 of domain1 to another face s2 of domain2 corresponding to orientation
/// orient is a 0-vector if domain* are 1D
///           a 1-vector if domain* are 2D
///           a 3-vector if domain* are 3D
///           ---cannot be generalized and we are forced to use orient this way because it is so defined
///           ---in boundaryInterface
/// the returned function must be deleted by the user
template <typename T>
gsFunction<T>* makeInterfaceMap ( const boundary::side s1, const boundary::side s2, const gsVector<bool> &orient, const  gsMatrix<T> &domain1, const  gsMatrix<T> &domain2)
{
    GISMO_ASSERT(domain1.rows()==domain2.rows(), "Impossible to map faces of diffeerent size");

    const index_t dir1=direction(s1);
    const index_t dir2=direction(s2);

    gsVector<T> shift1 = domain1.col(0);
    gsVector<T> shift2 = domain2.col(0);
    gsMatrix<T> map;
    gsVector<T> tran;
    switch (orient.rows())
    {
    case 0: // interfaces of a 1D domain
        GISMO_ASSERT(domain1.rows()==1, "Incompatible orientation and domains." );
        if( s2%2)
            return new gsConstantFunction<T>( domain2(0,0) );
        else
            return new gsConstantFunction<T>( domain2(0,1) );
        break;
    case 1: // interfaces of a 2D domain
        GISMO_ASSERT(domain1.rows()==2, "Incompatible orientation and domains." );
        map.setZero(2,2);
        tran.setZero(2);
        if (dir1==dir2)
        {
            map(0,0) = (domain2(0,1)-domain2(0,0))/(domain1(0,1)-domain1(0,0));
            map(1,1) = (domain2(1,1)-domain2(1,0))/(domain1(1,1)-domain1(1,0));
            if(orient(0))
            {
                map.col(1-dir1)*=-1;
                shift2(1-dir1)=domain2(1-dir1,1);
            }
        }
        else
        {
            map(1,0) = (domain2(1,1)-domain2(1,0))/(domain1(0,1)-domain1(0,0));
            map(0,1) = (domain2(0,1)-domain2(0,0))/(domain1(1,1)-domain1(1,0));
            if(orient(0))
            {
                map.col(1-dir1)*=-1;
                shift2(dir1)=domain2(dir1,1);
            }
        }
        break;
    case 2:
        GISMO_ERROR("can not associate an orientation to a vector of 2 bools!!!!! bad request.");
        break;
    case 3: // interfaces of a 3D domain
    {
        GISMO_ASSERT(domain1.rows()==3, "Incompatible orientation and domains." );
        map.setZero(3,3);
        tran.setZero(3);

        index_t dest[3]={0,1,2};
        if(orient(0))
        {
            dest[1]=2;
            dest[2]=1;
        }

        map(dir2,dir1) = (domain2(dir2,1)-domain2(dir2,0))/(domain1(dir1,1)-domain1(dir1,0));
        for (int i=1; i<3; ++i)
        {
            index_t io=(dir1+i)%3;
            index_t id=(dir2+dest[i])%3;
            if(!orient(i))
                map(id,io) = (domain2(id,1)-domain2(id,0))/(domain1(io,1)-domain1(io,0));
            else
            {
                map(id,io) = -(domain2(id,1)-domain2(id,0))/(domain1(io,1)-domain1(io,0));
                shift2(id) = domain2(id,1);
            }
        }
    }
        break;
    default:
        GISMO_ERROR("can not associate an orientation to a vector of more then 3 bools!!!!! bad request.");
        break;
    }
    if (!(s1%2))
        shift1(dir1)+=domain1(dir1,1)-domain1(dir1,0);
    if (!(s2%2))
        shift2(dir2)+=domain2(dir2,1)-domain2(dir2,0);
    tran += shift2-map*shift1;
    return new gsAffineFunction<T>(map,tran);
}


} // namespace gismo
