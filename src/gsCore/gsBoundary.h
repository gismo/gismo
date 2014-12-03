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


    /// Identifiers for topological corners.
    enum corner { southwestfront = 1, southeastfront = 2, northwestfront = 3, northeastfront = 4,
                  southwestback  = 5, southeastback  = 6, northwestback  = 7, northeastback  = 8,
                  southwest      = 1, southeast      = 2, northwest      = 3, northeast      = 4 };
};


// forward declarations of types
struct box_corner;
struct box_side;


/**
 * Struct side represent a side of a box
**/
struct box_side
{
public:
    index_t m_index;
public:
    box_side (index_t dir, bool par) : m_index(par?2*dir+2:2*dir+1) { GISMO_ASSERT(dir>=0,"invalid side");}
    box_side (index_t a=0) : m_index(a) { GISMO_ASSERT(a>=0,"invalid side");}
    box_side (boundary::side a) : m_index(a) { GISMO_ASSERT(a>=0,"invalid side");}

    // conversions
    operator boundary::side() const {return static_cast<boundary::side>(m_index);}
    box_side& operator= (boundary::side a) {m_index= a; return *this;}

    /**
     *\brief Returns the parametric direction orthogonal to this side.
     *  \return Integer which says which parameter has to be fixed
     *  in order to get the boundary.
     *
     *  <b>Example:</b>\n
     *  In 2D, let the parameter domain be defined by \f$(u,v)\in [0,1]^2\f$.
     *  Since the side with index \em 3 corresponds to "south", i.e. to \f$ \{ (u,v):\ v = 0 \} \f$,
     *  calling parameter(3) will return <em>1</em>, because \em v (i.e., parameter direction with index \em 1) is fixed/set to zero.\n
    **/
    index_t direction () const {return (m_index-1) / 2;}

    /**
     *  \brief Returns the parameter value (false=0=start, true=1=end) that corresponds to this side
     *
     *  \return \em false, if side \em s is defined by setting the
     *  corresponding parameter to 0, and\n
     *  \em true, if it is defined by setting the corresponding
     *  parameter to 1.
     *
     *  <b>Example:</b>\n
     *  In 2D, let the parameter domain be defined by \f$(u,v)\in [0,1]^2\f$.
     *  Since the side with index \em 3 corresponds to "south", i.e. to \f$ \{ (u,v):\ v = 0 \} \f$,
     *  calling parameter(3) will return <em>0=false</em>.
    **/
    bool    parameter () const {return (m_index-1)%2;}

    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the first valid side in an dim-dimensional box
    **/
    static  box_side    getFirst     (int dim) {return box_side(1);}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the last valid side in an dim-dimensional box
    **/
    static  box_side    getLast      (int dim) {return box_side(2*dim+1);}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the first valid side in an dim-dimensional box
    **/
    void    first     (int dim) {*this=getFirst(dim);}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the last valid side in an dim-dimensional box
    **/
    void    last      (int dim) {*this=getLast(dim);}
    /**
     * @brief set to next box_side
     */
    void    next      () { ++m_index;}
    /**
     * @brief set to previous box_side
     */
    void    previous  () { --m_index;}
    /**
     * @brief good
     * @param dim
     * @return true if the side is appropriate in dim dimensional boxes
     */
    bool good(int dim){return m_index>0 && m_index<(2*dim+1);}
};

/** 
    Struct patch_side represents a side of a patch
*/  
struct patch_side : public box_side
{
public:
    index_t patch;              ///< The index of the patch.
public:

    patch_side() : box_side(), patch(0) { }

    patch_side(int p, boundary::side s)
        : box_side(s), patch(p) { }

    patch_side(index_t p, box_side s)
        : box_side(s), patch(p) { }

    // getters
    box_side& side() {return *this;}
    const box_side& side() const {return *this;}

    operator box_side() const {return box_side(this->direction(), this->parameter());}

    //comparison
    // comparison
    bool operator== (const patch_side & other) const
    {
        return patch==other.patch && m_index == other.m_index;
    }


    // iteration
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the first valid side in an dim-dimensional box
    **/
    static  patch_side    getFirst     (int dim, index_t num_patches) { return patch_side(0,box_side::getFirst(dim));}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the last valid side in an dim-dimensional box
    **/
    static  patch_side    getLast      (int dim, index_t num_patches) { return patch_side(num_patches-1,box_side::getLast(dim));}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the first valid side in an dim-dimensional box
    **/
    void    first     (int dim, index_t num_patches) {*this=getFirst(dim,num_patches);}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the last valid side in an dim-dimensional box
    **/
    void    last      (int dim, index_t num_patches) {*this=getLast(dim,num_patches);}
    /**
     * @brief set to next box_side
     */
    void    next      (int dim) { box_side::next(); if (!box_side::good(dim)) {++patch; box_side::first(dim);} }
    /**
     * @brief set to previous box_side
     */
    void    previous  (int dim) { box_side::previous(); if (!box_side::good(dim)) {--patch; box_side::last(dim);} }
    /**
     * @brief good
     * @param dim
     * @return true if the patch side is appropriate dim dimensional boxes
     */
    bool good (int dim, index_t num_patches) {return box_side::good(dim) && 0<=patch && patch< num_patches;}


};

/// Print (as string) a patch side
inline std::ostream &operator<<(std::ostream &os, patch_side const & i)
{
    //os<<"Side: patch="<< i.patch<<", side="<< int(i.side)<<", ie. "<< i.side << ".\n";
    os<<i.patch<<" "<< int(i.side()) ;
    return os;
}



/**
 * Struct side represent a side of a box
**/
struct box_corner
{
public:
    index_t m_index;
public:
    box_corner (index_t a=0) : m_index(a) { GISMO_ASSERT(a>=0,"invalid side");}
    box_corner (boundary::corner a) : m_index(a) { GISMO_ASSERT(a>=0,"invalid side");}
    box_corner (gsVector<bool>  v) : m_index(1)
    {
        for (index_t p=0;p<v.rows();++p)
        {
            if(v(p))
                m_index+= 1<<p;
        }
    }

    // conversions
    operator boundary::corner() const {return static_cast<boundary::corner>(m_index);}
    box_corner& operator= (boundary::corner a) {m_index= a; return *this;}

    /**
     * @brief returns a vector of parameters describing the position
     *        of the corner
     * @param dim
     * @param param(i) is 1 if the corner is contained in box_face(i,1)
     *        and 0 if it is contained in box_face(i,0)
     */
    void parameters_into (int dim, gsVector<bool> &param) const
    {
        param.resize(dim);
        for (index_t i=0; i<dim; ++i)
            param(i)=((m_index-1)>>i)&1;
    }
    gsVector<bool> parameters(int dim) const
    {
        gsVector<bool> r;
        parameters_into(dim,r);
        return r;
    }


    /**
     * @brief returns a vector of sides containing the corner
     * @param dim is the ambient dimension
     * @param sides
     */
    void getContainingSides (int dim, std::vector<box_side> &sides) const
    {
        GISMO_ASSERT(dim>=0, "Dimension must be non negative");
        sides.resize(dim);
        gsVector<bool> param;
        parameters_into(dim,param);
        for (index_t i=0;i<dim;++i)
            sides[i]=box_side(i,param(i));
    }

    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the first valid side in an dim-dimensional box
    **/
    static  box_corner    getFirst     (int dim) {return box_corner(1);}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the last valid side in an dim-dimensional box
    **/
    static  box_corner    getLast      (int dim) {return box_corner((2<<dim)+1);}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the first valid side in an dim-dimensional box
    **/
    void    first     (int dim) {*this=getFirst(dim);}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the last valid side in an dim-dimensional box
    **/
    void    last      (int dim) {*this=getLast(dim);}
    /**
     * @brief set to next box_side
     */
    void    next      () { ++m_index;}
    /**
     * @brief set to previous box_side
     */
    void    previous  () { --m_index;}
    /**
     * @brief good
     * @param dim
     * @return true if the side is appropriate in dim dimensional boxes
     */
    bool good(int dim){return m_index>0 && m_index<((2<<dim)+2);}
};


/**
    Struct patch_corner represents a corner of a patch
*/
struct patch_corner : public box_corner
{
public:
    index_t patch;
public:
    patch_corner() : box_corner(0) { }
    patch_corner(int p,boundary::corner c)
        : box_corner(c), patch (p) { }

    patch_corner(int p, int c)
        : box_corner(c), patch (p) { }

    patch_corner(int p, box_corner c)
        : box_corner(c), patch (p) { }


    /**
     * @brief returns a vector of patch_sides that contain this corner
     * @param dim is the ambient dimension
     * @param sides
     */
    void getContainingSides (int dim, std::vector<patch_side> &sides) const
    {
        GISMO_ASSERT(dim>=0, "Dimension must be non negative");
        sides.resize(dim);
        gsVector<bool> param;
        parameters_into(dim,param);
        for (index_t i=0;i<dim;++i)
            sides[i]=patch_side(patch, box_side(i, param(i)));
    }

    // iteration
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the first valid side in an dim-dimensional box
    **/
    static  patch_corner    getFirst     (int dim, index_t num_patches) { return patch_corner(0,box_corner::getFirst(dim));}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the last valid side in an dim-dimensional box
    **/
    static  patch_corner    getLast      (int dim, index_t num_patches) { return patch_corner(num_patches-1,box_corner::getLast(dim));}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the first valid side in an dim-dimensional box
    **/
    void    first     (int dim, index_t num_patches) {*this=getFirst(dim,num_patches);}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the last valid side in an dim-dimensional box
    **/
    void    last      (int dim, index_t num_patches) {*this=getLast(dim,num_patches);}
    /**
     * @brief set to next box_side
     */
    void    next      (int dim) { box_corner::next(); if (!box_corner::good(dim)) {++patch; box_corner::first(dim);} }
    /**
     * @brief set to previous box_side
     */
    void    previous  (int dim) { box_corner::previous(); if (!box_corner::good(dim)) {--patch; box_corner::last(dim);} }
    /**
     * @brief good
     * @param dim
     * @return true if the patch side is appropriate dim dimensional boxes
     */
    bool good (int dim, index_t num_patches) {return box_corner::good(dim) && 0<=patch && patch< num_patches;}
};


/// temporary compatibility fix (to be removed): use patch_corner(ps.patch, box_corner(totParam)) instead
inline patch_corner getPatchCorner (const patch_side ps, const gsVector<bool> locParams)
{
    GISMO_ASSERT(ps.direction()<=locParams.cols(), "incompatible side and parameters, the dimension must agree");
    gsVector<bool> totParam;
    totParam.resize(locParams.cols()+1);
    index_t i=0;
    for (; i<ps.direction();++i)
        totParam(i)=locParams(i);
    totParam(i)=ps.parameter();
    for (++i; i<locParams.cols();++i)
        totParam(i+1)=locParams(i);
    return patch_corner(ps.patch, box_corner(totParam));
}
/// temporary compatibility fix (to be removed): use pc.getParameters(dim,result) instead
inline void getParsOnSide (const patch_corner pc, const patch_side ps, const int dim, gsVector<bool> &locParam)
{
    gsVector<bool> totParam;
    pc.parameters_into(dim,totParam);
    locParam.resize(dim-1);
    index_t i=0;
    for (; i<ps.direction();++i)
        locParam(i)=totParam(i);
    for (++i; i<locParam.cols();++i)
        locParam(i)=totParam(i+1);
}


/** 
    Struct boundaryInterface represents an interface between two patches
    by storing the two sides which are joined.
*/  
struct boundaryInterface
{
public:
    boundaryInterface() { }

    // special constructor for the 2d case
    boundaryInterface(patch_side const & _ps1,
                      patch_side const & _ps2,
                      bool o1)
        : ps1(_ps1), ps2(_ps2)
    {
        directionMap.resize(2);
        directionOrientation.resize(2);
        directionMap(ps1.direction())=ps2.direction();
        directionOrientation(ps1.direction())= (ps1.parameter()!=ps2.parameter());
        directionMap(1-ps1.direction())=1-ps2.direction();
        directionMap(1-ps1.direction())=o1;
    }


    //
    boundaryInterface(patch_side const & _ps1,
                      patch_side const & _ps2,
                      int dim)
        : ps1(_ps1), ps2(_ps2)
    {
        directionMap.resize(dim);
        directionOrientation.resize(dim);

        directionMap(ps1.direction())=ps2.direction();
        directionOrientation(ps1.direction())= (ps1.parameter()!=ps2.parameter());
        index_t i=1;
        for ( ;i<dim;++i)
        {
            const index_t o = (ps1.direction()+i)%dim;
            const index_t d = (ps2.direction()+i)%dim;

            directionMap(o)=d;
            directionOrientation=true;
            /// TODO: discuss and define default orientation
        }
    }

    boundaryInterface(patch_side const & _ps1,
                      patch_side const & _ps2,
                      gsVector<index_t> const & map_info,
                      gsVector<bool>    const & orient_flags)
        : ps1(_ps1), ps2(_ps2), directionMap(map_info), directionOrientation(orient_flags)
    {  }

    boundaryInterface(patch_side const & _ps1,
                      patch_side const & _ps2,
                      gsVector<bool>    const & orient_flags)
        : ps1(_ps1), ps2(_ps2)
    {
        const index_t dim = orient_flags.cols()+1;
        directionMap.resize(dim);
        directionOrientation.resize(dim);

        directionMap(ps1.direction())=ps2.direction();
        directionOrientation(ps1.direction())= (ps1.parameter()!=ps2.parameter());

        for (index_t i=1 ;i<dim;++i)
        {
            const index_t o = (ps1.direction()+i)%dim;
            const index_t d = (ps2.direction()+i)%dim;

            directionMap(o)=d;
            directionOrientation(o)=orient_flags(i-1);
        }
    }

    // DEPRECATED
    boundaryInterface(gsVector<int>     const & p,
                      gsVector<bool>    const & orient_flags)
        : ps1(patch_side(p(0),p(1))), ps2(patch_side(p(2),p(3)))
    {
        const index_t dim = orient_flags.cols()+1;
        directionMap.resize(dim);
        directionOrientation.resize(dim);

        directionMap(ps1.direction())=ps2.direction();
        directionOrientation(ps1.direction())= (ps1.parameter()!=ps2.parameter());

        for (index_t i=1 ;i<dim;++i)
        {
            const index_t o = (ps1.direction()+i)%dim;
            const index_t d = (ps2.direction()+i)%dim;

            directionMap(o)=d;
            directionOrientation(o)=orient_flags(i-1);
        }
    }

    bool operator== (const boundaryInterface & other) const
    {
        return ps1==other.ps1 && ps2==other.ps2
                && directionMap==other.directionMap
                && directionOrientation==other.directionOrientation;
    }

    /**
     * @brief first, returns the first patch_side of this interface
    **/
    patch_side & first  ()        {return ps1;}
    patch_side first    () const  {return ps1;}

    /**
     * @brief second, returns the second patch_side of this interface
    **/
    patch_side & second ()        {return ps2;}
    patch_side second   () const  {return ps2;}

    // use boundaryInterface.first() and boundaryInterface.second()
    // DEPRECATED
    patch_side  operator [] (size_t i) const
    {
        if (i==0)
            return ps1;
        else if (i==1)
            return ps2;
        else
            GISMO_ERROR("Invalid index "<<i<<": Interface has 2 elements(sides).");
    }
    // use boundaryInterface.first() and boundaryInterface.second()
    //DEPRECATED
    patch_side & operator [] (size_t i)
    {
        if (i==0)
            return ps1;
        else if (i==1)
            return ps2;
        else
            GISMO_ERROR("Invalid index "<<i<<": Interface has 2 elements(sides).");
    }


    // this is a work in progress, the old code was 2D specific, look at mapInfo and mapOrientation
    //DEPRECATED
    gsVector<bool> orient () const
    {
        GISMO_ASSERT(directionOrientation.size()==2, "This is deprecated and does not work if dim>2");
        gsVector<bool> result(1);
        result(0)=directionOrientation(1-ps1.direction());
        return result;
    }


    boundaryInterface getInverse() const
    {
        boundaryInterface result;
        result.directionMap.resize(directionMap.rows());
        result.directionOrientation.resize(directionOrientation.rows());
        result.ps1=ps2;
        result.ps2=ps1;
        for (index_t i=0;i<directionMap.rows();++i)
        {
            result.directionMap(directionMap(i))=i;
            result.directionOrientation(i)=directionOrientation(i);
        }
        return result;
    }

    patch_corner mapCorner ( const patch_corner c) const
    {
        gsVector<bool> par=c.parameters(directionMap.rows());
        const index_t dim=par.rows();
        gsVector<bool> new_par(dim);
        if (c.patch == ps1.patch && par(ps1.direction()) == ps1.parameter() )
        {
           index_t i=0;
          for (; i<ps1.direction();++i)
          {
            new_par(directionMap(i)) = directionOrientation(i) ?
                par(i) : !par(i);
          }
          new_par(i) = ps2.parameter();
          for (++i; i<dim;++i)
          {
            new_par(directionMap(i)) = directionOrientation(i) ?
                par(i) : !par(i);
          }
          return patch_corner(ps2.patch, box_corner(new_par));
        }
        else if (c.patch == ps2.patch && par(ps2.direction()) == ps2.parameter() )
        {
           index_t i=0;
          for (; i<ps1.direction();++i)
          {
            new_par(i) = directionOrientation(directionMap(i)) ?
                par(directionMap(i)) : !par(directionMap(i));
          }
          new_par(i) = ps1.parameter();
          for (++i; i<dim;++i)
          {
            new_par(i) = directionOrientation(directionMap(i)) ?
                par(directionMap(i)) : !par(directionMap(i));
          }
          return patch_corner(ps1.patch, box_corner(new_par));
        }
        else
        {
           gsWarn<<"cannot map corners that are not in the interface";
            return c;
        }
    }

public:

    patch_side ps1; ///< The first patch side.
    patch_side ps2; ///< The second patch side.


    /// We describe a permutation of the coordinates by storing
    /// a vector of integers:
    ///    - directionMap[i] stores the destination of coordinate i
    gsVector<index_t> directionMap;
    /// For each coordinate direction we save if the original coordinate and the destination one have the same orientation
    gsVector<bool>    directionOrientation;

    /// TODO: the information could be stored in a single vector of signed integers: the sign gives the orientation
protected:
    friend std::ostream &operator<<(std::ostream &os, const boundaryInterface i);
};


/// Print (as string) an interface
inline std::ostream &operator<<(std::ostream &os, const boundaryInterface i)
{
    os<<"interface between: "<<i.ps1.patch<<":"<< i.ps1.side()<<" and "<<i.ps2.patch<<":"<<i.ps2.side()<<" [ ";
    index_t j = 0;
    for ( ; j<i.directionMap.cols()-1; ++j)
    {
        os << j << "->" << (i.directionOrientation(j) ? "+" : "-") << i.directionMap(j)<<", ";
    }
    os << j << "->" << (i.directionOrientation(j) ? "+" : "-") << i.directionMap(j)<<"] ";
    return os;
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
// use box_side struct instead of enumerated values
//GS_DEPRECATED
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
// use box_side struct instead of enumerated values
// GS_DEPRECATED
inline bool parameter (int s)
{
    GISMO_ASSERT( s>0, "Requested parameter of none boundary.\n");
    return ( (s+1) % 2 == 0 ? false : true ) ;
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
}



/// Returns the face of the box corresponding to the side
template <typename T>
gsMatrix<T> getFace (const box_side side, const gsMatrix<T> box)
{
    GISMO_ASSERT(side.direction()< box.rows(), "the specified side is not appropriate for the domain dimension");
    gsMatrix<T> temp=box;
    const index_t dir=side.direction();
    if (side.parameter()) // bottom face
        temp(dir,0)=box(dir,1);
    else    // top face
        temp(dir,1)=box(dir,0);
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
GS_DEPRECATED gsFunction<T>* makeInterfaceMap ( const boundary::side s1, const boundary::side s2, const gsVector<bool> &orient, const  gsMatrix<T> &domain1, const  gsMatrix<T> &domain2)
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
