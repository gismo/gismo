/** @file gsBoundary.h

    @brief Provides structs and classes related to interfaces and boundaries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Bressan, F. Buchegger, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

#include <gsCore/gsExport.h>

namespace gismo
{
/** 
    @brief Struct that defines the boundary sides and corners and types of a geometric object.

    These definitions are used by, e.g., boxSide, boxCorner, etc.

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
    
    \ingroup Core
*/
struct boundary
{
    /// Identifiers for topological sides.
    enum side { west = 1, east = 2, south = 3, north= 4, front=5, back=6 ,
                left = 1, right= 2, down  = 3, up   = 4 , none = 0 };


    /// Identifiers for topological corners.
    // warning: naming southwest etc ambiguous for 3D (corresponds to an edge)
    enum corner { southwestfront = 1, southeastfront = 2, northwestfront = 3, northeastfront = 4,
                  southwestback  = 5, southeastback  = 6, northwestback  = 7, northeastback  = 8,
                  southwest      = 1, southeast      = 2, northwest      = 3, northeast      = 4
    };
};

struct boxCorner;// defined later

/**
   @brief Struct which represents a certain side of a box.
   
*/
class GISMO_EXPORT boxSide
{
public:
    enum side { none = 0,
                west = 1, east = 2, south = 3, north= 4, front=5, back=6 ,
                left = 1, right= 2, down  = 3, up   = 4   };

    /** \brief Index of the side
    *
    * ...stored as number, specified in boundary::side
    **/
    index_t m_index;
public:
    boxSide (): m_index(0) {}
    boxSide (index_t dir, bool par) : m_index(index(dir,par))
    { GISMO_ASSERT(dir>=0,"invalid side");}
    boxSide (index_t a) : m_index(a) { GISMO_ASSERT(a>=0,"invalid side");}
    boxSide (boundary::side a) : m_index(a) { GISMO_ASSERT(a>=0,"invalid side");}

    // conversions
    operator boundary::side() const {return static_cast<boundary::side>(m_index);}

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
    bool    parameter () const {return (m_index-1)%2 != 0;}

    /**
     * @brief returns the parallel opposite side
     * @return
     */
    boxSide opposite() const {return boxSide(((m_index-1)^1)+1);}


    /**
     *  \brief Returns the index (as specified in boundary::side) of the box side
    **/
    int    index () const {return m_index;}

    /**
     *  \brief Returns the index of the box side implied by input
     *  direction \a dir and parameter \a par
    **/
    static inline int index (index_t dir, bool par) {return par?2*dir+2:2*dir+1;}


    /** 
     * @brief returns the vector of the corners contained in the side 
     * @param dim is the ambient dimension 
     * @param corners 
     */ 
    void getContainedCorners (int dim, std::vector<boxCorner> &corners) const; 
    
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the first valid side in an dim-dimensional box
    **/
    static  boxSide    getFirst     (int dim) {return boxSide(1);}

    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the last valid side in an dim-dimensional box
    **/

    static  boxSide    getLast      (int dim) {return boxSide(2*dim);}
    /**
     * @brief helper for iterating on sides of an n-dimensional box
     * @param dim
     * @return the (invalid) side after the last one in dim-dimensional box
    **/

    static  boxSide    getEnd       (int dim) {return boxSide(2*dim+1);}

    /**
     * @brief Incrementset boxSide
     */
    boxSide& operator++ ()    { ++m_index; return *this;} //prefix
    //boxSide  operator++ (int) { boxSide temp(*this); ++m_index; return temp;} //postfix
 
    /**
     * @brief Decrement boxSide
     */
    boxSide& operator-- ()    { --m_index; return *this;} //prefix
    //boxSide  operator-- (int) { boxSide temp(*this); --m_index; return temp;} //postfix

};

/** 
    @brief  Struct which represents a certain side of a patch.

    Basically a boxSide with an additional index for the patch.
*/  
struct patchSide: public boxSide
{
public:
    index_t patch;              ///< The index of the patch.
public:

    patchSide() : boxSide(), patch(0) { }

    patchSide(index_t p, boxSide s)
        : boxSide(s), patch(p) { }

    patchSide(index_t p, boundary::side s)
        : boxSide(s), patch(p) { }

    // Accessors
    boxSide& side()       {return *this;}
    const boxSide& side() const {return *this;}

    bool operator== (const patchSide & other) const
    {
        return patch==other.patch && m_index==other.m_index;
    }
    bool operator!= (const patchSide& other) const
    {
        return patch==other.patch || m_index==other.m_index;
    }
    bool operator>  (const patchSide& other) const
    {
        return patch>other.patch || (patch==other.patch && m_index>other.m_index);
    }
    bool operator<  (const patchSide& other) const
    {
        return patch<other.patch || (patch==other.patch && m_index<other.m_index);
    }
    bool operator<= (const patchSide& other) const
    {
        return patch<other.patch || (patch==other.patch && m_index<=other.m_index);
    }
    bool operator>= (const patchSide& other) const
    {
        return patch>other.patch || (patch==other.patch && m_index>=other.m_index);
    }
};

/// Print (as string) a patch side
inline std::ostream &operator<<(std::ostream &os, patchSide const & i)
{
    //os<<"Side: patch="<< i.patch<<", side="<< int(i.side)<<", ie. "<< i.side << ".\n";
    os<<i.patch<<" "<< i.side();
    return os;
}



/**
   @brief Struct which represents a certain corner of a hyper-cube.
**/
struct boxCorner
{
public:
    /** \brief Index of the corner
    *
    * ...stored as number, specified in boundary::corner
    **/
    index_t m_index;
public:
    boxCorner (index_t a=0) : m_index(a) { GISMO_ASSERT(a>=0,"invalid side");}
    boxCorner (boundary::corner a) : m_index(a) { GISMO_ASSERT(a>=0,"invalid side");}
    boxCorner (gsVector<bool>  v) : m_index(1)
    {
        for (index_t p=0;p<v.rows();++p)
        {
            if(v(p))
                m_index+= 1<<p;
        }
    }

    // conversions
    operator boundary::corner() const {return static_cast<boundary::corner>(m_index);}
    boxCorner& operator= (boundary::corner a) {m_index= a; return *this;}

    /**
     * @brief returns a vector of parameters describing the position
     *        of the corner
     * @param dim
     * @param param is 1 if the corner is contained in box_face(i,1)
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
    void getContainingSides (int dim, std::vector<boxSide> &sides) const
    {
        GISMO_ASSERT(dim>=0, "Dimension must be non negative");
        sides.resize(dim);
        gsVector<bool> param;
        parameters_into(dim,param);
        for (index_t i=0;i<dim;++i)
            sides[i]=boxSide(i,param(i));
    }

    /**
     * @brief helper for iterating on corners of an n-dimensional box
     * @param dim
     * @return the first valid corners in an dim-dimensional box
    **/
    static  boxCorner    getFirst     (int dim) {return boxCorner(1);}

    /**
     * @brief helper for iterating on corners of an n-dimensional box
     * @param dim
     * @return the last valid corners in an dim-dimensional box
    **/

    static  boxCorner    getLast      (int dim) {return boxCorner((1<<dim));}
    /**
     * @brief helper for iterating on corners of an n-dimensional box
     * @param dim
     * @return the (invalid) corners after the last one in dim-dimensional box
    **/

    static  boxCorner    getEnd       (int dim) {return boxCorner((1<<dim)+1);}
    /**
     * @brief set to next boxCorner
     */
    boxCorner& operator++ () { ++m_index; return *this;} //prefix
    //boxCorner operator++ (int ) { boxCorner temp(*this); ++m_index; return temp;} //postfix

    /**
     * @brief set to previous boxCorner
     */
    boxCorner& operator-- () { --m_index; return *this;} //prefix
    //boxCorner operator-- (int ) { boxCorner temp(*this); --m_index; return temp;} //postfix

};


/**
    @brief Struct which represents a certain corner of a patch.

*/
struct patchCorner : public boxCorner
{
public:
    index_t patch;
public:
    patchCorner() : boxCorner(0) { }
    patchCorner(int p,boundary::corner c)
        : boxCorner(c), patch (p) { }

    patchCorner(int p, int c)
        : boxCorner(c), patch (p) { }

    patchCorner(int p, boxCorner c)
        : boxCorner(c), patch (p) { }


    bool operator== (const patchCorner& other) const
    {
        return patch==other.patch && m_index==other.m_index;
    }


    /**
     * @brief returns a vector of patchSides that contain this corner
     * @param dim is the ambient dimension
     * @param sides
     */
    void getContainingSides (int dim, std::vector<patchSide> &sides) const
    {
        GISMO_ASSERT(dim>=0, "Dimension must be non negative");
        sides.resize(dim);
        gsVector<bool> param;
        parameters_into(dim,param);
        for (index_t i=0;i<dim;++i)
            sides[i]=patchSide(patch, boxSide(i, param(i)));
    }
};


/** 
    @brief Struct which represents an interface between two patches.
    
    
*/  
struct GISMO_EXPORT boundaryInterface
{
public:
    boundaryInterface() { }

    // special constructor for the 2d case
    boundaryInterface(patchSide const & _ps1,
                      patchSide const & _ps2,
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
    boundaryInterface(patchSide const & _ps1,
                      patchSide const & _ps2,
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
            directionOrientation(o)=true;
            /// TODO: discuss and define default orientation
        }
    }

    boundaryInterface(gsVector<int>     const & p,
                      gsVector<index_t> const & map_info,
                      gsVector<bool>    const & orient_flags)
    : ps1(p(0),p(1)), ps2(p(2),p(3)),
      directionMap(map_info), 
      directionOrientation(orient_flags)
    {  
        GISMO_ASSERT(p.size() == 4, "Expecting four integers");
    }

    boundaryInterface(patchSide const & _ps1,
                      patchSide const & _ps2,
                      gsVector<index_t> const & map_info,
                      gsVector<bool>    const & orient_flags)
        : ps1(_ps1), ps2(_ps2), directionMap(map_info), directionOrientation(orient_flags)
    {  }

    GISMO_DEPRECATED boundaryInterface(patchSide const & _ps1,
                      patchSide const & _ps2,
                      gsVector<bool>    const & orient_flags)
    {
        init(_ps1,_ps2,orient_flags);
    }

    GISMO_DEPRECATED boundaryInterface(gsVector<int>     const & p,
                      gsVector<bool>    const & orient_flags)
    {
        init(patchSide(p(0),boxSide(p(1))),patchSide(p(2),boxSide(p(3))) ,orient_flags);
    }

    //DEPRECATED
    void init (patchSide const & _ps1,
                      patchSide const & _ps2,
                      gsVector<bool>    const & orient_flags)
    {
        ps1=_ps1;
        ps2=_ps2;

        const index_t dim = orient_flags.cols()+1;
        directionMap.resize(dim);
        directionOrientation.resize(dim);

        directionMap(ps1.direction())=ps2.direction();
        directionOrientation(ps1.direction())= (ps1.parameter()!=ps2.parameter());

        directionMap(1-ps1.direction())=1-ps2.direction();
        directionOrientation(1-ps1.direction())= orient_flags(0);

    }


    bool operator== (const boundaryInterface & other) const
    {
        return ps1==other.ps1 && ps2==other.ps2
                && directionMap==other.directionMap
                && directionOrientation==other.directionOrientation;
    }

    /**
     * @brief first, returns the first patchSide of this interface
    **/
          patchSide& first ()        {return ps1;}
    const patchSide& first () const  {return ps1;}

    /**
     * @brief second, returns the second patchSide of this interface
    **/
          patchSide& second ()        {return ps2;}
    const patchSide& second () const  {return ps2;}

    /**
     * @brief Returns the second side if ps is the first side,
     * otherwise returns the second side
    **/
          patchSide&   other (const patchSide & ps)        {return ps==ps1 ? ps2 : ps1;}
    const patchSide&   other (const patchSide & ps) const  {return ps==ps1 ? ps2 : ps1;}

    // use boundaryInterface.first() and boundaryInterface.second()
    GISMO_DEPRECATED patchSide  operator [] (size_t i) const
    {
        if (i==0)
            return ps1;
        else if (i==1)
            return ps2;
        else
            GISMO_ERROR("Invalid index "<<i<<": Interface has 2 elements(sides).");
    }
    // use boundaryInterface.first() and boundaryInterface.second()
    GISMO_DEPRECATED patchSide & operator [] (size_t i)
    {
        if (i==0)
            return ps1;
        else if (i==1)
            return ps2;
        else
            GISMO_ERROR("Invalid index "<<i<<": Interface has 2 elements(sides).");
    }


    // this is a work in progress, the old code was 2D specific, look at mapInfo and mapOrientation
    GISMO_DEPRECATED gsVector<bool> orient () const
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
            result.directionOrientation(directionMap(i))=directionOrientation(i);
        }
        return result;
    }

    patchCorner mapCorner ( const patchCorner c) const
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
            new_par(directionMap(i)) = ps2.parameter();
            for (++i; i<dim;++i)
            {
                new_par(directionMap(i)) = directionOrientation(i) ?
                            par(i) : !par(i);
            }
            return patchCorner(ps2.patch, boxCorner(new_par));
        }
        else if (c.patch == ps2.patch && par(ps2.direction()) == ps2.parameter() )
        {
            return getInverse().mapCorner(c);
        }
        else
        {
            gsWarn<<"cannot map corners that are not in the interface";
            return c;
        }
    }

    bool dirOrientation(const patchSide & ps,index_t dir) const
    {
        if(ps==ps1)
            return directionOrientation(dir);
        else
            return getInverse().dirOrientation(ps,dir);
    }

    index_t dirMap(const patchSide & ps,index_t dir) const
    {
        if(ps==ps1)
            return directionMap(dir);
        else
            return getInverse().dirMap(ps,dir);
    }
    
    /// Accessor for boundaryInterface::directionOrientation
    gsVector<bool> dirOrientation(const patchSide & ps) const
    {
        if(ps==ps1)
            return directionOrientation;
        else
            return getInverse().dirOrientation(ps);
    }

    /// Accessor for boundaryInterface::directionMap
    gsVector<index_t>  dirMap(const patchSide & ps) const
    {
        if(ps==ps1)
            return directionMap;
        else
            return getInverse().dirMap(ps);
    }

    const gsVector<index_t> & dirMap() const
    { return directionMap; }

    const gsVector<bool> & dirOrientation()  const
    { return directionOrientation; }

    
    void matchDofs(gsVector<int>    bSize, 
                   gsMatrix<unsigned> & b1, 
                   const gsMatrix<unsigned> & b2) const;

private:

    patchSide ps1; ///< The first patch side.
    patchSide ps2; ///< The second patch side.

    /**
     * @brief store the combinatorial data about the interface
     *
     *
     * Accessed by dirMap()
     *
     * we describe the permutation and orientation of the coordinate directions
     * through an affine map that puts ps1.patch next to ps2.patch in such a way that
     * ps1 coincide to ps2.
     *
     * \a boundaryInterface::directionMap stores the permutation of the coordinate directions, i.e. the rotation of the map.\n
     * \a boundaryInterface::directionOrientation stores the corresponding orientation.
     *
     *
     * <b>Example 1:</b>
     * \f[
     \mathrm{patch~0}~
     \begin{array}{|cc|cc|}
     \hline
     \uparrow x &     &  \uparrow x  &   &\\
      & \rightarrow y &  & \rightarrow y &\\\hline
     \end{array}
     ~\mathrm{patch~1}
     \f]
     * In Example 1, the image of
     * the x-axis of patch 0 is oriented such that it corresponds to the
     * image of the x-axis of patch 1. Hence, the coordinate direction "0" is
     * mapped to coordinate direction "0", and direction "1" is
     * mapped to direction "1". In this case, directionMap is stored
     * as the vector <em>[0,1]</em>.
     *
     * The orientations of both coordinate directions is the same, so
     * boundaryInterface::directionOrientation is <em>[ 1, 1 ]</em>
     *
     * <b>Example 2:</b>
     * \f[
     \mathrm{patch~0}~
     \begin{array}{|cc|cc|cc|}
     \hline
     \uparrow x &       &               & \rightarrow x \\
         & \rightarrow y & \downarrow y & \\\hline
     \end{array}
     ~\mathrm{patch~1}
     \f]
     *
     * In Example 2, the image of
     * the x-axis of patch 0 is oriented such that it corresponds to the
     * image of the y-asis of patch 1. This means that the coordinate direction "0" is
     * mapped to coordinate direction "1". Also, direction "1" is
     * mapped to direction "0". In this case, directionMap is <em>[1,0]</em>.
     *
     * The orientations of the x-axis of patch 0 and its corresponding
     * counterpart (i.e., the y-axis of patch 1) are reversed. The orientation of
     * the y-axis and its counterpart are the same. Hence,
     * boundaryInterface::directionOrientation is <em>[ 0, 1 ]</em>
     */
    gsVector<index_t> directionMap;
    /** For each coordinate direction we save if the original
    * coordinate and the destination one have the same orientation
    *
    * Accessed by dirOrientation()
    *
    * See boundaryInterface::directionMap for documentation.
    */
    gsVector<bool>    directionOrientation;

    /// TODO: the information could be stored in a single vector of signed integers: the sign gives the orientation
    /// the problem is that it is necessary to shift the indices as there is no -0
protected:
    friend std::ostream &operator<<(std::ostream &os, const boundaryInterface & i);
};


/// Print (as string) an interface
inline std::ostream &operator<<(std::ostream &os, const boundaryInterface & i)
{
    os << "interface between: "<<i.ps1.patch<<":"<< i.ps1.side()<<" and "
       << i.ps2.patch<<":"<<i.ps2.side()<<" [ ";
    for (index_t j = 0; j<i.directionMap.size(); ++j)
    {
        if ( i.ps1.direction() == j ) 
            continue;
        os << j << "~" << (i.directionOrientation(j) ? "(+" : "(-") << i.directionMap(j)<<") ";
    }
    os <<"]";
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
/// use boxSide struct instead of enumerated values
GISMO_DEPRECATED inline int direction (int s)
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
/// use boxSide struct instead of enumerated values
GISMO_DEPRECATED inline bool parameter (int s)
{
    GISMO_ASSERT( s>0, "Requested parameter of none boundary.\n");
    return ( (s+1) % 2 == 0 ? false : true ) ;
}



/// Print (as string) a boundary side
inline std::ostream &operator<<(std::ostream &os, const boxSide& o)
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



/**
 *  \brief get the matrix containing the lower and upper corner of
 *  the specified side of the given box
 */
template <typename T>
gsMatrix<T> getFace (const boxSide side, const gsMatrix<T> &box)
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


} // namespace gismo

