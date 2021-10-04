/** @file gsBoundaryConditions.h

    @brief Provides gsBoundaryConditions class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBoundary.h>


namespace gismo
{

/// @brief Specifies the type of boundary condition
///
/// \ingroup Pde
struct condition_type
{
    /// Specifies the type of boundary condition
    enum type
    {
        unknownType = -1,
        dirichlet = 0, ///< Dirichlet type
        weak_dirichlet = 10, ///< Dirichlet type
        neumann   = 1, ///< Neumann type
        robin     = 2, ///< Robin type
        clamped   = 3, ///< Robin type
        weak_clamped = 30,
        collapsed = 4  ///< Robin type
        //mixed BD means: there are both dirichlet and neumann sides
        //robin: a linear combination of value and derivative
        //cauchy: there are 2 conditions (value+deriv) defined on the same side
    };

};

// Print (as string) a boundary type
inline std::ostream &operator<<(std::ostream &os, const condition_type::type& o)
{
    switch (o)
    {
    case condition_type::dirichlet:
    {
        os<< "Dirichlet";
        break;
    }
    case condition_type::weak_dirichlet:
    {
        os<< "Weak Dirichlet";
        break;
    }
    case condition_type::neumann:
    {
        os<< "Neumann";
        break;
    }
    case condition_type::robin:
    {
        os<< "Robin";
        break;
    }
    case condition_type::clamped:
    {
        os<< "Clamped";
        break;
    }
    case condition_type::weak_clamped:
    {
        os<< "Weak Clamped";
        break;
    }
    case condition_type::collapsed:
    {
        os<< "Collapsed";
        break;
    }
    default:
        os<< "condition type not known.\n";
    };
    return os;
}

/**
    @brief Class that defines a boundary condition for a side of a
    patch for some unknown variable of a PDE.

    \todo rename to boundaryCondition

    \ingroup Pde
*/
template<class T>
struct boundary_condition
{
    typedef typename gsFunction<T>::Ptr function_ptr;

    boundary_condition( int p, boxSide s, const function_ptr & f_shptr,
                        const std::string & label, short_t unknown,
                        short_t unkcomp, bool parametric)
    : ps(p, s),
      m_function(f_shptr),
      m_label(label),
      m_unknown(unknown),
      m_unkcomp(unkcomp),
      m_parametric(parametric)
    {
        if (m_label == "Dirichlet") m_type = condition_type::dirichlet;
        else if (m_label == "Weak Dirichlet") m_type = condition_type::weak_dirichlet;
        else if (m_label == "Neumann")   m_type = condition_type::neumann;
        else if (m_label == "Robin")     m_type = condition_type::robin;
        else if (m_label == "Clamped")   m_type = condition_type::clamped;
        else if (m_label == "Weak Clamped")   m_type = condition_type::weak_clamped;
        else if (m_label == "Collapsed") m_type = condition_type::collapsed;
        else m_type = condition_type::unknownType;
    }

    boundary_condition( int p, boxSide s, const function_ptr & f_shptr,
                        condition_type::type t, short_t unknown, bool parametric)
    : ps(p, s),
      m_function(f_shptr),
      m_type(t),
      m_unknown(unknown),
      m_unkcomp(-1),
      m_parametric(parametric)
    {
        switch (t)
        {
        case condition_type::dirichlet:
        {
            m_label = "Dirichlet";
            break;
        }
        case condition_type::weak_dirichlet:
        {
            m_label = "Weak Dirichlet";
            break;
        }
        case condition_type::neumann:
        {
            m_label = "Neumann";
            break;
        }
        case condition_type::robin:
        {
            m_label = "Robin";
            break;
        }
        case condition_type::clamped:
        {
            m_label = "Clamped";
            break;
        }
        case condition_type::weak_clamped:
        {
            m_label = "weak Clamped";
            break;
        }
        case condition_type::collapsed:
        {
            m_label = "Collapsed";
            break;
        }
        default:
            m_label = "Unknown";
            break;
        };
    }

    boundary_condition( int p, boxSide s, const function_ptr & f_shptr,
                        condition_type::type t, int unknown, int unkcomp, bool parametric)
    : ps(p, s),
      m_function(f_shptr),
      m_type(t),
      m_unknown(unknown),
      m_unkcomp(unkcomp),
      m_parametric(parametric)
    {
        switch (t)
        {
        case condition_type::dirichlet:
        {
            m_label = "Dirichlet";
            break;
        }
        case condition_type::weak_dirichlet:
        {
            m_label = "Weak Dirichlet";
            break;
        }
        case condition_type::neumann:
        {
            m_label = "Neumann";
            break;
        }
        case condition_type::robin:
        {
            m_label = "Robin";
            break;
        }
        case condition_type::clamped:
        {
            m_label = "Clamped";
            break;
        }
        case condition_type::weak_clamped:
        {
            m_label = "Weak Clamped";
            break;
        }
        case condition_type::collapsed:
        {
            m_label = "Collapsed";
            break;
        }
        default:
            m_label = "Unknown";
            break;
        };
    }

    /// Returns true if there is no function data (homogeneous condition)
    bool isHomogeneous() const { return m_function.get() == NULL; }

    /// Returns the function data pointer of the boundary condition
    function_ptr function() const { return m_function; }

    // Returns a reference to the function data
    //const gsFunction<T> & function() const { return *m_function; }

    /// Returns the type of the boundary condition
    condition_type::type  type() const { return m_type; }

    /// Returns the type of the boundary condition
    const std::string & ctype() const { return m_label; }

    /// Returns the patch to which this boundary condition refers to
    index_t patch()    const { return ps.patch; }

    /// Returns the side to which this boundary condition refers to
    boxSide side()     const { return ps.side(); }

    /// Returns the unknown to which this boundary condition refers to
    short_t     unknown()  const { return m_unknown; }

    /// Returns the component of the unknown which this boundary condition refers to
    short_t     unkComponent()  const { return m_unkcomp; }

    /// Returns true if the function data for this boundary condition
    /// is defined in parametric coordinates
    bool    parametric()  const { return m_parametric; }


    patchSide ps;                ///< Side of a patch for this boundary condition

    function_ptr m_function;     ///< Function data for this boundary condition

    // TO DO : robin coefficients?

    condition_type::type m_type;// todo: remove

    std::string m_label;         ///< Description of type of the boundary condition

    short_t m_unknown;               ///< Unknown to which this boundary condition refers to

    short_t m_unkcomp;               ///< Component of unknown to which this boundary condition refers to

    bool m_parametric;
};

/**
    @brief Class prescribing a value related to a corner of a patch
*/
template<class T>
struct corner_value
{
    corner_value(index_t p, boxCorner c, T v, short_t unk = 0)
        : patch(p), corner(c), value(v), unknown(unk) { }

    index_t patch;     ///< The index of the patch.
    boxCorner corner; ///< The corner
    T value;          ///< The value
    short_t   unknown;    ///< Unknown to which this boundary condition refers to
};

/** @brief
    Class containing a set of  boundary conditions.

    The boundary conditions are stored in the form of a list of boundary_condition
    instances.

    \ingroup Pde
*/
template<class T>
class GISMO_EXPORT gsBoundaryConditions
{

public:

    typedef typename std::deque<boundary_condition<T> > bcContainer;
    typedef typename bcContainer::iterator iterator;
    typedef typename bcContainer::const_iterator const_iterator;

    typedef typename std::deque<corner_value<T> >       cornerContainer;
    typedef typename cornerContainer::iterator citerator;
    typedef typename cornerContainer::const_iterator const_citerator;

    typedef typename std::deque<boundaryInterface> ppContainer;
    typedef typename ppContainer::iterator ppiterator;
    typedef typename ppContainer::const_iterator const_ppiterator;

    // Format: std::pair<type,bcContainer>
    typedef std::map<std::string,bcContainer> bcData;
    typedef typename bcData::iterator bciterator;
    typedef typename bcData::const_iterator const_bciterator;

    typedef memory::shared_ptr<gsBoundaryConditions> Ptr;
    typedef memory::unique_ptr<gsBoundaryConditions> uPtr;

    typedef typename boundary_condition<T>::function_ptr function_ptr;

    typedef std::list<util::reference_wrapper<const boundary_condition<T> > > bcRefList;
public:

    /*
    gsBoundaryConditions & operator= (uPtr other)
    {
        if ( other.get() != NULL )
        {
            this->swap(*other);
            other.reset();
        }
        return *this;
    }
    */

    void swap(gsBoundaryConditions & other)
    {
        m_bc.swap(other.m_bc);
        corner_values.swap(other.corner_values);
        m_periodicPairs.swap(other.m_periodicPairs);
    }

public:

    void clear()
    {
        m_bc.clear();
        corner_values.clear();
        m_periodicPairs.clear();
    }

    size_t size() const
    {
        size_t sz = 0;
        for (typename bcData::const_iterator it = m_bc.begin(); it != m_bc.end(); ++it)
            sz += it->second.size();
        return sz + corner_values.size();
    }

    /// Return a reference to boundary conditions of certain type
    const bcContainer & container(const std::string & label) const {return m_bc[label]; }

    /// Return a reference to boundary conditions of certain type for
    /// unknown \a unk
    bcRefList get(const std::string & label, const short_t unk = 0) const
    {
        bcRefList result;
        const const_bciterator it = m_bc.find(label);
        if ( it != m_bc.end() )
            for (const_iterator c = it->second.begin(); c!= it->second.end(); ++c)
                if ( c->m_unknown == unk )
                    result.push_back(*c);
        return result;
    }

    /// Return a reference to the Dirichlet sides
    const bcContainer & dirichletSides() const {return m_bc["Dirichlet"]; }

    /// Return a reference to the Weak Dirichlet sides
    const bcContainer & weakDirichletSides() const {return m_bc["Weak Dirichlet"]; }

    /// Return a reference to the Neumann sides
    const bcContainer & neumannSides()   const {return m_bc["Neumann"]; }

    /// Return a reference to the Robin sides
    const bcContainer & robinSides()     const {return m_bc["Robin"]; }

    const cornerContainer & cornerValues() const  {return corner_values;  }

    /// Extracts the BC, comming from a certain component.
    bcContainer reducedContainer(const bcContainer & container, short_t unknown) const
    {
        bcContainer red;
        //red.reserve(container.size());
        for(typename bcContainer::const_iterator iter=container.begin(); iter!=container.end();++iter)
        {
            if(iter->unknown()==unknown)
                red.push_back(*iter);
        }
        return red;
    }

    bcContainer allConditions() const
    {
        bcContainer all;
        //all.reserve( size() - corner_values.size() );
        for (typename bcData::const_iterator it = m_bc.begin(); it != m_bc.end(); ++it)
            all.insert( all.end(), it->second.begin(), it->second.end() );
        return all;
    }

    /// Returns a const-iterator to the beginning of the Bc container of type \a label
    const_iterator begin(const std::string & label) const {return m_bc[label].begin(); }

    /// Returns an iterator to the beginning of the Bc container of type \a label
    iterator begin(const std::string & label) { return m_bc[label].begin(); }

    /// Returns a const-iterator to the end of the Bc container of type \a label
    const_iterator end(const std::string & label) const {return m_bc[label].end(); }

    /// Returns an iterator to the end of the Bc container of type \a label
    iterator end(const std::string & label) { return m_bc[label].end(); }

    const_bciterator beginAll() const {return m_bc.begin(); }
    bciterator beginAll() {return m_bc.begin(); }

    const_bciterator endAll() const {return m_bc.end(); }
    bciterator endAll() {return m_bc.end(); }

    /// Get a const-iterator to the beginning of the Dirichlet sides
    /// \return an iterator to the beginning of the Dirichlet sides
    const_iterator dirichletBegin() const
    { return m_bc["Dirichlet"].begin(); }

    /// Get a const-iterator to the end of the Dirichlet sides
    /// \return an iterator to the end of the Dirichlet sides
    const_iterator dirichletEnd() const
    { return m_bc["Dirichlet"].end(); }

    /// Get an iterator to the beginning of the Dirichlet sides
    /// \return an iterator to the beginning of the Dirichlet sides
    iterator dirichletBegin()
    { return m_bc["Dirichlet"].begin(); }

    /// Get an iterator to the end of the Dirichlet sides
    /// \return an iterator to the end of the Dirichlet sides
    iterator dirichletEnd()
    { return m_bc["Dirichlet"].end(); }

    /// Get a const-iterator to the beginning of the Weak Dirichlet sides
    /// \return an iterator to the beginning of the Weak Dirichlet sides
    const_iterator weakDirichletBegin() const
    { return m_bc["Weak Dirichlet"].begin(); }

    /// Get a const-iterator to the end of the Weak Dirichlet sides
    /// \return an iterator to the end of the Weak Dirichlet sides
    const_iterator weakDirichletEnd() const
    { return m_bc["Weak Dirichlet"].end(); }

    /// Get an iterator to the beginning of the Weak Dirichlet sides
    /// \return an iterator to the beginning of the Weak Dirichlet sides
    iterator weakDirichletBegin()
    { return m_bc["Weak Dirichlet"].begin(); }

    /// Get an iterator to the end of the Weak Dirichlet sides
    /// \return an iterator to the end of the Weak Dirichlet sides
    iterator weakDirichletEnd()
    { return m_bc["Weak Dirichlet"].end(); }


    /// Get a const-iterator to the beginning of the Neumann sides
    /// \return an iterator to the beginning of the Neumann sides
    const_iterator neumannBegin() const
    { return m_bc["Neumann"].begin(); }

    /// Get a const-iterator to the end of the Neumann sides
    /// \return an iterator to the end of the Neumann sides
    const_iterator neumannEnd() const
    { return m_bc["Neumann"].end(); }

    /// Get an iterator to the beginning of the Neumann sides
    /// \return an iterator to the beginning of the Neumann sides
    iterator neumannBegin()
    { return m_bc["Neumann"].begin(); }

    /// Get an iterator to the end of the Neumann sides
    /// \return an iterator to the end of the Neumann sides
    iterator neumannEnd()
    { return m_bc["Neumann"].end(); }

    /// Get a const-iterator to the beginning of the Robin sides
    /// \return an iterator to the beginning of the Robin sides
    const_iterator robinBegin() const
    { return m_bc["Robin"].begin(); }

    /// Get a const-iterator to the end of the Robin sides
    /// \return an iterator to the end of the Robin sides
    const_iterator robinEnd() const
    { return m_bc["Robin"].end(); }

    /// Get an iterator to the beginning of the corner values
    /// \return an iterator to the beginning of the corner values
    const_citerator cornerBegin() const
    { return corner_values.begin(); }

    /// Get an iterator to the end of corner values
    /// \return an iterator to the end of the corner values
    const_citerator cornerEnd() const
    { return corner_values.end(); }

    /// Get an iterator to the beginning of the Robin sides
    /// \return an iterator to the beginning of the Robin sides
    iterator robinBegin()
    { return m_bc["Robin"].begin(); }

    /// Get an iterator to the end of the Robin sides
    /// \return an iterator to the end of the Robin sides
    iterator robinEnd()
    { return m_bc["Robin"].end(); }

    /// Get an iterator to the beginning of the corner values
    /// \return an iterator to the beginning of the corner values
    citerator cornerBegin()
    { return corner_values.begin(); }

    /// Get an iterator to the end of corner values
    /// \return an iterator to the end of the corner values
    citerator cornerEnd()
    { return corner_values.end(); }

    void add(int p, boxSide s, const std::string & label,
             const function_ptr & f_ptr, short_t unknown = 0,
             int comp = -1, bool parametric = false)
    {
        m_bc[label].push_back(
            boundary_condition<T>(p, s, f_ptr, label, unknown, comp, parametric) );
    }

    void add(int p, boxSide s, const std::string & label,
             gsFunction<T> * f, short_t unknown = 0,
             int comp = -1, bool parametric = false)
    {
        function_ptr f_ptr = memory::make_shared_not_owned(f);
        m_bc[label].push_back(
            boundary_condition<T>(p, s, f_ptr, label, unknown, comp, parametric) );
    }

    void add(int p, boxSide s, const std::string & label,
             const gsFunction<T> & f, short_t unknown = 0,
             int comp = -1, bool parametric = false)
    {
        function_ptr fun = memory::make_shared(f.clone().release());
        add(p,s,label,fun,unknown,comp,parametric);
    }

    /** \brief Adds another boundary condition
     *
     * Creates an object of type boundary_condition and adds is to
     * the list of corresponding boundary conditions.
     *
     * \param p Index of the patch
     * \param s Side of the patch
     * \param t Type of boundary condition (see condition_type::type)
     * \param f Function defining the boundary condition
     * \param unknown Specifies which unknown variable the boundary condition
     * refers to (to be used if more than one variable is unknown, e.g.,
     * velocity and pressure)
     * \param parametric True if the function data for this boundary condition
     * is defined in parametric coordinates.
     */
    void addCondition(int p, boxSide s, condition_type::type t,
                      gsFunction<T> * f, short_t unknown = 0, bool parametric = false, int comp = -1)
    {
        function_ptr fun = memory::make_shared_not_owned(f);
        addCondition(p,s,t,fun,unknown,parametric,comp);
    }

    void addCondition(int p, boxSide s, condition_type::type t,
                      const function_ptr & f_shptr, short_t unknown = 0,
                      bool parametric = false, int comp = -1)
    {
        switch (t)
        {
        case condition_type::dirichlet :
            // this->add(p,s,f_shptr,"Dirichlet",unknown,comp,parametric);
            m_bc["Dirichlet"].push_back( boundary_condition<T>(p,s,f_shptr,t,unknown,comp,parametric) );
            break;
        case condition_type::weak_dirichlet :
            // this->add(p,s,f_shptr,"Dirichlet",unknown,comp,parametric);
            m_bc["Weak Dirichlet"].push_back( boundary_condition<T>(p,s,f_shptr,t,unknown,comp,parametric) );
            break;
        case condition_type::neumann :
            m_bc["Neumann"].push_back( boundary_condition<T>(p,s,f_shptr,t,unknown,comp,parametric) );
            break;
        case condition_type::robin :
            m_bc["Robin"].push_back( boundary_condition<T>(p,s,f_shptr,t,unknown,comp,parametric) );
            break;
        case condition_type::clamped :
            m_bc["Clamped"].push_back( boundary_condition<T>(p,s,f_shptr,t,unknown,comp,parametric) );
            break;
        case condition_type::weak_clamped :
            m_bc["Weak Clamped"].push_back( boundary_condition<T>(p,s,f_shptr,t,unknown,comp,parametric) );
            break;
        case condition_type::collapsed :
            m_bc["Collapsed"].push_back( boundary_condition<T>(p,s,f_shptr,t,unknown,comp,parametric) );
            break;
        default:
            gsWarn<<"gsBoundaryConditions: Unknown boundary condition.\n";
        }
    }

    void addCondition(int p, boxSide s, condition_type::type t,
                      const gsFunction<T> & func, short_t unknown = 0,
                      bool parametric = false, int comp = -1)
    {
        function_ptr fun(func.clone().release());
        addCondition(p,s,t,fun,unknown,parametric,comp);
    }

    void addCondition( boxSide s, condition_type::type t,
                       gsFunction<T> * f, short_t unknown = 0, bool parametric = false, int comp = -1)
    {
        // for single-patch only
        addCondition(0,s,t,f,unknown,parametric,comp);
    }

    void addCondition(const patchSide& ps, condition_type::type t,
                      gsFunction<T> * f, short_t unknown = 0, bool parametric = false, int comp = -1)
    {
        addCondition(ps.patch, ps.side(), t, f, unknown,parametric,comp);
    }

    void addCondition(const patchSide& ps, condition_type::type t,
                      const function_ptr & f_shptr, short_t unknown = 0, bool parametric = false, int comp = -1)
    {
        addCondition(ps.patch, ps.side(), t, f_shptr, unknown,parametric,comp);
    }

    void addCondition(const patchSide& ps, condition_type::type t,
                      const gsFunction<T> & func, short_t unknown = 0, bool parametric = false, int comp = -1)
    {
        addCondition(ps.patch, ps.side(), t, func, unknown,parametric,comp);
    }

    void addCornerValue(boxCorner c, T value, int p = 0, short_t unknown = 0)
    {
        corner_values.push_back( corner_value<T>(p,c,value,unknown) );
    }

    /// Prints the object as a string.
    std::ostream & print(std::ostream &os) const
    {
        //os << "gsBoundaryConditions :\n";
        for (typename bcData::const_iterator it = m_bc.begin(); it != m_bc.end(); ++it)
            os << "* "<<std::setw(13)<<std::left<<it->first<<" : "<< it->second.size() <<"\n";

        if (!corner_values.empty())
            os << "* Corner values : "<< corner_values.size() <<"\n";

        return os;
    }

    /**
     * @brief   getSideCondition
     * @param   ps the patch side

     * @return the boundary condition associated to ps or NULL if no
     * condition is associated to ps

     It is the task of the user of this function to check if the
     returned pointer is NULL.

     Do not use this function if you want to apply boundary conditions during matrix assembly.
     Instead, iterate over all conditions of the type you need (eg. Neumann, Dirichlet)

     */
    const boundary_condition<T>* getConditionFromSide (patchSide ps) const
    {
        const_iterator beg, end, cur;
        patchSideComparison psRef(ps);
        beg = dirichletBegin();
        end = dirichletEnd();
        cur=std::find_if(beg,end,psRef);
        if (cur != end)
            return &(*cur);
        beg = neumannBegin();
        end = neumannEnd();
        cur=std::find_if(beg,end,psRef);
        if (cur != end)
            return &(*cur);
        beg = robinBegin();
        end = robinEnd();
        cur = std::find_if(beg,end,psRef);
        if (cur != end)
            return &(*cur);

        return NULL;
    }

    /**
     * @brief   getConditionFromSide returns the boundary conditions associated to the given patch side
     * @param[in] ps patch side
     * @param[out] result bcContainer containing the boundary conditions associated to ps
     */

    void getConditionFromSide (patchSide ps, bcContainer& result) const
    {
        result.clear();
        const_iterator beg, end, cur;

        beg = dirichletBegin();
        end = dirichletEnd();
        for(cur=beg; cur!=end; cur++)
            if(cur->ps == ps)
                result.push_back(*cur);

        beg = neumannBegin();
        end = neumannEnd();
        for(cur=beg; cur!=end; cur++)
            if(cur->ps == ps)
                result.push_back(*cur);

        beg = robinBegin();
        end = robinEnd();
        for(cur=beg; cur!=end; cur++)
            if(cur->ps == ps)
                result.push_back(*cur);
    }


    /**
     * @brief returns the set of all boundary conditions which refer to patch \a np
     * @param[in] np the patch index
     * @param[out] result the new set of boundary conditions
     */
    void getConditionsForPatch(const index_t np, gsBoundaryConditions& result) const
    {
        result.clear();
        bcContainer bc_all = allConditions(); //inefficient, but fewer code
        for(const_iterator it = bc_all.begin(); it!= bc_all.end();it++)
        {
            if((*it).patch()==np)
            {
                if(it->type() == condition_type::dirichlet || it->type() == condition_type::neumann || it->type() == condition_type::robin)
                    result.addCondition(0,(*it).side(),(*it).type(),(*it).function(),(*it).unknown());
                else
                   result.add(0,(*it).side(),it->ctype(),(*it).function(),(*it).unknown());
            }
        }

        for(const_citerator it = cornerBegin(); it!= cornerEnd();it++)
        {
            if((*it).patch==np)
                result.addCornerValue( (*it).corner, (*it).value,  0, (*it).unknown);
        }
    }

    // Periodic conditions

    /// Get number of periodic pairs
    size_t numPeriodic() const { return m_periodicPairs.size(); }

    /// Return a reference to the periodic sides
    const ppContainer & periodicPairs() const {return m_periodicPairs; }

    /// Get a const-iterator to the beginning of the periodic sides
    /// \return an iterator to the beginning of the periodic sides
    const_ppiterator periodicBegin() const
    { return m_periodicPairs.begin(); }

    /// Get a const-iterator to the end of the periodic sides
    /// \return an iterator to the end of the periodic sides
    const_ppiterator periodicEnd() const
    { return m_periodicPairs.end(); }

    /// Get an iterator to the beginning of the periodic sides
    /// \return an iterator to the beginning of the periodic sides
    ppiterator periodicBegin()
    { return m_periodicPairs.begin(); }

    /// Get an iterator to the end of the periodic sides
    /// \return an iterator to the end of the periodic sides
    ppiterator periodicEnd()
    { return m_periodicPairs.end(); }

    /// Add a periodic condition between side \a s1 of box \a p1 and side \a s2 of box \a p2.
    void addPeriodic(int p1, boxSide s1, int p2, boxSide s2, short_t dim)
    { m_periodicPairs.push_back( boundaryInterface(patchSide(p1,s1), patchSide(p2,s2), dim) ); }

    /// Removes all periodic pairs
    void clearPeriodicPairs() { m_periodicPairs.clear(); }

    /// Set transformation matrix for the periodic pairs of sides
    void setTransformMatrix(gsMatrix<T> trMatrix)
    { m_trMatrix = trMatrix; }

    /// Set identity transformation matrix for the periodic pairs of sides
    void setIdentityMatrix(short_t dim)
    { m_trMatrix = gsMatrix<T>::Identity(dim, dim); }

    /// Get transformation matrix for the periodic pairs of sides
    gsMatrix<T> getTransformMatrix() const
    {
        GISMO_ASSERT(m_trMatrix.rows() > 0, "Transformation matrix for periodic conditions not set!");
        return m_trMatrix;
    }

    /// Set the geometry map to evaluate boundary conditions.
    void setGeoMap(const gsFunctionSet<T> & gm)
    {
      m_patches = &gm;
    }

    /// Returns the geometry map
    const gsFunctionSet<T> & geoMap() const
    {
        GISMO_ASSERT(nullptr!=m_patches, "Geometry map was not provided in BC.");
        return *m_patches;
    }

private: // Data members
    struct patchSideComparison
    {
        const patchSide m_ps;
        patchSideComparison(patchSide ps)
            : m_ps(ps)
        {}

        bool operator() (const boundary_condition<T> &bc) const
        {
            return bc.ps==m_ps;
        }
    };

    cornerContainer corner_values; ///< List of corners with fixed value

    mutable bcData m_bc;  ///< Containers for BCs of various types

    ppContainer m_periodicPairs; // TODO: add read from xml
    gsMatrix<T> m_trMatrix;

    // Pointer to associated multipatch domain
    const gsFunctionSet<T> * m_patches;

}; // class gsBoundaryConditions

/// Print (as string)
template<class T>
std::ostream &operator<<(std::ostream &os, const gsBoundaryConditions<T>& bvp)
{return bvp.print(os); }

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBoundaryConditions.hpp)
#endif
