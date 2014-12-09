
#pragma once

#include <iostream>
#include <map>

#include <gsCore/gsMultiPatch.h>
#include <gsPde/gsPde.h>
#include <gsPde/gsBoundaryConditions.h>


namespace gismo
{

    template <class T> class gsFunction;
    
/** @brief
    Full specification of a boundary value problem including the PDE,
    computational domain, and boundary conditions.
    
    The PDE is stored in the form of an instance of gsPde.
    The computational domain may be given either as a single gsGeometry patch
    or as a collection of patches in the form of a gsMultiPatch.
    (Internally, it is always a gsMultiPatch that is stored.)
    The boundary conditions are stored in the form of a list of boundary_condition
    instances.
*/

template<class T>
class gsBVProblem 
{
    
public:

    typedef typename std::vector<boundary_condition<T> >::iterator iterator;
    typedef typename std::vector<boundary_condition<T> >::const_iterator const_iterator;

public:

    /// Default empty constructor
    gsBVProblem() { }

    /// Constructor for single-patch PDE (without boundary conditions)
    gsBVProblem( gsGeometry<T> * geo, gsPde<T> * pde ) : m_patches(*geo), m_pde(pde)
        { }

    /// Constructor for multi-patch PDE (without boundary conditions)
    gsBVProblem( const gsMultiPatch<T>& mp, gsPde<T> * pde ) : m_patches(mp), m_pde(pde)
        { }

    /// Constructor for multi-patch PDE (without boundary conditions)
    gsBVProblem( gsMovable< gsMultiPatch<T> > mp, gsPde<T> * pde ) : m_pde(pde)
    {
        m_patches.swap( mp.ref() );
    }

    ~gsBVProblem() // Destructor
    {
        delete m_pde;
    }

public:

    gsMultiPatch<T> & patches() { return m_patches; }

    const gsMultiPatch<T> & patches() const { return m_patches; }

    gsPde<T> & pde() const      { return *m_pde; }

    gsBoundaryConditions<T> & boundaryConditions() { return m_BCs; }
    const gsBoundaryConditions<T> & boundaryConditions() const { return m_BCs; }


    int nPatches() const {return m_patches.nPatches(); }

    gsGeometry<T> & patch(unsigned i) const {return m_patches.patch(i); }
    
    std::vector<boundary_condition<T> > dirichletSides()
    {return m_BCs.dirichletSides(); }

    std::vector<boundary_condition<T> > neumannSides()  
    {return m_BCs.neumannSides();   }

    std::vector<boundary_condition<T> > robinSides()    
    {return m_BCs.robinSides();  }
    
    std::vector<boundary_condition<T> > allConditions()
    { return m_BCs.allConditions(); }
    
    /// Get a const-iterator to the beginning of the Dirichlet sides
    /// \return an iterator to the beginning of the Dirichlet sides
    const_iterator dirichletBegin() const
	{ return m_BCs.dirichletBegin(); }
    
    /// Get a const-iterator to the end of the Dirichlet sides
    /// \return an iterator to the end of the Dirichlet sides
    const_iterator dirichletEnd() const
	{ return m_BCs.dirichletEnd(); }
    
    /// Get an iterator to the beginning of the Dirichlet sides
    /// \return an iterator to the beginning of the Dirichlet sides
    iterator dirichletBegin()
	{ return m_BCs.dirichletBegin(); }
    
    /// Get an iterator to the end of the Dirichlet sides
    /// \return an iterator to the end of the Dirichlet sides
    iterator dirichletEnd()
    { return m_BCs.dirichletEnd(); }

    /// Get a const-iterator to the beginning of the Neumann sides
    /// \return an iterator to the beginning of the Neumann sides
    const_iterator neumannBegin() const
	{ return m_BCs.neumannBegin(); }
    
    /// Get a const-iterator to the end of the Neumann sides
    /// \return an iterator to the end of the Neumann sides
    const_iterator neumannEnd() const
	{ return m_BCs.neumannEnd(); }
    
    /// Get an iterator to the beginning of the Neumann sides
    /// \return an iterator to the beginning of the Neumann sides
    iterator neumannBegin()
	{ return m_BCs.neumannBegin(); }
    
    /// Get an iterator to the end of the Neumann sides
    /// \return an iterator to the end of the Neumann sides
    iterator neumannEnd()
	{ return m_BCs.neumannEnd(); }

    /// Get a const-iterator to the beginning of the Robin sides
    /// \return an iterator to the beginning of the Robin sides
    const_iterator robinBegin() const
    { return m_BCs.robinBegin(); }
    
    /// Get a const-iterator to the end of the Robin sides
    /// \return an iterator to the end of the Robin sides
    const_iterator robinEnd() const
	{ return m_BCs.robinEnd(); }
    
    /// Get an iterator to the beginning of the Robin sides
    /// \return an iterator to the beginning of the Robin sides
    iterator robinBegin()
	{ return m_BCs.robinBegin(); }
    
    /// Get an iterator to the end of the Robin sides
    /// \return an iterator to the end of the Robin sides
    iterator robinEnd()
	{ return m_BCs.robinEnd(); }

    // Returns the type of condition for the boundary side s of patch p
    condition_type::type typeOf(boxSide const &s, const int& p = 0) const
    {
        return m_BCs.typeOf(s,p);
    }
  
    void addCondition(int p, boxSide s, condition_type::type t, gsFunction<T> * f)
    {
        m_BCs.addCondition(p,s,t,f);
    }

    void addCondition( boxSide s, condition_type::type t, gsFunction<T> * f)
    {
        // for single-patch only
        GISMO_ASSERT( m_patches.size() == 1, 
        "Tried to use addCondition without indicating the patch index." );
        addCondition(0,s,t,f);
    }

    void addCondition(const patchSide& ps, condition_type::type t, gsFunction<T> * f)
    {
        addCondition(ps.patch, ps.side(), t, f);
    }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { 
        os << "gsBVProblem (Boundary Value Problem):\n";
        os << "* PDE: "    << *m_pde    ;
        os << "* Domain: "<< m_patches ;
        os << m_BCs ;
        return os; 
    }

// Data members
private:

    gsMultiPatch<T> m_patches;

    gsPde<T> * m_pde;

    gsBoundaryConditions<T> m_BCs;

}; // class gsBVProblem


//////////////////////////////////////////////////
//////////////////////////////////////////////////

/// Print (as string)
template<class T>
std::ostream &operator<<(std::ostream &os, const gsBVProblem<T>& bvp)
{return bvp.print(os); };
    
}; // namespace gismo
