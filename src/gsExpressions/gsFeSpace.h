/** @file gsFeSpace.h

    @brief Defines a space expression

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/

#pragma once

#include <gsAssembler/gsDirichletValues.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsExpressions/gsFeSpaceData.h>

namespace gismo
{
namespace expr
{

/**
   Expression for finite element variable in an isogeometric function
   space

   This corresponds to an FE variable that
   belongs to an isogeometric function space
*/
template<class T>
class gsFeSpace :public symbol_expr< gsFeSpace<T> >
{
    friend class gsNullExpr<T>;

protected:
    typedef symbol_expr< gsFeSpace<T> > Base;

    // contains id, mapper, fixedDofs, etc
    gsFeSpaceData<T> * m_sd;

public:
    enum{Space = 1, ScalarValued=0, ColBlocks=0};// test space

    typedef const gsFeSpace Nested_t; //no ref

    typedef T Scalar;

    const gsFeSpace<T> & rowVar() const {return *this;}

    gsDofMapper & mapper()
    {
        GISMO_ASSERT(NULL!=m_sd, "Space/mapper not properly initialized.");
        return m_sd->mapper;
    }

    const gsDofMapper & mapper() const
    {return const_cast<gsFeSpace*>(this)->mapper();}

    inline const gsMatrix<T> & fixedPart() const {return m_sd->fixedDofs;}
    gsMatrix<T> & fixedPart() {return m_sd->fixedDofs;}

    index_t   id() const { return (m_sd ? m_sd->id : -101); }
    void setSpaceData(gsFeSpaceData<T>& sd) {m_sd = &sd;}

    index_t   interfaceCont() const {return m_sd->cont;}
    index_t & setInterfaceCont(const index_t _r) const
    {
        GISMO_ASSERT(_r>-2 && _r<1, "Invalid or not implemented (r="<<_r<<").");
        return m_sd->cont = _r;
    }

    gsFeSolution<T> function(const gsMatrix<T>& solVector) const
    { return gsFeSolution<T>(*this); }

    void getCoeffs(const gsMatrix<T>& solVector, gsMatrix<T> & result,
                   const index_t p = 0) const
    {
        const index_t dim = this->dim();

        const gsMultiBasis<T> & mb = static_cast<const gsMultiBasis<T>&>(this->source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<T>*>(&this->source()), "error");

        // Reconstruct solution coefficients on patch p
        const index_t sz  = mb[p].size();
        result.resize(sz, dim!=1 ? dim : solVector.cols()); // (!)

        for (index_t c = 0; c!=dim; c++) // for all components
        {
            // loop over all basis functions (even the eliminated ones)
            for (index_t i = 0; i < sz; ++i)
            {
                const int ii = m_sd->mapper.index(i, p, c);
                if ( m_sd->mapper.is_free_index(ii) ) // DoF value is in the solVector
                  for(index_t s = 0; s != solVector.cols(); ++s )
                    result(i,c+s) = solVector(ii,s); //assume dim==1 xor solVector.cols()==1
                else // eliminated DoF: fill with Dirichlet data
                    result(i,c) =  m_sd->fixedDofs.at( m_sd->mapper.global_to_bindex(ii) );
            }
        }
    }

    // space restrictTo(boundaries);
    // space restrictTo(bcRefList domain);

    void setupMapper(gsDofMapper dofsMapper) const
    {
        GISMO_ASSERT( dofsMapper.isFinalized(), "The provided dof-mapper is not finalized.");
        GISMO_ASSERT( dofsMapper.mapSize()==static_cast<size_t>(this->source().size()*dofsMapper.numComponents()), "The dof-mapper is not consistent: mapSize()="<<dofsMapper.mapSize()<<"!="<<static_cast<size_t>(this->source().size())<<"=this->source().size()");
        m_sd->mapper = give(dofsMapper);
    }

    void setup(const index_t _icont = -1) const
    {
        this->setInterfaceCont(_icont);
        m_sd->mapper = gsDofMapper();

        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(&this->source()) )
        {
            m_sd->mapper = gsDofMapper(*mb, this->dim() );
            //m_mapper.init(*mb, this->dim()); //bug
            if ( 0==this->interfaceCont() ) // Conforming boundaries ?
            {
                for ( gsBoxTopology::const_iiterator it = mb->topology().iBegin();
                      it != mb->topology().iEnd(); ++it )
                {
                    mb->matchInterface(*it, m_sd->mapper);
                }
            }
        }

        if (const gsMappedBasis<2,T> * mb =
            dynamic_cast<const gsMappedBasis<2,T>*>(&this->source()) )
        {
            m_sd->mapper.setIdentity(mb->nPatches(), mb->size() , this->dim());
        }

        m_sd->mapper.finalize();
    }

    void setup(const gsBoundaryConditions<T> & bc, const index_t dir_values,
               const index_t _icont = -1) const
    {
        this->setInterfaceCont(_icont);
        m_sd->mapper = gsDofMapper();
        const index_t dim = this->dim();
        const gsMultiBasis<T> *mb = dynamic_cast<const gsMultiBasis<T> *>(&this->source());
        if (mb != nullptr)
        {
            m_sd->mapper = gsDofMapper(*mb,bc, dim,this->id(),(0 == this->interfaceCont())); // second to last index assumes the unknown in the BCs is the same as the space ID
        }
        else if (const gsBasis<T> *b =
                   dynamic_cast<const gsBasis<T> *>(&this->source()))
        {
            m_sd->mapper = gsDofMapper(*b,bc, dim,this->id()); // last index assumes the unknown in the BCs is the same as the space ID
        }
        else if (const gsMappedBasis<2, T> *mapb =
                   dynamic_cast<const gsMappedBasis<2, T> *>(&this->source()))
        {
            m_sd->mapper = gsDofMapper(*mapb, bc, dim, this->id(),false);

            /* 
             * NOTE (from @hverhelst): In earlier code (before 18-12-2025)
             * the clamped BCs were applied here for MappedBasis.
             * However, the following line contained a bug, leading to
             * slightly different results compared to the one using 
             * the current gsDofMapper initialization with BCs above.
                    
                    for (index_t c = 0; c!=dim; c++) // for all components
                    {
                        if (c==cc || cc==-1 )
                            for (index_t k = 0; k < bnd.size() - 1; ++k)  //<<-- the -1 here should not be there
                                m_sd->mapper.matchDof(  it->ps.patch, (bnd)(k, 0),
                                                        it->ps.patch, (bnd1)(k, 0), c);
                    }
            */
        }
            
        m_sd->mapper.finalize();
        // Compute Dirichlet node values
        gsDirichletValues(bc, dir_values, *this);
    }

    void print(std::ostream &os) const { os << "u"; }

protected:
    friend class gismo::gsExprHelper<Scalar>;
    friend class symbol_expr<gsFeSpace>;
    explicit gsFeSpace(index_t _d = 1) : Base(_d), m_sd(nullptr) { }
};

template<class T> inline bool
operator== (const gsFeSpace<T> & a, const gsFeSpace<T> & b)
{ return a.id()== b.id() && a.isAcross()==b.isAcross(); }

}// namespace expr
}// namespace gismo