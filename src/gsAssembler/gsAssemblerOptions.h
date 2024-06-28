/** @file gsAssemblerOptions.h

    @brief Provides assembler and solver options.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, C. Hofer
*/

#pragma once

namespace gismo
{

struct dirichlet
{	
    enum strategy
    {
        elimination  = 11, ///< Enforce Dirichlet BCs by eliminating them from the system

        penalize     = 13, ///< Penalize the diagonal at the position of Dirichlet DoFs,

        nitsche      = 12, ///< Enforce the boundary condition weakly by a penalty term
        
        /// Compute Dirichlet DoFs in the normal direction (for a vector valued function),
        /// The tangential component are handled with the Nitsche method.
        eliminatNormal = 14,

        none         = 0 ///<< Do absolutely nothing for Dirichlet boundary conditions.
    };

    enum values
    {
        homogeneous   = 100, ///< Assume homogeneous Dirichlet conditions

        interpolation = 101, ///< Compute Dirichlet DoFs by using interpolation on the boundary
        
        l2Projection  = 102, ///< Compute Dirichlet DoFs by using L2 projection on the boundary
        
        user          = 103 ///< User will provide values of the Dirichlet dofs
    };
};

struct iFace
{	
    enum strategy
    {
        /// Glue patches together by merging DoFs across an
        /// interface into one. This only works for conforming
        /// interfaces.
        conforming = 1,
        glue       = 1,

        /// Use discontinuous Galerkin-like coupling between
        /// adjacent patches.
        dg = 2,

        /// Use enhanced smoothness splines between interfaces of adjacent patches.
        smooth = 3,
        
        /// Do absolutely nothing for coupling the interfaces.
        none = 0
    };

};

/*
    enum iFaceTopology
    {
        nested   = 1,

        clamped  = 2,
    }
*/

struct transform
{
    enum type
    {
        Hgrad = 1, // covariant, inverse_composition
        Hdiv  = 2, // Piola
        Hcurl = 3
    };
};

// for mixed formulations
struct discreteSpace
{	
    enum type
    {
        taylorHood    = 1,
        //instead of raviartThomas, there should be nested_Space and BubbleElement
        //raviartThomas should go away here.
        raviartThomas = 2,

        none          = 0
    };
};

/**
 * Info Stucture, which should keep track of the structure of the considered PDE and Discretization. This
 * class brings together information about the different unknown and eachs components, which basis
 * corresponds to which component and the transformation for the several unknowns.
 *
 * It is assumed that we order the components consecutive, e.g. we have a PDE with 2 unknown, the first has 1 component
 * and the second one 2. This means, we have in total 3 components, the first corresponds to unk 1 the second and the third
 * to unk 2.
 * If we want to have different solution and test spaces,
 * we have the double ammount of components, where the second half corresponds to the trial variables. These have the same
 * structure as the solution variables (same dimension of each variable).  Each component (including the trial space)
 * is assigned a basis, which might be equal. The purpose is to keep track, which basis belongs to which component, to prevent
 * evaluation the same basis twice, which is a very costly part of IgA.
 *
 * Furthermore, the class stores, which transformation is used for which unknown (all components of this unknown must have
 * the same transformation). Only in the case of HGrad conforming transformation, the transformation applies to all components
 * of the unknown in the same way. Hence, the transformation is connected to the whole unknown.
 *
class gsAssemblerSetup
{
public:

    //Standard case is Poisson Case
    gsAssemblerSetup()
    {
        m_comp2basis.setOnes(1);
        m_unkDims.setOnes(1);

        nUnk= 1;
        nComp = 1;

        m_transforms.push_back (transform::Hgrad);

        index_t c=0;
        m_unk2comp.resize(nUnk);
        m_comp2unk.resize(nComp);
        for(index_t unk=0; unk<nUnk;++unk)
        {
            m_unk2comp[unk].resize(m_unkDims[unk]);
            for(unsigned comp=0; comp<m_unkDims[unk];++comp,++c)
            {
                m_unk2comp[unk][comp]= c;
                m_comp2unk[c]=unk;
            }
        }
    }

     // initialize the info structure with your settings of the
     //  discretization and the PDE.  comp2basis describes a map from
     //  every component (sequentially numbered) to each basis see
     //  also SparseSystem. This information is contained in both
     //  clases.
     
    template<typename T>
    void init(const gsPde<T>& pde, const gsVector<index_t>& comp2basis)
    {
        m_comp2basis = comp2basis;
        m_unkDims= pde.unknownDim();

        nUnk= m_unkDims.size();
        nComp = m_unkDims.sum();

        usingDifferentTrialSpaces = false;
        if(nComp==2*m_comp2basis.size())
            usingDifferentTrialSpaces=true;
        else if(nComp!=m_comp2basis.size())
            GISMO_ERROR("You messed up something, number of components does not match the ones from the PDE");

        m_transforms.resize(nUnk);
        for(index_t unk=0; unk<nUnk;++unk)
            m_transforms[unk] = transform::Hgrad;

        index_t c=0;
        m_unk2comp.resize(nUnk);
        m_comp2unk.resize(nComp);
        for(index_t unk=0; unk<nUnk;++unk)
        {
            m_unk2comp[unk].resize(m_unkDims[unk]);
            for(unsigned comp=0; comp<m_unkDims[unk];++comp,++c)
            {
                m_unk2comp[unk][comp]= c;
                m_comp2unk[c]=unk;
            }
        }
    }


public:

//Accessor funtions, to be completed.

    inline bool trialSpaceNotIdentical() const {return usingDifferentTrialSpaces;}
    inline bool isTrialSpace(index_t compWithTrial)   const {return compWithTrial>= nComp ? true: false;}

    inline index_t getUnk(index_t compWithTrial) const {return m_comp2unk[compWithTrial%nComp];}

    inline const gsVector<index_t>& getAllBasis() const {return m_comp2basis;}
    inline index_t getBasis(index_t compWithTrial) const {return m_comp2basis[compWithTrial];}

    inline index_t getBasis(index_t comp, bool testSpace) const
    {
        return m_comp2basis[usingDifferentTrialSpaces && testSpace ? comp+nComp : comp];
    }

    inline index_t getBasis(index_t unk, index_t comp_unk, bool testSpace=false) const
    {
        return getBasis(m_unk2comp[unk][comp_unk], testSpace);
    }

    inline void setTransform(index_t unk, transform::type trans) { m_transforms[unk]=trans;}
    inline const transform::type& getTransform(index_t unk)const {return m_transforms[unk];}

    inline const gsVector<index_t> & getComp(index_t unk)const {return m_unk2comp[unk];}

    inline const gsVector<unsigned>& unknownDims() const { return m_unkDims; }
    inline unsigned  unknownDims(index_t unk) const { return m_unkDims[unk]; }

    inline index_t  nUnknown() const { return nUnk; }
    inline index_t  nComponents(bool withTrialSpace=false) const { return withTrialSpace ? 2*nComp : nComp ; }

private:

    gsVector<index_t> m_comp2basis;
    gsVector<unsigned> m_unkDims;
    std::vector<transform::type> m_transforms;

    std::vector<gsVector< index_t> > m_unk2comp;
    std::vector< index_t> m_comp2unk;

    index_t nUnk;
    index_t nComp;

    bool usingDifferentTrialSpaces;
};
*/



struct gsAssemblerOptions
{
public:
    // Default constructor
    gsAssemblerOptions()
        : dirValues    (dirichlet::l2Projection ),
          dirStrategy  (dirichlet::elimination  ),
          intStrategy  (iFace    ::conforming   ),
          transformType(transform::Hgrad        ),
          spaceType    (discreteSpace::taylorHood   ),

          bdA(2.0),
          bdB(1  ),
          memOverhead(0.33334),
          quA(1.0),
          quB(1  )
    { }

public:
    //gsAssemblerSetup      info;

    dirichlet::values    dirValues;

    dirichlet::strategy  dirStrategy;

    iFace::strategy      intStrategy;

    transform::type      transformType;
    discreteSpace::type   spaceType;

    // If set to a value different than zero, it controls the
    // allocation of the sparse matrix, ie. the maximum number of
    // non-zero entries per column (set to: A * p + B)
    double  bdA;
    int bdB;

    // more memory is allocated then required for efficency reasons,
    // more precise, (1+memOverhead) times the original memory is allocated
    // default value is 0.33334 -> 75% of the allocated memory is used.
    double memOverhead;

    // The formula for the number of quadrature points for all
    // integral computations will be set to the integer which is
    // closest to (A * p + B), where \a p is the (coordinate-wise)
    // degree of the basis
    double  quA;
    int quB;

public: // Utility functions that return values implied by the settings


    index_t numQuNodes(const gsBasis<real_t> & b) const
    {
        return numQuNodes(b,quA,quB);
    }

    static index_t numQuNodes(const gsBasis<real_t> & b,
                              double _quA, int _quB)
    {
        index_t res = 1;
        for(short_t i=0; i<b.domainDim(); ++i )
        {
            res *= static_cast<index_t>(_quA * b.degree(i) + _quB + 0.5);
        }
        
        return res;
    }


    index_t numColNz(const gsBasis<real_t> & b) const
    {
        return numColNz(b,bdA,bdB,memOverhead);
    }

    static index_t numColNz(const gsBasis<real_t> & b,
                            double _bdA, int _bdB, double _mem)
    {
        index_t nz = 1;
        for (short_t i = 0; i != b.dim(); ++i)
            nz *= static_cast<index_t>(_bdA * b.degree(i) + _bdB + 0.5);
        return static_cast<index_t>(nz*(1.0+_mem));
    }

};



}
