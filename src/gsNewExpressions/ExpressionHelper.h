/** @file ExpressionHelper.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsNewExpressions/ExpressionForwardDeclarations.h>

namespace gismo
{

template <class T>
class ExpressionHelper
{
    typedef T Scalar;
private:

    // NEW
    typedef util::gsThreaded<gsFuncData<T> > thFuncData;
    typedef util::gsThreaded<gsMapData<T> >  thMapData;
    typedef std::map<const gsFunctionSet<T>*,thFuncData>  FuncData;
    typedef std::map<const gsFunctionSet<T>*,thMapData>  MapData;
    // typedef std::pair<const gsFunctionSet<T>*,thMapData*> CFuncKey;
    // typedef std::map<CFuncKey,thFuncData>  CFuncData;

    typedef typename FuncData::iterator FuncDataIt;
    typedef typename MapData ::iterator MapDataIt;
    // typedef typename CFuncData ::iterator CFuncDataIt;

    util::gsThreaded<gsMatrix<T> > m_points;
    util::gsThreaded<gsVector<T> > m_weights;
    FuncData  m_fdata;///< functions
    MapData   m_mdata;///< maps
    // CFuncData m_cdata;///< compositions

    // memory::shared_ptr<gsExprHelper> m_mirror;

    typename gsDomain<T>::Ptr m_domain;

    // mutable pair of VariableObject and data,
    // ie. not uniquely assigned to a gsFunctionSet
    const gsFunctionSet<T> * mutSrc;
    const gsFunctionSet<T> * mutMap;
    thFuncData               mutData;

    // Represents the current element
    // expr::gsFeElement<T> m_element; //sharedby all threads

public:

    /// @brief @todo
    gsMatrix<T> & points()    { return m_points; }
    // /// @brief @todo
    // gsMatrix<T> & pointsIfc() { return this->iface().m_points; }


    /// @brief @todo
    gsVector<T> & weights()    { return m_weights; }
    /// @brief @todo
    const gsVector<T> & weights() const { return m_weights; }

    // /// @brief @todo
    // gsVector<T> & weightsIfc() { return this->iface().m_weights; }
    // /// @brief @todo
    // const gsVector<T> & weightsIfc() const { return this->iface().m_weights; }

    // /// @brief @todo
    // bool isMirrored() const { return nullptr!=m_mirror; }

    // /// @brief @todo
    // static uPtr make() { return uPtr(new gsExprHelper()); }

    /// @brief @todo
    void reset()
    {
        points().clear();
        //mapVar.reset();
    }

    // /// @brief @todo
    // void cleanUp()
    // {
    //     #pragma omp single
    //     {
    //         m_mdata.clear();
    //         m_fdata.clear();
    //         m_cdata.clear();
    //         //mutSrc = nullptr;
    //         mutMap = nullptr;
    //         mutData.mine().clear();
    //         if (isMirrored())
    //         {
    //             m_mirror->m_mdata.clear();
    //             m_mirror->m_fdata.clear();
    //             m_mirror->m_cdata.clear();
    //             //m_mirror->mutSrc = nullptr;
    //             m_mirror->mutMap = nullptr;
    //             m_mirror->mutData.mine().clear();
    //         }
    //     }//implicit barrier
    // }


public:

    Expr::ConstantObject<T,0>      getConstant(const T s, std::string label="c")
    {
        // Create a new constant scalar expression
        Expr::ConstantObject<T,0> expr(std::array<size_t, 0>{}, label);
        gsMatrix<T> val(1,1);
        val<<s;
        expr.setValue(val);
        return expr;
    }

    Expr::ConstantObject<T,1>      getConstant(const gsVector<T> & v, std::string label="C")
    {
        // Create a new constant vector expression
        Expr::ConstantObject<T,1> expr(std::array<size_t,1>{(size_t)v.rows()}, label);
        expr.setValue(v);
        return expr;
    }

    Expr::ConstantObject<T,2>      getConstant(const gsMatrix<T> & m, std::string label="C")
    {
        // Create a new constant matrix expression
        Expr::ConstantObject<T,2> expr(std::array<size_t,2>{(size_t)m.rows(),(size_t)m.cols()}, label);
        expr.setValue(m);
        return expr;
    }

    Expr::VariableObject<T,0,true> getScalarFunction(const gsConstantFunction<T> & cfunc, std::string label="f")
    {
        GISMO_ASSERT(cfunc.targetDim()==1,"Function is not scalar");
        // Create a new VariableObject scalar expression
        Expr::VariableObject<T,0,true> expr(cfunc.domainDim(),{}, label);
        expr.setSource(cfunc);
        return expr;
    }

    Expr::VariableObject<T,0,false> getScalarFunction(const gsFunction<T> & func, std::string label="f")
    {
        GISMO_ASSERT(func.targetDim()==1,"Function is not scalar");
        // Create a new VariableObject scalar expression
        Expr::VariableObject<T,0,false> expr(func.domainDim(),{}, label);
        expr.setSource(func);
        return expr;
    }

    Expr::VariableObject<T,1,true> getVectorFunction(const gsConstantFunction<T> & cfunc, std::string label="F")
    {
        GISMO_ASSERT(cfunc.targetDim()!=1,"Function is not scalar");
        // Create a new VariableObject scalar expression
        Expr::VariableObject<T,1,true> expr(cfunc.domainDim(),std::array<size_t,1>{(size_t)cfunc.targetDim()}, label);
        expr.setSource(cfunc);
        return expr;
    }

    Expr::VariableObject<T,1,false> getVectorFunction(const gsFunction<T> & func, std::string label="F")
    {
        GISMO_ASSERT(func.targetDim()!=1,"Function is not scalar");
        // Create a new VariableObject scalar expression
        Expr::VariableObject<T,1,false> expr(func.domainDim(),std::array<size_t,1>{(size_t)func.targetDim()}, label);
        expr.setSource(func);
        return expr;
    }

    Expr::SpaceObject<T,Expr::SpaceType::Test,0> getScalarTestSpace(const gsFunctionSet<T> & space, size_t id = 0, std::string label="φ")
    {
        // TODO: Assert if ID exists already
        GISMO_ASSERT(space.targetDim()==1,"Function is not scalar");
        // Create a new VariableObject scalar expression
        Expr::SpaceObject<T,Expr::SpaceType::Test,0> expr(space.domainDim(),std::array<size_t,0>{}, id, label);
        expr.setSource(space);
        return expr;
    }

    // TODO: replace 3rd template with Expr::SpaceType::Trial
    Expr::SpaceObject<T,Expr::SpaceType::Trial,0> getScalarTrialSpace(const gsFunctionSet<T> & space, size_t id = 0, std::string label="ψ")
    {
        // TODO: Assert if ID exists already
        GISMO_ASSERT(space.targetDim()==1,"Function is not scalar");
        // Create a new VariableObject scalar expression
        Expr::SpaceObject<T,Expr::SpaceType::Trial,0> expr(space.domainDim(),std::array<size_t,0>{}, id, label);
        expr.setSource(space);
        return expr;
    }

    Expr::SpaceObject<T,Expr::SpaceType::Test,1> getVectorTestSpace(const gsFunctionSet<T> & space, size_t dim, size_t id = 0, std::string label="φ")
    {
        GISMO_ASSERT(dim!=1,"Function is not a vector");
        // Create a new VariableObject scalar expression
        Expr::SpaceObject<T,Expr::SpaceType::Test,1> expr(space.domainDim(),std::array<size_t,1>{(size_t)dim}, id, label);
        expr.setSource(space);
        return expr;
    }

    Expr::SpaceObject<T,Expr::SpaceType::Trial,1> getVectorTrialSpace(const gsFunctionSet<T> & space, size_t dim, size_t id = 0, std::string label="ψ")
    {
        GISMO_ASSERT(dim!=1,"Function is not a vector");
        // Create a new VariableObject scalar expression
        Expr::SpaceObject<T,Expr::SpaceType::Trial,1> expr(space.domainDim(),std::array<size_t,1>{(size_t)dim}, id, label);
        expr.setSource(space);
        return expr;
    }

    template <size_t _Order>
    Expr::SolutionObject<T,Expr::SpaceType::Trial,_Order> getSolution(const Expr::SpaceObject<T,Expr::SpaceType::Trial,_Order> & space,
                                                      gsMatrix<T> & solVector, std::string label="u")
    {
        return Expr::SolutionObject<T,Expr::SpaceType::Trial,_Order>(space, solVector, label);
    }

    template <size_t _Order>
    Expr::SolutionObject<T,Expr::SpaceType::Trial,_Order> getSolution(const Expr::SpaceObject<T,Expr::SpaceType::Trial,_Order> & space,
                                                      gsVector<T> & solVector, std::string label="u")
    {
        return Expr::SolutionObject<T,Expr::SpaceType::Trial,_Order>(space, solVector, label);
    }

    template <size_t _Order>
    Expr::SolutionObject<T,Expr::SpaceType::Test,_Order> getSolution(const Expr::SpaceObject<T,Expr::SpaceType::Test,_Order> & space,
                                                      gsMatrix<T> & solVector, std::string label)
    {
        GISMO_ERROR("Solution can only be constructed from Trial space");
    }

    template <size_t Order, bool _IsConstant>
    void add(const Expr::VariableObject<T,Order,_IsConstant> & VariableObject)
    {
        const_cast<Expr::VariableObject<T,Order,_IsConstant>&>(VariableObject)
            .setData(this->m_fdata[&VariableObject.source()]);
    }

    template <size_t _Space, size_t order>
    void add(const Expr::SpaceObject<T,_Space,order> & space)
    {
        const_cast<Expr::SpaceObject<T,_Space,order>&>(space)
            .setData(this->m_fdata[&space.source()]);
    }

    void precompute(const index_t patchIndex = 0,
                    boundary::side bs = boundary::none)
    {
        //First compute the maps
        for (MapDataIt it = m_mdata.begin(); it != m_mdata.end(); ++it)
        {
            it->second.mine().points.swap(m_points.mine());//swap
            it->second.mine().side    = bs;
            it->second.mine().patchId = patchIndex;
            it->first->function(patchIndex).computeMap(it->second.mine());
            it->second.mine().points.swap(m_points.mine());
        }

        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
        {
            it->second.mine().patchId = patchIndex;
            it->first->piece(patchIndex)
                .compute(m_points, it->second.mine());
        }

        // for (CFuncDataIt it = m_cdata.begin(); it != m_cdata.end(); ++it)
        // {
        //     it->first.first->piece(patchIndex)
        //         .compute(it->first.second->mine().values[0], it->second.mine());
        //     it->second.mine().patchId = patchIndex;
        // }

        // // Mutable VariableObject to treat BCs
        // if (nullptr!=mutSrc && 0!=mutData.mine().flags)
        // {
        //     mutSrc->piece(patchIndex)
        //         .compute( mutMap ? m_mdata[mutMap].mine().values[0]
        //                   : m_points, mutData.mine() );
        // }
    }

};
}//namespace gismo