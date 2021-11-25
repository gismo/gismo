/** @file gsExprHelper.h

    @brief Generic expressions helper

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsExpressions.h>
#include <gsUtils/gsThreaded.h>

namespace gismo
{

/**
   Class holding an expression environment
 */
template<class T>
class gsExprHelper
{
private:
    gsExprHelper(const gsExprHelper &);

    gsExprHelper() : m_mirror(nullptr), mesh_ptr(nullptr)
    { mutVar.setData(mutData); }

private:
    typedef std::map<const gsFunctionSet<T>*,util::gsThreaded<gsFuncData<T> > >  FuncData;
    typedef std::map<const gsFunctionSet<T>*,util::gsThreaded<gsMapData<T> > >  MapData;
    typedef std::map<const expr::gsFeVariable<T>*,util::gsThreaded<gsFuncData<T> > >  VarData;

    typedef std::pair<const gsFunctionSet<T>*,util::gsThreaded<gsMapData<T> >*> CFuncKey;
    typedef std::map<CFuncKey,util::gsThreaded<gsFuncData<T> > >  CFuncData;

    typedef typename FuncData::iterator FuncDataIt;
    typedef typename MapData ::iterator MapDataIt;
    typedef typename VarData ::iterator VarDataIt;
    typedef typename CFuncData ::iterator CFuncDataIt;

    util::gsThreaded<gsMatrix<T> > m_points;
    FuncData  m_fdata;
    MapData   m_mdata;
    CFuncData m_cdata;
    VarData   m_vdata;// for BCs etc. //do we need it interfaced?

    memory::shared_ptr<gsExprHelper> m_mirror;
private:

    // mutable pair of variable and data,
    // ie. not uniquely assigned to a gsFunctionSet
    expr::gsFeVariable<T> mutVar ; //ADD: symbol.setSource(.);
    gsFuncData<T>         mutData;
    bool mutParametric;
    const gsMultiBasis<T> * mesh_ptr;

public:
    typedef memory::unique_ptr<gsExprHelper> uPtr;
    typedef memory::shared_ptr<gsExprHelper>  Ptr;

    typedef const expr::gsGeometryMap<T>   geometryMap;
    typedef const expr::gsFeElement<T>     element;
    typedef const expr::gsFeVariable<T>    variable;
    typedef const expr::gsFeSpace<T>       space;
    typedef const expr::gsComposition<T>   composition;
    typedef const expr::gsNullExpr<T>      nullExpr;

private:
    explicit gsExprHelper(gsExprHelper * m)
    : mesh_ptr(m->mesh_ptr) { m_mirror = memory::make_shared_not_owned(m); }

public:

    ~gsExprHelper() { }

    gsMatrix<T> & points() { return m_points; }
    inline gsExprHelper & iface()
    {
        if (nullptr==m_mirror )
            m_mirror = memory::make_shared(new gsExprHelper(this));
        return *m_mirror;
    }

    bool isMirrored() const { return nullptr!=m_mirror; }

    static uPtr make() { return uPtr(new gsExprHelper()); }

    void reset()
    {
        points().clear();
        //mapVar.reset();
    }

    void cleanUp()
    {
        mutData.clear();
    }

    void setMultiBasis(const gsMultiBasis<T> & mesh) { mesh_ptr = &mesh; }

    bool multiBasisSet() { return NULL!=mesh_ptr;}

    const gsMultiBasis<T> & multiBasis()
    {
        GISMO_ASSERT(multiBasisSet(), "Integration elements not set.");
        return *mesh_ptr;
    }

    const gsMultiPatch<T> & multiPatch() const
    {
        GISMO_ASSERT(!m_mdata.empty(), "Geometry map not set.");
        GISMO_ASSERT(nullptr!=dynamic_cast<const gsMultiPatch<T>*>(m_mdata.begin()->first), "Multipatch geometry map not set.");
        return *static_cast<const gsMultiPatch<T>*>(m_mdata.begin()->first);
    }

    const gsMapData<T> & multiPatchData() const
    {
        GISMO_ASSERT(!m_mdata.empty(), "Geometry map not set.");
        return m_mdata.begin()->second;
    }

    geometryMap getMap(const gsFunction<T> & mp)
    {
        expr::gsGeometryMap<T> gm;
        gm.setSource(mp);
        gsDebugVar(gm.m_fs);
        return gm;
    }

    geometryMap getMap(const gsMultiPatch<T> & mp)
    {
        expr::gsGeometryMap<T> gm;
        gm.setSource(mp);
        return gm;
    }

    expr::gsFeVariable<T> getVar(const gsFunctionSet<T> & mp, index_t dim = 1)
    {
        expr::gsFeVariable<T> var;
        var.setSource(mp);
        var.setDim(dim);
        return var;
    }

    composition getVar(const gsFunctionSet<T> & mp, geometryMap & G)
    {
        expr::gsComposition<T> var(G);
        var.setSource(mp);
        return var;
    }

    expr::gsFeSpace<T> getSpace(const gsFunctionSet<T> & mp, index_t dim = 1)
    {
        expr::gsFeSpace<T> var;
        var.setSource(mp);
        var.setDim(dim);
        return var;
    }

    variable getMutVar() const { return mutVar; }

    void setMutSource(const gsFunction<T> & func, bool param)
    {
        mutVar.setSource(func);
        mutParametric = param;
    }

private:
    template <class E1>
    void _parse(const expr::_expr<E1> & a1)
    {
        a1.parse(*this);
        //a1.print(gsInfo);
    }

    template <class E1, class... Rest>
    void _parse(const expr::_expr<E1> & a1, Rest &... restArgs)
    {
        _parse(a1);
        _parse(restArgs...);
    }

    template<size_t I, typename... Ts>
    void _parse_tuple_i (const std::tuple<Ts...> &tuple)
    {
        std::get<I>(tuple).parse(*this);
        if (I + 1 < sizeof... (Ts))
            _parse_tuple_i<(I+1 < sizeof... (Ts) ? I+1 : I)> (tuple);
    }

    template<typename... Ts>
    void _parse_tuple (const std::tuple<Ts...> &tuple) {_parse_tuple_i<0>(tuple);}
    //template<> void _parse_tuple<> (const std::tuple<> &) { }

public:

    template<class... Ts>
    void parse(const std::tuple<Ts...> &tuple)
    {
#pragma omp single
        {
            m_mdata.clear();
            m_fdata.clear();
            m_vdata.clear();
            m_cdata.clear();

            if (isMirrored())
            {
                m_mirror->m_mdata.clear();
                m_mirror->m_fdata.clear();
                m_mirror->m_vdata.clear();
                m_mirror->m_cdata.clear();
            }
        }//implicit barrier

        _parse_tuple(tuple);

        // Additional evaluation flags
        for (MapDataIt it  = m_mdata.begin(); it != m_mdata.end(); ++it)
            it->second.mine().flags |= SAME_ELEMENT|NEED_ACTIVE;
        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
            it->second.mine().flags |= SAME_ELEMENT|NEED_ACTIVE;
        for (VarDataIt it  = m_vdata.begin(); it != m_vdata.end(); ++it)
            it->second.mine().flags |= SAME_ELEMENT|NEED_ACTIVE;
        for (CFuncDataIt it  = m_cdata.begin(); it != m_cdata.end(); ++it)
            it->second.mine().flags |= SAME_ELEMENT|NEED_ACTIVE;

        // gsInfo<< "\nfdata: "<< m_fdata.size()<<"\n";
        // gsInfo<< "mdata: "<< m_mdata.size()<<"\n";
        // gsInfo<< "vdata: "<< m_vdata.size()<<"\n";
        // gsInfo<< "cdata: "<< m_cdata.size()<<std::endl;
    }

//    #ifndef _OPENMP
    template<class... expr>
    void parse(const expr &... args)
    {
#pragma omp single
        {
        m_mdata.clear();
        m_fdata.clear();
        m_vdata.clear();
        m_cdata.clear();
        }//implicit barrier

        #pragma omp single
        {
        m_mdata.clear();
        m_fdata.clear();
        m_vdata.clear();
        m_cdata.clear();
        }//implicit barrier

        _parse(args...);

        // Add initial evaluation flags
        for (MapDataIt it  = m_mdata.begin(); it != m_mdata.end(); ++it)
            it->second.mine().flags |= SAME_ELEMENT|NEED_ACTIVE;
        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
            it->second.mine().flags |= SAME_ELEMENT|NEED_ACTIVE;
        for (VarDataIt it  = m_vdata.begin(); it != m_vdata.end(); ++it)
            it->second.mine().flags |= SAME_ELEMENT|NEED_ACTIVE;
        for (CFuncDataIt it  = m_cdata.begin(); it != m_cdata.end(); ++it)
            it->second.mine().flags |= SAME_ELEMENT|NEED_ACTIVE;

        /*
        gsInfo<< "\nfdata: "<< m_fdata.size()<<"\n";
        gsInfo<< "mdata: "<< m_mdata.size()<<"\n";
        gsInfo<< "vdata: "<< m_vdata.size()<<"\n";
        gsInfo<< "cdata: "<< m_cdata.size()<<"\n";
        if (m_mirror)
        {
            gsInfo<< "\nfdata2: "<< iface().m_fdata.size()<<"\n";
            gsInfo<< "mdata2: "  << iface().m_mdata.size()<<"\n";
            gsInfo<< "vdata2: "  << iface().m_vdata.size()<<"\n";
            gsInfo<< "cdata2: "  << iface().m_cdata.size()<<"\n";
        }
        */
    }
//#endif

    void add(const expr::gsGeometryMap<T> & sym)
    {//TODO: inherit gsGeomatryMap from symbol_expr
        GISMO_ASSERT(NULL!=sym.m_fs, "Geometry map "<<&sym<<" is invalid");
        gsExprHelper & eh = (sym.isAcross() ? iface() : *this);
#       pragma omp critical (m_mdata_first_touch)
            const_cast<expr::gsGeometryMap<T>&>(sym)
                .setData(eh.m_mdata[sym.m_fs]);
    }

    void add(const expr::gsComposition<T> & sym)
    {
        GISMO_ASSERT(NULL!=sym.m_fs, "Composition "<<&sym<<" is invalid");
        add(sym.inner());
        sym.inner().data().flags |= NEED_VALUE;
        auto k = std::make_pair(sym.m_fs,&m_mdata[sym.inner().m_fs]);
        auto it = m_cdata.find(k);
        gsExprHelper & eh = (sym.isAcross() ? iface() : *this);
        if (m_cdata.end()==it)
            // when the variable is added for the first time,
            // we have to be thread-safe (atomic).
#           pragma omp critical (m_cdata_first_touch)
            const_cast<expr::gsComposition<T>&>(sym)
                .setData(eh.m_cdata[ give(k) ]);
        else
            const_cast<expr::gsComposition<T>&>(sym)
                .setData(eh.m_cdata[ give(k) ]);
    }

    template <class E>
    void add(const expr::symbol_expr<E> & sym)
    {
        //parallel: variables become thread-local
        // for each variable we provide a gsFuncData pointer
        // in the same thread this can be the same ptr (as done now)
        gsExprHelper & eh = (sym.isAcross() ? iface() : *this);

        if (NULL!=sym.m_fs)
        {
            if ( 1==sym.m_fs->size() &&
                 sym.m_fs->domainDim()<=sym.m_fs->targetDim() )// map?
            {
                //gsDebug<<"+ Map "<< sym.m_fs <<"\n";
#               pragma omp critical (m_mdata_first_touch)
                const_cast<expr::symbol_expr<E>&>(sym)
                    .setData( eh.m_mdata[sym.m_fs] );
            }
            else
            {
                //gsDebug<<"+ Func "<< sym.m_fs <<"\n";
#               pragma omp critical (m_fdata_first_touch)
                const_cast<expr::symbol_expr<E>&>(sym)
                    .setData( eh.m_fdata[sym.m_fs] );
            }
        }
        else
        {
            gsDebug<<"\ngsExprHelper: Boundary condition is NOT working yet for ["<< sym <<"] \n";
            // eg. mutable var?
        }
    }

    // template <class E>
    // void add(expr::avg_expr<E> & avg)
    // {
    //     // avg(u) * avg(g) // specialize * op ???
    //     // how about nv(.)

    //     // ??mode B11,B12,B21,B22
    //     avg.parse(*this); // pass a state: now you are left/right
    //     // OR
    //     avg.parse(this->iface());
    // }

    void precompute(const index_t patchIndex = 0,
                    boundary::side bs = boundary::none)
    {
        //First compute the maps
        for (MapDataIt it = m_mdata.begin(); it != m_mdata.end(); ++it)
        {
            it->second.mine().points.swap(m_points.mine());//swap
            it->second.mine().side    = bs;
            it->first->function(patchIndex).computeMap(it->second);
            it->second.mine().patchId = patchIndex;
            it->second.mine().points.swap(m_points.mine());
        }

        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
        {
            it->first->piece(patchIndex)
                .compute(m_points, it->second);
            it->second.mine().patchId = patchIndex;
        }

        for (VarDataIt it = m_vdata.begin(); it != m_vdata.end(); ++it)
        {
            it->first->source().piece(patchIndex)
                .compute(m_points, it->second);
            it->second.mine().patchId = patchIndex;
        }

        for (CFuncDataIt it = m_cdata.begin(); it != m_cdata.end(); ++it)
        {
            it->first.first->piece(patchIndex)
                .compute(it->first.second->mine().values[0], it->second);
            it->second.mine().patchId = patchIndex;
        }
    }

};//class


} //namespace gismo
