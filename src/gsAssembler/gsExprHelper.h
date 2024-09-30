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

    gsExprHelper() : m_mirror(nullptr), mesh_ptr(nullptr),
                     mutSrc(nullptr), mutMap(nullptr)
    { }

    explicit gsExprHelper(gsExprHelper * m)
    : m_mirror(memory::make_shared_not_owned(m)),
      mesh_ptr(m->mesh_ptr), mutSrc(nullptr), mutMap(nullptr)
    { }

private:
    typedef util::gsThreaded<gsFuncData<T> > thFuncData;
    typedef util::gsThreaded<gsMapData<T> >  thMapData;
    typedef std::map<const gsFunctionSet<T>*,thFuncData>  FuncData;
    typedef std::map<const gsFunctionSet<T>*,thMapData>  MapData;
    typedef std::pair<const gsFunctionSet<T>*,thMapData*> CFuncKey;
    typedef std::map<CFuncKey,thFuncData>  CFuncData;

    typedef typename FuncData::iterator FuncDataIt;
    typedef typename MapData ::iterator MapDataIt;
    typedef typename CFuncData ::iterator CFuncDataIt;

    util::gsThreaded<gsMatrix<T> > m_points;
    FuncData  m_fdata;///< functions
    MapData   m_mdata;///< maps
    CFuncData m_cdata;///< compositions

    memory::shared_ptr<gsExprHelper> m_mirror;

    const gsMultiBasis<T> * mesh_ptr;

    // mutable pair of variable and data,
    // ie. not uniquely assigned to a gsFunctionSet
    const gsFunctionSet<T> * mutSrc;
    const gsFunctionSet<T> * mutMap;
    thFuncData               mutData;

    // Represents the current element
    expr::gsFeElement<T> m_element;

public:
    typedef memory::unique_ptr<gsExprHelper> uPtr;
    typedef memory::shared_ptr<gsExprHelper>  Ptr;

    typedef const expr::gsGeometryMap<T>   geometryMap;
    typedef       expr::gsFeElement<T> &   element;
    typedef const expr::gsFeVariable<T>    variable;
    typedef const expr::gsFeSpace<T>       space;
    typedef const expr::gsComposition<T>   composition;
    typedef const expr::gsNullExpr<T>      nullExpr;

public:

    ~gsExprHelper() { }

    gsMatrix<T> & points()    { return m_points; }
    gsMatrix<T> & pointsIfc() { return this->iface().m_points; }

    bool isMirrored() const { return nullptr!=m_mirror; }

    static uPtr make() { return uPtr(new gsExprHelper()); }

    void reset()
    {
        points().clear();
        //mapVar.reset();
    }

    void cleanUp()
    {
        #pragma omp single
        {
            m_mdata.clear();
            m_fdata.clear();
            m_cdata.clear();
            //mutSrc = nullptr;
            mutMap = nullptr;
            mutData.mine().clear();
            if (isMirrored())
            {
                m_mirror->m_mdata.clear();
                m_mirror->m_fdata.clear();
                m_mirror->m_cdata.clear();
                //m_mirror->mutSrc = nullptr;
                m_mirror->mutMap = nullptr;
                m_mirror->mutData.mine().clear();
            }
        }//implicit barrier
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
        if ( !m_mdata.empty() )
        {
        GISMO_ASSERT(nullptr!=dynamic_cast<const gsMultiPatch<T>*>(m_mdata.begin()->first),
                     "Multipatch geometry map not set.");
            return *static_cast<const gsMultiPatch<T>*>(m_mdata.begin()->first);
        }
        if (isMirrored() && !m_mirror->m_mdata.empty() )
        {
            GISMO_ASSERT(nullptr!=dynamic_cast<const gsMultiPatch<T>*>(m_mirror->m_mdata.begin()->first),
                         "Multipatch geometry map not set.");
            return *static_cast<const gsMultiPatch<T>*>(m_mirror->m_mdata.begin()->first);
        }
        GISMO_ERROR("Geometry map not set.");
    }

    const gsMapData<T> & multiPatchData() const
    {
        GISMO_ASSERT(!m_mdata.empty(), "Geometry map not set.");
        return m_mdata.begin()->second;
    }

    geometryMap getMap(const gsFunctionSet<T> & mp)
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

    variable getMutVar() const
    {
        expr::gsFeVariable<T> var;
        return var;
    }

    element getElement() { return m_element; }

    composition getMutVar(geometryMap & G)
    {
        expr::gsComposition<T> var(G);
        //mutMap = &G.source();
        return var;
    }

    void setMutSource(const gsFunctionSet<T> & func)
    {
        mutSrc = &func;
    }

    //void clearMutSource() ?

    void activateFlags(unsigned flg)
    {
        // Additional evaluation flags
        for (MapDataIt it  = m_mdata.begin(); it != m_mdata.end(); ++it)
            it->second.mine().flags |= flg;
        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
            it->second.mine().flags |= flg;
        for (CFuncDataIt it  = m_cdata.begin(); it != m_cdata.end(); ++it)
            it->second.mine().flags |= flg;
        // gsInfo<< "\n-fdata: "<< m_fdata.size()<<"\n";
        // gsInfo<< "-mdata: "<< m_mdata.size()<<"\n";
        // gsInfo<< "-cdata: "<< m_cdata.size()<<std::endl;
    }

private:

    inline gsExprHelper & iface()
    {
        if (nullptr==m_mirror )
            m_mirror = memory::make_shared(new gsExprHelper(this));
        return *m_mirror;
    }

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

    void setInitialFlags()
    {
        // Additional evaluation flags
        for (MapDataIt it  = m_mdata.begin(); it != m_mdata.end(); ++it)
            it->second.mine().flags |= NEED_ACTIVE;
        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
            it->second.mine().flags |= NEED_ACTIVE;
        for (CFuncDataIt it  = m_cdata.begin(); it != m_cdata.end(); ++it)
        it->second.mine().flags |= NEED_ACTIVE;
        //gsInfo<< "\n-fdata: "<< m_fdata.size()<<"\n";
        //gsInfo<< "-mdata: "<< m_mdata.size()<<"\n";
        //gsInfo<< "-cdata: "<< m_cdata.size()<<std::endl;

        if (isMirrored())
        {
            for (MapDataIt it  = m_mirror->m_mdata.begin(); it != m_mirror->m_mdata.end(); ++it)
                it->second.mine().flags |= NEED_ACTIVE;
            for (FuncDataIt it = m_mirror->m_fdata.begin(); it != m_mirror->m_fdata.end(); ++it)
                it->second.mine().flags |= NEED_ACTIVE;
            for (CFuncDataIt it  = m_mirror->m_cdata.begin(); it != m_mirror->m_cdata.end(); ++it)
                it->second.mine().flags |= NEED_ACTIVE;
            // gsInfo<< "+fdata: "<< m_mirror->m_fdata.size()<<"\n";
            // gsInfo<< "+mdata: "<< m_mirror->m_mdata.size()<<"\n";
            // gsInfo<< "+cdata: "<< m_mirror->m_cdata.size()<<std::endl;
        }
    }

public:

    template<class... Ts>
    void parse(const std::tuple<Ts...> &tuple)
    {
        cleanUp(); //assumes parse is called once.
        _parse_tuple(tuple);
        setInitialFlags();
    }

    template<class... expr>
    void parse(const expr &... args)
    {
        cleanUp(); //assumes parse is called once.
        _parse(args...);
        setInitialFlags();
    }

    void add(const expr::gsGeometryMap<T> & sym)
    {
        GISMO_ASSERT(NULL!=sym.m_fs, "Geometry map "<<&sym<<" is invalid");
        gsExprHelper & eh = (sym.isAcross() ? iface() : *this);
#       pragma omp critical (m_mdata_first_touch)
            const_cast<expr::gsGeometryMap<T>&>(sym)
                .setData(eh.m_mdata[sym.m_fs]);
    }

    void add(const expr::gsComposition<T> & sym)
    {
        //GISMO_ASSERT(NULL!=sym.m_fs, "Composition "<<&sym<<" is invalid");
        add(sym.inner());//the map
        sym.inner().data().flags |= NEED_VALUE;
        if (nullptr==sym.m_fs)
        {
            //gsInfo<<"\nGot BC composition\n";
            mutMap = &sym.inner().source();
            if (nullptr!=mutSrc)
            {
#               pragma omp critical (m_fdata_first_touch)
                const_cast<expr::gsComposition<T>&>(sym)
                    .setData( mutData );

                const_cast<expr::gsComposition<T>&>(sym)
                    .setSource(*mutSrc);
            }
            else
                gsWarn<<"\nSomething went terribly wrong here (add gsComposition).\n";
            return;
        }

        //register the function //if !=nullptr?
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
            /*
            if ( 1==sym.m_fs->size() &&
                 sym.m_fs->domainDim()<=sym.m_fs->targetDim() )// map?
            {
                //gsDebug<<"+ Map "<< sym.m_fs <<"\n";
#               pragma omp critical (m_mdata_first_touch)
                const_cast<expr::symbol_expr<E>&>(sym)
                    .setData( eh.m_mdata[sym.m_fs] );
            }
            else
            */
            {
                //gsDebug<<"+ Func "<< sym.m_fs <<"\n";
#               pragma omp critical (m_fdata_first_touch)
                const_cast<expr::symbol_expr<E>&>(sym)
                    .setData( eh.m_fdata[sym.m_fs] );
            }
        }
        else
        {
            //gsDebug<<"\nGot a mutable variable.\n";
            if (nullptr!=mutSrc)
            {
#               pragma omp critical (m_fdata_first_touch)
                const_cast<expr::symbol_expr<E>&>(sym)
                    .setData( mutData );

                const_cast<expr::symbol_expr<E>&>(sym)
                    .setSource(*mutSrc);
            }
            else
                gsWarn<<"\nSomething went wrong here (add symbol_expr).\n";
        }
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

        for (CFuncDataIt it = m_cdata.begin(); it != m_cdata.end(); ++it)
        {
            it->first.first->piece(patchIndex)
                .compute(it->first.second->mine().values[0], it->second.mine());
            it->second.mine().patchId = patchIndex;
        }

        // Mutable variable to treat BCs
        if (nullptr!=mutSrc && 0!=mutData.mine().flags)
        {
            mutSrc->piece(patchIndex)
                .compute( mutMap ? m_mdata[mutMap].mine().values[0]
                          : m_points, mutData.mine() );
        }
    }

    void precompute(const boundaryInterface & iFace)
    {
        this->precompute( iFace.first ().patch, iFace.first().side() );
        if ( isMirrored() )
            m_mirror->precompute(iFace.second().patch, iFace.second().side());
    }

};//class gsExprHelper


} //namespace gismo
