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

    gsExprHelper() : mesh_ptr(NULL)
    { mutVar.setData(mutData); }

private:
    typedef std::map<const gsFunctionSet<T>*,gsFuncData<T> >  FuncData;
    typedef std::map<const gsFunctionSet<T>*,gsMapData<T> >  MapData;
    typedef std::map<const expr::gsFeVariable<T>*,gsFuncData<T> >  VarData;

    typedef std::pair<const gsFunctionSet<T>*,gsMapData<T>*> CFuncKey;
    typedef std::map<CFuncKey,gsMapData<T> >  CFuncData;

    typedef typename FuncData::iterator FuncDataIt;
    typedef typename MapData ::iterator MapDataIt;
    typedef typename VarData ::iterator VarDataIt;
    typedef typename CFuncData ::iterator CFuncDataIt;

    CFuncData m_cdata;
    FuncData m_fdata, m_fdata2;
    MapData m_mdata, m_mdata2;
    VarData m_vdata;// for BCs etc. //do we need it interfaced?
    //gsFunctionSet* : interface, thread

    gsMapData<T> mapData;
private:

    // mutable pair of variable and data,
    // ie. not uniquely assigned to a gsFunctionSet
    expr::gsFeVariable<T> mutVar ; //ADD: symbol.setSource(.);
    gsFuncData<T>         mutData;
    bool mutParametric;

    const gsMultiBasis<T> * mesh_ptr;

public:
    typedef const expr::gsGeometryMap<T>   geometryMap;
    typedef const expr::gsFeElement<T>   & element;
    typedef const expr::gsFeVariable<T>    variable;
    typedef const expr::gsFeSpace<T>     & space; // reference needed for init.
    typedef const expr::gsNullExpr<T>      nullExpr;

    typedef expr::gsFeVariable<T>    nonConstVariable;
    typedef expr::gsFeSpace<T>       nonConstSpace;

    typedef memory::unique_ptr<gsExprHelper> uPtr;
    typedef memory::shared_ptr<gsExprHelper>  Ptr;
public:

    gsMatrix<T> & points() { return mapData.points; }

    static uPtr make() { return uPtr(new gsExprHelper()); }

    void reset()
    {
        points().clear();
        //mapVar.reset();
    }

    void print() const
    {

    }

    void cleanUp()
    {
        mapData.clear();
        mutData.clear();
    }

    void clean(bool space=false)
    {
    }

    void setMultiBasis(const gsMultiBasis<T> & mesh) { mesh_ptr = &mesh; }

    bool multiBasisSet() { return NULL!=mesh_ptr;}

    const gsMultiBasis<T> & multiBasis()
    {
        GISMO_ASSERT(multiBasisSet(), "Integration elements not set.");
        return *mesh_ptr;
    }

    // gsFuncData<T> * getData()
    // {
    //     return
    // }

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

    //geometryMap getMap() const

    //geometryMap getMap2() const

    expr::gsFeVariable<T> getVar(const gsFunctionSet<T> & mp, index_t dim = 1)
    {
        expr::gsFeVariable<T> var;
        var.setSource(mp);
        var.setDim(dim);
        //var.setData( m_fdata[&mp] );//NO
        return var;
    }

    expr::gsComposition<T> getVar(const gsFunctionSet<T> & mp, geometryMap G)
    {
        //GISMO_UNUSED(G);
        //GISMO_ASSERT(&G==&mapVar, "geometry map not known");
        //m_vlist.push_back( expr::gsFeVariable<T>() );
        expr::gsComposition<T> var(G);
        //gsFuncData<T> & fd = m_itable[&mp];
        //fd.dim = mp.dimensions();
        //gsDebugVar(&fd);
        //var.registerData(mp, fd, 1, mapData);
        var.setSource(mp);

        // const_cast<expr::gsGeometryMap<T>&>(sym)
        //     .m_fd = & m_mdata[sym.m_fs];

        //var.setMap( m_mdata[&G.source()] );
        //gsInfo<<"aaaaaaaaaaaaaaaaa Reg.  symb "<< &G.data() <<"\n";
        return var;
    }

    /*
    nonConstVariable searchVar(const gsFunctionSet<T> & mp)
    {
    // Search in the variable list and get back the first object which source \a mp
        GISMO_UNUSED(G);
        GISMO_ASSERT(&G==&mapVar, "geometry map not known");
        m_vlist.push_back( expr::gsFeVariable<T>() );
        expr::gsFeVariable<T> & var = m_vlist.back();
        gsFuncData<T> & fd = m_itable[&mp];
        //fd.dim = mp.dimensions();
        //gsDebugVar(&fd);
        var.registerData(mp, fd, 1, mapData);
        return var;
    }
    */

    std::list<expr::gsFeSpace<T> > m_slist;//to be moved in ExprAssembler
    expr::gsFeSpace<T> & getSpace(const gsFunctionSet<T> & mp, index_t dim = 1)
    {
        m_slist.push_back( expr::gsFeSpace<T>() );
        expr::gsFeSpace<T> & var = m_slist.back();
        var.setSource(mp);
        var.setData( m_fdata[&mp] );
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
        a1.print(gsInfo);
    }

    // void _parse(const expr::gsFeVariable<T> & a1)
    // {
    //     a1.parse(*this);
    //     gsDebug<< "got gsFeVariable "<< &a1 <<"\n";
    // }

    template <class E1, class... Rest>
    void _parse(const expr::_expr<E1> & a1, Rest &... restArgs)
    {
        _parse(a1);
        _parse(restArgs...);
    }
public:

    template<class... expr>
    void parse(const expr &... args)
    {
        m_mdata.clear();
        m_fdata.clear();
        m_vdata.clear();
        m_cdata.clear();

        _parse(args...);

        // Add initial evaluation flags
        for (MapDataIt it  = m_mdata.begin(); it != m_mdata.end(); ++it)
            it->second.flags |= SAME_ELEMENT|NEED_ACTIVE;
        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
            it->second.flags |= SAME_ELEMENT|NEED_ACTIVE;
        for (VarDataIt it  = m_vdata.begin(); it != m_vdata.end(); ++it)
            it->second.flags |= SAME_ELEMENT|NEED_ACTIVE;
        for (CFuncDataIt it  = m_cdata.begin(); it != m_cdata.end(); ++it)
            it->second.flags |= SAME_ELEMENT|NEED_ACTIVE;

        // gsInfo<< "\nfdata: "<< m_fdata.size()<<"\n";
        // gsInfo<< "mdata: "<< m_mdata.size()<<"\n";
        // gsInfo<< "vdata: "<< m_vdata.size()<<"\n";
        // gsInfo<< "cdata: "<< m_cdata.size()<<"\n";
    }

    void add(const expr::gsGeometryMap<T> & sym)
    {//TODO: inherit gsGeomatryMap from symbol_expr
        GISMO_ASSERT(NULL!=sym.m_fs, "Geometry map "<<&sym<<" is invalid");
        const_cast<expr::gsGeometryMap<T>&>(sym)
            .setData(m_mdata[sym.m_fs]);
        //gsDebug<<"+ gMap "<< sym.m_fs <<" ("<<sym.m_fd<<")\n";
    }

    void add(const expr::gsComposition<T> & sym)
    {
        GISMO_ASSERT(NULL!=sym.m_fs, "Composition "<<&sym<<" is invalid");
        add(sym.inner());
        sym.inner().data().flags |= NEED_VALUE;
        const_cast<expr::gsComposition<T>&>(sym)
            .setData(m_cdata[ std::make_pair(sym.m_fs,
            const_cast<gsMapData<T>*>(sym.inner().m_fd))]);
        //gsDebug<<"+ Comp "<< sym.m_fs <<" ("<<sym.m_fd<<")\n";
    }

    template <class E>
    void add(const expr::symbol_expr<E> & sym)
    {
        //parallel: variables become thread-local
        // for each variable we provide a gsFuncData pointer
        // in the same thread this can be the same ptr (as done now)

        //grab &sym.source();

        if (NULL!=sym.m_fs)
        {
            if ( 1==sym.m_fs->size() &&
                 sym.m_fs->domainDim()<=sym.m_fs->targetDim() )// map?
            {
                //gsDebug<<"+ Map "<< sym.m_fs <<"\n";
                const_cast<expr::symbol_expr<E>&>(sym)
                    .setData( m_mdata[sym.m_fs] );
            }
            else
            {
                //gsDebug<<"+ Func "<< sym.m_fs <<"\n";
                const_cast<expr::symbol_expr<E>&>(sym)
                    .setData( m_fdata[sym.m_fs] );
            }
        }
        else
        {
            gsDebug<<"- No source for "<< sym.m_fs <<"\n";
            // eg. mutable var?
        }

        // get variable (all): create and set m_fs and m_d [var is thread-local]
        //parse (all): set m_fd and flags [stored in m_fdata, th-local]
    }

    void precompute(const index_t patchIndex = 0)
    {
        //TO DO! : mapData.points  --> gsMatrix points.

        //First compute the maps
        for (MapDataIt it = m_mdata.begin(); it != m_mdata.end(); ++it)
        {
            it->second.points.swap(mapData.points);//swap
            it->first->function(patchIndex).computeMap(it->second);
            it->second.patchId = patchIndex;
            it->second.points.swap(mapData.points);
        }

        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
        {
            it->first->piece(patchIndex)
                .compute(mapData.points, it->second);
            it->second.patchId = patchIndex;
        }

        for (VarDataIt it = m_vdata.begin(); it != m_vdata.end(); ++it)
        {
            it->first->source().piece(patchIndex)
                .compute(mapData.points, it->second);
            it->second.patchId = patchIndex;
        }

        for (CFuncDataIt it = m_cdata.begin(); it != m_cdata.end(); ++it)
        {
            it->first.first->piece(patchIndex)
                .compute(it->first.second->values[0], it->second);
            it->second.patchId = patchIndex;
        }
    }

/*
    void precompute(const index_t patch1, const index_t patch2);
*/


};//class


} //namespace gismo
