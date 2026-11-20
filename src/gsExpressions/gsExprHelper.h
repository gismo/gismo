/** @file gsExprHelper.h

    @brief Generic expressions helper

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsExpressions/gsExpressions.h>
#include <gsExpressions/gsFeElement.h>
#include <gsUtils/gsThreaded.h>
#include <gsAssembler/gsQuadrature.h>

namespace gismo
{

/**
 * @brief Helper class for expressions, holding an expression environment
 */
template<class T>
class gsExprHelper
{
private:
    gsExprHelper(const gsExprHelper &);

    gsExprHelper() : m_mirror(nullptr), m_domain(nullptr),
                     mutSrc(nullptr), mutMap(nullptr), m_element(*this)
    { }

    explicit gsExprHelper(gsExprHelper * m)
    : m_mirror(memory::make_shared_not_owned(m)),
      m_domain(m->m_domain), mutSrc(nullptr), mutMap(nullptr), m_element(*this)
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

    /// Key of a stacked two-sided entry: (source, symbolSide::mode stored as short_t).
    typedef std::pair<const gsFunctionSet<T>*,short_t> TwoSidedKey;
    typedef std::map<TwoSidedKey,thFuncData>           TwoSidedData;
    typedef typename TwoSidedData::iterator            TwoSidedDataIt;

    util::gsThreaded<gsMatrix<T> > m_points;
    util::gsThreaded<gsVector<T> > m_weights;
    FuncData     m_fdata;///< functions
    MapData      m_mdata;///< maps
    CFuncData    m_cdata;///< compositions
    TwoSidedData m_tsdata;///< stacked two-sided (jump/avg) function data

    memory::shared_ptr<gsExprHelper> m_mirror;

    typename gsDomain<T>::Ptr m_domain;

    // mutable pair of variable and data,
    // ie. not uniquely assigned to a gsFunctionSet
    const gsFunctionSet<T> * mutSrc;
    const gsFunctionSet<T> * mutMap;
    thFuncData               mutData;

    // Represents the current element
    expr::gsFeElement<T> m_element; //sharedby all threads

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

    /// @brief @todo
    gsMatrix<T> & points()    { return m_points; }
    /// @brief @todo
    gsMatrix<T> & pointsIfc() { return this->iface().m_points; }

    /// @brief @todo
    gsVector<T> & weights()    { return m_weights; }
    /// @brief @todo
    const gsVector<T> & weights() const { return m_weights; }

    /// @brief @todo
    gsVector<T> & weightsIfc() { return this->iface().m_weights; }
    /// @brief @todo
    const gsVector<T> & weightsIfc() const { return this->iface().m_weights; }

    /// @brief @todo
    bool isMirrored() const { return nullptr!=m_mirror; }

    /// @brief @todo
    static uPtr make() { return uPtr(new gsExprHelper()); }

    /// @brief @todo
    void reset()
    {
        points().clear();
        //mapVar.reset();
    }

    /// @brief @todo
    void cleanUp()
    {
        #pragma omp single
        {
            m_mdata.clear();
            m_fdata.clear();
            m_cdata.clear();
            m_tsdata.clear();
            //mutSrc = nullptr;
            mutMap = nullptr;
            mutData.mine().clear();
            if (isMirrored())
            {
                m_mirror->m_mdata.clear();
                m_mirror->m_fdata.clear();
                m_mirror->m_cdata.clear();
                m_mirror->m_tsdata.clear();
                //m_mirror->mutSrc = nullptr;
                m_mirror->mutMap = nullptr;
                m_mirror->mutData.mine().clear();
            }
        }//implicit barrier
    }

    /// @brief @todo
    void setDomain(typename gsDomain<T>::Ptr domain) { m_domain = give(domain); }

    /// @brief @todo
    bool domainSet() { return NULL!=m_domain;}

    /// @brief @todo
    const gsDomain<T> & domain() const { return *m_domain; }

    /// @brief @todo
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

    /// @brief @todo
    const gsMapData<T> & multiPatchData() const
    {
        GISMO_ASSERT(!m_mdata.empty(), "Geometry map not set.");
        return m_mdata.begin()->second;
    }

    /// @brief @todo
    geometryMap getMap(const gsFunctionSet<T> & mp)
    {
        expr::gsGeometryMap<T> gm;
        gm.setSource(mp);
        return gm;
    }

    /// @brief @todo
    expr::gsFeVariable<T> getVar(const gsFunctionSet<T> & mp, index_t dim = 1)
    {
        expr::gsFeVariable<T> var;
        var.setSource(mp);
        var.setDim(dim);
        return var;
    }

    /// @brief @todo
    composition getVar(const gsFunctionSet<T> & mp, geometryMap & G)
    {
        expr::gsComposition<T> var(G);
        var.setSource(mp);
        return var;
    }

    /// @brief @todo
    expr::gsFeSpace<T> getSpace(const gsFunctionSet<T> & mp, index_t dim = 1)
    {
        expr::gsFeSpace<T> var;
        var.setSource(mp);
        var.setDim(dim);
        return var;
    }

    /// @brief @todo
    variable getMutVar() const
    {
        expr::gsFeVariable<T> var;
        return var;
    }

    /// @brief @todo
    element getElement() { return m_element; }

    /// @brief @todo
    composition getMutVar(geometryMap & G)
    {
        expr::gsComposition<T> var(G);
        //mutMap = &G.source();
        return var;
    }

    /// @brief @todo
    void setMutSource(const gsFunctionSet<T> & func)
    {
        mutSrc = &func;
    }

    //void clearMutSource() ?

    /// @brief @todo
    void activateFlags(unsigned flg)
    {
        // Additional evaluation flags
        for (MapDataIt it  = m_mdata.begin(); it != m_mdata.end(); ++it)
            it->second.mine().flags |= flg;
        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
            it->second.mine().flags |= flg;
        for (CFuncDataIt it  = m_cdata.begin(); it != m_cdata.end(); ++it)
            it->second.mine().flags |= flg;
        for (TwoSidedDataIt it = m_tsdata.begin(); it != m_tsdata.end(); ++it)
            it->second.mine().flags |= flg;
        // gsInfo<< "\n-fdata: "<< m_fdata.size()<<"\n";
        // gsInfo<< "-mdata: "<< m_mdata.size()<<"\n";
        // gsInfo<< "-cdata: "<< m_cdata.size()<<std::endl;
    }

    /// @brief Initializes the mirror helper, if it is not already initialized. This is needed for interface evaluation.
    void initializeIface() { iface(); }

private:

    inline gsExprHelper & iface()
    {
        gsExprHelper * mirror;
#       pragma omp critical (m_mirror_init)
        {
            if (nullptr==m_mirror )
                m_mirror = memory::make_shared(new gsExprHelper(this));
            mirror = m_mirror.get();
        }
        return *mirror;
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

    template<size_t I, typename... Ts>
    void _parse_tt_tuple_i (const std::tuple<Ts...> &tuple)
    {
        if (std::get<I>(tuple).isMatrix())
        {
            std::get<I>(tuple).rowVar().parse(*this);
            std::get<I>(tuple).colVar().parse(*this);
        }
        if (I + 1 < sizeof... (Ts))
            _parse_tt_tuple_i<(I+1 < sizeof... (Ts) ? I+1 : I)> (tuple);
    }

    template<typename... Ts>
    void _parse_tt_tuple (const std::tuple<Ts...> &tuple) {_parse_tt_tuple_i<0>(tuple);}

    void setInitialFlags()
    {
        // Additional evaluation flags
        for (MapDataIt it  = m_mdata.begin(); it != m_mdata.end(); ++it)
            it->second.mine().flags |= NEED_ACTIVE;
        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
            it->second.mine().flags |= NEED_ACTIVE;
        for (CFuncDataIt it  = m_cdata.begin(); it != m_cdata.end(); ++it)
        it->second.mine().flags |= NEED_ACTIVE;
        for (TwoSidedDataIt it = m_tsdata.begin(); it != m_tsdata.end(); ++it)
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

    template<class... Ts>
    void parsePattern(const std::tuple<Ts...> &tuple)
    {
        cleanUp(); //assumes parse is called once.
        _parse_tt_tuple(tuple);
        for (FuncDataIt it = m_fdata.begin(); it != m_fdata.end(); ++it)
            it->second.mine().flags = NEED_ACTIVE;
        if (isMirrored())
            for (FuncDataIt it = m_mirror->m_fdata.begin(); it != m_mirror->m_fdata.end(); ++it)
                it->second.mine().flags = NEED_ACTIVE;
        for (TwoSidedDataIt it = m_tsdata.begin(); it != m_tsdata.end(); ++it)
            it->second.mine().flags = NEED_ACTIVE;
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
        GISMO_ENSURE(!sym.isTwoSided(),
                     "jump()/avg() on a composition is not supported.");
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
        if (sym.isTwoSided())
        {
            // Row-stacking assumes actives.rows()>1 and symbol_expr::rows()
            // == m_fs->targetDim(): true for FE spaces (and, through their
            // space, gsFeSolution), but a function-backed variable's
            // gsFunctionSet::compute gives actives.rows()==1 with values[0]
            // rows indexing target components instead -- stacking that would
            // silently produce a mis-shaped matrix, not a jump.
            GISMO_ENSURE(0 != E::Space,
                         "jump()/avg() are supported on FE spaces (and "
                         "gsFeSolution) only; a function-backed variable's "
                         "gsFuncData rows index target components, not "
                         "actives, so row-stacking cannot express its jump.");
            GISMO_ENSURE(NULL != sym.m_fs, "Symbol is invalid");
#           pragma omp critical (m_fdata_first_touch)
            {
                this->m_fdata[sym.m_fs];          // left  source entry (this helper)
                this->iface().m_fdata[sym.m_fs];  // right source entry (mirror helper)
                const_cast<expr::symbol_expr<E>&>(sym)
                    .setData( m_tsdata[TwoSidedKey(sym.m_fs,(short_t)sym.sideMode())] );
            }
            return;
        }

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

    /// True while a jump()/avg() symbol is registered (after parse()).
    bool hasTwoSided() const { return !m_tsdata.empty(); }

    /// One-sided evaluation. Not valid while a jump()/avg() symbol is
    /// registered: those need both sides of a face.
    void precompute(const index_t patchIndex = 0,
                    boundary::side bs = boundary::none)
    {
        GISMO_ENSURE(m_tsdata.empty(), "A jump()/avg() symbol is registered: "
                     "two-sided symbols are only valid in a skeleton/ghost face loop "
                     "(assembleSkeleton/assembleGhost, integralSkeleton/integralGhost).");
        _precompute(patchIndex, bs);
    }

    void precompute(const boundaryInterface & iFace)
    {
        if (!m_tsdata.empty())
        {
            GISMO_ENSURE(isMirrored(), "A jump()/avg() symbol is registered "
                         "but this helper has no mirror to evaluate the other "
                         "side of the face.");
            // gsFuncData carries a single patchId and _eval::push maps the
            // stacked actives with it: both sides must be the same patch.
            GISMO_ENSURE(iFace.first().patch == iFace.second().patch,
                         "jump()/avg() require both sides of the face on the "
                         "same patch (multipatch DG via jump/avg is not "
                         "implemented).");
            // Requests made through the stacked entry (NEED_*, derivOrder and
            // SAME_ELEMENT) must reach the two entries that are actually
            // computed -- in particular SAME_ELEMENT, without which the mirror
            // would produce one actives column per point while this side
            // produces one in total, and the two could not be stacked.
            // The union is symmetrised (ld/rd/sd all merged into both ld and
            // rd): a one-sided symbol sharing the same source (e.g.
            // dnk(u.left(),G,2) alongside dnk(u.jump(),G,1)) may have written
            // extra flags/derivOrder into only one of ld/rd, and both sides
            // of a stacked entry must be evaluated to the SAME order for
            // _stackFuncData's per-order layout to line up (its
            // L.values.size()==R.values.size() check).
            for (TwoSidedDataIt it = m_tsdata.begin(); it != m_tsdata.end(); ++it)
            {
                const gsFunctionSet<T> * fs = it->first.first;
                gsFuncData<T> & sd = it->second.mine();
                gsFuncData<T> & ld = m_fdata[fs].mine();
                gsFuncData<T> & rd = m_mirror->m_fdata[fs].mine();
                const unsigned f = ld.flags | rd.flags | sd.flags;
                ld.flags = rd.flags = f;
                const index_t dOrd = math::max(math::max(ld.derivOrder, rd.derivOrder), sd.derivOrder);
                ld.derivOrder = rd.derivOrder = dOrd;
            }
        }

        this->_precompute( iFace.first ().patch, iFace.first().side() );
        if ( isMirrored() )
            m_mirror->_precompute(iFace.second().patch, iFace.second().side());

        for (TwoSidedDataIt it = m_tsdata.begin(); it != m_tsdata.end(); ++it)
            _stackFuncData( m_fdata[it->first.first].mine(),
                            m_mirror->m_fdata[it->first.first].mine(),
                            (expr::symbolSide::mode)it->first.second,
                            it->second.mine() );
    }

private:

    void _precompute(const index_t patchIndex, boundary::side bs)
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

    /// Fills the two-sided (jump/avg) entry \a S by stacking the already
    /// computed one-sided entries \a L (left) and \a R (right) row-wise:
    /// actives are concatenated (no sign, no scale) so that both sides
    /// scatter through the ordinary one-sided assembly machinery
    /// (gsExprAssembler::_eval::push, _patternFace); each order's values
    /// block is concatenated per-active with sign/scale so that an
    /// expression evaluating row i*bsz+off (i indexing actives) picks up the
    /// left contribution for i<nL and the (signed, scaled) right one for
    /// i>=nL -- exactly the layout dnk_expr and symbol_expr::eval read.
    /// \a mode is expr::symbolSide::jump or expr::symbolSide::avg.
    static void _stackFuncData(const gsFuncData<T> & L, const gsFuncData<T> & R,
                                expr::symbolSide::mode mode, gsFuncData<T> & S)
    {
        const index_t nL = L.actives.rows();
        const index_t nR = R.actives.rows();
        GISMO_ENSURE(0 != nL && 0 != nR,
                     "jump()/avg(): actives were not computed on one of the "
                     "two sides of the face.");
        GISMO_ENSURE(L.actives.cols() == R.actives.cols(),
                     "jump()/avg(): the two sides of the face do not agree "
                     "on the number of quadrature points (SAME_ELEMENT "
                     "was not propagated to both sides).");
        GISMO_ENSURE(L.values.size() == R.values.size(),
                     "jump()/avg(): the two sides of the face do not agree "
                     "on the number of computed derivative orders.");
        GISMO_ENSURE(L.patchId == R.patchId,
                     "jump()/avg() require both sides of the face on the "
                     "same patch.");

        const T sgnR  = (expr::symbolSide::jump == mode) ? T(-1) : T(1);
        const T scale = (expr::symbolSide::avg  == mode) ? T(0.5) : T(1);

        S.flags      = L.flags;
        S.derivOrder = L.derivOrder;
        S.patchId    = L.patchId;
        S.dim        = L.dim;

        S.actives.resize(nL+nR, L.actives.cols());
        S.actives.topRows(nL)    = L.actives;
        S.actives.bottomRows(nR) = R.actives;

        S.values.resize(L.values.size());
        for (size_t n = 0; n != L.values.size(); ++n)
        {
            const gsMatrix<T> & Ln = L.values[n];
            const gsMatrix<T> & Rn = R.values[n];
            if (0 == Ln.rows() && 0 == Rn.rows())
                continue; // order n not evaluated on either side: leave empty
            GISMO_ENSURE(Ln.cols() == Rn.cols(),
                         "jump()/avg(): the two sides of the face do not "
                         "agree on the number of quadrature points for "
                         "derivative order "<<n<<".");
            GISMO_ENSURE(Ln.rows()*nR == Rn.rows()*nL,
                         "jump()/avg(): the two sides of the face do not "
                         "agree on the per-active block size for derivative "
                         "order "<<n<<".");
            const index_t bsL = Ln.rows()/nL;
            const index_t bsR = Rn.rows()/nR;
            gsMatrix<T> & Sn = S.values[n];
            Sn.resize(nL*bsL + nR*bsR, Ln.cols());
            Sn.topRows(nL*bsL)    = scale*Ln;
            Sn.bottomRows(nR*bsR) = scale*sgnR*Rn;
        }

        const bool curlsOk = (0!=L.curls.rows() && 0!=R.curls.rows());
        if (curlsOk)
        {
            S.curls.resize(L.curls.rows()+R.curls.rows(), L.curls.cols());
            S.curls.topRows(L.curls.rows())    = scale*L.curls;
            S.curls.bottomRows(R.curls.rows()) = scale*sgnR*R.curls;
        }
        else
            S.curls.resize(0,0);

        const bool divsOk = (0!=L.divs.rows() && 0!=R.divs.rows());
        if (divsOk)
        {
            S.divs.resize(L.divs.rows()+R.divs.rows(), L.divs.cols());
            S.divs.topRows(L.divs.rows())    = scale*L.divs;
            S.divs.bottomRows(R.divs.rows()) = scale*sgnR*R.divs;
        }
        else
            S.divs.resize(0,0);

        const bool laplOk = (0!=L.laplacians.rows() && 0!=R.laplacians.rows());
        if (laplOk)
        {
            S.laplacians.resize(L.laplacians.rows()+R.laplacians.rows(), L.laplacians.cols());
            S.laplacians.topRows(L.laplacians.rows())    = scale*L.laplacians;
            S.laplacians.bottomRows(R.laplacians.rows()) = scale*sgnR*R.laplacians;
        }
        else
            S.laplacians.resize(0,0);
    }

};//class gsExprHelper

namespace internal
{

/// \brief Returns the "faceShift" option of a skeleton/ghost face loop
/// (gsExprAssembler::assembleSkeleton, gsExprEvaluator::integralSkeleton).
///
/// The face rule nodes sit on the knot line; the two sides of the face are
/// evaluated at points moved by faceShift*h into the left/right neighbouring
/// element. Values and lower derivatives are then O(faceShift*h*|grad u|)
/// accurate, while the top derivative of a degree-k spline is exact (it is
/// elementwise constant). The shift must lie in (0,1): 0 makes both sides
/// coincide (every jump vanishes), a negative shift swaps the two sides, and
/// a shift >= 1 leaves the neighbouring element (outside the parameter domain
/// at a boundary element). shiftFacePoints() additionally rejects, per face,
/// shifts lost to rounding and shifts so close to 1 that a point reaches the
/// opposite knot line up to rounding (a conservative check, see there).
template<class T>
T faceShift(const gsOptionList & opt)
{
    const T shift = (T)opt.askReal("faceShift", 1e-6);
    GISMO_ENSURE(shift > 0 && shift < 1, "faceShift must lie in (0,1): with shift 0 both sides of a face land on the knot line and every jump term silently vanishes, a negative shift swaps the two sides, and a shift >= 1 leaves the neighbouring element.");
    return shift;
}

/// \brief For every fixed direction dir of a patch: the plain Gauss face rule
/// (one node across the face, quA*deg+quB along it, degrees read from \a sizing)
/// and the same-patch interface (boxSide(dir,true) | boxSide(dir,false)) that
/// precompute(boundaryInterface) evaluates both sides on.
///
/// Face integrals always use a plain Gauss rule and ignore any custom
/// quadrature factory.
template<class T>
void makeFaceRules(const gsDomain<T> & sizing, const gsOptionList & opt,
                   const index_t patch,
                   std::vector<typename gsQuadRule<T>::uPtr> & rules,
                   std::vector<boundaryInterface> & faceIfc)
{
    const short_t d = sizing.dim();
    rules.resize(d);
    faceIfc.resize(d);
    for (short_t dir = 0; dir != d; ++dir)
    {
        rules[dir] = gsGaussRule<T>::make(
            gsQuadrature::numNodes(sizing, opt.getReal("quA"), opt.getInt("quB"), dir));
        faceIfc[dir] = boundaryInterface(patchSide(patch, boxSide(dir,true)),
                                         patchSide(patch, boxSide(dir,false)), d);
    }
}

/// \brief Copies points() into pointsIfc(), then moves points() into the left
/// element by shift*h_left and pointsIfc() into the right element by
/// shift*h_right, along the face direction.
///
/// This order is the invariant precompute(boundaryInterface) relies on.
///
/// All face nodes share the face-direction coordinate, so column 0 suffices
/// for two O(1) checks. First, the left points must lie strictly below and
/// the right points strictly above the knot line x0; otherwise shift*h is
/// below the floating-point resolution of x0 and every jump would vanish
/// silently. Second, the left points must lie strictly above fl(x0 - h_left)
/// and the right points strictly below fl(x0 + h_right); otherwise a shift
/// close to 1 may have rounded a point onto or past the opposite knot line,
/// into the wrong element. The second check is sufficient, not necessary:
/// fl(x0 +- h) can be one ulp off the true knot, and a point can round
/// exactly onto the patch end or the lower knot of the left element, so it
/// also rejects some shifts whose points are evaluated in the correct
/// element. This happens only for 1 - shift below about spacing(x)/(2h), with
/// x the opposite knot and spacing(x) the gap from x to the next smaller
/// double (e.g. 500*2^-53 for h = 1e-3 next to x = 1).
template<class T>
void shiftFacePoints(gsExprHelper<T> & eh,
                     const typename gsDomain<T>::iterator & face, const T shift)
{
    const short_t dir = face.side().direction();
    const T x0 = eh.points().cols() > 0 ? eh.points()(dir,0) : (T)0;
    const T hL = face.getPerpendicularCellSize();
    const T hR = face.getPerpendicularCellSizeRight();
    eh.pointsIfc() = eh.points();
    eh.points()   .row(dir).array() -= shift * hL;
    eh.pointsIfc().row(dir).array() += shift * hR;
    GISMO_ENSURE(eh.points().cols() == 0 ||
                 (eh.points()(dir,0) < x0 && eh.pointsIfc()(dir,0) > x0),
                 "faceShift*h is below the floating-point resolution of the face coordinate: both sides of the face would stay on the knot line. Increase faceShift.");
    GISMO_ENSURE(eh.points().cols() == 0 ||
                 (eh.points()(dir,0) > x0 - hL && eh.pointsIfc()(dir,0) < x0 + hR),
                 "faceShift is too close to 1: a shifted face point reaches the opposite knot line of its neighbouring element up to rounding (the check is conservative and may reject a shift whose points still lie in the correct element). Decrease faceShift.");
}

} // namespace internal


} //namespace gismo
