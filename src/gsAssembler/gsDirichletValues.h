
namespace gismo {

namespace expr
{
template<class T> class gsFeSpace;
};

template<class T>
void gsDirichletValues(
    const gsBoundaryConditions<T> & bc,
    const index_t dir_values,
    const expr::gsFeSpace<T> & u)
{
    //if ( m_options.getInt("DirichletStrategy") == dirichlet::nitsche)
    //    return; // Nothing to compute

    //const gsMultiBasis<T> & mbasis = dynamic_cast<const gsMultiBasis<T>&>(u.source());

    // eg. not penalize
    const gsDofMapper & mapper = u.mapper();

    switch ( dir_values )
    {
    case dirichlet::homogeneous :
    case dirichlet::user :
        gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
        // If we have a homogeneous problem then fill with zeros
        fixedDofs.setZero( mapper.boundarySize() );
        break;
    case dirichlet::interpolation:
        gsDirichletValuesByInterpolation(u,bc);
        break;
    case dirichlet::l2Projection:

        break;
    default:
        GISMO_ERROR("Something went wrong with Dirichlet values: "<< dir_values);
    }

    /* Corner values -- todo
       for ( typename gsBoundaryConditions<T>::const_citerator
       it = bbc.cornerBegin();
       it != bbc.cornerEnd(); ++it )
       {
       if(it->unknown == unk)
       {
       const int i  = mbasis[it->patch].functionAtCorner(it->corner);
       const int ii = mapper.bindex( i , it->patch );
       u.fixedPart().row(ii).setConstant(it->value);
       }
       else
       continue;
       }
    */
}

template<class T>
void gsDirichletValuesByInterpolation(const expr::gsFeSpace<T> & u,
                                      const gsBoundaryConditions<T> & bc)
{
    const index_t parDim = u.source().domainDim();
    const gsMultiBasis<T> & mbasis =
        *dynamic_cast<const gsMultiBasis<T>*>(&u.source());

    const gsFunctionSet<T> & gmap = bc.geoMap();

    std::vector< gsVector<T> > rr;
    gsMatrix<unsigned> boundary;
    gsVector<T> b(1);
    gsMatrix<T> fpts, tmp;

    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();

    // Iterate over all patch-sides with Boundary conditions
    typedef gsBoundaryConditions<T> bcList;
    for ( typename bcList::const_iterator it =  bc.begin("Dirichlet");
          it != bc.end("Dirichlet") ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;

        const index_t com = it->unkComponent();
        //
        for (index_t r = 0; r!=u.dim(); ++r)
        {
            if (com!=-1 && r!=com) continue;
            fixedDofs.resize(u.mapper().boundarySize(), u.dim() );

            const int k = it->patch();
            const gsBasis<T> & basis = mbasis[k];

            // Get dofs on this boundary
            boundary = basis.boundary(it->side());

            // If the condition is homogeneous then fill with zeros
            if ( it->isHomogeneous() )
            {
                for (index_t i=0; i!= boundary.size(); ++i)
                {
                    const int ii = u.mapper().bindex( boundary.at(i) , k, com );
                    fixedDofs.at(ii) = 0;
                }
                continue;
            }

            // Get the side information
            int dir = it->side().direction( );
            index_t param = (it->side().parameter() ? 1 : 0);

            // Compute grid of points on the face ("face anchors")
            rr.clear();
            rr.reserve( parDim );

            for ( int i=0; i < parDim; ++i)
            {
                if ( i==dir )
                {
                    b[0] = ( basis.component(i).support() ) (0, param);
                    rr.push_back(b);
                }
                else
                {
                    rr.push_back( basis.component(i).anchors().transpose() );
                }
            }

            // GISMO_ASSERT(it->function()->targetDim() == u.dim(),
            //              "Given Dirichlet boundary function does not match problem dimension."
            //              <<it->function()->targetDim()<<" != "<<u.dim()<<"\n");

            // Compute dirichlet values
            if ( it->parametric() )
                fpts = it->function()->eval( gsPointGrid<T>( rr ) );
            else
                fpts = it->function()->eval(  gmap.piece(it->patch()).eval(  gsPointGrid<T>( rr ) )  );

            if ( fpts.rows() != u.dim() )
            {
                // assume scalar
                tmp.resize(u.dim(), fpts.cols());
                tmp.setZero();
                gsDebugVar(!dir);
                tmp.row(!dir) = (param ? 1 : -1) * fpts; // normal !
                fpts.swap(tmp);
            }

            // Interpolate dirichlet boundary
            typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());
            typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
            const gsMatrix<T> & dVals =  geo->coefs();

            // Save corresponding boundary dofs
            for (index_t l=0; l!= boundary.size(); ++l)
            {
                const int ii = u.mapper().bindex( boundary.at(l) , k, com );
                fixedDofs.at(ii) = dVals.at(l);
            }
        }
    }
}

}; // namespace gismo
