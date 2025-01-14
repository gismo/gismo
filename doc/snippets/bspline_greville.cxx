    gsKnotVector<> kv (-1, 0, 3, 3, 1 );
    gsBSplineBasis<> basis(kv);

    gsMatrix<> greville  = basis.anchors();
    gsInfo << greville <<  "\n\n";

    gsMatrix<index_t> actmat;
    basis.active_into( greville, actmat );
    gsInfo << actmat <<  "\n\n";

    gsMatrix<> evaluate  = basis.eval( greville );
    gsInfo << evaluate <<  "\n\n";

    
