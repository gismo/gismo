    gsKnotVector<> kv (-1, 0, 3,3, 1 );
    gsBSplineBasis<> bsp(kv);

    gsMatrix<> greville  = bsp.anchors();
    gsInfo << greville <<  "\n";

    gsMatrix<> evaluate  = bsp.eval( greville );
    gsInfo << evaluate <<  "\n";
