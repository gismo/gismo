    gsKnotVector<>   kv (-1, 0, 3, 3, 1 );
    gsBSplineBasis<> basis(kv);
    gsInfo << basis.detail() << "\n";

    gsMatrix<> pt(1,1);
    pt << 0;

    gsMatrix<index_t> act;
    gsMatrix<> val;
    basis.active_into(pt, act);
    basis.eval_into(pt, val);
    gsInfo << "Active basis functions on " << pt << ":\n" << act << "\n";
    gsInfo << "And their values:\n" << val << "\n";