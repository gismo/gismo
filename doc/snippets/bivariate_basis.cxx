    gsKnotVector<>   ku (-1, 1, 3, 4); // fist, last, interior, multEnd
    gsKnotVector<>   kv (-1, 1, 1, 4); // fist, last, interior, multEnd
    
    gsTensorBSplineBasis<2> basis( ku, kv );
    
    gsMesh<> mesh(basis);

    gsWriteParaview( basis, "basis", 1000); // object, name, number of plotting samples
    gsWriteParaview( mesh, "mesh");