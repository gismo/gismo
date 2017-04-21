
    // Task: multiply the rows of M by the corresponding scalars in v
    gsVector<> v(4);
    gsMatrix<> M(4,6);
    v.setRandom();
    M.setRandom();

    // Manually
    gsMatrix<> R1(4,6);
    for (index_t i = 0; i<4; ++i) // for all rows
        for (index_t j = 0; j<6; ++j) // for all columns
            R1(i,j) = M(i,j) * v[i];

    // Using row-blocks
    gsMatrix<> R2(4,6);
    for (index_t i = 0; i<4; ++i) // for all rows
        R2.row(i).noalias() = M.row(i) * v[i]; // matrix x scalar

    // Using only one line
    gsMatrix<> R3;
    R3.noalias() = v.asDiagonal() * M; // matrix x matrix

    gsInfo << "v is: \n" << v << "\n\n";
    gsInfo << "M is: \n" << M << "\n\n";

    gsInfo << "R1 is: \n" << R1 << "\n";
    gsInfo << "R2 is: \n" << R2 << "\n";
    gsInfo << "R3 is: \n" << R3 << "\n";
