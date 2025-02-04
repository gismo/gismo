/** @file

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/

#include <gismo.h>
#include <gsJSON/gsJSON.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    gsJSON j;

    j["a"] = 2;
    j["b"] = 3.0;
    j["c"]["d"] = "four";
    j["e"] = {"e1", "e2", "e3"};

    gsMatrix<> mat(2,2);
    mat << 1, 2, 3, 4;
    j["A"] = mat;

    gsVector<index_t> vec(3);
    vec << 1, 2, 3;
    j["v"] = vec;

    gsKnotVector<> kv(0,1,0,2);
    j["kv"] = kv;

    gsBSplineBasis<> basis(kv);
    j["basis"] = basis;

    gsTensorBSplineBasis<2> tbasis(kv,kv);
    j["tbasis"] = tbasis;

    gsOptionList opt;
    opt.addInt("a", "", 2);
    opt.addReal("b", "", 3.0);
    gsJSON j2(opt);
    gsInfo<<j2<<"\n";
    gsOptionList opt2 = j2.get<gsOptionList>();
    gsInfo<<opt2;

    for (auto it = j.begin(); it != j.end(); ++it)
    {
        gsInfo<<it.key()<<": "<<it.value()<<"\n";
    }

    j.save("test.json");

    return EXIT_SUCCESS;

}//main
