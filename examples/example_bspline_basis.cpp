
# include <gismo.h>

using namespace gismo;


int main(int, char**)
{

    gsKnotVector<>   kv (-1, 0, 3,3, 1 );
    gsBSplineBasis<> bsp(kv);
    gsInfo << bsp.detail() << "\n";

    bsp.uniformRefine();
    gsInfo << bsp.detail() << "\n";


    return 0;
}
