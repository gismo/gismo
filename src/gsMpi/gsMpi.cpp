

#include <gsMpi/gsMpi.h>

namespace gismo
{


gsMpi & gsMpiSingleton(const int& argc, char** argv)
{
    // create singleton instance
    static gsMpi singleton (argc, argv);
    return singleton;
}

};
