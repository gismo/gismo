

#include <gsMpi/gsMpi.h>

namespace gismo
{


gsMpi & gsMpiSingleton(const int& argc, char** argv)
{
    // create singleton instance
    static gsMpi singleton (argc, argv);
    return singleton;
}

gsNoMpi & gsNoMpiSingleton(const int& argc, char** argv)
{
    // create singleton instance
    static gsNoMpi singleton (argc, argv);
    return singleton;
}

};
