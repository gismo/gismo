

#include <gsMpi/gsMpi.h>

namespace gismo
{


gsMpi & gsMpiSingleton(const int& argc, char** argv)
{
    // create singleton instance
    static gsMpi singleton (argc, argv);
    return singleton;
}

#if !defined(NDEBUG) && defined(GISMO_WITH_MPI)
MPI_Errhandler gsMpiComm::ErrHandler = MPI_ERRORS_ARE_FATAL;
#endif

};
