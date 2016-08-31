

#include <gsMpi/gsMpiHelper.h>

namespace gismo
{

// The singleton function
gsMpiComm & gsMpiCommSingleton(const int& argc, char** argv)
{
    // create singleton instance
    static gsMpiComm singleton (argc, argv);
    return singleton;
}

};

