

#include <gsMpi/gsMpiHelper.h>

namespace gismo
{

// The singleton function
gsMPIHelper & gsMpiHelperSingleton(const int& argc, char** argv)
{
    // create singleton instance
    static gsMPIHelper singleton (argc, argv);
    return singleton;
}

};

