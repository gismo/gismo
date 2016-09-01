

#include <gsMpi/gsMpiHelper.h>

namespace gismo
{


gsMpiComm & gsMpiCommSingleton(const int& argc, char** argv)
{
    // create singleton instance
    static gsMpiComm singleton (argc, argv);
    return singleton;
}

gsSerialComm & gsSerialCommSingleton(const int& argc, char** argv)
{
    // create singleton instance
    static gsSerialComm singleton (argc, argv);
    return singleton;
}

};
