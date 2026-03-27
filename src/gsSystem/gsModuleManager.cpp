
#include <gsSystem/gsModuleManager.h>

namespace gismo
{

gsModuleManager & gsModuleManagerSingleton()
{
    // create singleton instance
    static gsModuleManager singleton;
    return singleton;
}

};
