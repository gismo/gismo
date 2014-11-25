

#include <gsIO/gsCmdLineArgs.h>


namespace gismo
{



gsArgSwitch::gsArgSwitch(const std::string& flag, 
                         const std::string& name, 
                         const std::string& desc,
                         bool def ,
        TCLAP::Visitor* v) : Base(flag,name,desc,def,v) {}

gsArgSwitch::gsArgSwitch(const std::string& flag, 
                const std::string& name, 
                const std::string& desc,
                TCLAP::CmdLineInterface& parser,
                bool def,
        TCLAP::Visitor* v ) : Base(flag,name,desc,parser,def,v) {};





} //namespace gismo

