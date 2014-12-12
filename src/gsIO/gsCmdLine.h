
#pragma once

#include <gsCore/gsExport.h>
#include <gsCore/gsDebug.h>

#include <tclap/CmdLine.h>   // External file

#include <gsIO/gsCmdLineArgs.h>

namespace gismo
{

  /** 
      Class for command-line argument parsing
      add(..) and parse(..) are the main functions

      For arguments, getValue() is the main member
  */

typedef TCLAP::ArgException gsArgException; // For catching exceptions
  
class GISMO_EXPORT gsCmdLine :  public TCLAP::CmdLine
{
public:
    typedef TCLAP::CmdLine Base;
public:
    /**
     * Command line constructor. Defines how the arguments will be
     * parsed.
     * \param message - The message to be used in the usage
     * output.
     * \param delimiter - The character that is used to separate
     * the argument flag/name from the value.  Defaults to ' ' (space).
     * \param version - The version number to be used in the
     * --version switch.
     * \param helpAndVersion - Whether or not to create the Help and
     * Version switches. Defaults to true.
     */
    gsCmdLine(const std::string& message,	const char delimiter = ' ',
              const std::string& version = "none",
              bool helpAndVersion = true);
    
    ~gsCmdLine();

 private:
    class GismoCmdOut : public TCLAP::StdOutput
    {
    public:
	void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e);
	void usage(TCLAP::CmdLineInterface& c);
	void version(TCLAP::CmdLineInterface& c);
    };

    GismoCmdOut cmdout;

}; // class gsCmdLine


}; // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCmdLine.cpp)
#endif
