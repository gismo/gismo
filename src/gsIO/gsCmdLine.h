/** @file gsCmdLine.h

    @brief Provides input command line arguments.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

namespace TCLAP
{
//forward declarations
class CmdLineInterface;
class Arg;
}


namespace gismo
{

class gsCmdLinePrivate;

/* 
    @brief Type used for exceptions while reading command-line arguments
*/
//typedef TCLAP::ArgException gsArgException;

/** 
    @brief Class for command-line argument parsing

    The members add(..) and parse(..) are the main functions.
    
    For arguments, getValue() is the main member.
    
    \ingroup IO
*/
class GISMO_EXPORT gsCmdLine
{
public:
    typedef std::map<std::string,TCLAP::Arg*> ArgTable;

public:
    /**
     * Command line constructor. Defines how the arguments will be
     * parsed.
     * \param message - The message to be used in the usage
     * output.
     * \param delimiter - The character that is used to separate
     * the argument flag/name from the value.  Defaults to ' ' (space).
     * \param helpAndVersion - Whether or not to create the Help and
     * Version switches. Defaults to true.
     */
    gsCmdLine(const std::string& message,	const char delimiter = ' ',
              bool helpAndVersion = true);

    gsCmdLine(int argc, char *argv[], const std::string& message);
    
    ~gsCmdLine();

    operator TCLAP::CmdLineInterface &();

    static void printVersion();
    
public:

    bool getValues(int argc, char *argv[]);

    std::string & getMessage();

    void addInt(const std::string& flag, 
                const std::string& name, 
                const std::string& desc, 
                int              & value);

    void addReal(const std::string& flag, 
                 const std::string& name, 
                 const std::string& desc, 
                 real_t & value);

    void addString(const std::string& flag, 
                   const std::string& name, 
                   const std::string& desc, 
                   std::string & value);

    void addSwitch(const std::string& name, 
                   const std::string& desc, 
                   bool & value);

    void addSwitch(const std::string& flag, 
                   const std::string& name, 
                   const std::string& desc, 
                   bool & value);

    void addPlainString(const std::string& name, 
                        const std::string& desc, 
                        std::string & value);


    // -----------------
    
    int getInt(const std::string& flag, 
               const std::string& name, 
               const std::string& desc, 
               const int        & value);
    
    real_t getReal(const std::string& flag, 
                   const std::string& name, 
                   const std::string& desc, 
                   const real_t & value);

    std::string getString(const std::string& flag, 
                          const std::string& name, 
                          const std::string& desc, 
                          const std::string & value);

    bool getSwitch(const std::string& flag, 
                   const std::string& name, 
                   const std::string& desc, 
                   const bool & value);

    bool getSwitch(const std::string& name, 
                   const std::string& desc, 
                   const bool & value);

    std::string getPlainString(const std::string& name, 
                               const std::string& desc, 
                               const std::string & value);

    bool valid() const;

private:

    gsCmdLinePrivate * my;
    
}; // class gsCmdLine


}; // namespace gismo
