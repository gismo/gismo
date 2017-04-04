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

/** 
 *  @brief Class for command-line argument parsing
 *
 *  A typical setup looks as follows:
 *
 *  \code{.cpp}
 *  int main( int argc, char** argv )
 *  {
 *      std::string name("G+SMO");
 *
 *      gsCmdLine cmd("This is my hello world example."); // Constructor
 *
 *      cmd.addString("n", "name", "Your name", name);
 *      // Here more of the addXxxx member functions might follow...
 *
 *      // AFTER all options have been registered with the addXxxx functions,
 *      // we parse the command line. If the return value is false, we exit.
 *      // At this point, the variable "name" is updated to the value given
 *      // by the user.
 *      if (!cmd.getValues(argc, argv)) return 1;
 *
 *      // Here, no more of the addXxxx function are allowed to follow.
 *
 *      // If everything was fine, we proceed:
 *      gsInfo << "Hello " << name << "!\n";
 *      return 0;
 *  }
 *  \endcode
 *
 *  \ingroup IO
 */
class GISMO_EXPORT gsCmdLine
{
public:
    typedef std::map<std::string,TCLAP::Arg*> ArgTable;

public:
    ///
    /// @brief Command line constructor. Defines how the arguments will be
    /// parsed.
    ///
    /// @param message        The message to be used in the usage output.
    /// @param delimiter      The character that is used to separate
    ///                       the argument flag/name from the value. 
    ///                       Defaults to ' ' (space).
    /// @param helpAndVersion Whether or not to create the Help and
    ///                       Version switches. Defaults to true.
    ///
    /// The options that are allowed in the command line have to be defined
    /// using the member functions \a addInt, \a add Real, \a addString,
    /// \a addSwitch and \a addPlainString. Finally, the parsing is invoked
    /// by calling \a getValues
    gsCmdLine(const std::string& message,
              const char delimiter = ' ',
              bool helpAndVersion = true);

    /// @brief Parses the command line based on the specified parameters
    ///
    /// The specification has to be done using \a addInt, \a add Real
    /// \a addString, \a addSwitch and \a addPlainString before calling
    /// this member function.
    ///
    /// The parameters \a argc and \a argv are those of the main function.
    ///
    /// The function returns true if the parsing succeeded. Oterwise,
    /// it prints an error message and returns false.
    bool getValues(int argc, char *argv[]);

    /// @brief Register an int option for the command line
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used
    /// @param name       Long form of the flag
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const variable, initialized
    ///                   with the default value. If \a getValues is
    ///                   invoked, the user-provided value is written
    ///                   to that variable.
    ///
    /// If the flag is "n", the user might call "-n 10" at the command line.
    /// It the name is "size", the user might call "--size 10" at the command line.
    void addInt(const std::string& flag, 
                const std::string& name, 
                const std::string& desc, 
                int              & value);

    /// @brief Register a real option for the command line
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used
    /// @param name       Long form of the flag
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const variable, initialized
    ///                   with the default value. If \a getValues is
    ///                   invoked, the user-provided value is written
    ///                   to that variable.
    ///
    /// If the flag is "t", the user might call "-t .5" at the command line.
    /// It the name is "tau", the user might call "--tau .5" at the command line.
     void addReal(const std::string& flag, 
                  const std::string& name, 
                  const std::string& desc, 
                  real_t           & value);

    /// @brief Register a string option for the command line
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used
    /// @param name       Long form of the flag
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const variable, initialized
    ///                   with the default value. If \a getValues is
    ///                   invoked, the user-provided value is written
    ///                   to that variable.
    ///
    /// If the flag is "f", the user might call "-f foo.xml" at the command line.
    /// If the name is "file", the user might call "--file foo.xml" at the command line.
    void addString(const std::string& flag, 
                   const std::string& name, 
                   const std::string& desc, 
                   std::string      & value);

    // TODO: do we need this variant?
    void addSwitch(const std::string& name, 
                   const std::string& desc, 
                   bool             & value);

    /// @brief Register a switch option for the command line
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used
    /// @param name       Long form of the flag
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const bool variable with
    ///                   value "false". If \a getValues is invoked and
    ///                   the user has called the swich, the variable is
    ///                   set to true
    ///
    /// If the flag is "l", the user might call "-l" at the command line.
    /// If the name is "log", the user might call "--log" at the command line.
    void addSwitch(const std::string& flag, 
                   const std::string& name, 
                   const std::string& desc, 
                   bool             & value);

    /// @brief Register a string parameter that has to be given directly (not
    /// as an option, i.e., not after a flag starting with "-" or "--")
    void addPlainString(const std::string& name, 
                        const std::string& desc, 
                        std::string      & value);


    /// Writes all given options (as specified by \a addInt
    /// \a addReal, \a addString or \a addSwitch) into a
    /// gsOptionList object. The values from \a addPlainString
    /// are not copied.
    ///
    /// Must be invoked after \a getValues
    gsOptionList getOptionList();

    ~gsCmdLine();

    operator TCLAP::CmdLineInterface &();

    static void printVersion();

    std::string& getMessage();

    // ----------------- TODO: do we really need these functions?
   
    gsCmdLine(int argc, char *argv[], const std::string& message);

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
