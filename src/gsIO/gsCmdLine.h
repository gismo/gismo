/** @file gsCmdLine.h

    @brief Provides input command line arguments.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, S. Takacs
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

// Forward declaration
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
 *      try { cmd.getValues(argc, argv); } catch(int rv) { return rv; }
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
    typedef index_t intVal_t;

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
    /// using the member functions \a addInt, \a addMultiInt, \a addReal,
    /// \a addMultiReal, \a addString, \a addMultiString, \a addSwitch and
    /// \a addPlainString. Finally, the parsing is invoked
    /// by calling \a getValues.
    gsCmdLine(const std::string& message,
              const char delimiter = ' ',
              bool helpAndVersion = true);

    /// @brief Register an int option for the command line
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used.
    /// @param name       Long form of the flag.
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const variable, initialized
    ///                   with the default value. When \a getValues is
    ///                   invoked and the user provided a value, it is
    ///                   written to that variable.
    ///
    /// If the flag is "n", the user might call "-n 10" at the command line.
    /// It the name is "size", the user might call "--size 10" at the command line.
    void addInt(const std::string& flag,
                const std::string& name,
                const std::string& desc,
                intVal_t         & value);

    /// @brief Register an int option for the command line, which can be assigned more than once
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used.
    /// @param name       Long form of the flag.
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const vector. When \a getValues is
    ///                   invoked, the vector is filled with that values. Pre-existing
    ///                   values are kept only if the option has been used zero times.
    ///
    /// If the flag is "n", the user might call "-n 10" at the command line.
    /// It the name is "size", the user might call "--size 10" at the command line.
    void addMultiInt(const std::string    & flag,
                     const std::string    & name,
                     const std::string    & desc,
                     std::vector<intVal_t>& value);

    /// @brief Register a real option for the command line
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used
    /// @param name       Long form of the flag
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const variable, initialized
    ///                   with the default value. When \a getValues is
    ///                   invoked and the user provided a value, it is
    ///                   written to that variable.
    ///
    /// If the flag is "t", the user might call "-t .5" at the command line.
    /// It the name is "tau", the user might call "--tau .5" at the command line.
     void addReal(const std::string& flag,
                  const std::string& name,
                  const std::string& desc,
                  real_t           & value);

    /// @brief Register a real option for the command line, which can be assigned more than once
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used.
    /// @param name       Long form of the flag.
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const vector. When \a getValues is
    ///                   invoked, the vector is filled with that values. Pre-existing
    ///                   values are kept only if the option has been used zero times.
    ///
    /// If the flag is "t", the user might call "-t .5" at the command line.
    /// It the name is "tau", the user might call "--tau .5" at the command line.
    void addMultiReal(const std::string  & flag,
                      const std::string  & name,
                      const std::string  & desc,
                      std::vector<real_t>& value);

    /// @brief Register a string option for the command line
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used.
    /// @param name       Long form of the flag.
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const variable, initialized
    ///                   with the default value. When \a getValues is
    ///                   invoked and the user provided a value, it is
    ///                   written to that variable.
    ///
    /// If the flag is "f", the user might call "-f foo.xml" at the command line.
    /// If the name is "file", the user might call "--file foo.xml" at the command line.
    void addString(const std::string& flag,
                   const std::string& name,
                   const std::string& desc,
                   std::string      & value);

    /// @brief Register a string option for the command line, which can be assigned more than once
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used.
    /// @param name       Long form of the flag.
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const vector. When \a getValues is
    ///                   invoked, the vector is filled with that values. Pre-existing
    ///                   values are kept only if the option has been used zero times.
    ///
    /// If the flag is "f", the user might call "-f foo.xml" at the command line.
    /// If the name is "file", the user might call "--file foo.xml" at the command line.
    void addMultiString(const std::string      & flag,
                       const std::string       & name,
                       const std::string       & desc,
                       std::vector<std::string>& value);

    /// @brief Register a switch option for the command line
    ///
    /// @param flag       One character flag for using the option.
    ///                   If empty, no such flag can be used.
    /// @param name       Long form of the flag.
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const bool variable.
    ///                   When \ref getValues is invoked and the user
    ///                   has added the switch on the command line,
    ///                   the \a value is toggled (i.e. if false it
    ///                   becomes true, if true it becomes false)
    ///
    /// If the flag is "l", the user might call "-l" at the command line.
    /// If the name is "log", the user might call "--log" at the command line.
    void addSwitch(const std::string& flag,
                   const std::string& name,
                   const std::string& desc,
                   bool             & value);

    /// @brief Register a switch option for the command line without flag
    ///
    /// \see gsCmdLine::addSwitch
    ///
    void addSwitch(const std::string& name,
                   const std::string& desc,
                   bool             & value) { addSwitch("",name,desc,value); }

    /// @brief Register a string parameter that has to be given directly (not
    /// as an option, i.e., not after a flag starting with "-" or "--")
    ///
    /// @param name       Name of the option.
    /// @param desc       Description (printed if --help is invoked)
    /// @param value      This should be a non-const variable. When
    ///                   \a getValues is invoked, the user-provided value
    ///                   is written to that variable.
    ///
    /// You must not declare more than one plain string.
    void addPlainString(const std::string& name,
                        const std::string& desc,
                        std::string      & value);

    /// @brief Parses the command line based on the specified parameters
    ///
    /// The specification has to be done using \a addInt, \a addMultiInt,
    /// \a addReal, \a addMultiReal, \a addString, \a addMultiString,
    /// \a addSwitch and \a addPlainString before calling
    /// this member function.
    ///
    /// The parameters \a argc and \a argv are those of the main function.
    ///
    /// If the parsing did non succeed, the function throws.
    void getValues(int argc, char *argv[]);

    /// Writes all given options (as specified by \a addInt, \a addReal,
    /// \a addString or \a addSwitch or \a addPlainString) into a
    /// gsOptionList object.
    ///
    /// Must be invoked after \a getValues. This function takes its values
    /// from the registered variables, so changes in thoes are taken into
    /// account.
    gsOptionList getOptionList();

    /// Prints the version information
    static void printVersion();

    /// Returns the program's description (as specified in the constructor)
    std::string& getMessage();

    /// Parses the command line and returns true iff everything is okey.
    /// This function should be called *after* the parameters have been registered.
    ///
    /// If the user has invoked --help or --version, the result is true.
    bool valid(int argc, char *argv[]) const;

    /// Sets exception handling (true/false)
    void setExceptionHandling(const bool state);

    /// Gets the exception handling status (true/false)
    bool getExceptionHandling() const;

    // Destructor
    ~gsCmdLine();

private:

    gsCmdLinePrivate * my;

}; // class gsCmdLine

} // namespace gismo
