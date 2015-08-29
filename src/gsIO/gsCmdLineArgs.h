/** @file gsCmdLineArgs.h
    
    @brief Provides input command line arguments.
    
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

namespace TCLAP
{
//forward declarations
class CmdLineInterface;
class Visitor;
template<class C> class MultiArg;
}

namespace gismo
{

template<class C> class gsArgMultiValPrivate;

template<class C> 
class gsArgMultiVal
{

public:

    /**
     * Constructor.
     * \param flag - The one character flag that identifies this
     * argument on the command line.
     * \param name - A one word name for the argument.  Can be
     * used as a long flag on the command line.
     * \param desc - A description of what the argument is for or
     * does.
     * \param req - Whether the argument is required on the command
     * line.
     * \param typeDesc - A short, human readable description of the
     * type that this object expects.  This is used in the generation
     * of the USAGE statement.  The goal is to be helpful to the end user
     * of the program.
     * \param parser - A CmdLine parser object to add this Arg to
     * \param v - An optional visitor.  You probably should not
     * use this unless you have a very good reason.
     */
    gsArgMultiVal( const std::string& flag, 
                   const std::string& name,
                   const std::string& desc,
                   bool req,
                   const std::string& typeDesc,
                   TCLAP::CmdLineInterface& parser,
                   TCLAP::Visitor* v = NULL );

    /*
     * Constructor.
     * \param flag - The one character flag that identifies this
     * argument on the command line.
     * \param name - A one word name for the argument.  Can be
     * used as a long flag on the command line.
     * \param desc - A description of what the argument is for or
     * does.
     * \param req - Whether the argument is required on the command
     * line.
     * \param typeDesc - A short, human readable description of the
     * type that this object expects.  This is used in the generation
     * of the USAGE statement.  The goal is to be helpful to the end user
     * of the program.
     * \param v - An optional visitor.  You probably should not
     * use this unless you have a very good reason.

    gsArgMultiVal( const std::string& flag,
                   const std::string& name,
                   const std::string& desc,
                   bool req,
                   const std::string& typeDesc,
                   TCLAP::Visitor* v = NULL); 
    */
        
    /*
     * Constructor.
     * \param flag - The one character flag that identifies this
     * argument on the command line.
     * \param name - A one word name for the argument.  Can be
     * used as a long flag on the command line.
     * \param desc - A description of what the argument is for or
     * does.
     * \param req - Whether the argument is required on the command
     * line.
     * \param constraint - A pointer to a Constraint object used
     * to constrain this Arg.
     * \param v - An optional visitor.  You probably should not
     * use this unless you have a very good reason.

    gsArgMultiVal( const std::string& flag,
                   const std::string& name,
                   const std::string& desc,
                   bool req,
                   TCLAP::Constraint<C>* constraint,
                   TCLAP::Visitor* v = NULL );
    */
    
    /*
     * Constructor.
     * \param flag - The one character flag that identifies this
     * argument on the command line.
     * \param name - A one word name for the argument.  Can be
     * used as a long flag on the command line.
     * \param desc - A description of what the argument is for or
     * does.
     * \param req - Whether the argument is required on the command
     * line.
     * \param constraint - A pointer to a Constraint object used
     * to constrain this Arg.
     * \param parser - A CmdLine parser object to add this Arg to
     * \param v - An optional visitor.  You probably should not
     * use this unless you have a very good reason.

    gsArgMultiVal( const std::string& flag, 
                   const std::string& name,
                   const std::string& desc,
                   bool req,
                   TCLAP::Constraint<C>* constraint,
                   TCLAP::CmdLineInterface& parser,
                   TCLAP::Visitor* v = NULL );
    */

	const std::vector<C>& getValue();

    operator TCLAP::MultiArg<C> &();

private:

    gsArgMultiValPrivate<C> * my;

}; // class gsArgMultiVal

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCmdLineArgs.hpp)
#endif
