/** @file gsCmdLineArgs.hpp
    
    @brief Provides implementation of input command line arguments.
    
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <tclap/MultiArg.h>

namespace gismo
{

template<class C> class gsArgMultiValPrivate
{
public:

    gsArgMultiValPrivate( const std::string& flag, 
                          const std::string& name,
                          const std::string& desc,
                          bool req,
                          const std::string& typeDesc,
                          TCLAP::CmdLineInterface& parser,
                          TCLAP::Visitor* v )
    : arg(flag,name,desc,req,typeDesc,parser,v)
    { }
    
    TCLAP::MultiArg<C> arg;  
};

template<class C>
gsArgMultiVal<C>::gsArgMultiVal( const std::string& flag, 
                                 const std::string& name,
                                 const std::string& desc,
                                 bool req,
                                 const std::string& typeDesc,
                                 TCLAP::CmdLineInterface& parser,
                                 TCLAP::Visitor* v )
: my(new gsArgMultiValPrivate<C>(flag,name,desc,req,typeDesc,parser,v))
{ } 

template<class C>
const std::vector<C> & gsArgMultiVal<C>::getValue()
{
    return my->arg.getValue();
}

template<class C>
gsArgMultiVal<C>::operator TCLAP::MultiArg<C> &()
{
    return my->arg;
}


/*
template<class C>
gsArgMultiVal<C>::gsArgMultiVal( const std::string& flag, 
                                 const std::string& name,
                                 const std::string& desc,
                                 bool req,
                                 const std::string& typeDesc,
                                 TCLAP::Visitor* v ): Base(flag,name,desc,req,typeDesc,v) {}; 

template<class C>
gsArgMultiVal<C>::gsArgMultiVal( const std::string& flag,
                                 const std::string& name,
                                 const std::string& desc,
                                 bool req,
                                 TCLAP::Constraint<C>* constraint,
                                 TCLAP::Visitor* v ) : Base(flag,name,desc,req,constraint,v) {} 
    
template<class C>
gsArgMultiVal<C>::gsArgMultiVal( const std::string& flag, 
                                 const std::string& name,
                                 const std::string& desc,
                                 bool req,
                                 TCLAP::Constraint<C>* constraint,
                                 TCLAP::CmdLineInterface& parser,
                                 TCLAP::Visitor* v ) : Base(flag,name,desc,req,constraint,parser,v) {} 
*/
   
} // namespace gismo
