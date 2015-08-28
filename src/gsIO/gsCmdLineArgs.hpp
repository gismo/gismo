/** @file gsCmdLineArgs.hpp
    
    @brief Provides implementation of input command line arguments.
    
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsCore/gsDebug.h>

namespace gismo
{

template<class C>
gsArgVal<C>::gsArgVal( const std::string& flag, 
                       const std::string& name, 
                       const std::string& desc, 
                       bool req, 
                       C value,
                       const std::string& typeDesc,
                       TCLAP::Visitor* v) 
    : Base(flag,name,desc,req,value,typeDesc){ }

template<class C>
gsArgVal<C>::gsArgVal( const std::string& flag, 
                           const std::string& name, 
                           const std::string& desc, 
                           bool req, 
                           C value,
                           const std::string& typeDesc,
                           TCLAP::CmdLineInterface& parser,
                           TCLAP::Visitor* v ) 
    : Base(flag,name,desc,req,value,typeDesc,parser) { }

template<class C>
gsArgVal<C>::gsArgVal( const std::string& flag, 
                           const std::string& name, 
                           const std::string& desc, 
                           bool req, 
                           C value,
                           TCLAP::Constraint<C>* constraint,
                           TCLAP::CmdLineInterface& parser,
        TCLAP::Visitor* v )  : Base(flag,name,desc,req,value,constraint,parser,v){ }

template<class C>
gsArgVal<C>::gsArgVal(const std::string& flag, 
                      const std::string& name, 
                      const std::string& desc, 
                      bool req, 
                      C value,
                      TCLAP::Constraint<C>* constraint,
        TCLAP::Visitor* v ): Base(flag,name,desc,req,value,constraint) { }
 
template<class C>
gsArgValPlain<C>::gsArgValPlain( const std::string& name, 
                           const std::string& desc, 
                           bool req,
                           C value,
                           const std::string& typeDesc,
                           bool ignoreable ,
        TCLAP::Visitor* v ) : Base(name,desc,req,value,typeDesc,ignoreable,v) {} 
    
template<class C>
gsArgValPlain<C>::gsArgValPlain( const std::string& name, 
                           const std::string& desc, 
                           bool req,
                           C value,
                           const std::string& typeDesc,
                           TCLAP::CmdLineInterface& parser,
                           bool ignoreable ,
                           TCLAP::Visitor* v ): Base(name,desc,req,value,typeDesc,parser,ignoreable,v) {} 
    
template<class C>
gsArgValPlain<C>::gsArgValPlain( const std::string& name, 
               const std::string& desc, 
               bool req,
               C value,
               TCLAP::Constraint<C>* constraint,
               bool ignoreable ,
               TCLAP::Visitor* v  ): Base(name,desc,req,value,constraint,ignoreable,v) {}; 

template<class C>
    gsArgValPlain<C>::gsArgValPlain( const std::string& name, 
                              const std::string& desc, 
                              bool req,
                              C value,
                              TCLAP::Constraint<C>* constraint,
                              TCLAP::CmdLineInterface& parser,
                              bool ignoreable ,
                              TCLAP::Visitor* v ): Base(name,desc,req,value,constraint,parser,ignoreable,v) {} 


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
                                 const std::string& typeDesc,
                                 TCLAP::CmdLineInterface& parser,
                                 TCLAP::Visitor* v ): Base(flag,name,desc,req,typeDesc,parser,v) {} 
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
 
   
} // namespace gismo
