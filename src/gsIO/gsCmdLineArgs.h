
#pragma once

#include <gsCore/gsConfig.h>
#include <gsCore/gsExport.h>
#include <gsCore/gsDebug.h>

#include <tclap/SwitchArg.h>
#include <tclap/MultiSwitchArg.h>
#include <tclap/UnlabeledValueArg.h>
#include <tclap/UnlabeledMultiArg.h>


namespace gismo
{


template<class C = real_t>
class gsArgVal :  public TCLAP::ValueArg<C>
{
public:
    typedef typename TCLAP::ValueArg<C> Base;
public:

        /**
         * Labeled ValueArg constructor.
         * You could conceivably call this constructor with a blank flag, 
         * but that would make you a bad person.  It would also cause
         * an exception to be thrown.   If you want an unlabeled argument, 
         * use the other constructor.
         * \param flag - The one character flag that identifies this
         * argument on the command line.
         * \param name - A one word name for the argument.  Can be
         * used as a long flag on the command line.
         * \param desc - A description of what the argument is for or
         * does.
         * \param req - Whether the argument is required on the command
         * line.
         * \param value - The default value assigned to this argument if it
         * is not present on the command line.
         * \param typeDesc - A short, human readable description of the
         * type that this object expects.  This is used in the generation
         * of the USAGE statement.  The goal is to be helpful to the end user
         * of the program.
         * \param v - An optional visitor.  You probably should not
         * use this unless you have a very good reason.
         */
    gsArgVal( const std::string& flag, 
              const std::string& name, 
              const std::string& desc, 
              bool req, 
              C value,
              const std::string& typeDesc,
              TCLAP::Visitor* v = NULL) ;
				 
				 
        /**
         * Labeled ValueArg constructor.
         * You could conceivably call this constructor with a blank flag, 
         * but that would make you a bad person.  It would also cause
         * an exception to be thrown.   If you want an unlabeled argument, 
         * use the other constructor.
         * \param flag - The one character flag that identifies this
         * argument on the command line.
         * \param name - A one word name for the argument.  Can be
         * used as a long flag on the command line.
         * \param desc - A description of what the argument is for or
         * does.
         * \param req - Whether the argument is required on the command
         * line.
         * \param value - The default value assigned to this argument if it
         * is not present on the command line.
         * \param typeDesc - A short, human readable description of the
         * type that this object expects.  This is used in the generation
         * of the USAGE statement.  The goal is to be helpful to the end user
         * of the program.
         * \param parser - A CmdLine parser object to add this Arg to
         * \param v - An optional visitor.  You probably should not
         * use this unless you have a very good reason.
         */
        gsArgVal( const std::string& flag, 
                  const std::string& name, 
                  const std::string& desc, 
                  bool req, 
                  C value,
                  const std::string& typeDesc,
                  TCLAP::CmdLineInterface& parser,
                  TCLAP::Visitor* v = NULL );
 
        /**
         * Labeled ValueArg constructor.
         * You could conceivably call this constructor with a blank flag, 
         * but that would make you a bad person.  It would also cause
         * an exception to be thrown.   If you want an unlabeled argument, 
         * use the other constructor.
         * \param flag - The one character flag that identifies this
         * argument on the command line.
         * \param name - A one word name for the argument.  Can be
         * used as a long flag on the command line.
         * \param desc - A description of what the argument is for or
         * does.
         * \param req - Whether the argument is required on the command
         * line.
         * \param value - The default value assigned to this argument if it
         * is not present on the command line.
         * \param constraint - A pointer to a Constraint object used
		 * to constrain this Arg.
         * \param parser - A CmdLine parser object to add this Arg to.
         * \param v - An optional visitor.  You probably should not
         * use this unless you have a very good reason.
         */
        gsArgVal( const std::string& flag, 
                  const std::string& name, 
                  const std::string& desc, 
                  bool req, 
                  C value,
                  TCLAP::Constraint<C>* constraint,
                  TCLAP::CmdLineInterface& parser,
                  TCLAP::Visitor* v = NULL );
	  
        /**
         * Labeled ValueArg constructor.
         * You could conceivably call this constructor with a blank flag, 
         * but that would make you a bad person.  It would also cause
         * an exception to be thrown.   If you want an unlabeled argument, 
         * use the other constructor.
         * \param flag - The one character flag that identifies this
         * argument on the command line.
         * \param name - A one word name for the argument.  Can be
         * used as a long flag on the command line.
         * \param desc - A description of what the argument is for or
         * does.
         * \param req - Whether the argument is required on the command
         * line.
         * \param value - The default value assigned to this argument if it
         * is not present on the command line.
         * \param constraint - A pointer to a Constraint object used
		 * to constrain this Arg.
         * \param v - An optional visitor.  You probably should not
         * use this unless you have a very good reason.
         */
        gsArgVal( const std::string& flag, 
                  const std::string& name, 
                  const std::string& desc, 
                  bool req, 
                  C value,
                  TCLAP::Constraint<C>* constraint,
                  TCLAP::Visitor* v = NULL );


}; // class gsArgVal

class GISMO_EXPORT gsArgSwitch :  public TCLAP::SwitchArg
{
public:
    typedef TCLAP::SwitchArg Base;
public:

        /**
         * SwitchArg constructor.
         * \param flag - The one character flag that identifies this
         * argument on the command line.
         * \param name - A one word name for the argument.  Can be
         * used as a long flag on the command line.
         * \param desc - A description of what the argument is for or
         * does.
         * \param def - The default value for this Switch. 
         * \param v - An optional visitor.  You probably should not
         * use this unless you have a very good reason.
         */
    gsArgSwitch(const std::string& flag, 
                const std::string& name, 
                const std::string& desc,
                bool def = false,
                TCLAP::Visitor* v = NULL);
    
				  
    /**
     * SwitchArg constructor.
     * \param flag - The one character flag that identifies this
     * argument on the command line.
     * \param name - A one word name for the argument.  Can be
     * used as a long flag on the command line.
     * \param desc - A description of what the argument is for or
     * does.
     * \param parser - A CmdLine parser object to add this Arg to
     * \param def - The default value for this Switch.
     * \param v - An optional visitor.  You probably should not
     * use this unless you have a very good reason.
     */
    gsArgSwitch(const std::string& flag, 
                const std::string& name, 
                const std::string& desc,
                TCLAP::CmdLineInterface& parser,
                bool def = false,
                TCLAP::Visitor* v = NULL);



}; // class gsArgSwitch


template<class C>
class gsArgValPlain : public TCLAP::UnlabeledValueArg<C>
{
public:
    typedef typename TCLAP::UnlabeledValueArg<C> Base;

public:

    /**
     * UnlabeledValueArg constructor.
     * \param name - A one word name for the argument.  Note that this is used for
     * identification, not as a long flag.
     * \param desc - A description of what the argument is for or
     * does.
     * \param req - Whether the argument is required on the command
     * line.
     * \param value - The default value assigned to this argument if it
     * is not present on the command line.
     * \param typeDesc - A short, human readable description of the
     * type that this object expects.  This is used in the generation
     * of the USAGE statement.  The goal is to be helpful to the end user
     * of the program.
     * \param ignoreable - Allows you to specify that this argument can be
     * ignored if the '--' flag is set.  This defaults to false (cannot
     * be ignored) and should  generally stay that way unless you have 
     * some special need for certain arguments to be ignored.
     * \param v - Optional Vistor.  You should leave this blank unless
     * you have a very good reason.
     */
    gsArgValPlain( const std::string& name, 
                   const std::string& desc, 
                   bool req,
                   C value,
                   const std::string& typeDesc,
                   bool ignoreable = false,
                   TCLAP::Visitor* v = NULL); 
    
		/**
		 * UnlabeledValueArg constructor.
		 * \param name - A one word name for the argument.  Note that this is used for
		 * identification, not as a long flag.
		 * \param desc - A description of what the argument is for or
		 * does.
		 * \param req - Whether the argument is required on the command
		 * line.
		 * \param value - The default value assigned to this argument if it
		 * is not present on the command line.
		 * \param typeDesc - A short, human readable description of the
		 * type that this object expects.  This is used in the generation
		 * of the USAGE statement.  The goal is to be helpful to the end user
		 * of the program.
		 * \param parser - A CmdLine parser object to add this Arg to
		 * \param ignoreable - Allows you to specify that this argument can be
		 * ignored if the '--' flag is set.  This defaults to false (cannot
		 * be ignored) and should  generally stay that way unless you have 
		 * some special need for certain arguments to be ignored.
		 * \param v - Optional Vistor.  You should leave this blank unless
		 * you have a very good reason.
		 */
		gsArgValPlain( const std::string& name, 
                               const std::string& desc, 
                               bool req,
                               C value,
                               const std::string& typeDesc,
                               TCLAP::CmdLineInterface& parser,
                               bool ignoreable = false,
                               TCLAP::Visitor* v = NULL ); 					
						
    /**
     * UnlabeledValueArg constructor.
     * \param name - A one word name for the argument.  Note that this is used for
     * identification, not as a long flag.
     * \param desc - A description of what the argument is for or
     * does.
     * \param req - Whether the argument is required on the command
     * line.
     * \param value - The default value assigned to this argument if it
     * is not present on the command line.
     * \param constraint - A pointer to a Constraint object used
     * to constrain this Arg.
     * \param ignoreable - Allows you to specify that this argument can be
     * ignored if the '--' flag is set.  This defaults to false (cannot
     * be ignored) and should  generally stay that way unless you have 
     * some special need for certain arguments to be ignored.
     * \param v - Optional Vistor.  You should leave this blank unless
     * you have a very good reason.
     */
    gsArgValPlain( const std::string& name, 
                   const std::string& desc, 
                   bool req,
                   C value,
                   TCLAP::Constraint<C>* constraint,
                   bool ignoreable = false,
                   TCLAP::Visitor* v = NULL ); 
    
    
    /**
     * UnlabeledValueArg constructor.
     * \param name - A one word name for the argument.  Note that this is used for
     * identification, not as a long flag.
     * \param desc - A description of what the argument is for or
     * does.
     * \param req - Whether the argument is required on the command
     * line.
     * \param value - The default value assigned to this argument if it
     * is not present on the command line.
     * \param constraint - A pointer to a Constraint object used
     * to constrain this Arg.
     * \param parser - A CmdLine parser object to add this Arg to
     * \param ignoreable - Allows you to specify that this argument can be
     * ignored if the '--' flag is set.  This defaults to false (cannot
     * be ignored) and should  generally stay that way unless you have 
     * some special need for certain arguments to be ignored.
     * \param v - Optional Vistor.  You should leave this blank unless
     * you have a very good reason.
     */
    gsArgValPlain( const std::string& name, 
                   const std::string& desc, 
                   bool req,
                   C value,
                   TCLAP::Constraint<C>* constraint,
                   TCLAP::CmdLineInterface& parser,
                   bool ignoreable = false,
                   TCLAP::Visitor* v = NULL);
    

}; // class gsArgValPlain

template<class C>
class gsArgMultiVal :  public TCLAP::MultiArg<C>
{
public:
    typedef typename TCLAP::MultiArg<C> Base;

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
     * \param v - An optional visitor.  You probably should not
     * use this unless you have a very good reason.
     */
    gsArgMultiVal( const std::string& flag,
                   const std::string& name,
                   const std::string& desc,
                   bool req,
                   const std::string& typeDesc,
                   TCLAP::Visitor* v = NULL); 
    
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
     * \param constraint - A pointer to a Constraint object used
     * to constrain this Arg.
     * \param v - An optional visitor.  You probably should not
     * use this unless you have a very good reason.
     */
    gsArgMultiVal( const std::string& flag,
                   const std::string& name,
                   const std::string& desc,
                   bool req,
                   TCLAP::Constraint<C>* constraint,
                   TCLAP::Visitor* v = NULL );
    
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
     * \param constraint - A pointer to a Constraint object used
     * to constrain this Arg.
     * \param parser - A CmdLine parser object to add this Arg to
     * \param v - An optional visitor.  You probably should not
     * use this unless you have a very good reason.
     */
    gsArgMultiVal( const std::string& flag, 
                   const std::string& name,
                   const std::string& desc,
                   bool req,
                   TCLAP::Constraint<C>* constraint,
                   TCLAP::CmdLineInterface& parser,
                   TCLAP::Visitor* v = NULL );


}; // class gsArgMultiVal

} // namespace gismo

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCmdLineArgs.hpp)
#include GISMO_HPP_HEADER(gsCmdLineArgs.cpp)
#endif
