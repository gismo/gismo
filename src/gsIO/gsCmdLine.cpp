/** @file gsCmdLine.cpp
    
    @brief Provides implemementation of input command line arguments.
    
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsIO/gsCmdLine.h>
//#include <gsIO/gsCmdLineArgs.h>

// --- start External files
#include <tclap/CmdLine.h>   
#include <tclap/ValueArg.h>
//#include <tclap/UnlabeledValueArg.h>
//#include <tclap/MultiArg.h>
//#include <tclap/UnlabeledMultiArg.h>
#include <tclap/SwitchArg.h>
//#include <tclap/MultiSwitchArg.h>
// --- end External files



namespace gismo
{

class gsCmdLinePrivate
{
public:
    /*
    typedef TCLAP::Arg                             Arg;
    typedef TCLAP::ValueArg<int>                   IntArg;
    typedef TCLAP::ValueArg<real_t>                RealArg;
    typedef TCLAP::ValueArg<std::string>           StrArg;
    typedef TCLAP::SwitchArg                       SwitchArg;
    typedef TCLAP::UnlabeledValueArg<std::string>  PlainStrArg;
    */
    
    gsCmdLinePrivate(const std::string& message,	
                     const char delimiter = ' ',
                     bool helpAndVersion = true)
        : cmd(message,delimiter,GISMO_VERSION,helpAndVersion), plainString(NULL)
    {
        cmd.setOutput( &cmdout );
    }

	~gsCmdLinePrivate()  
	{  
        freeAll( intVals   ); 
        freeAll( realVals  ); 
        freeAll( stringVals); 
        freeAll( switches  ); 
        delete plainString; 
    }

public:
    
    TCLAP::CmdLine cmd;

    std::map<std::string,TCLAP::Arg*> args;

    std::vector<std::string>        argstr;
    
    // Stores integer arguments
    std::vector<TCLAP::ValueArg<int>*    > intVals ;
    std::vector<int*>                      intRes ;

    // Stores real_t arguments
    std::vector<TCLAP::ValueArg<real_t>*> realVals;
    std::vector<real_t*>                  realRes ;

    // Stores string arguments
    std::vector<TCLAP::ValueArg<std::string>*> stringVals;
    std::vector<std::string*>                     strRes ;

    // Stores switch arguments
    std::vector<TCLAP::SwitchArg*      > switches;
    std::vector<bool*>                     swRes ;

    // Stores plain string argument
    TCLAP::UnlabeledValueArg<std::string> * plainString;
    std::string *                               pstrRes;

    // Stores config filename
    //std::string config;

    class GismoCmdOut : public TCLAP::StdOutput
    {
    public:
        void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e);
        void usage(TCLAP::CmdLineInterface& c);
        void version(TCLAP::CmdLineInterface& c);
    };
    
    GismoCmdOut cmdout;
};

gsCmdLine::gsCmdLine( const std::string& message,	
                      const char delimiter,
                      bool helpAndVersion)
: my(new gsCmdLinePrivate(message,delimiter,helpAndVersion))
{
/*
    // Config file
    my->stringVals.push_back(
        new TCLAP::ValueArg<std::string>("","config",
        "File containing configuration options",false,my->config,"string",my->cmd) );
    my->strRes.push_back(&config);
    //my->parsed = false;
    */
}

gsCmdLine::gsCmdLine(int argc, char *argv[], const std::string& message)
: my(new gsCmdLinePrivate(message,' ', true))
{
    my->argstr.assign(argv, argv + argc);
}

gsCmdLine::operator TCLAP::CmdLineInterface &()
{
    return my->cmd;
}

std::string & gsCmdLine::getMessage()
{
    return my->cmd.getMessage();
}

void gsCmdLine::addInt( const std::string& flag, 
                        const std::string& name, 
                        const std::string& desc, 
                        int              & value)
{
    //value = getInt(flag,name,desc,value);
    my->intVals.push_back(new TCLAP::ValueArg<int>(flag,name,desc,false,value,"int",my->cmd) );
    my->intRes.push_back(&value);
}

int gsCmdLine::getInt( const std::string& flag,
                       const std::string& name, 
                       const std::string& desc, 
                       const int        & value)
{
    ArgTable::const_iterator it = my->args.find(name);
    if ( it != my->args.end() )
    {
        gsWarn<< "Value already defined.\n";
        return value;
    }
    TCLAP::ValueArg<int> * a = 
        new TCLAP::ValueArg<int>(flag,name,desc,false,value,"int",my->cmd);
    my->intVals.push_back( a );
    my->args[name] = a;
    try 
    {
        for (int i = 1; static_cast<std::size_t>(i) < my->argstr.size(); i++) 
            if ( a->processArg( &i, my->argstr ) )
                break;
    }
    catch ( TCLAP::ArgException& e )
    { 
        gsWarn << "\nSomething went wrong when reading the command line.\n";
        gsWarn << "Error: " << e.error() << " " << e.argId() << "\n";
    }

    a->TCLAP::Arg::reset();
    return a->getValue();
}

void gsCmdLine::addReal( const std::string & flag, 
                         const std::string& name, 
                         const std::string& desc, 
                         real_t & value)
{
    value = getReal(flag,name,desc,value);
    my->realRes.push_back(&value);
}

real_t gsCmdLine::getReal( const std::string& flag,
                           const std::string& name, 
                           const std::string& desc, 
                           const real_t & value)
{
    ArgTable::const_iterator it = my->args.find(name);
    if ( it != my->args.end() )
    {
        gsWarn<< "Value already defined.\n";
        return value;
    }
    TCLAP::ValueArg<real_t> * a = 
        new TCLAP::ValueArg<real_t>(flag,name,desc,false,value,"float",my->cmd);
    my->realVals.push_back( a );
    my->args[name] = my->realVals.back();
    try 
    {
        for (int i = 1; static_cast<std::size_t>(i) < my->argstr.size(); i++) 
            if ( a->processArg( &i, my->argstr ) )
                break;
    }
    catch ( TCLAP::ArgException& e )
    { 
        gsWarn << "\nSomething went wrong when reading the command line.\n";
        gsWarn << "Error: " << e.error() << " " << e.argId() << "\n";
    }

    a->TCLAP::Arg::reset();
    return my->realVals.back()->getValue();
}

void gsCmdLine::addString( const std::string& flag, 
                           const std::string& name, 
                           const std::string& desc, 
                           std::string & value)
{
    my->stringVals.push_back(
    new TCLAP::ValueArg<std::string>(flag,name,desc,false,value,"string",my->cmd) );
    //value = getString(flag,name,desc,value);
    my->strRes.push_back(&value);
}

std::string gsCmdLine::getString( const std::string& flag,
                                  const std::string& name, 
                                  const std::string& desc, 
                                  const std::string & value)
{
    ArgTable::const_iterator it = my->args.find(name);
    if ( it != my->args.end() )
    {
        gsWarn<< "Value already defined.\n";
        return value;
    }
    TCLAP::ValueArg<std::string> * a =
        new TCLAP::ValueArg<std::string>(flag,name,desc,false,value,"string",my->cmd);
    my->stringVals.push_back( a );
    my->args[name] = a;
    try 
    {
        for (int i = 1; static_cast<std::size_t>(i) < my->argstr.size(); i++) 
            if ( a->processArg( &i, my->argstr ) )
                break;
    }
    catch ( TCLAP::ArgException& e )
    { 
        gsWarn << "\nSomething went wrong when reading the command line.\n";
        gsWarn << "Error: " << e.error() << " " << e.argId() << "\n";
    }

    a->TCLAP::Arg::reset();
    return my->stringVals.back()->getValue();
}
    
void gsCmdLine::addSwitch( const std::string& name, 
                           const std::string& desc, 
                           bool & value)
{
    addSwitch("",name,desc,value);
}

void gsCmdLine::addSwitch( const std::string& flag, 
                           const std::string& name, 
                           const std::string& desc, 
                           bool & value)
{
    my->switches.push_back(new TCLAP::SwitchArg("",name,desc,my->cmd) );
    //value = getSwitch(flag,name,desc,value);
    my->swRes.push_back(&value);
}

bool gsCmdLine::getSwitch( const std::string& flag,
                           const std::string& name, 
                           const std::string& desc, 
                           const bool & value)
{
    ArgTable::const_iterator it = my->args.find(name);
    if ( it != my->args.end() )
    {
        gsWarn<< "Value already defined.\n";
        return value;
    }
    TCLAP::SwitchArg * a = new TCLAP::SwitchArg(flag,name,desc,my->cmd);
    my->switches.push_back(a);
    my->args[name] = a;
    try 
    {
        for (int i = 1; static_cast<std::size_t>(i) < my->argstr.size(); i++) 
            if ( a->processArg( &i, my->argstr ) )
                break;
    }
    catch ( TCLAP::ArgException& e )
    { 
        gsWarn << "\nSomething went wrong when reading the command line.\n";
        gsWarn << "Error: " << e.error() << " " << e.argId() << "\n";
    }

    a->TCLAP::Arg::reset();
    return value | a->getValue();
}

bool gsCmdLine::getSwitch( const std::string& name, 
                           const std::string& desc, 
                           const bool & value)
{
    return getSwitch("",name,desc,value);
}

void gsCmdLine::addPlainString( const std::string& name, 
                                const std::string& desc, 
                                std::string & value)
{
    //value = getPlainString(name,desc,value);

    if ( my->plainString )
    {
        gsWarn<<"Plain string already added.\n";
        return;
    }
    else
    {
        my->plainString = 
            new TCLAP::UnlabeledValueArg<std::string>(name,desc,false,value,"string",my->cmd);
        my->pstrRes     = &value;
    }
    my->pstrRes = &value;
}

std::string gsCmdLine::getPlainString( const std::string& name, 
                                       const std::string& desc, 
                                       const std::string & value)
{
    if ( my->plainString && my->plainString->getValue() !=  name)
    {
        gsWarn<<"Plain string "<<my->plainString->getValue()<<" pre-exists.\n";
        return value;
    }

    my->plainString = 
        new TCLAP::UnlabeledValueArg<std::string>(name,desc,false,value,"string",my->cmd);
    my->args[name] = my->plainString;

    try 
    {
        for (int i = 1; static_cast<std::size_t>(i) < my->argstr.size(); i++) 
            if ( my->plainString->processArg( &i, my->argstr ) )
                break;
    }
    catch ( TCLAP::ArgException& e )
    { 
        gsWarn << "\nSomething went wrong when reading the command line.\n";
        gsWarn << "Error: " << e.error() << " " << e.argId() << "\n";
    }

    my->plainString->TCLAP::Arg::reset();
    return my->plainString->getValue();
}




bool gsCmdLine::valid() const 
{
    try 
    {
        my->cmd.parse(my->argstr);
    }
    catch ( TCLAP::ArgException& e )
    { 
        gsWarn << "\nSomething went wrong when reading the command line.\n";
        gsWarn << "Error: " << e.error() << " " << e.argId() << "\n";
        return false;
    }
    return true;
}

bool gsCmdLine::getValues(int argc, char *argv[])
{
    if (argc == 0)
        gsInfo << "Add \"-h\" to see a list of available command-line arguments.\n\n";

    try 
    {
        my->cmd.parse(argc,argv);
        
        for( std::size_t i=0; i!=my->intVals.size(); ++i)
            *my->intRes[i] = my->intVals[i]->getValue();

        for( std::size_t i=0; i!=my->realVals.size(); ++i)
            *my->realRes[i] = my->realVals[i]->getValue();

        for( std::size_t i=0; i!=my->stringVals.size(); ++i)
            *my->strRes[i] = my->stringVals[i]->getValue();

        for( std::size_t i=0; i!=my->switches.size(); ++i)
            *my->swRes[i] |= my->switches[i]->getValue();
        
        if ( my->plainString )
            *my->pstrRes = my->plainString->getValue();
    }
    catch ( TCLAP::ArgException& e )
    { 
        gsWarn << "\nSomething went wrong when reading the command line.\n";
        gsWarn << "Error: " << e.error() << " " << e.argId() << "\n"; 
        return false; 
    }
    
    return true;
}


gsCmdLine::~gsCmdLine() 
{ 
    delete my;
}

void gsCmdLinePrivate::GismoCmdOut::failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e)
{
    gsInfo << e.what() << "\n";
    gsInfo <<"\n USAGE: \n";
    _longUsage( c, gsInfo );
    throw;
}

void gsCmdLinePrivate::GismoCmdOut::usage(TCLAP::CmdLineInterface& c)
{
    std::string head = c.getProgramName();
    head += " is part of G+Smo.";
    spacePrint( gsInfo, head , 75, 3, 0 );
    gsInfo << "\n";
    spacePrint( gsInfo, c.getMessage() , 75, 3, 0 );
    
    gsInfo <<"\n Usage: \n";
    std::list<TCLAP::Arg*> argList = c.getArgList();
    TCLAP::XorHandler xorHandler   = c.getXorHandler();
    std::vector< std::vector<TCLAP::Arg*> > xorList = xorHandler.getXorList();
    
    // first the xor 
    for ( int i = 0; static_cast<unsigned int>(i) < xorList.size(); i++ )
	{
	    for ( TCLAP::ArgVectorIterator it = xorList[i].begin(); 
		  it != xorList[i].end(); 
		  it++ )
		{
		    spacePrint( gsInfo, (*it)->longID(), 75, 3, 3 );
		    spacePrint( gsInfo, (*it)->getDescription(), 75, 5, 0 );
		    
		    if ( it+1 != xorList[i].end() )
			spacePrint(gsInfo, "-- OR --", 75, 9, 0);
		}
	    gsInfo << "\n\n";
	}
    
    // then the rest
    for (TCLAP::ArgListIterator it = argList.begin(); it != argList.end(); it++)
	if ( !xorHandler.contains( (*it) ) )
	    {
		spacePrint( gsInfo, (*it)->longID(), 75, 3, 3 ); 
		spacePrint( gsInfo, (*it)->getDescription(), 75, 5, 0 ); 
		gsInfo << std::endl;
	    }    
}


void gsCmdLinePrivate::GismoCmdOut::version(TCLAP::CmdLineInterface& c)
{
    GISMO_UNUSED(c);
    //gsInfo <<"Executable: "<< c.getProgramName() <<", part of \n\n";
    gsInfo << "\n";
    gsInfo << "                   G+Smo \n";
    gsInfo << "      Geometry plus Simulation modules\n";
    gsInfo << "               version "<< GISMO_VERSION<<"\n";
    gsInfo << "   Compiled by ";
//https://sourceforge.net/p/predef/wiki/Compilers, see also boost/predef.h
#if defined(_MSC_VER)
    gsInfo << "MSVC "<<_MSC_FULL_VER <<" ("<<__cplusplus <<", ";
#elif defined(__clang__ )
    gsInfo << "Clang "<<"XX" <<" ("<<__cplusplus <<", ";
#elif defined(_INTEL_COMPILER)
    gsInfo << "Intel C++ "<<__INTEL_COMPILER<<" ("<<__cplusplus <<", ";
#elif defined(__MINGW64__)
    gsInfo << "MinGW "<<__MINGW64_VERSION_MAJOR<<"."<<__MINGW64_VERSION_MINOR<<" ("<<__cplusplus <<", ";
#elif defined(__SUNPRO_CC)
    gsInfo << "Solaris Studio "<<__SUNPRO_CC<<" ("<<__cplusplus <<", ";
#elif defined(__GNUG__)
    gsInfo << "GNU GCC "<<__GNUC__<<"."<<__GNUC_MINOR__<<" ("<<__cplusplus <<", ";
#else
    gsInfo << "C++ ("<<__cplusplus <<", ";
#endif

#ifdef __INTEL_MKL__
    gsInfo << "MKL "<<INTEL_MKL_VERSION<<", ";
#endif

#ifdef _LIBCPP_VERSION
    gsInfo << "libc++ "<<_LIBCPP_VERSION <<")\n";
#  elif defined(__GLIBCXX__)
    gsInfo << "glibstdc++ "<< __GLIBCXX__ <<")\n";
#  elif defined(__GLIBCPP__)
    gsInfo << "glibstdc++ "<< __GLIBCPP__ <<")\n";
#elif defined(__LIBCOMO__)
    gsInfo << "Comeau STL "<< __LIBCOMO__ <<")\n";
#  elif defined(__STL_CONFIG_H)
    gsInfo << "SGI STL)\n";
#  elif defined(__MSL_CPP__)
    gsInfo << "MSL standard lib)\n";
#  elif defined(__IBMCPP__)
    gsInfo << "VACPP STL)\n";
#  elif defined(MSIPL_COMPILE_H)
    gsInfo << "Modena C++ STL)\n";
#  elif (defined(_YVALS) && !defined(__IBMCPP__)) || defined(_CPPLIB_VER)
    gsInfo << "Dinkumware STL)\n";
#  elif defined(__STD_RWCOMPILER_H__) || defined(_RWSTD_VER)
    gsInfo << "Rogue Wave lib)\n";
#else
    gsInfo << "Unknown-STD)\n";
#endif
    gsInfo << "   RICAM-Linz 2012 - 2017, http://gs.jku.at/gismo\n";
}


} //namespace gismo
