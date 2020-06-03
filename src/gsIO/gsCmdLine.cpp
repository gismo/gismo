/** @file gsCmdLine.cpp

    @brief Provides implemementation of input command line arguments.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, S. Takacs
*/

#include <gsIO/gsCmdLine.h>

// --- start External files
#include <tclap/CmdLine.h>
#include <tclap/ValueArg.h>
//#include <tclap/UnlabeledValueArg.h>
//#include <tclap/MultiArg.h>
//#include <tclap/UnlabeledMultiArg.h>
#include <tclap/SwitchArg.h>
//#include <tclap/MultiSwitchArg.h>
// --- end External files

#include <gsIO/gsOptionList.h>

namespace gismo
{

class gsCmdLinePrivate
{
    typedef index_t intVal_t;
public:
    /*
    typedef TCLAP::Arg                             Arg;
    typedef TCLAP::ValueArg<intVal_t>              IntArg;
    typedef TCLAP::ValueArg<real_t>                RealArg;
    typedef TCLAP::ValueArg<std::string>           StrArg;
    typedef TCLAP::SwitchArg                       SwitchArg;
    typedef TCLAP::UnlabeledValueArg<std::string>  PlainStrArg;
    */

    gsCmdLinePrivate(const std::string& message,
                     const char delimiter = ' ',
                     bool helpAndVersion = true)
        : cmd(message,delimiter,GISMO_VERSION,helpAndVersion), plainStringVal(NULL)
#ifndef NDEBUG
          , didParseCmdLine(false)
#endif
    {
        cmd.setOutput( &cmdout );
    }

    ~gsCmdLinePrivate()
    {
        freeAll( intVals        );
        freeAll( multiIntVals   );
        freeAll( realVals       );
        freeAll( multiRealVals  );
        freeAll( stringVals     );
        freeAll( multiStringVals);
        freeAll( switchVals     );
        delete plainStringVal;
#ifndef NDEBUG
        if (!didParseCmdLine)
        {
            gsWarn<< "gsCmdLine::getValues was never called. "
                "Arguments were not parsed from command line.\n"<<std::flush;
            abort();
        }
#endif
    }

public:

    TCLAP::CmdLine cmd;

    // Stores integer arguments
    std::vector<TCLAP::ValueArg<intVal_t>*>    intVals;
    std::vector<intVal_t*>                     intRes;

    // Stores multi integer arguments
    std::vector<TCLAP::MultiArg<intVal_t>*>    multiIntVals;
    std::vector<std::vector<intVal_t>*>        multiIntRes;

    // Stores real_t arguments
    std::vector<TCLAP::ValueArg<real_t>*>      realVals;
    std::vector<real_t*>                       realRes;

    // Stores multi real_t arguments
    std::vector<TCLAP::MultiArg<real_t>*>      multiRealVals;
    std::vector<std::vector<real_t>*>          multiRealRes;

    // Stores string arguments
    std::vector<TCLAP::ValueArg<std::string>*> stringVals;
    std::vector<std::string*>                  stringRes;

    // Stores multi string arguments
    std::vector<TCLAP::MultiArg<std::string>*> multiStringVals;
    std::vector<std::vector<std::string>*>     multiStringRes;

    // Stores switch arguments
    std::vector<TCLAP::SwitchArg*>             switchVals;
    std::vector<bool*>                         switchRes;

    // Stores plain string argument
    TCLAP::UnlabeledValueArg<std::string> *    plainStringVal;
    std::string *                              plainStringRes;

    // Stores config filename
    //std::string config;

    class GismoCmdOut : public TCLAP::StdOutput
    {
    public:
        void failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e);
        void usage(TCLAP::CmdLineInterface& c);
        void version(TCLAP::CmdLineInterface& c);
    };

    class GismoNullOut : public TCLAP::CmdLineOutput
    {
    public:
        void failure(TCLAP::CmdLineInterface& , TCLAP::ArgException& ) { }
        void usage(TCLAP::CmdLineInterface& )  { }
        void version(TCLAP::CmdLineInterface& ) { }
    };


    static GismoCmdOut cmdout;
    static GismoNullOut cmdNullOut;

#ifndef NDEBUG
    bool didParseCmdLine;
#endif
};

gsCmdLinePrivate::GismoCmdOut gsCmdLinePrivate::cmdout;
gsCmdLinePrivate::GismoNullOut gsCmdLinePrivate::cmdNullOut;

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
    my->stringRes.push_back(&config);
    //my->parsed = false;
    */
}

void gsCmdLine::addInt( const std::string& flag,
                        const std::string& name,
                        const std::string& desc,
                        intVal_t         & value)
{
    GISMO_ASSERT( flag.size() < 2, "The flag (short form) must be empty or one character only." );
    GISMO_ASSERT( name.size() > 1, "The name (long form of the flag) must be at least two characters." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->intVals.push_back(new TCLAP::ValueArg<intVal_t>(flag,name,desc,false,value,"int",my->cmd) );
    my->intRes.push_back(&value);
}

void gsCmdLine::addMultiInt( const std::string    & flag,
                             const std::string    & name,
                             const std::string    & desc,
                             std::vector<intVal_t>& value)
{
    GISMO_ASSERT( flag.size() < 2, "The flag (short form) must be empty or one character only." );
    GISMO_ASSERT( name.size() > 1, "The name (long form of the flag) must be at least two characters." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->multiIntVals.push_back(new TCLAP::MultiArg<intVal_t>(flag,name,desc,false,"int",my->cmd) );
    my->multiIntRes.push_back(&value);
}

void gsCmdLine::addReal( const std::string& flag,
                         const std::string& name,
                         const std::string& desc,
                         real_t           & value)
{
    GISMO_ASSERT( flag.size() < 2, "The flag (short form) must be empty or one character only." );
    GISMO_ASSERT( name.size() > 1, "The name (long form of the flag) must be at least two characters." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->realVals.push_back(new TCLAP::ValueArg<real_t>(flag,name,desc,false,value,"float",my->cmd) );
    my->realRes.push_back(&value);
}

void gsCmdLine::addMultiReal( const std::string  & flag,
                              const std::string  & name,
                              const std::string  & desc,
                              std::vector<real_t>& value)
{
    GISMO_ASSERT( flag.size() < 2, "The flag (short form) must be empty or one character only." );
    GISMO_ASSERT( name.size() > 1, "The name (long form of the flag) must be at least two characters." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->multiRealVals.push_back(new TCLAP::MultiArg<real_t>(flag,name,desc,false,"float",my->cmd) );
    my->multiRealRes.push_back(&value);
}

void gsCmdLine::addString( const std::string& flag,
                           const std::string& name,
                           const std::string& desc,
                           std::string      & value)
{
    GISMO_ASSERT( flag.size() < 2, "The flag (short form) must be empty or one character only." );
    GISMO_ASSERT( name.size() > 1, "The name (long form of the flag) must be at least two characters." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->stringVals.push_back(new TCLAP::ValueArg<std::string>(flag,name,desc,false,value,"string",my->cmd));
    my->stringRes.push_back(&value);
}

void gsCmdLine::addMultiString( const std::string       & flag,
                                const std::string       & name,
                                const std::string       & desc,
                                std::vector<std::string>& value)
{
    GISMO_ASSERT( flag.size() < 2, "The flag (short form) must be empty or one character only." );
    GISMO_ASSERT( name.size() > 1, "The name (long form of the flag) must be at least two characters." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->multiStringVals.push_back(new TCLAP::MultiArg<std::string>(flag,name,desc,false,"string",my->cmd) );
    my->multiStringRes.push_back(&value);
}

void gsCmdLine::addSwitch( const std::string& flag,
                           const std::string& name,
                           const std::string& desc,
                           bool             & value)
{
    GISMO_ASSERT( flag.size() < 2, "The flag (short form) must be empty or one character only." );
    GISMO_ASSERT( name.size() > 1, "The name (long form of the flag) must be at least two characters." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->switchVals.push_back(new TCLAP::SwitchArg(flag,name,desc,my->cmd) );
    my->switchRes.push_back(&value);
}

void gsCmdLine::addPlainString( const std::string& name,
                                const std::string& desc,
                                std::string      & value)
{
    GISMO_ASSERT( !name.empty(), "The name (long form of the flag) must not be empty." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );

    GISMO_ENSURE( !my->plainStringVal, "Plain string already added." );

    my->plainStringVal = new TCLAP::UnlabeledValueArg<std::string>(name,desc,false,value,"string",my->cmd);
    my->plainStringRes = &value;
}


bool gsCmdLine::valid(int argc, char *argv[]) const
{
    const bool eh = my->cmd.getExceptionHandling();
    TCLAP::CmdLineOutput * o = my->cmd.getOutput();
    my->cmd.setExceptionHandling(false);
    my->cmd.setOutput( &my->cmdNullOut );
    bool result = true;
    try
    {
        my->cmd.parse(argc,argv);
    }
    catch ( TCLAP::ExitException& )   { /*result = true;*/  }
    //catch ( TCLAP::ArgException&  ) { result = false;     }
    catch (...)                       { result = false;     }

    my->cmd.reset();
    my->cmd.setExceptionHandling(eh);
    my->cmd.setOutput(o);
    return result;
}


void gsCmdLine::getValues(int argc, char *argv[])
{
    GISMO_ASSERT( !my->didParseCmdLine, "gsCmdLine::getValues must not be called twice." );
#ifndef NDEBUG
    my->didParseCmdLine = true;
#endif

    my->cmd.parse(argc,argv);

    for( size_t i=0; i!=my->intVals.size(); ++i)
        *my->intRes[i] = my->intVals[i]->getValue();

    for( size_t i=0; i!=my->realVals.size(); ++i)
        *my->realRes[i] = my->realVals[i]->getValue();

    for( size_t i=0; i!=my->stringVals.size(); ++i)
        *my->stringRes[i] = my->stringVals[i]->getValue();

    for( size_t i=0; i!=my->switchVals.size(); ++i)
        // Toggle switch-result if switch is present
        *my->switchRes[i] ^= my->switchVals[i]->getValue();

    if ( my->plainStringVal )
        *my->plainStringRes = my->plainStringVal->getValue();

    for( size_t i=0; i!=my->multiIntVals.size(); ++i)
        if( my->multiIntVals[i]->isSet() )
            *my->multiIntRes[i] = my->multiIntVals[i]->getValue();

    for( size_t i=0; i!=my->multiRealVals.size(); ++i)
        if( my->multiRealVals[i]->isSet() )
            *my->multiRealRes[i] = my->multiRealVals[i]->getValue();

    for( size_t i=0; i!=my->multiStringVals.size(); ++i)
        if( my->multiStringVals[i]->isSet() )
            *my->multiStringRes[i] = my->multiStringVals[i]->getValue();
}

void gsCmdLine::setExceptionHandling(const bool state)
{
    my->cmd.setExceptionHandling(state);
}

bool gsCmdLine::getExceptionHandling() const
{
    return my->cmd.getExceptionHandling();
}

#define ADD_OPTION_LIST_ENTRY(res,vals,addFct)                                      \
{                                                                                   \
    std::string nm = (vals)->getName() + ".";                                       \
    index_t sz = static_cast<index_t>((res).size());                                \
    for ( index_t j=0; j<sz; ++j )                                                  \
    { result.addFct( nm+util::to_string(j), (vals)->getDescription(), (res)[j] ); } \
    result.addInt( nm+"Size", (vals)->getDescription(), sz );                       \
}

gsOptionList gsCmdLine::getOptionList()
{
    GISMO_ASSERT( my->didParseCmdLine, "gsCmdLine::getOptionList can be called only after gsCmdLine::getValues." );

    gsOptionList result;
    for( size_t i=0; i!=my->intVals.size(); ++i)
        result.addInt( my->intVals[i]->getName(), my->intVals[i]->getDescription(), *my->intRes[i] );
    for( size_t i=0; i!=my->realVals.size(); ++i)
        result.addReal( my->realVals[i]->getName(), my->realVals[i]->getDescription(), *my->realRes[i] );
    for( size_t i=0; i!=my->stringVals.size(); ++i)
        result.addString( my->stringVals[i]->getName(), my->stringVals[i]->getDescription(), *my->stringRes[i] );
    for( size_t i=0; i!=my->switchVals.size(); ++i)
        result.addSwitch( my->switchVals[i]->getName(), my->switchVals[i]->getDescription(), *my->switchRes[i] );
    for( size_t i=0; i!=my->multiIntVals.size(); ++i)
        ADD_OPTION_LIST_ENTRY(*my->multiIntRes[i],my->multiIntVals[i],addInt)
    for( size_t i=0; i!=my->multiRealVals.size(); ++i)
        ADD_OPTION_LIST_ENTRY(*my->multiRealRes[i],my->multiRealVals[i],addReal)
    for( size_t i=0; i!=my->multiStringVals.size(); ++i)
        ADD_OPTION_LIST_ENTRY(*my->multiStringRes[i],my->multiStringVals[i],addString)
    if ( my->plainStringVal )
        result.addString( my->plainStringVal->getName(), my->plainStringVal->getDescription(), *my->plainStringRes );
    return result;
}

#undef ADD_OPTION_LIST_ENTRY

gsCmdLine::~gsCmdLine()
{
    delete my;
}

void gsCmdLinePrivate::GismoCmdOut::failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e)
{
    gsInfo << " ERROR: " << e.what() << "\n";
    gsInfo <<"\n USAGE: \n";
    //_longUsage( c, gsInfo );
    this->usage(c);
    throw TCLAP::ExitException(1);
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
    for ( size_t i = 0; static_cast<unsigned int>(i) < xorList.size(); i++ )
    {
        for ( TCLAP::ArgVectorIterator it = xorList[i].begin();
                it != xorList[i].end();
                it++
            )
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


void gsCmdLinePrivate::GismoCmdOut::version(TCLAP::CmdLineInterface&)
{
    gsCmdLine::printVersion();
}

void gsCmdLine::printVersion()
{
    //gsInfo <<"Executable: "<< c.getProgramName() <<", part of \n\n";
    gsInfo << "\n";
    gsInfo << "                   G+Smo \n";
    gsInfo << "      Geometry plus Simulation modules\n";
    gsInfo << "               version "<< GISMO_VERSION<<"\n";
    gsInfo << "Compiled by ";
//https://sourceforge.net/p/predef/wiki/Compilers, see also boost/predef.h
#if defined(_MSC_VER) && _MSC_VER < 1600
    gsInfo << "MSVC "<<_MSC_FULL_VER <<" ("<<"199711L" <<", ";
#elsif _MSC_VER >= 1900
    gsInfo << "MSVC "<<_MSC_FULL_VER <<" ("<<_MSVC_LANG <<", ";
#elsif _MSC_VER >= 1600
    gsInfo << "MSVC "<<_MSC_FULL_VER <<" ("<<"201103L" <<", ";
#elif defined(__clang__ )
    gsInfo << "Clang "<<__clang_version__<<" ("<<__cplusplus <<", ";
#elif defined(_INTEL_COMPILER)
    gsInfo << "Intel C++ "<<__INTEL_COMPILER<<" ("<<__cplusplus <<", ";
#elif defined(__MINGW64__)
    gsInfo << "MinGW "<<__MINGW64_VERSION_MAJOR<<"."<<__MINGW64_VERSION_MINOR<<" ("<<__cplusplus <<", ";
#elif defined(__SUNPRO_CC)
    gsInfo << "Solaris Studio "<<__SUNPRO_CC<<" ("<<__cplusplus <<", ";
#elif defined(__GNUG__)
    gsInfo << "GNU GCC "<<__GNUC__<<"."<<__GNUC_MINOR__<<"."<<__GNUC_PATCHLEVEL__<<" ("<<__cplusplus <<", ";
#else
    gsInfo << "C++ ("<<__cplusplus <<", ";
#endif

#ifdef __INTEL_MKL__
    gsInfo << "MKL "<<INTEL_MKL_VERSION<<", ";
#endif

#ifdef _LIBCPP_VERSION
    gsInfo << "libc++ "<<_LIBCPP_VERSION <<")\n";
#  elif defined(__GLIBCXX__)
    gsInfo << "glibc++ "<< __GLIBCXX__ <<")\n";
#  elif defined(__GLIBCPP__)
    gsInfo << "glibc++ "<< __GLIBCPP__ <<")\n";
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
    gsInfo << "Dinkumware STL "<< _CPPLIB_VER<<")\n";
#  elif defined(__STD_RWCOMPILER_H__) || defined(_RWSTD_VER)
    gsInfo << "Rogue Wave lib "<<_RWSTD_VER<<")\n";
#else
    gsInfo << "Unknown-STD)\n";
#endif
    //gsInfo << "Eigen "<< EIGEN_WORLD_VERSION<<"."<<EIGEN_MAJOR_VERSION<<"."<<EIGEN_MINOR_VERSION<<"\n";
    gsInfo << "web: http://github.com/gismo\n";
}

std::string & gsCmdLine::getMessage()
{
    return my->cmd.getMessage();
}

} //namespace gismo
