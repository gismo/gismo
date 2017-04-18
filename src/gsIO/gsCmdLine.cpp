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
        freeAll( switches       );
        delete plainString;
    }

public:

    TCLAP::CmdLine cmd;

    std::map<std::string,TCLAP::Arg*> args;
    std::vector<std::string>        argstr;

    // Stores integer arguments
    std::vector<TCLAP::ValueArg<int>*>         intVals;
    std::vector<int*>                          intRes;

    // Stores multi integer arguments
    std::vector<TCLAP::MultiArg<int>*>         multiIntVals;
    std::vector<std::vector<int>*>             multiIntRes;

    // Stores real_t arguments
    std::vector<TCLAP::ValueArg<real_t>*>      realVals;
    std::vector<real_t*>                       realRes;

    // Stores multi real_t arguments
    std::vector<TCLAP::MultiArg<real_t>*>      multiRealVals;
    std::vector<std::vector<real_t>*>          multiRealRes;

    // Stores string arguments
    std::vector<TCLAP::ValueArg<std::string>*> stringVals;
    std::vector<std::string*>                  strRes;

    // Stores multi string arguments
    std::vector<TCLAP::MultiArg<std::string>*> multiStringVals;
    std::vector<std::vector<std::string>*>     multiStrRes;

    // Stores switch arguments
    std::vector<TCLAP::SwitchArg*>             switches;
    std::vector<bool*>                         swRes;

    // Stores plain string argument
    TCLAP::UnlabeledValueArg<std::string> *    plainString;
    std::string *                              pstrRes;

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
#ifndef NDEBUG
    bool didParseCmdLine;
#endif
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

void gsCmdLine::addInt( const std::string& flag, 
                        const std::string& name, 
                        const std::string& desc, 
                        int              & value)
{
    GISMO_ASSERT( !name.empty(), "The name (long form of the flag) must not be empty." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->intVals.push_back(new TCLAP::ValueArg<int>(flag,name,desc,false,value,"int",my->cmd) );
    my->intRes.push_back(&value);
}

void gsCmdLine::addMultiInt( const std::string& flag, 
                             const std::string& name, 
                             const std::string& desc, 
                             std::vector<int> & value)
{
    GISMO_ASSERT( !name.empty(), "The name (long form of the flag) must not be empty." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->multiIntVals.push_back(new TCLAP::MultiArg<int>(flag,name,desc,false,"int",my->cmd) );
    my->multiIntRes.push_back(&value);
}

void gsCmdLine::addReal( const std::string& flag, 
                         const std::string& name, 
                         const std::string& desc, 
                         real_t           & value)
{
    GISMO_ASSERT( !name.empty(), "The name (long form of the flag) must not be empty." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->realVals.push_back(new TCLAP::ValueArg<real_t>(flag,name,desc,false,value,"float",my->cmd) );
    my->realRes.push_back(&value);
}

void gsCmdLine::addMultiReal( const std::string  & flag, 
                              const std::string  & name, 
                              const std::string  & desc, 
                              std::vector<real_t>& value)
{
    GISMO_ASSERT( !name.empty(), "The name (long form of the flag) must not be empty." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->multiRealVals.push_back(new TCLAP::MultiArg<real_t>(flag,name,desc,false,"float",my->cmd) );
    my->multiRealRes.push_back(&value);
}

void gsCmdLine::addString( const std::string& flag, 
                           const std::string& name, 
                           const std::string& desc, 
                           std::string      & value)
{
    GISMO_ASSERT( !name.empty(), "The name (long form of the flag) must not be empty." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->stringVals.push_back(new TCLAP::ValueArg<std::string>(flag,name,desc,false,value,"string",my->cmd));
    my->strRes.push_back(&value);
}

void gsCmdLine::addMultiString( const std::string       & flag, 
                                const std::string       & name, 
                                const std::string       & desc, 
                                std::vector<std::string>& value)
{
    GISMO_ASSERT( !name.empty(), "The name (long form of the flag) must not be empty." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->multiStringVals.push_back(new TCLAP::MultiArg<std::string>(flag,name,desc,false,"string",my->cmd) );
    my->multiStrRes.push_back(&value);
}

void gsCmdLine::addSwitch( const std::string& flag, 
                           const std::string& name, 
                           const std::string& desc, 
                           bool             & value)
{
    GISMO_ASSERT( !name.empty(), "The name (long form of the flag) must not be empty." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );
    my->switches.push_back(new TCLAP::SwitchArg(flag,name,desc,my->cmd) );
    my->swRes.push_back(&value);
}

void gsCmdLine::addPlainString( const std::string& name, 
                                const std::string& desc, 
                                std::string & value)
{
    GISMO_ASSERT( !name.empty(), "The name (long form of the flag) must not be empty." );
    GISMO_ASSERT( !my->didParseCmdLine, "Variables must not be registered after calling gsCmdLine::getValues." );

    GISMO_ENSURE( !my->plainString, "Plain string already added." );

    my->plainString =
        new TCLAP::UnlabeledValueArg<std::string>(name,desc,false,value,"string",my->cmd);

    my->pstrRes = &value;
}

bool gsCmdLine::getValues(int argc, char *argv[])
{
    GISMO_ASSERT( !my->didParseCmdLine, "gsCmdLine::getValues must not be called twice." );
#ifndef NDEBUG
    my->didParseCmdLine = true;
#endif

    try
    {
        my->cmd.parse(argc,argv);

        for( std::size_t i=0; i!=my->intVals.size(); ++i)
            *my->intRes[i] = my->intVals[i]->getValue();

        for( std::size_t i=0; i!=my->multiIntVals.size(); ++i)
            *my->multiIntRes[i] = my->multiIntVals[i]->getValue();

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
        GISMO_UNUSED(e);
        return false; 
    }

    return true;
}


#define ADD_OPTION_LIST_ENTRY(type,iterate,addFct)                                        \
{                                                                                         \
    std::string nm = iterate->getName() + ".";                                            \
    const std::vector<type>& vec = iterate->getValue();                                   \
    std::size_t sz = vec.size();                                                          \
    for ( std::size_t j=0; j<sz; ++j )                                                    \
        result.addFct( nm+util::to_string(j), iterate->getDescription(), vec[j] );        \
    result.addInt( nm+"Size", iterate->getDescription(), sz );                            \
}

gsOptionList gsCmdLine::getOptionList()
{
    GISMO_ASSERT( my->didParseCmdLine, "gsCmdLine::getOptionList can be called only after gsCmdLine::getValues." );

    gsOptionList result;
    for( std::size_t i=0; i!=my->intVals.size(); ++i)
        result.addInt( my->intVals[i]->getName(), my->intVals[i]->getDescription(), my->intVals[i]->getValue() );
    for( std::size_t i=0; i!=my->realVals.size(); ++i)
        result.addReal( my->realVals[i]->getName(), my->realVals[i]->getDescription(), my->realVals[i]->getValue() );
    for( std::size_t i=0; i!=my->stringVals.size(); ++i)
        result.addString( my->stringVals[i]->getName(), my->stringVals[i]->getDescription(), my->stringVals[i]->getValue() );
    for( std::size_t i=0; i!=my->switches.size(); ++i)
        result.addSwitch( my->switches[i]->getName(), my->switches[i]->getDescription(), my->switches[i]->getValue() );
    for( std::size_t i=0; i!=my->multiIntVals.size(); ++i)
        ADD_OPTION_LIST_ENTRY(index_t,my->multiIntVals[i],addInt)
    for( std::size_t i=0; i!=my->multiRealVals.size(); ++i)
        ADD_OPTION_LIST_ENTRY(real_t,my->multiRealVals[i],addReal)
    for( std::size_t i=0; i!=my->multiStringVals.size(); ++i)
        ADD_OPTION_LIST_ENTRY(std::string,my->multiStringVals[i],addString)
    if ( my->plainString )
        result.addString( my->plainString->getName(), my->plainString->getDescription(), my->plainString->getValue() );
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


void gsCmdLinePrivate::GismoCmdOut::version(TCLAP::CmdLineInterface& c)
{
    GISMO_UNUSED(c);
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
    gsInfo << "RICAM-Linz 2012 - 2017, http://gs.jku.at/gismo\n";
}

std::string & gsCmdLine::getMessage()
{
    return my->cmd.getMessage();
}

} //namespace gismo
