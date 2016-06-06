/** @file gsCmdLine.cpp
    
    @brief Provides implemementation of input command line arguments.
    
    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <vector> // temp

#include <gsCore/gsConfig.h>
#include <gsCore/gsDebug.h>
#include <gsCore/gsMemory.h>

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

    // Stores integer arguments
    std::vector<TCLAP::ValueArg<int>*    > intVals ;
    std::vector<int*>               intRes ;

    // Stores real_t arguments
    std::vector<TCLAP::ValueArg<real_t>*> realVals;
    std::vector<real_t*>           realRes ;

    // Stores string arguments
    std::vector<TCLAP::ValueArg<std::string>*> stringVals;
    std::vector<std::string*>           strRes ;

    // Stores switch arguments
    std::vector<TCLAP::SwitchArg*      > switches;
    std::vector<bool*>              swRes ;

    // Stores plain string argument
    TCLAP::UnlabeledValueArg<std::string> * plainString;
    std::string *                pstrRes;

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
{ }

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
                             int & value)
{
    my->intVals.push_back(new TCLAP::ValueArg<int>(flag,name,desc,false,value,"int",my->cmd) );
    my->intRes.push_back(&value);
}

void gsCmdLine::addReal( const std::string & flag, 
                         const std::string& name, 
                             const std::string& desc, 
                             real_t & value)
{

    my->realVals.push_back(new TCLAP::ValueArg<real_t>(flag,name,desc,false,value,"float",my->cmd) );
    my->realRes.push_back(&value);
}

void gsCmdLine::addString( const std::string& flag, 
                           const std::string& name, 
                           const std::string& desc, 
                           std::string & value)
{

    my->stringVals.push_back(
    new TCLAP::ValueArg<std::string>(flag,name,desc,false,value,"string",my->cmd) );
    my->strRes.push_back(&value);
}

void gsCmdLine::addSwitch( const std::string& name, 
                           const std::string& desc, 
                           bool & value)
{
    my->switches.push_back(new TCLAP::SwitchArg("",name,desc,my->cmd) );
    my->swRes.push_back(&value);
}

void gsCmdLine::addSwitch( const std::string& flag, 
                           const std::string& name, 
                           const std::string& desc, 
                           bool & value)
{
    my->switches.push_back(new TCLAP::SwitchArg(flag,name,desc,my->cmd) );
    my->swRes.push_back(&value);
}

void gsCmdLine::addPlainString( const std::string& name, 
                                const std::string& desc, 
                                std::string & value)
{
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
    gsInfo << "             version "<< GISMO_VERSION <<"\n";
    gsInfo << "   Copyright (C) JKU-RICAM-Linz, 2012 - 2016\n";
}


} //namespace gismo
