
#include <gsIO/gsCmdLine.h>

#include <tclap/CmdLine.cpp> // External file

namespace gismo
{

gsCmdLine::gsCmdLine( const std::string& message,	const char delimiter,
                         const std::string& version,
                         bool helpAndVersion)  
    : Base(message,delimiter,version,helpAndVersion) 
{ 
    this->setOutput( &cmdout );
}

gsCmdLine::~gsCmdLine() { }

void gsCmdLine::GismoCmdOut::failure(TCLAP::CmdLineInterface& c, TCLAP::ArgException& e)
{
    gsInfo << e.what() << "\n";
    gsInfo <<"\n USAGE: \n";
    _longUsage( c, gsInfo );
    throw;
}

void gsCmdLine::GismoCmdOut::usage(TCLAP::CmdLineInterface& c)
{
    std::string head= c.getProgramName();
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


void gsCmdLine::GismoCmdOut::version(TCLAP::CmdLineInterface& c)
{
    //gsInfo <<"Executable: "<< c.getProgramName() <<", part of \n\n";
    gsInfo << "\n";
    gsInfo << "                   G+Smo \n";
    gsInfo << "      Geometry plus Simulation modules\n";
    gsInfo << "             version "<< GISMO_VERSION <<"\n";
    gsInfo << "   Copyright (C) RICAM-Linz, 2012 - 2014\n";
}


} //namespace gismo
