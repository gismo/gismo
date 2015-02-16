
#include <gsIO/gsCmdLine.h>

#include <gsCore/gsMemory.h>

#include <tclap/CmdLine.cpp> // External file

namespace gismo
{

gsCmdLine::gsCmdLine( const std::string& message,	
                      const char delimiter,
                      bool helpAndVersion)  
    : Base(message,delimiter,GISMO_VERSION,helpAndVersion) 
{ 
    this->setOutput( &cmdout );
}

void gsCmdLine::addInt( const std::string& flag, 
                             const std::string& name, 
                             const std::string& desc, 
                             int & value)
{

    m_intVals.push_back(new gsArgVal<int>(flag,name,desc,false,value,"int",*this) );
    m_intRes.push_back(&value);
}

void gsCmdLine::addReal( const std::string & flag, 
                         const std::string& name, 
                             const std::string& desc, 
                             real_t & value)
{

    m_realVals.push_back(new gsArgVal<real_t>(flag,name,desc,false,value,"float",*this) );
    m_realRes.push_back(&value);
}

void gsCmdLine::addString( const std::string& flag, 
                           const std::string& name, 
                           const std::string& desc, 
                           std::string & value)
{

    m_stringVals.push_back(
    new gsArgVal<std::string>(flag,name,desc,false,value,"string",*this) );
    m_strRes.push_back(&value);
}

void gsCmdLine::addSwitch( const std::string& name, 
                           const std::string& desc, 
                           bool & value)
{
    m_switches.push_back(new gsArgSwitch("",name,desc,*this) );
    m_swRes.push_back(&value);
}


bool gsCmdLine::getValues(int argc, char *argv[])
{
    try 
    {
        this->parse(argc,argv);
        
        for( std::size_t i=0; i!=m_intVals.size(); ++i)
            *m_intRes[i] = m_intVals[i]->getValue();

        for( std::size_t i=0; i!=m_realVals.size(); ++i)
            *m_realRes[i] = m_realVals[i]->getValue();

        for( std::size_t i=0; i!=m_stringVals.size(); ++i)
            *m_strRes[i] = m_stringVals[i]->getValue();

        for( std::size_t i=0; i!=m_switches.size(); ++i)
            *m_swRes[i] = m_switches[i]->getValue();

    }
    catch ( gsArgException& e )
    { gsWarn << "Error: " << e.error() << " " << e.argId() << "\n"; return false; }
    return true;
}


gsCmdLine::~gsCmdLine() 
{ 
    freeAll( m_intVals   );
    freeAll( m_realVals  );
    freeAll( m_stringVals);
    freeAll( m_switches  );
}

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
    gsInfo << "   Copyright (C) JKU-RICAM-Linz, 2012 - 2015\n";
}


} //namespace gismo
