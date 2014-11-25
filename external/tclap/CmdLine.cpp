///////////////////////////////////////////////////////////////////////////////
//Begin CmdLine.cpp
///////////////////////////////////////////////////////////////////////////////

#include <tclap/CmdLine.h>

namespace TCLAP {

inline 
CmdLine::CmdLine(const std::string& m,
                 const char delim,
                 const std::string& v,
		 bool help )
    :
  _argList(std::list<Arg*>()),
  _progName("not_set_yet"),
  _message(m),
  _version(v),
  _numRequired(0),
  _delimiter(delim),
  _xorHandler(XorHandler()),
  _argDeleteOnExitList(std::list<Arg*>()),
  _visitorDeleteOnExitList(std::list<Visitor*>()),
  _output(0),
  _handleExceptions(true),
  _userSetOutput(false),
  _helpAndVersion(help)
{
	_constructor();
}

inline CmdLine::~CmdLine()
{
	ClearContainer(_argDeleteOnExitList);
	ClearContainer(_visitorDeleteOnExitList);

	if ( !_userSetOutput ) {
		delete _output;
		_output = 0;
	}
}

inline void CmdLine::_constructor()
{
	_output = new StdOutput;

	Arg::setDelimiter( _delimiter );

	Visitor* v;

	if ( _helpAndVersion )
	{
		v = new HelpVisitor( this, &_output );
		SwitchArg* help = new SwitchArg("h","help",
		                      "Displays usage information and exits.",
		                      false, v);
		add( help );
		deleteOnExit(help);
		deleteOnExit(v);

		v = new VersionVisitor( this, &_output );
		SwitchArg* vers = new SwitchArg("","version",
		                      "Displays version information and exits.",
		                      false, v);
		add( vers );
		deleteOnExit(vers);
		deleteOnExit(v);
	}

	v = new IgnoreRestVisitor();
	SwitchArg* ignore  = new SwitchArg(Arg::flagStartString(),
	          Arg::ignoreNameString(),
	          "Ignores the rest of the labeled arguments following this flag.",
	          false, v);
	add( ignore );
	deleteOnExit(ignore);
	deleteOnExit(v);
}

inline void CmdLine::xorAdd( std::vector<Arg*>& ors )
{
	_xorHandler.add( ors );

	for (ArgVectorIterator it = ors.begin(); it != ors.end(); it++)
	{
		(*it)->forceRequired();
		(*it)->setRequireLabel( "OR required" );
		add( *it );
	}
}

inline void CmdLine::xorAdd( Arg& a, Arg& b )
{
	std::vector<Arg*> ors;
	ors.push_back( &a );
	ors.push_back( &b );
	xorAdd( ors );
}

inline void CmdLine::add( Arg& a )
{
	add( &a );
}

inline void CmdLine::add( Arg* a )
{
	for( ArgListIterator it = _argList.begin(); it != _argList.end(); it++ )
		if ( *a == *(*it) )
			throw( SpecificationException(
			        "Argument with same flag/name already exists!",
			        a->longID() ) );

	a->addToList( _argList );

	if ( a->isRequired() )
		_numRequired++;
}


inline void CmdLine::parse(int argc, const char * const * argv)
{
		// this step is necessary so that we have easy access to
		// mutable strings.
		std::vector<std::string> args;
		for (int i = 0; i < argc; i++)
			args.push_back(argv[i]);

		parse(args);
}

inline void CmdLine::parse(std::vector<std::string>& args)
{
	bool shouldExit = false;
	int estat = 0;

	try {
		_progName = args.front();
		args.erase(args.begin());

		int requiredCount = 0;

		for (int i = 0; static_cast<unsigned int>(i) < args.size(); i++) 
		{
			bool matched = false;
			for (ArgListIterator it = _argList.begin();
			     it != _argList.end(); it++) {
				if ( (*it)->processArg( &i, args ) )
				{
					requiredCount += _xorHandler.check( *it );
					matched = true;
					break;
				}
			}

			// checks to see if the argument is an empty combined
			// switch and if so, then we've actually matched it
			if ( !matched && _emptyCombined( args[i] ) )
				matched = true;

			if ( !matched && !Arg::ignoreRest() )
				throw(CmdLineParseException("Couldn't find match "
				                            "for argument",
				                            args[i]));
		}

		if ( requiredCount < _numRequired )
			missingArgsException();

		if ( requiredCount > _numRequired )
			throw(CmdLineParseException("Too many arguments!"));

	} catch ( ArgException& e ) {
		// If we're not handling the exceptions, rethrow.
		if ( !_handleExceptions) {
			throw;
		}

		try {
			_output->failure(*this,e);
		} catch ( ExitException &ee ) {
			estat = ee.getExitStatus();
			shouldExit = true;
		}
	} catch (ExitException &ee) {
		// If we're not handling the exceptions, rethrow.
		if ( !_handleExceptions) {
			throw;
		}

		estat = ee.getExitStatus();
		shouldExit = true;
	}

	if (shouldExit)
		exit(estat);
}

inline bool CmdLine::_emptyCombined(const std::string& s)
{
	if ( s.length() > 0 && s[0] != Arg::flagStartChar() )
		return false;

	for ( int i = 1; static_cast<unsigned int>(i) < s.length(); i++ )
		if ( s[i] != Arg::blankChar() )
			return false;

	return true;
}

inline void CmdLine::missingArgsException()
{
		int count = 0;

		std::string missingArgList;
		for (ArgListIterator it = _argList.begin(); it != _argList.end(); it++)
		{
			if ( (*it)->isRequired() && !(*it)->isSet() )
			{
				missingArgList += (*it)->getName();
				missingArgList += ", ";
				count++;
			}
		}
		missingArgList = missingArgList.substr(0,missingArgList.length()-2);

		std::string msg;
		if ( count > 1 )
			msg = "Required arguments missing: ";
		else
			msg = "Required argument missing: ";

		msg += missingArgList;

		throw(CmdLineParseException(msg));
}

inline void CmdLine::deleteOnExit(Arg* ptr)
{
	_argDeleteOnExitList.push_back(ptr);
}

inline void CmdLine::deleteOnExit(Visitor* ptr)
{
	_visitorDeleteOnExitList.push_back(ptr);
}

inline CmdLineOutput* CmdLine::getOutput()
{
	return _output;
}

inline void CmdLine::setOutput(CmdLineOutput* co)
{
	if ( !_userSetOutput )
		delete _output;
	_userSetOutput = true;
	_output = co;
}

inline std::string& CmdLine::getVersion()
{
	return _version;
}

inline std::string& CmdLine::getProgramName()
{
	return _progName;
}

inline std::list<Arg*>& CmdLine::getArgList()
{
	return _argList;
}

inline XorHandler& CmdLine::getXorHandler()
{
	return _xorHandler;
}

inline char CmdLine::getDelimiter()
{
	return _delimiter;
}

inline std::string& CmdLine::getMessage()
{
	return _message;
}

inline bool CmdLine::hasHelpAndVersion()
{
	return _helpAndVersion;
}

inline void CmdLine::setExceptionHandling(const bool state)
{
	_handleExceptions = state;
}

inline bool CmdLine::getExceptionHandling() const
{
	return _handleExceptions;
}

inline void CmdLine::reset()
{
	for( ArgListIterator it = _argList.begin(); it != _argList.end(); it++ )
		(*it)->reset();
	
	_progName.clear();
}

///////////////////////////////////////////////////////////////////////////////
//End CmdLine.cpp
///////////////////////////////////////////////////////////////////////////////

} //namespace TCLAP
