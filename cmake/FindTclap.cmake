######################################################################
## FindTclap.cmake
## This file is part of the G+Smo library. 
##
## Author: Harald Weiner
######################################################################


find_path(Tclap_DIR
			NAMES 
				Arg.h ArgException.h ArgTraits.h
				CmdLine.h CmdLineInterface.h CmdLineOutput.h
				Constraint.h DocBookOutput.h HelpVisitor.h
				IgnoreRestVisitor.h MultiArg.h MultiSwitchArg.h
				OptionalUnlabeledTracker.h StandardTraits.h
				StdOutput.h SwitchArg.h UnlabeledMultiArg.h
				UnlabeledValueArg.h ValueArg.h 
				ValuesConstraint.h VersionVisitor.h
				Visitor.h XorHandler.h ZshCompletionOutput.h
			HINTS /usr/include/tclap
			PATHS
					/usr/include/tclap/
					/usr/local/include/tclap/
					/opt/tclap/
					${Tclap_DIR}/
)

if(Tclap_DIR)
	set(Tclap_FOUND TRUE)
	if(NOT Tclap_FIND_QUIETLY)
		message("Tclap_FOUND")
	endif(NOT Tclap_FIND_QUIETLY)
else(Tclap_DIR)
	set(Tclap_DIR "Tclap_DIR-NOTFOUND " CACHE PATH "The path to the Tclap header files")
	if(Tclap_FIND_REQUIRED)
      message("Could not find Tclap header files.")
      message(FATAL_ERROR "Set the variable Tclap_DIR and try again.")
   endif()
endif(Tclap_DIR)
