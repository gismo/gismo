/** @file main.cpp

    @brief Main file for unittests

    This file includes the magic necessary in order to get your unit tests
    that you create with UnitTest++ to automatically run. There should
    never be a reason to modify this file.
    For a reference on creating tests, see gsTutorial.cpp

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    To run with gdb:
    gdb --args ./unittests -R [testname]

    (gdb) catch throw
    (gdb) run
    ...
    (gdb) bt
    ...     see backtrack to find error
    (gdb) quit

    Author(s): A. Mantzaflaris,  H. Weiner, J. Vogl
 **/

#include "gismo_unittest.h"
#include "TestReporterStdout.h"

// Tolerance for approximate comparisons
const real_t EPSILON = std::pow(10.0, - REAL_DIG * 0.75);

// Macros for excluding some test name
#define EXCLUDE_SUITE(name) if ( !strcmp(testCase->m_details.suiteName,#name) ) return false


/// Selects tests by matching their names to the input prefix given by
/// command line argument
namespace gismo {

template<class T>
struct gsUTSelector
{
    bool operator()(const UnitTest::Test * const testCase) const
    { return true; }
};

#ifdef GISMO_WITH_GMP
template<>
struct gsUTSelector<mpq_class>
{
    bool operator()(const UnitTest::Test * const testCase) const
    {
        EXCLUDE_SUITE(gsPoissonSolver_test);
        EXCLUDE_SUITE(gsIterativeSolvers_test);
        EXCLUDE_SUITE(gsPreconditioner_test);
        return true;
    }
};
#endif

class gsUTSelectorCmd
{
private:
    int            m_argc;
    char        ** m_argv;

public:
    gsUTSelectorCmd(int argc, char* argv[])
    : m_argc(argc), m_argv(argv)
    { }

    bool operator()(const UnitTest::Test * const testCase) const
    {
        bool toRun = false;
        for (int i=1; i<m_argc; ++i)
        {
            const size_t n = strlen(m_argv[i]);
            toRun |= !strncmp(testCase->m_details.suiteName, m_argv[i], n);// prefix match
            toRun |= !strncmp(testCase->m_details.testName , m_argv[i], n);// prefix match
            toRun |= gsFileManager::pathEqual(testCase->m_details.filename, m_argv[i]);// exact match up to path sep.
        }
        return toRun;
    }
};
} // namespace gismo

int main(int argc, char* argv[])
{
    gsCmdLine::printVersion();

    UnitTest::TestReporterStdout reporter;
    UnitTest::TestRunner runner(reporter);

    if (argc > 1)
    {
        gsUTSelectorCmd sel(argc,argv);
        int result = runner.RunTestsIf(UnitTest::Test::GetTestList(), NULL, sel, 0);
        if ( 0 == runner.GetTestResults()->GetTotalTestCount() )
        {
            gsInfo << "Did not find any matching test.\n";
            return EXIT_FAILURE;
        }
        return result;
    }
    else
    {
        
        gsUTSelector<real_t> sel;
        int result = runner.RunTestsIf(UnitTest::Test::GetTestList(), NULL, sel, 0);
        return result;
        //return UnitTest::RunAllTests();
    }
}
