/** @file common/CheckMacros.h

    @brief common UnitTest++ macro definitions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris,  H. Weiner
**/

#ifndef COMMON_CHECKSMACROS_H
#define COMMON_CHECKSMACROS_H

#include "../UnitTest++/Config.h"
#include "../UnitTest++/TestResults.h"
#include "../UnitTest++/MemoryOutStream.h"
#include "CheckMatrix.h"
#include "CheckFiles.h"

namespace UnitTest {

#define CHECK_MATRIX_CLOSE(expected, actual, tolerance) \
	UNITTEST_MULTILINE_MACRO_BEGIN \
        UT_TRY \
		({ \
            UnitTest::CheckMatrixClose(*UnitTest::CurrentTest::Results(), expected, actual, \
              tolerance, UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__)); \
        }) \
 		UT_CATCH (std::exception, e, \
		{ \
			UnitTest::MemoryOutStream message; \
			message << "Unhandled exception (" << e.what() << ") in CHECK_MATRIX_CLOSE(" #expected ", " #actual ")"; \
			UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
				message.GetText()); \
		}) \
        UT_CATCH_ALL \
		({ \
            UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
                    "Unhandled exception in CHECK_MATRIX_CLOSE(" #expected ", " #actual ")"); \
        }) \
	UNITTEST_MULTILINE_MACRO_END


// all the colum sums are equal to one (1)
#define CHECK_MATRIX_PARTITION_OF_UNIT_CLOSE(actual, tolerance) \
	UNITTEST_MULTILINE_MACRO_BEGIN \
        UT_TRY \
		({ \
            UnitTest::CheckMatrixPartitionOfUnitClose(*UnitTest::CurrentTest::Results(), actual, \
              tolerance, UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__)); \
        }) \
 		UT_CATCH (std::exception, e, \
		{ \
			UnitTest::MemoryOutStream message; \
			message << "Unhandled exception (" << e.what() << ") in CHECK_MATRIX_PARTITION_OF_UNIT_CLOSE(" #actual ")"; \
			UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
				message.GetText()); \
		}) \
        UT_CATCH_ALL \
		({ \
            UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
                    "Unhandled exception in CHECK_MATRIX_PARTITION_OF_UNIT_CLOSE(" #actual ")"); \
        }) \
	UNITTEST_MULTILINE_MACRO_END

#define CHECK_FILES_EQUAL(actual, expected) \
	UNITTEST_MULTILINE_MACRO_BEGIN \
        UT_TRY \
		({ \
			UnitTest::CheckFilesEqual(*UnitTest::CurrentTest::Results(), actual, \
			expected, UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__)); \
        }) \
 		UT_CATCH (std::exception, e, \
		{ \
			UnitTest::MemoryOutStream message; \
			message << "Unhandled exception (" << e.what() << ") in CHECK_FILES_EQUAL(" #actual ")"; \
			UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
				message.GetText()); \
		}) \
        UT_CATCH_ALL \
		({ \
            UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
                    "Unhandled exception in CHECK_FILES_EQUAL(" #actual ")"); \
        }) \
	UNITTEST_MULTILINE_MACRO_END

#define CHECK_FILES_EQUAL_IGNORE(actual, expected, ignoreLines) \
	UNITTEST_MULTILINE_MACRO_BEGIN \
        UT_TRY \
		({ \
			UnitTest::CheckFilesEqual(*UnitTest::CurrentTest::Results(), actual, \
			expected, ignoreLines, UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__)); \
        }) \
 		UT_CATCH (std::exception, e, \
		{ \
			UnitTest::MemoryOutStream message; \
			message << "Unhandled exception (" << e.what() << ") in CHECK_FILES_EQUAL(" #actual ")"; \
			UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
				message.GetText()); \
		}) \
        UT_CATCH_ALL \
		({ \
            UnitTest::CurrentTest::Results()->OnTestFailure(UnitTest::TestDetails(*UnitTest::CurrentTest::Details(), __LINE__), \
                    "Unhandled exception in CHECK_FILES_EQUAL(" #actual ")"); \
        }) \
	UNITTEST_MULTILINE_MACRO_END

#ifndef NDEBUG
#define CHECK_THROW_IN_DEBUG(condition, expectedException) \
		CHECK_THROW(condition, expectedException)
#else
#define CHECK_THROW_IN_DEBUG(condition, expectedException)
#endif

} // end namespace UnitTest

#endif // COMMON_CHECKSMACROS_H

