/** @file common/CheckFiles.h

    @brief common UnitTest++ macro definitions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):H. Weiner
**/

#ifndef COMMONS_CHECKFILES_H
#define COMMONS_CHECKFILES_H

#include <string>

#include "../UnitTest++/TestResults.h"

namespace UnitTest {

void CheckFilesEqual(TestResults& results, const std::string & actual,
		const std::string & expected, TestDetails const& details);

void CheckFilesEqual(TestResults& results, const std::string & actual,
		const std::string & expected, std::vector<int> ignoreLines,
		TestDetails const& details);

} // namespace UnitTest

#endif // COMMONS_CHECKFILES_H
