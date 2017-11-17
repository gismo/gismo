/** @file common/CheckFiles.cpp

    @brief common UnitTest++ macro definitions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):H. Weiner
**/

#include <algorithm>
#include <iterator>
#include <string>
#include <fstream>
#include <vector>

#include "CheckFiles.h"
#include "../UnitTest++/MemoryOutStream.h"

namespace UnitTest {

bool contains(std::vector<int>& vector, int value) {
	bool result = std::find(vector.begin(), vector.end(), value) != vector.end();
	return result;
}

void CheckFilesEqual(TestResults& results, std::string const& actual,
		std::string const& expected, TestDetails const& details)
{
	std::vector<int> ignoreLines(0);
	CheckFilesEqual(results, actual, expected, ignoreLines, details);
}

void CheckFilesEqual(TestResults& results, std::string const& actual,
		std::string const& expected, std::vector<int> ignoreLines,
		TestDetails const& details)
{
	std::ifstream file1(expected.c_str());
	std::ifstream file2(actual.c_str());
	if(!file1.is_open()) {
		UnitTest::MemoryOutStream stream;
		stream << "could not open expected file with filename ['" << expected << "']";
		results.OnTestFailure(details, stream.GetText());
		return;
	}
	if (!file2.is_open())
	{
		UnitTest::MemoryOutStream stream;
		stream << "could not open actual file with filename ['" << actual << "']";
		results.OnTestFailure(details, stream.GetText());
		file1.close();
		return;
	}
	bool result = true;
	std::string line1="";
	std::string line2="";
	unsigned lineNo = 0;
	while(!file1.eof() && !(file1.eof()) && result)
	{
		getline(file1,line1);
		getline(file2,line2);
		if (!contains(ignoreLines, lineNo)) {
			result = (line1 == line2);
		}
		lineNo++;
	}
	if (!result) {
		UnitTest::MemoryOutStream stream;
		stream << "expected file ['" << expected << "'] and actual file ['";
		stream << actual << "']  do not have the  same content (lineNo='";
		stream << lineNo << "' is different)!";
		results.OnTestFailure(details, stream.GetText());
	}
	file1.close();
	file2.close();
}

} // namespace UnitTest
