/** @file common/Suppress.h

    @brief common UnitTest++ macro definitions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):H. Weiner
**/

#ifndef COMMONS_SUPPRESSOUTPUT_H
#define COMMONS_SUPPRESSOUTPUT_H

namespace UnitTest {

void deactivate_output();

void reactivate_output();

} // namespace UnitTest

#endif // #ifndef COMMONS_SUPPRESSOUTPUT_H
