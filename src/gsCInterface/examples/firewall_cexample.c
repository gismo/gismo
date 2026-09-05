/** @file firewall_cexample.c

    @brief Verifies that C++ exceptions do not cross the C boundary:
    a failing call returns a sentinel and sets gsCLastError().

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <stdio.h>
#include <string.h>

#include <gsCInterface/Cgismo.h>

int main(void)
{
    /* Setting an option that was never added throws a GISMO error
       inside the library; the firewall must convert that into a plain
       return plus a per-thread message - not a crash. */
    gsCOptionList * ol = gsOptionList_create();
    gsOptionList_setInt(ol, "DoesNotExist", 1);
    if (strlen(gsCLastError()) == 0)
    {
        printf("FAIL: expected a non-empty gsCLastError()\n");
        gsOptionList_delete(ol);
        return 1;
    }
    printf("caught: '%.60s...'\n", gsCLastError());
    gsOptionList_delete(ol);

    /* A successful call clears the error state. */
    gsCKnotVector * kv;
    double knots[] = {0.,0.,1.,1.};
    kv = gsKnotVector_create(knots, 4);
    if (kv == NULL || strlen(gsCLastError()) != 0)
    {
        printf("FAIL: valid call left error state\n");
        return 1;
    }
    gsKnotVector_delete(kv);

    printf("firewall OK\n");
    return 0;
}
