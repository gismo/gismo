/** @file modules_cexample.c

    @brief Module availability from C: gsC_hasModule / gsC_modules.

    Doubles as the test of the feature-detection functions: exits with
    a nonzero code if a compiled-in module is reported missing or an
    unknown one is reported present.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gsCInterface/Cgismo.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[])
{
    int errors = 0;
    printf("modules compiled into gismo: %s\n", gsC_modules());

    if (!gsC_hasModule("gsNurbs"))   { printf("gsNurbs missing\n");   ++errors; }
    if (!gsC_hasModule("gsCore"))    { printf("gsCore missing\n");    ++errors; }
    if ( gsC_hasModule("gsNoSuchModule")) { printf("gsNoSuchModule reported present\n"); ++errors; }
    if ( gsC_hasModule("gsNurb"))    { printf("prefix match reported present\n"); ++errors; }
    if ( gsC_hasModule(NULL))        { printf("NULL reported present\n"); ++errors; }
    if (0 != strcmp(gsCLastError(), "")) { printf("unexpected error: %s\n", gsCLastError()); ++errors; }

    printf("gsC_hasModule(\"gsNurbs\") = %d, gsC_hasModule(\"gsNoSuchModule\") = %d -> %s\n",
           gsC_hasModule("gsNurbs"), gsC_hasModule("gsNoSuchModule"), errors ? "FAILED" : "ok");
    return errors;
}
