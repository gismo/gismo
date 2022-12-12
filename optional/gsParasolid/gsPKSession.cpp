/** @file gsPKSession.cpp

    @brief Manages starting and stopping Parasolid session

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsParasolid/gsPKSession.h>

#include <gsParasolid/gsFrustrum.h>

#include <string>


namespace gismo {

namespace extensions {


bool gsPKSession::start()
{
    // Set the schema path
    if (!getenv("P_SCHEMA")) 
    {
        std::string env_schema("P_SCHEMA=" GISMO_P_SCHEMA);
        putenv( const_cast<char*>(env_schema.c_str()) );
        //gsWarn<<"Setting P_SCHEMA to"<<GISMO_P_SCHEMA<<"\n";
    }

    PK_SESSION_start_o_t start_options;

    register_frustrum();

    // Set initialization. options to default values
    PK_SESSION_start_o_m(start_options);

    // Start session
    PK_SESSION_start(&start_options);

    return true;
}


bool gsPKSession::stop()
{
    // Stop Parasolid session
    PK_SESSION_stop();

    return true;
}

}//extensions

}//gismo

