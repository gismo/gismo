
#include <gsParasolid/gsPKSession.h>

#include <string>

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
