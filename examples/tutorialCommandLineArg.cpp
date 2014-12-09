/** @file tutorialCommandLineArg.cpp

    @brief Tutorial on how to use command line parser in G+Smo.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

// TO DO: Check if this is appropriate tutorial...

#include <iostream>
#include <string>
#include <gismo.h>

using namespace gismo;

int main(int argc, char* argv[])
{

    std::string plainString("");
    std::string string("");
    int number = 0;
    bool boolean = false;

    try
    {
        gsCmdLine cmd("Tutorial Command Line Arguments");

        // =================================================================
        // Defining command line arguments
        // =================================================================


        // -----------------------------------------------------------------
        // plain argument
        // -----------------------------------------------------------------
        
        std::string name = "plain";
        std::string desc =  "Description of the plain command line argument.";
        bool req = false; // whether the argument is required
        std::string value = "default_plain_value"; 
        std::string typeDesc = "string"; // type description
        
        gsArgValPlain<std::string> plainArg(name, desc, req, value, typeDesc, cmd);
        
        // -----------------------------------------------------------------
        // argument value
        // -----------------------------------------------------------------
        
        // example 1

        std::string flag = "s";
        name = "string";
        desc = "Description of string command line argument.";
        req = false;
        value = "default_string";
        typeDesc = "string";
        
        gsArgVal<std::string> strArg(flag, name, desc, req, value, typeDesc, cmd);
        
        // example 2
        
        gsArgVal<int> intArg("n", "num", "Description of int command line argument", 
                             false, 0, "int", cmd);
        
        // -----------------------------------------------------------------
        // argument switch
        // -----------------------------------------------------------------
        
        flag = "b";
        name = "bool";
        desc = "Description of the switch argument.";
        gsArgSwitch switchArg(flag, name, desc, cmd);
        
        
        // =================================================================
        // parsing command line arguments
        // =================================================================
        
        cmd.parse(argc, argv);
        plainString = plainArg.getValue();
        string = strArg.getValue();
        number = intArg.getValue();
        boolean = switchArg.getValue();
    }
    catch (gsArgException& e)
    {
        std::cout << "Error: " << e.error() << " " << e.argId() << std::endl;
        return -1;
    }
    
    std::cout << "Printing command line arguments:\n\n\n"
              << "Plain string: " << plainString << "\n\n"
              << "String:       " << string << "\n\n"
              << "Number:       " << number << "\n\n"
              << "Switch:       " << boolean << "\n" << std::endl;

    return 0;
}







