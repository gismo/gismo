/** @file gsG1OptionList.h
 *
    @brief Option list for the G1 Basis

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/

#pragma once

#include <gismo.h>


namespace gismo
{

struct gluingData
{
    enum strategy
    {
        approximate  = 0,
        l2projection = 1, // global L2-projection
        local = 2 // local L2-projection
    };
};

struct g1BasisEdge
{
    enum strategy
    {
        l2projection = 0, // global L2-projection
        local = 1 // local L2-projection
    };
};

struct g1BasisVertex
{
    enum strategy
    {
        local = 0, // local transversal vector
        global = 1 // global transversal vetor
    };
};

struct user
{
    enum name
    {
        pascal = 0,
        andrea = 1,
    };

};

class gsG1OptionList
{

public:
    gsG1OptionList()
    {

    }

    gsG1OptionList(gsOptionList & list)
    {
        optionList = list;
    }

    index_t getInt(const std::string label) { return optionList.getInt(label); };
    bool getSwitch(const std::string label) { return optionList.getSwitch(label); };
    real_t getReal(const std::string label) { return optionList.getReal(label); };
    std::string getString(const std::string label) { return optionList.getString(label); };

    void addInt(const std::string & label, const std::string & desc, const index_t & value) { optionList.addInt(label, desc, value);  };
    void addSwitch(const std::string & label, const std::string & desc, const bool & value) { optionList.addSwitch(label, desc, value);  };
    void addReal(const std::string & label, const std::string & desc, const real_t & value) { optionList.addReal(label, desc, value);  };
    void addString(const std::string & label, const std::string & desc, const std::string & value) { optionList.addString(label, desc, value);  };

    void setInt(const std::string & label, const index_t & value) { optionList.setInt(label, value);  };
    void setSwitch(const std::string & label, const bool & value) { optionList.setSwitch(label, value);  };
    void setReal(const std::string & label, const real_t & value) { optionList.setReal(label, value);  };
    void setString(const std::string & label, const std::string & value) { optionList.setString(label, value);  };

protected:

    gsOptionList optionList;

};

} // namespace gismo