/** @file gsParaviewDataSet.h

    @brief Provides a helper class to write Paraview (.vts) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Karampatzakis, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsMSplines/gsMappedBasis.h>   // Only to make linker happy
#include <gsCore/gsDofMapper.h>         // Only to make linker happy
#include <gsAssembler/gsExprHelper.h>  
#include <gsAssembler/gsExprEvaluator.h>

#include<fstream>

namespace gismo 
{

class gsParaviewDataSet // a collection of .vts files 
{

public:
    int m_numPatches;
private:
    index_t plotPrec;
    int part;
    std::string m_basename;
    gsExprHelper<real_t>::geometryMap * m_geoMap;
    gsExprEvaluator<real_t> * m_evaltr;

    
    // Will only be called by save() function, to make sure it is the last element in the file
    void outputGeometry();

    
public:
    gsParaviewDataSet(std::string basename, 
                      gsExprHelper<real_t>::geometryMap * geoMap,
                      gsExprEvaluator<real_t> * eval);

                       


    template<class E>
    void addField(const expr::_expr<E> & expr,
                  std::string label);

    // void addFields( expr... )
    // {
    //     // maybe I can use the variable names as the labels by
    //     // using the # "stringify" 
    //     // see: https://stackoverflow.com/questions/3386861/converting-a-variable-name-to-a-string-in-c
    // }

    std::vector<std::string> filenames();

    void save();


};
} // End namespace gismo
