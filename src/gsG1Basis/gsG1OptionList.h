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
        exact  = 0,
        global = 1, // global L2-projection
        local = 2 // local L2-projection
    };
};

struct g1BasisEdge
{
    enum strategy
    {
        global = 0, // global L2-projection
        local = 1 // local L2-projection
    };
};

struct g1BasisVertex
{
    enum strategy
    {
        global = 0, // global L2-projection
        local = 1 // local L2-projection
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

    void initialize(int argc, char *argv[]);

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

void gsG1OptionList::initialize(int argc, char *argv[])
{

    // Geometry data
    index_t loop = 1; // Number of refinement steps
    index_t geometry = 0; // Which geometry

    index_t P_geo = 0; // geometry degree elevate

    index_t numRefine = 4; // Initial refinement
    index_t R_geo = 1; // Regularity

    // For the spline space of the gluing data
    index_t p_gd = 1;
    index_t r_gd = 0;

    index_t threads = 12; // For parallel computing



    real_t threshold = 1e-5; // For computing the kernel
    real_t zero = 1e-12; // For setting the matrix for the kernel

    real_t lambda = 1e-12; // lambda value

    bool plot = false;
    bool latex = false;
    bool latex_plot = false;

    bool neumann_bdy = false;

    bool localGd = false;
    bool exactGd = false;
    bool localEdge = false;
    bool localVertex = false;
    bool isogeometric = false;
    bool h1projection = false;
    bool h2projection = false;
    bool h1projectionProof = false;

    bool info = false;

    gluingData::strategy gluingData_strategy = gluingData::global;
    g1BasisEdge::strategy g1BasisEdge_strategy = g1BasisEdge::global;
    g1BasisVertex::strategy g1BasisVertex_strategy = g1BasisVertex::global;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("k", "refine", "Number of refinement steps", numRefine);

    cmd.addInt("p", "p_tilde", "Polynomial degree for tilde{p}", p_gd);
    cmd.addInt("P", "P_geo", "Polynomial degree for geometry", P_geo);
    cmd.addInt("r", "r_tilde", "Regularity for tilde{r}", r_gd);
    cmd.addInt("R", "R_geo", "Regularity for geometry", R_geo);

    cmd.addInt("g", "geometry", "Geometry", geometry);
    cmd.addInt("t", "threads", "Threads", threads);
    cmd.addInt( "l", "loop", "The number of refinement steps", loop);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addSwitch( "localGd", "To compute the gluing data with local support", localGd );
    cmd.addSwitch( "exactGd", "To compute the gluing data exact", exactGd );
    cmd.addSwitch( "localEdge", "To compute the G1 edge basis functions with local support", localEdge );
    cmd.addSwitch( "localVertex", "To compute the G1 vertex basis functions with the average dd_ik", localVertex );
    cmd.addSwitch("latex","Print the rate and error latex-ready",latex);
    cmd.addSwitch("latex_plot","Print the rate and error latex-ready",latex_plot);
    cmd.addSwitch("neumann","Compute the biharmonic with neumann bdy",neumann_bdy);
    cmd.addSwitch( "isogeometric", "Project the basis in isogeometric concept", isogeometric );
    cmd.addSwitch( "h1projection", "Project the basis in H1 norm", h1projection );
    cmd.addSwitch( "h2projection", "Project the basis in H2 norm", h2projection );
    cmd.addSwitch( "h1projectionProof", "Project the basis in H1 norm", h1projectionProof );
    cmd.addSwitch( "info", "Print information", info );
    cmd.addReal("e","threshold", "The threshold for computing the kernel", threshold);
    cmd.addReal("z","zero", "When the value should be set to zero", zero);
    try { cmd.getValues(argc,argv); } catch (int rv) {  }

    optionList.addInt("loop","Loop", loop);
    optionList.addInt("geometry","Geometry", geometry);
    optionList.addInt("threads","Threads", threads);
    optionList.addInt("numRefine","Number of refinement", numRefine);

    optionList.addInt("p_tilde","Grad",p_gd);
    optionList.addInt("r_tilde","Reg",r_gd);
    optionList.addInt("regularity","Regularity of the initial geometry",R_geo);

    optionList.addInt("P_geo","Grad", P_geo);

    optionList.addSwitch("plot","Plot in Paraview",plot);
    optionList.addSwitch("latex","Latex output",latex);
    optionList.addSwitch("latex_plot","Latex output",latex_plot);

    optionList.addSwitch("twoPatch", "For the two-patch paper",false);

    optionList.addSwitch("neumann", "For computing the neumann bdy",neumann_bdy);

    optionList.addSwitch("isogeometric", "Project the basis in isogeometric concept",isogeometric);
    optionList.addSwitch("h1projection", "Project the basis in H1 norm",h1projection);
    optionList.addSwitch("h2projection", "Project the basis in H2 norm",h2projection);
    optionList.addSwitch("h1projectionProof", "Project the basis in H1 norm",h1projectionProof);

    optionList.addSwitch("info", "Print information!",info);

    optionList.addReal("threshold","Threshold",threshold);
    optionList.addReal("zero","Zero",zero);

    optionList.addReal("lambda","lambda for two Patch", lambda);
    optionList.addReal("lambda2","lambda for two Patch", lambda);

    optionList.addReal("factor","factor for two Patch", lambda);

    if (localGd)
        gluingData_strategy = gluingData::local;
    else if (exactGd)
        gluingData_strategy = gluingData::exact;
    if (localEdge)
        g1BasisEdge_strategy = g1BasisEdge::local;
    if (localVertex)
        g1BasisVertex_strategy = g1BasisVertex::local;


    optionList.addInt("gluingData","The strategy for the gluing data",gluingData_strategy);
    optionList.addInt("g1BasisEdge","The strategy for the g1 basis edge",g1BasisEdge_strategy);
    optionList.addInt("g1BasisVertex","The strategy for the g1 basis vertex",g1BasisVertex_strategy);

}

} // namespace gismo