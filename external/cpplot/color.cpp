/**
 * Copyright 2011, 2012 Jonatan Olofsson
 *
 * This file is part of cpplot.
 *
 * cpplot is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * cpplot is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with cpplot.  If not, see <http://www.gnu.org/licenses/>.
 */

/****************************************************************************
License: Gnu Public license (GPL) v3
* Author: Jonatan Olofsson (jonatan.olofsson@gmail.com)
* Version: 0.1
* Based on
Author: Yuichi Katori (yuichi.katori@gmail.com)
Project:MATPLOT++ (MATLAB-like plotting tool in C++).
Version:0.3.13
****************************************************************************/

#include "cpplot_common.hpp"

namespace cpplot {
    std::vector<float> ColorSpec2RGB(const std::string& c) {
        float r,g,b;
        float h=0.6;
        float f=0.5;
        //line
             if( c == "k" ) { r = 0;g = 0;b = 0; }// black
        else if( c == "r" ) { r = 1;g = 0;b = 0; }// red
        else if( c == "b" ) { r = 0;g = 0;b = 1; }// blue
        else if( c == "g" ) { r = 0;g = 1;b = 0; }// green
        else if( c == "c" ) { r = 0;g = 1;b = 1; }// cyan
        else if( c == "m" ) { r = 1;g = 0;b = 1; }// magenta
        else if( c == "y" ) { r = 1;g = 1;b = 0; }// yellow
        else if( c == "w" ) { r = 1;g = 1;b = 1; }// white

        //dark color
        else if( c == "dr" ) { r = h;g = 0;b = 0; }// red
        else if( c == "db" ) { r = 0;g = 0;b = h; }// blue
        else if( c == "dg" ) { r = 0;g = h;b = 0; }// green
        else if( c == "dc" ) { r = 0;g = h;b = h; }// cyan
        else if( c == "dm" ) { r = h;g = 0;b = h; }// magenta
        else if( c == "dy" ) { r = h;g = h;b = 0; }// yellow

        //light color
        else if( c == "lr" ) { r = 1;g = f;b = f; }// red
        else if( c == "lb" ) { r = f;g = f;b = 1; }// blue
        else if( c == "lg" ) { r = f;g = 1;b = f; }// green
        else if( c == "lc" ) { r = f;g = 1;b = 1; }// cyan
        else if( c == "lm" ) { r = 1;g = f;b = 1; }// magenta
        else if( c == "ly" ) { r = 1;g = 1;b = f; }// yellow

        //universal color
        else if( c == "ur" ) { r = 1;   g = 0.2; b = 0; }//red
        else if( c == "ub" ) { r = 0;   g = 0.25;b = 1; }//blue
        else if( c == "ug" ) { r = 0.2; g = 0.6; b = 0.4; }//green
        else if( c == "uy" ) { r = 1;   g = 1;   b = 1; }//yellow
        else if( c == "uc" ) { r = 0.4; g = 0.8; b = 1; }//sky blue
        else if( c == "up" ) { r = 1;   g = 0.6; b = 0.6; }//pink
        else if( c == "uo" ) { r = 1;   g = 0.6; b = 0; }//orange
        else if( c == "um" ) { r = 0.6; g = 0;   b = 0.4; }//perple
        else if( c == "ubr") { r = 0.4; g = 0.2; b = 0; }//brown

        std::vector<float> out(3);
        out[0] = r;
        out[1] = g;
        out[2] = b;
        return out;

    }

    std::string rgb2colorspec(std::vector<float> rgb) {
        char c[100];
        sprintf(c,"[%f %f %f]",rgb[0],rgb[1],rgb[2]);
        std::string s = c;
        return s;
    }
}
