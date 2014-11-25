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
    text_t Text::text(double x, double y, const std::string s) {
        // text on current axes
        //if(is_debug1) { cout<<"mode:"<<mode<<" Text (text):"<<i_text<<endl; }
        position[0] = x;
        position[1] = y;
        String = s;
        return shared_from_this();
    }

    void Text::clear() {
        String.clear();
    }
    //~ void Text::set_font(char font_[], int size) {
        //font=font_;
        //font_size=size;
    //~ }
    void Text::draw() {
        glColor3f(0,0,0);
        glRasterPos2f( ctx(position[0]), cty(position[1]) );
        gl2psText(String.c_str(), "Arial", 12);
        for(int i = 0; i< (int)String.size(); ++i) {
            glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, String[i] );
        }
    }
}
