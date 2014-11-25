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
    /// subplot
    axes_t layer_t_t::subplot(const int m, const int n, const int p) {
        int pp = p-1;
        ca = axes[pp];
        if(ca == NULL) {
            axes_t a(new axes_t_t(shared_from_this()));
            axes[pp] = ca = a;
            int ix = pp/m;
            int iy = (m-1) - pp%m;

            //~ ix = (p-1)%n;
            //~ iy = (m-1)-(p-1)/n;
            ca->position[0] = (ix + 0.13)/n;
            ca->position[1] = (iy + 0.11)/m;
            ca->position[2] = 0.775/n;
            ca->position[3] = 0.815/m;

            ca->viewport3d[0] = 1.0*ix/n;
            ca->viewport3d[1] = 1.0*iy/m;
            ca->viewport3d[2] = 1.0/n;
            ca->viewport3d[3] = 1.0/m;
        }

        return ca;
    }

    bool layer_t_t::mouse(const int button, const int state, const int x, const int y) {
        float X,Y;
        double rx,ry,mx,my;//mouse
        X = (float)                  x  / figure->window_w;
        Y = (float)(figure->window_h-y) / figure->window_h;
        float l,b,w,h;

        // mouse capture axes //
        if(selected_axes && selected_axes->visible) {
            l = selected_axes->position[0];
            b = selected_axes->position[1];
            w = selected_axes->position[2];
            h = selected_axes->position[3];

            if( ( selected_axes->Mouse ) && ( (l<=X)&&(X<=l+w)&&(b<=Y)&&(Y<=b+h) ) ) {
                rx = (X-l)/w;
                ry = (Y-b)/h;
                mx = rx*( selected_axes->XLim[1] - selected_axes->XLim[0] ) + selected_axes->XLim[0];
                my = ry*( selected_axes->YLim[1] - selected_axes->YLim[0] ) + selected_axes->YLim[0];
                if(selected_axes->mouse(button, state, mx, my)) return true;
            }
        }

            // axes select //
        if ( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ) {// Left Click


            if(selected_axes) selected_axes->selected = false;
            for(axess_t::reverse_iterator a = axes.rbegin(); a != axes.rend(); ++a) {
                if(ca->visible) {
                    l = a->second->position[0];
                    b = a->second->position[1];
                    w = a->second->position[2];
                    h = a->second->position[3];

                    if( (l<=X)&&(X<=l+w)&&(b<=Y)&&(Y<=b+h) ) {
                        // Click was inside axes
                        selected_axes = a->second;
                        a->second->selected = true;
                        a->second->xButtonDown = x;
                        a->second->yButtonDown = y;
                        return true; // Stop propagation
                    }
                }
            }
        }// left click

        return false; // Continue propagation
    }

    // events (mouse, motion)
    bool layer_t_t::motion(const int x, const int y) {
        if(selected_axes) return selected_axes->motion(x,y);
        else return false;
    }

    void layer_t_t::toggle_visibility() {
        set_visibility(!visible);
    }

    void layer_t_t::set_visibility(bool v) {
        visible = v;
        for(axess_t::iterator a = axes.begin(); a != axes.end(); ++a) {
            a->second->visible = v;
        }
    }

    void layer_t_t::draw() {
        if(!visible) return;
        for(axess_t::iterator a = axes.begin(); a != axes.end(); ++a) {
            a->second->draw();
        }
    }
}
