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

#include "figure.hpp"
#include "glut.hpp"

namespace cpplot {
    figure_t_t::figure_t_t(const std::string name, const bool viz)
                :   keyboard_callback(NULL),
                    xPassive(100),
                    yPassive(0),
                    window_name(name),
                    window_number(0),
                    visible(viz)
    {
        position[0] = position[1] = 0;
        position[2] = position[3] = 500;

        if ( name == " ")
        {
            char c[30] ;
            sprintf( c, "Gismo Plot %d", getpid() );
            window_name = c;
        }
        
    }

    figure_t_t::~figure_t_t() {
        glutDestroyWindow(window_number);
        layers.clear();
    }

    layer_t figure_t_t::layer(const std::string name, const bool visible) {
        cl = layers[name];
        if(cl == NULL) {
            boost::shared_ptr<layer_t_t> p(new layer_t_t(shared_from_this(), name, visible));
            layers[name] = cl = p;
        }
        return cl;
    }

    void figure_t_t::draw() {
        if(!visible) return;

        glutSetWindow(window_number);

        glEnable(GL_DEPTH_TEST);
        //glDepthFunc(GL_NEVER);

        glClearColor(1, 1, 1, 0.);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glViewport(0,0, (int)(window_w), (int)(window_h));
        glLoadIdentity();
        //gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );
        glOrtho(0, 1, 0, 1, -1, 3);

        for(layers_t::iterator l = layers.begin(); l != layers.end(); ++l) {
            l->second->draw();
        }

        draw_layer_list();

        glFlush();
        glViewport(0,0, window_w, window_h);
        glutSwapBuffers();
    }

    void figure_t_t::draw_layer_list() {
        if(layers.size() <= 1) return;
        int l,w,h,r;
        std::string s;
        l = 1;
        //t = 1;
        w = 20;//button_width;
        h = 20;//button_height;
        r = 3;

        if(xPassive < 25) {
            glViewport(0,0, (int)(window_w),(int)(window_h));
            glLoadIdentity();
            gluOrtho2D( 0.0, (int)(window_w), (int)(window_h),0 );

            glDisable(GL_LINE_STIPPLE);
            gl2psDisable(GL2PS_LINE_STIPPLE);

            glLineWidth(2);
            glColor3d(0,0,1);

            int j = 0;
            for(layers_t::iterator lit = layers.begin(); lit != layers.end(); ++j, ++lit) {
                glBegin(GL_LINE_STRIP);// Draw the box
                    glVertex2d(l+r   ,h*j+r );
                    glVertex2d(l+r   ,h*j+h-r );
                    glVertex2d(l+w-r ,h*j+h-r );
                    glVertex2d(l+w-r ,h*j+r );
                    glVertex2d(l+r   ,h*j+r );
                glEnd();

                if(lit->second->is_visible()) {// Tick the box
                    glBegin(GL_LINE_STRIP);
                        glVertex2d(l+9 ,h*j+5 );
                        glVertex2d(l+8 ,h*j+15 );
                        glVertex2d(l+15  ,h*j+7 );
                    glEnd();
                }

                glColor3f(0,0,1);
                glRasterPos2f( 22, h*j+h-6);
                s = lit->second->layername;
                gl2psText(s.c_str(), "Arial", 12);
                for(unsigned int i = 0; i < s.size(); ++i) {
                    glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, s[i] );
                }
            }
        }
    }


    // events (mouse)
    void figure_t_t::mouse(const int button, const int state, const int x, const int y) {
        // Layer list click
        int l,w,h;
        l = 1;
        //t = 1;
        w = 20;//button_width;
        h = 20;//button_height;

        if ( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ) {// Left Click
            int j = 0;
            for(layers_t::iterator p = layers.begin(); p != layers.end(); ++j, ++p) {
                if((l<x)&&(x<w)&&(h*j<y)&&(y<h*j+h)) {
                    // Layer was clicked!
                    //Toggle visibility
                    p->second->toggle_visibility();
                    //double click to focus layer
                    if(p->second->time_clicked.elapsed() < 0.2) {
                        for(layers_t::iterator k = layers.begin(); k != layers.end(); ++k) { k->second->set_visibility(false); }
                        p->second->set_visibility(true);
                        //cout <<"!!"<<endl;
                    }
                    p->second->time_clicked.restart();

                    return;
                }
            }
        }

        for(layers_t::reverse_iterator l = layers.rbegin(); l != layers.rend(); ++l) {
            if(l->second->mouse(button, state, x, y)) {
                selected_layer = l->second;
            }
        }
    }



    // Events ///
    void figure_t_t::reshape(int w, int h) {
        window_w = w;
        window_h = h;
        glViewport(0,0,w,h);
        //std::cout <<"window size: "<< w <<" "<<h<<std::endl;
    }
    void figure_t_t::motion(int x, int y) {
        if(selected_layer) {
            selected_layer->motion(x, y);
        }
    }
    void figure_t_t::passivemotion(const int x, const int y) {
        xPassive = x;
        yPassive = y;

        //~ for(layers_t::reverse_iterator l = layers.rbegin(); l != layers.rend(); ++l) {
            //~ if(l->second->passivemotion(x, y)) return;
        //~ }
        //cout <<"Passive: "<<x<<" "<<y<<endl;
    }
    void figure_t_t::keyboard(unsigned char key, int x, int y) {
        switch(key) {
            case 'q':
                glutDestroyWindow(window_number);
                break;
            case 'p':
                print();
                break;
        }
        if(keyboard_callback) keyboard_callback(key, x, y);
    }
    // print ///
    void figure_t_t::print(const std::string filename) {
        FILE *fp;
        int state = GL2PS_OVERFLOW, buffsize = 0;
        
        char s[30] ;
        if ( filename == "gismo_plot.eps")
            sprintf( s, "gismo_plot_%d.eps", getpid() );
        else
            sprintf( s, "%s", filename.c_str() );
            
            fp = fopen( s, "wb");
            printf("\nWriting %s ... ",s);
            //std::cout << "Writing '" << filename << "'... ";
            while(state == GL2PS_OVERFLOW) {
                buffsize += 2024*2024;
                gl2psBeginPage("test", "gl2ps", NULL, GL2PS_EPS, GL2PS_SIMPLE_SORT,
                               GL2PS_USE_CURRENT_VIEWPORT,
                               GL_RGBA, 0, NULL, 0, 0, 0, buffsize, fp, s );
                draw();
                state = gl2psEndPage();
    }
            fclose(fp);
            std::cout << "Done!" << std::endl;
    }
}
