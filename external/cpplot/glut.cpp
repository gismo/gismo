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
#include "figure.hpp"
#include <iostream>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/bind.hpp>

namespace cpplot {
    namespace glut {

        void idle();
        void display();
        void reshape(int w, int h);
        void mouse(int button, int state, int x, int y);
        void motion(int x, int y);
        void passivemotion(int x,int y);
        void keyboard(unsigned char key, int x, int y);


        int tool_window_number;
        typedef std::map<int, figure_t> glutmap;
        glutmap windows;
        typedef std::list<figure_t> window_queue_t;
        window_queue_t window_queue;
        boost::mutex wq_mutex;

        void set_window_title(const int wn, const std::string name) {
            assert(wn);
            glutSetWindow(wn);
            glutSetWindowTitle(name.c_str());
        }

        void create_window(const figure_t fig) {
            glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
            glutInitWindowPosition( fig->position[0], fig->position[1] );
            glutInitWindowSize( fig->position[2], fig->position[3] );

            fig->window_number = glutCreateWindow(fig->window_name.c_str());
            glutIdleFunc(glut::idle); // eats the CPU !
            glutReshapeFunc(glut::reshape);
            glutDisplayFunc(glut::display);
            glutMotionFunc(glut::motion);
            glutMouseFunc(glut::mouse);
            glutPassiveMotionFunc(glut::passivemotion);
            glutKeyboardFunc(glut::keyboard);
            glutShowWindow();

            windows.insert(glutmap::value_type(fig->window_number, fig));
        }

        void register_figure(const figure_t fig) {
            //std::cout << "Registering new figure " << std::endl;            
            boost::mutex::scoped_lock l(wq_mutex);
            window_queue.push_back(fig);
        }

        void tool() {
            boost::mutex::scoped_lock l(wq_mutex);
            window_queue_t::iterator it = window_queue.begin();
            while(it != window_queue.end()) {
                create_window(*it);
                window_queue.erase(it++);
            }
            boost::this_thread::sleep(boost::posix_time::milliseconds(200));
        }

        void run(int& argc, char* argv[]) {
            glutInit(&argc, argv);
            glutInitWindowPosition( -50, -50 );
            glutInitWindowSize( 0, 0  );
            tool_window_number = glutCreateWindow("");
            glutIdleFunc(tool);
            glutHideWindow();
            glutMainLoop();
        }

        boost::thread glut_thread;
        void init(int& argc, char* argv[]) {
            glut_thread = boost::thread(boost::bind(run, argc, argv));
        }

        void idle() { glutPostRedisplay(); usleep(1000); }
        void display() {
            //~ std::cout << "Display..."<< std::endl;
            int window = glutGetWindow();
            if(window == 0) return;

            //~ glutmap::iterator w = windows.find(window);
            //~ if(w == windows.end()) return;
            for(glutmap::iterator w = windows.begin(); w != windows.end(); ++w)
                w->second->draw();
        }
        void reshape(int w, int h) {
            //~ std::cout << "reshape" << std::endl;
            int window = glutGetWindow();
            if(window == 0) return;

            glutmap::iterator win = windows.find(window);
            if(win == windows.end()) return;
            win->second->reshape(w,h);
        }
        void mouse(int button, int state, int x, int y) {
            int window = glutGetWindow();
            if(window == 0) return;

            glutmap::iterator w = windows.find(window);
            if(w == windows.end()) return;
            w->second->mouse(button,state,x,y);
        }
        void motion(int x, int y) {
            int window = glutGetWindow();
            glutPostRedisplay();
            if(window == 0) return;

            glutmap::iterator w = windows.find(window);
            if(w == windows.end()) return;
            w->second->motion(x,y);
        }
        void passivemotion(int x,int y) {
            int window = glutGetWindow();
            if(window == 0) return;

            glutmap::iterator w = windows.find(window);
            if(w == windows.end()) return;
            w->second->passivemotion(x,y);
        }
        void keyboard(unsigned char key, int x, int y) {
            int window = glutGetWindow();
            if(window == 0) return;

            glutmap::iterator w = windows.find(window);
            if(w == windows.end()) return;
            w->second->keyboard(key,x,y);
        }
    }
}
