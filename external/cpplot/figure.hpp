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

#ifndef _CPPLOT_FIGURE_HPP_
#define _CPPLOT_FIGURE_HPP_

#include "cpplot_common.hpp"

namespace cpplot {
    /**
     * The figure class represents a single window, and may contain
     * one or more layers (which in turn contain the axes, see layer.hpp)
     *
     * When several layers are used, a list of the layers will be drawn
     * in the upper left corner with tickboxes available to control the
     * visibility of each layer.
     */
    class figure_t_t : public boost::enable_shared_from_this<figure_t_t>, boost::noncopyable {
        private:
            layers_t layers; ///< List of all associated layers in figure

        public:
            // User defined event callbacks
            void(*keyboard_callback)(unsigned char, int, int); ///< User defined keyboard callback

            // Variables changed by events //
            int window_w, window_h; ///< GL window dimensions
            int xPassive, yPassive; ///< Passive mouse motion

            std::string window_name; ///< Name of the window title
            int window_number; ///< GL window number
            int position[4]; /// Window position and size; left top width height
            bool visible; ///< Decides if the window should be visible or hidden
            layer_t cl, selected_layer; ///< Pointers to current (in program) and selected (by user) layers

            /**
             * The window is initiated with a window title and a boolean indicating its visibility
             */
            figure_t_t(const std::string name = " ", const bool viz = true);
            /**
             * The destructor makes sure the GL window is closed along with the object
             */
            ~figure_t_t();

            void set_window_name(const std::string name) {
                window_name = name;
                if(window_number) glut::set_window_title(window_number, name);
            }

            /**
             * Create a new layer in the figure
             */
            layer_t layer(const std::string name = "default", const bool viz = true);
            layer_t gcl() { return cl ? cl : layer(); } ///< Get the current layer, or create a new one if needed
            figure_t clear() { cl.reset(); layers.clear(); return shared_from_this(); } ///< Clear the figure
            void draw(); ///< Draw the layers
            void draw_layer_list(); ///< Draw the layer list

            // GLUT Callback Functions ///
            void reshape(const int w, const int h); ///< Callback on window resize
            void mouse(const int button, const int state, const int x, const int y); ///< Callback on mouse click
            void motion(const int x, const int y); ///< Callback on mouse motion when button is held down
            void passivemotion(const int x, const int y); ///< Callback on mouse motion when NO button is held down
            void keyboard(const unsigned char key, const int x, const int y); ///< Callback on key press

            figure_t gcf() { return shared_from_this(); } ///< Get pointer to self
            // interface ///
            axes_t gca() { return gcl()->gca(); } ///< Get the currently active axes object (created if nescessary)
            axes_t subplot(const int m, const int n, const int p) { return gcl()->subplot(m,n,p); }

            template<typename T>
            void set(const std::string v) { gca()->gco<T>()->set(v); } ///< Set property
            void set(const std::string v) { gca()->gco<Line>()->set(v); } ///< Set property
            template<typename T>
            void set(const float v) { gca()->gco<T>()->set(v); } ///< Set property
            void set(const float v) { gca()->gco<Line>()->set(v); } ///< Set property
            template<typename T>
            void set(const std::string p, const std::string v) { gca()->gco<T>()->set(p,v); } ///< Set property
            void set(const std::string p, const std::string v) { gca()->gco<Line>()->set(p,v); } ///< Set property
            template<typename T>
            void set(const std::string p, float v) { gca()->gco<T>()->set(p,v); } ///< Set property
            void set(const std::string p, float v) { gca()->gco<Line>()->set(p,v); } ///< Set property

            // Axes
            void axis(const double xMin, const double xMax, const double yMin, const double yMax) {
                gca()->axis(xMin, xMax, yMin, yMax);
            } ///< See axes.hpp
            void axis(const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax) {
                gca()->axis(xMin, xMax, yMin, yMax, zMin, zMax);
            } ///< See axes.hpp
            void axis(const std::string s) { gca()->axis(s); } ///< See axes.hpp
            void axis(const bool s) { gca()->axis(s); } ///< See axes.hpp
            void grid(const std::string s) { gca()->grid(s); } ///< See axes.hpp
            void grid(const bool s) { gca()->grid(s); } ///< See axes.hpp
            void ticklabel(const bool s) { gca()->ticklabel(s); } ///< See axes.hpp
            void title(const std::string s) { gca()->title(s); } ///< See axes.hpp
            void xlabel(const std::string s) { gca()->xlabel(s); } ///< See axes.hpp
            void ylabel(const std::string s) { gca()->ylabel(s); } ///< See axes.hpp
            void mouse_capture(const bool y) { gca()->mouse_capture(y); } ///< See axes.hpp

            axes_t colorbar() { return gca()->colorbar(); } ///< See axes.hpp
            void gray() { gca()->gray(); }; ///< See axes.hpp
            void jet() { gca()->jet(); } ///< See axes.hpp
            void hsv() { gca()->hsv(); } ///< See axes.hpp
            void hot() { gca()->hot(); } ///< See axes.hpp
            void cool() { gca()->hot(); } ///< See axes.hpp
            void spring() { gca()->spring(); } ///< See axes.hpp
            void summer() { gca()->summer(); } ///< See axes.hpp
            void autumn() { gca()->autumn(); } ///< See axes.hpp
            void winter() { gca()->winter(); } ///< See axes.hpp

            // Line
            void vertex(const double x, const double y) { gca()->gco<Line>()->vertex(x,y); } ///< See line.hpp
            void vertex(const double x, const double y, const double z) { gca()->gco<Line>()->vertex(x,y,z); } ///< See line.hpp

            line_t plot(const dvec& y) { return gca()->add<Line>()->plot(y); } ///< See line.hpp
            line_t plot(const dvec& x,const dvec& y) { return gca()->add<Line>()->plot(x,y); } ///< See line.hpp
            line_t plot(const dvec& x, const dvec& y, const dvec& z) { return gca()->add<Line>()->plot(x,y,z); } ///< See line.hpp

            line_t semilogx(const dvec& x, const dvec& y) { return gca()->add<Line>()->semilogx(x,y); } ///< See line.hpp
            line_t semilogy(const dvec& x, const dvec& y) { return gca()->add<Line>()->semilogy(x,y); } ///< See line.hpp
            line_t loglog(const dvec& x, const dvec& y)   { return gca()->add<Line>()->loglog(x,y); } ///< See line.hpp

            void vertex(const double x, const double y, const double ep, const double em)
                { gca()->gco<Line>()->vertex(x,y,ep,em); } ///< See line.hpp
            void errorbar(const dvec& x, const dvec& y, const dvec& e)
                { gca()->gco<Line>()->errorbar(x,y,e); } ///< See line.hpp
            void errorbar(const dvec& x, const dvec& y, const dvec& ep, const dvec& em)
                { gca()->gco<Line>()->errorbar(x,y,ep, em); } ///< See line.hpp


            // Surface, Contour
            surface_t surface(const dmat& Z) { return gca()->add<Surface>()->surface(Z); } ///< See surface.hpp
            surface_t surface(const dmat& Z, const dmat& C) { return gca()->add<Surface>()->surface(Z, C); } ///< See surface.hpp
            surface_t surface(const dmat& Z, const tcmat& C) { return gca()->add<Surface>()->surface(Z, C); } ///< See surface.hpp
            surface_t surface(const dvec& x, const dvec& y, const dmat& Z) { return gca()->add<Surface>()->surface(x,y,Z); } ///< See surface.hpp
            surface_t surface(const dvec& x, const dvec& y, const dmat& Z, const dmat& C) { return gca()->add<Surface>()->surface(x,y,Z,C); } ///< See surface.hpp
            surface_t surface(const dvec& x, const dvec& y, const dmat& Z, const tcmat& C) { return gca()->add<Surface>()->surface(x,y,Z,C); } ///< See surface.hpp
            surface_t surface(const dmat& X, const dmat& Y, const dmat& Z) { return gca()->add<Surface>()->surface(X,Y,Z); } ///< See surface.hpp
            surface_t surface(const dmat& X, const dmat& Y, const dmat& Z, const dmat& C) { return gca()->add<Surface>()->surface(X,Y,Z,C); } ///< See surface.hpp
            surface_t surface(const dmat& X, const dmat& Y, const dmat& Z, const tcmat& C) { return gca()->add<Surface>()->surface(X,Y,Z,C); } ///< See surface.hpp

            surface_t pcolor(const dmat& C) { return gca()->add<Surface>()->pcolor(C); } ///< See surface.hpp
            surface_t pcolor(const tcmat& C) { return gca()->add<Surface>()->pcolor(C); } ///< See surface.hpp
            surface_t pcolor(const dvec& x, const dvec& y, const dmat& C) { return gca()->add<Surface>()->pcolor(x,y,C); } ///< See surface.hpp
            surface_t pcolor(const dvec& x, const dvec& y, const tcmat& C) { return gca()->add<Surface>()->pcolor(x,y,C); } ///< See surface.hpp
            surface_t pcolor(const dmat& X, const dmat& Y, const dmat& C) { return gca()->add<Surface>()->pcolor(X,Y,C); } ///< See surface.hpp
            surface_t pcolor(const dmat& X, const dmat& Y, const tcmat& C) { return gca()->add<Surface>()->pcolor(X,Y,C); } ///< See surface.hpp

            surface_t contour(const dmat& Z) { return gca()->add<Surface>()->contour(Z); } ///< See surface.hpp
            surface_t contour(const dmat& Z,int n) { return gca()->add<Surface>()->contour(Z, n); } ///< See surface.hpp
            surface_t contour(const dmat& Z, const dvec& v) { return gca()->add<Surface>()->contour(Z,v); } ///< See surface.hpp
            surface_t contour(const dvec& x, const dvec& y, const dmat& Z) { return gca()->add<Surface>()->contour(x,y,Z); } ///< See surface.hpp
            surface_t contour(const dvec& x, const dvec& y, const dmat& Z, const int n) { return gca()->add<Surface>()->contour(x,y,Z,n); } ///< See surface.hpp
            surface_t contour(const dvec& x, const dvec& y, const dmat& Z, const dvec& v) { return gca()->add<Surface>()->contour(x,y,Z,v); } ///< See surface.hpp

            surface_t mesh(const dvec& x, const dvec& y, const dmat& Z) { return gca()->add<Surface>()->mesh(x,y,Z); } ///< See surface.hpp
            surface_t surf(const dvec& x, const dvec& y, const dmat& Z) { return gca()->add<Surface>()->surf(x,y,Z); } ///< See surface.hpp

            void shading(const std::string c) { gca()->gco<Surface>()->shading(c); } ///< See surface.hpp

            // Patch
            patch_t patch(const dmat& X, const dmat& Y) { return gca()->add<Patch>()->patch(X,Y); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const dvec& C) { return gca()->add<Patch>()->patch(X,Y,C); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const tcvec& C) { return gca()->add<Patch>()->patch(X,Y,C); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const dmat& Z) { return gca()->add<Patch>()->patch(X,Y,Z); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const dmat& Z, const dvec& C) { return gca()->add<Patch>()->patch(X,Y,Z,C); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const dmat& Z, const tcvec& C) { return gca()->add<Patch>()->patch(X,Y,Z,C); } ///< See patch.hpp

            patch_t bar(const dvec& y) { return gca()->add<Patch>()->bar(y); } ///< See patch.hpp
            patch_t bar(const dvec& y, const float width) { return gca()->add<Patch>()->bar(y, width); } ///< See patch.hpp
            patch_t bar(const dvec& x, const dvec& y) { return gca()->add<Patch>()->bar(x,y); } ///< See patch.hpp
            patch_t bar(const dvec& x, const dvec& y, const float width) { return gca()->add<Patch>()->bar(x,y,width); } ///< See patch.hpp

            // Text ///
            //TODO: more fonts
            text_t text(const double x, const double y, const std::string s) { return gca()->add<Text>()->text(x,y,s); } ///< See text.hpp
            //~ void set_font(const std::string font_, const int size);

            void print(const std::string filename = "gismo_plot.eps"); ///< Print the figure contents to an eps file
    };
}
#endif
