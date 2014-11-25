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

#ifndef _CPPLOT_LAYER_HPP_
#define _CPPLOT_LAYER_HPP_

#include "cpplot_common.hpp"
#include <boost/timer.hpp>

namespace cpplot {
    /**
     * Each layer is contained in a parent figure, which may toggle each
     * layers visibility. The layers in turn contain the axes in which the
     * plots are drawn, and is thus responsible for keeping track of
     * which axes objects that belongs to the layer (and create new when needed).
     */
    class layer_t_t : public boost::enable_shared_from_this<layer_t_t>, public boost::noncopyable {
        private:
            axes_t ca, selected_axes; ///< Pointers to current (in program) and selected (by user) axes
            bool visible; ///< Decides if the layer is drawn or not
            float xButtonDown, yButtonDown; /// Last clicked mouse position

        public:
            boost::timer time_clicked; ///< Time that the layer was last clicked. Used to detect double click
            std::string layername; ///< Name of layer
            figure_t figure; ///< Pointer to the figure that that the layer is drawn in

            axess_t axes; ///< Axes in the layer

            /**
             * The layer constructor connects the layer to its figure, assigns a name, and its visibility
             */
            layer_t_t(const figure_t fig, const std::string name, const bool viz)
                :   visible(viz),
                    layername(name),
                    figure(fig)
                {}

            void draw(); ///< Draw the layer to the figure
            axes_t subplot(const int m, const int n, const int p); ///< Select (and create if needed) an axes places as a subplot in the layer
            axes_t gca() { return ca ? ca : subplot(1,1,1); } ///< Get the current axes from the layer
            figure_t gcf() { return figure; } ///< Get the layer's figure
            layer_t clear() { ca.reset(); axes.clear(); return shared_from_this(); } ///< Clear the layer

            void toggle_visibility(); ///< Toggle the visibility of the layer
            void set_visibility(bool v); ///< Set the visibility of the layer
            bool is_visible() { return visible; } ///< Check if the layer is visible
            bool mouse(const int button, const int state, const int x, const int y); ///< Callback when a mouse button is clicked.
            bool motion(const int x, const int y); ///< Callback when mouse is moved with button down
            //~ const bool passivemotion(const int x, const int y);

            // interface ///
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


            // Surface, Contour ///
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

            // Patch ///
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
    };
}
#endif
