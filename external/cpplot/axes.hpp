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

#ifndef _CPPLOT_AXES_HPP_
#define _CPPLOT_AXES_HPP_

#include "cpplot_common.hpp"

namespace cpplot {
    /**
     * The axes objects are containers (visually and conceptually) to the
     * lines, patches, surfaces etc that constitutes the plotted data.
     * Axes are kept by layers, and in turn keep track of, and create, all
     * its childred so call drawings (lines, surfaces, etc..)
     */
    class axes_t_t : public boost::enable_shared_from_this<axes_t_t>, public boost::noncopyable {
        public:
            std::vector<std::vector<float> > cmap; ///< Colormap data container
        private:
            boost::mutex children_mutex; ///< Protects the list of drawings from collisions on multithreaded read/write
            layer_t layer; ///< Pointer to the layer to which the axes belongs
            axes_t color_bar_axes; ///< Pointer to an optional bar with color legend
            int window_h(); ///< Obtain the figure height
            int window_w(); ///< Obtain the figure width
            void draw2d(); ///< Draw 2D axes
            void draw3d(); ///< Draw 3D axes
            void draw_colorbar(); ///< Draw color legend bar
            float ctx(const double x); ///< Transform coordinate to local coordinate system (2D, x)
            float cty(const double y); ///< Transform coordinate to local coordinate system (2D, y)
            float ct3x(const double x); ///< Transform coordinate to local coordinate system (3D, x)
            float ct3y(const double y); ///< Transform coordinate to local coordinate system (3D, y)
            float ct3z(const double z); ///< Transform coordinate to local coordinate system (3D, z)

        public:
            /**
             * The constructor accepts as argument the layer in which
             * the axes object is placed, and initializes the member variables
             * with their default values.
             */
            axes_t_t(layer_t);
            ~axes_t_t(){children.clear();} ///< The destructor makes sure it retains no hooks on its children when its dead
            layer_t gcl() { return layer; } ///< Get the parent layer
            float cta0,phi0;///< default value or specified by command line
            float cta,phi;  ///< controlled by mouse
            // Mouse
            double XMouse,YMouse; ///< Last clicked position
            bool Mouse; ///< Capture the mouse actions
            int xButtonDown, yButtonDown;///< last clicked mouse position
            float ctaButtonDown, phiButtonDown; ///< View angles when the button is pressed
            bool(*mouse_callback)(int button, int state, int x, int); ///< Optional user-defined callback when a mouse button is pressed

            double xmin,xmax,ymin,ymax,zmin,zmax; ///< Contains the limits of the data contained in the child drawings

            void reset_limits(); ///< Reset the axes limits for recalculation
            void config(); ///< Configure the axes object to fit the data
            axes_t set(const std::string p, const std::string v); ///< Set property

            axes_t colorbar(); ///< Add legend colorbar to the axes


            void ptext(const float x, const float y, const std::string s); ///< Draw text to the axes
            void ptext3(const float x, const float y, const float z, const std::string s); ///< Draw text to the axes in 3D, left adjusted at the coordinates
            void ptext3c(const float x, const float y, const float z, const std::string s); ///< Draw text to the axes in 3D, centered at the coordinates

            bool mouse(const int button, const int state, const int x, const int y); ///< Callback when mouse button clicked
            bool motion(const int x, const int y); ///< Callback when mouse moved with button pressed

            dvec make_tick(const double min, const double max); ///< Draw ticks to the axes

            enum types {_2D, _3D, color_bar} type; ///< Indicates the type of plotted data

            // Matlab variables //
            // styles
            bool Box;///< Axes on/off
            std::string GridLineStyle; ///< Style of plotted grid
            float LineWidth; ///< Width of grid lines
            std::string TickDir; ///< Direction of axes ticks; {in} | out
            bool visible; ///< Axes drawn or not. This does not affect the visibility of the children, only the axes.
            bool XGrid,YGrid,ZGrid; ///< Draw grid

            // General Information
            drawings_t children; ///< Container for children drawings
            drawing_t co; ///< Currently active drawing object
            bool selected; ///< Axes have been clicked by the user
            float position[4]; ///< Position and size of the axes (left bottom width height)
            float viewport3d[4];///< Position and size of the axes used by 3D-view (left bottom width height)

            // Scale and axes
            std::string XAxisLocation; ///< Location where to draw the x-axis; top | {bottom}
            std::string YAxisLocation; ///< Location where to draw the x-axis; {left} | right

            double XLim[2],YLim[2],ZLim[2]; ///< Plot range
            enum modes { automatic, manual } XLimMode,YLimMode,ZLimMode; ///< Mode of axes adjustment to data; {automatic} | manual
            math::scale XScale,YScale,ZScale; ///< Scale of axes; {linear_scale} | {logarithmic_scale}

            dvec XTick,YTick,ZTick; ///< Storage of tick locations
            bool TickLabel; ///< Label the ticks

            //View
            float camera_position[3]; ///< Position of the 3D-view camera
            float camera_target[3]; ///< Position of the 3D-view target
            float camera_up_vector[3]; ///< Orientation of the 3D-view camera

            // Label
            std::string Title; ///< Axes title
            std::string XLabel,YLabel,ZLabel; ///< Axes labels

            double CLim[2]; ///< Color limits

            void draw(); ///< Draw the axes and its children


            /// interface
            /**
             * Set the axis limits (2D)
             */
            axes_t axis(const double xMin, const double xMax, const double yMin, const double yMax);
            /**
             * Set the axis limits (3D)
             */
            axes_t axis(const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax);
            axes_t axis(const std::string s); ///< Set the visibility of the axis "on" or "off"
            axes_t axis(const bool s = true); ///< Set the visibility of the axis
            axes_t grid(const std::string s); ///< Set the visibility of the grid "on" or "off"
            axes_t grid(bool s = true); ///< Set the visibility of the grid
            axes_t ticklabel(const bool s = true); ///< Set the visibility of the ticklabels
            axes_t title(const std::string s); ///< Set the title of the axes
            axes_t xlabel(const std::string s); ///< Set the label for the x-axis
            axes_t ylabel(const std::string s); ///< Set the label for the y-axis
            axes_t mouse_capture(const bool y = true); ///< Capture mouse events

            /**
             * Add a new child drawing to the axes. The added type must
             * inherit drawing_t_t.
             */
            template<typename T>
            boost::shared_ptr<T> add() {
                boost::shared_ptr<T> p(new T(shared_from_this()));
                co = boost::dynamic_pointer_cast<drawing_t_t, T>(p);
                assert(co);
                boost::mutex::scoped_lock l(children_mutex);
                children.push_back(co);
                return p;
            }

            /**
             * Get the current drawing object. If the type of object requested is
             * incompatible with the currently selected object, or if no
             * object has yet been created, a new one will be created, added and returned.
             */
            template<typename T>
            boost::shared_ptr<T> gco() {
                if(!co) return add<T>();
                boost::shared_ptr<T> ptr = boost::dynamic_pointer_cast<T, drawing_t_t>(co);
                return ptr ? ptr : add<T>();
            }

            /// Clear the figure
            axes_t clear() {
                boost::mutex::scoped_lock l(children_mutex);
                co.reset(); children.clear(); return shared_from_this();
            }

            // Colors ///
            void color(const float r, const float g, const float b); ///< Set color of axes
            std::vector<float> colormap(const std::string c, const float t); ///< Set the colormap of the axes children
            void colormap(const std::string c); ///< Set the colormap of the axes children
            void colormap(const std::vector<std::vector<float> >& c); ///< Set the colormap of the axes children


            void gray() { colormap("Gray"); }; ///< Set the colormap of the axes children
            void jet() { colormap("Jet"); } ///< Set the colormap of the axes children
            void hsv() { colormap("HSV"); } ///< Set the colormap of the axes children
            void hot() { colormap("Hot"); } ///< Set the colormap of the axes children
            void cool() { colormap("Cool"); } ///< Set the colormap of the axes children
            void spring() { colormap("Spring"); } ///< Set the colormap of the axes children
            void summer() { colormap("Summer"); } ///< Set the colormap of the axes children
            void autumn() { colormap("Autumn"); } ///< Set the colormap of the axes children
            void winter() { colormap("Winter"); } ///< Set the colormap of the axes children

            std::vector<float> map2color(const double x); ///< Convert color value to its mapped corresponding color

            void vertex(const double x, const double y) { gco<Line>()->vertex(x,y); } ///< See line.hpp
            void vertex(const double x, const double y, const double z) { gco<Line>()->vertex(x,y,z); } ///< See line.hpp

            line_t plot(const dvec& y) { return add<Line>()->plot(y); } ///< See line.hpp
            line_t plot(const dvec& x,const dvec& y) { return add<Line>()->plot(x,y); } ///< See line.hpp
            line_t plot(const dvec& x, const dvec& y, const dvec& z) { return add<Line>()->plot(x,y,z); } ///< See line.hpp

            line_t semilogx(const dvec& x, const dvec& y) { return add<Line>()->semilogx(x,y); } ///< See line.hpp
            line_t semilogy(const dvec& x, const dvec& y) { return add<Line>()->semilogy(x,y); } ///< See line.hpp
            line_t loglog(const dvec& x, const dvec& y)   { return add<Line>()->loglog(x,y); } ///< See line.hpp

            void vertex(const double x, const double y, const double ep, const double em)
                { gco<Line>()->vertex(x,y,ep,em); } ///< See line.hpp
            void errorbar(const dvec& x, const dvec& y, const dvec& e)
                { gco<Line>()->errorbar(x,y,e); } ///< See line.hpp
            void errorbar(const dvec& x, const dvec& y, const dvec& ep, const dvec& em)
                { gco<Line>()->errorbar(x,y,ep,em); } ///< See line.hpp


            // Surface, Contour ///
            surface_t surface(const dmat& Z) { return add<Surface>()->surface(Z); } ///< See surface.hpp
            surface_t surface(const dmat& Z, const dmat& C) { return add<Surface>()->surface(Z, C); } ///< See surface.hpp
            surface_t surface(const dmat& Z, const tcmat& C) { return add<Surface>()->surface(Z, C); } ///< See surface.hpp
            surface_t surface(const dvec& x, const dvec& y, const dmat& Z) { return add<Surface>()->surface(x,y,Z); } ///< See surface.hpp
            surface_t surface(const dvec& x, const dvec& y, const dmat& Z, const dmat& C) { return add<Surface>()->surface(x,y,Z,C); } ///< See surface.hpp
            surface_t surface(const dvec& x, const dvec& y, const dmat& Z, const tcmat& C) { return add<Surface>()->surface(x,y,Z,C); } ///< See surface.hpp
            surface_t surface(const dmat& X, const dmat& Y, const dmat& Z) { return add<Surface>()->surface(X,Y,Z); } ///< See surface.hpp
            surface_t surface(const dmat& X, const dmat& Y, const dmat& Z, const dmat& C) { return add<Surface>()->surface(X,Y,Z,C); } ///< See surface.hpp
            surface_t surface(const dmat& X, const dmat& Y, const dmat& Z, const tcmat& C) { return add<Surface>()->surface(X,Y,Z,C); } ///< See surface.hpp

            surface_t pcolor(const dmat& C) { return add<Surface>()->pcolor(C); } ///< See surface.hpp
            surface_t pcolor(const tcmat& C) { return add<Surface>()->pcolor(C); } ///< See surface.hpp
            surface_t pcolor(const dvec& x, const dvec& y, const dmat& C) { return add<Surface>()->pcolor(x,y,C); } ///< See surface.hpp
            surface_t pcolor(const dvec& x, const dvec& y, const tcmat& C) { return add<Surface>()->pcolor(x,y,C); } ///< See surface.hpp
            surface_t pcolor(const dmat& X, const dmat& Y, const dmat& C) { return add<Surface>()->pcolor(X,Y,C); } ///< See surface.hpp
            surface_t pcolor(const dmat& X, const dmat& Y, const tcmat& C) { return add<Surface>()->pcolor(X,Y,C); } ///< See surface.hpp

            surface_t contour(const dmat& Z) { return add<Surface>()->contour(Z); } ///< See surface.hpp
            surface_t contour(const dmat& Z,int n) { return add<Surface>()->contour(Z, n); } ///< See surface.hpp
            surface_t contour(const dmat& Z, const dvec& v) { return add<Surface>()->contour(Z,v); } ///< See surface.hpp
            surface_t contour(const dvec& x, const dvec& y, const dmat& Z) { return add<Surface>()->contour(x,y,Z); } ///< See surface.hpp
            surface_t contour(const dvec& x, const dvec& y, const dmat& Z, const int n) { return add<Surface>()->contour(x,y,Z,n); } ///< See surface.hpp
            surface_t contour(const dvec& x, const dvec& y, const dmat& Z, const dvec& v) { return add<Surface>()->contour(x,y,Z,v); } ///< See surface.hpp

            surface_t mesh(const dvec& x, const dvec& y, const dmat& Z) { return add<Surface>()->mesh(x,y,Z); } ///< See surface.hpp
            surface_t surf(const dvec& x, const dvec& y, const dmat& Z) { return add<Surface>()->surf(x,y,Z); } ///< See surface.hpp

            void shading(const std::string c) { gco<Surface>()->shading(c); } ///< See surface.hpp

            // Patch ///
            patch_t patch(const dmat& X, const dmat& Y) { return add<Patch>()->patch(X,Y); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const dvec& C) { return add<Patch>()->patch(X,Y,C); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const tcvec& C) { return add<Patch>()->patch(X,Y,C); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const dmat& Z) { return add<Patch>()->patch(X,Y,Z); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const dmat& Z, const dvec& C) { return add<Patch>()->patch(X,Y,Z,C); } ///< See patch.hpp
            patch_t patch(const dmat& X, const dmat& Y, const dmat& Z, const tcvec& C) { return add<Patch>()->patch(X,Y,Z,C); } ///< See patch.hpp

            patch_t bar(const dvec& y) { return add<Patch>()->bar(y); } ///< See patch.hpp
            patch_t bar(const dvec& y, const float width) { return add<Patch>()->bar(y, width); } ///< See patch.hpp
            patch_t bar(const dvec& x, const dvec& y) { return add<Patch>()->bar(x,y); } ///< See patch.hpp
            patch_t bar(const dvec& x, const dvec& y, const float width) { return add<Patch>()->bar(x,y,width); } ///< See patch.hpp

            // Text ///
            //TODO: more fonts
            text_t text(const double x, const double y, const std::string s) { return add<Text>()->text(x,y,s); } ///< See text.hpp
    };
}
#endif
