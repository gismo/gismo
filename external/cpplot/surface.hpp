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

#ifndef _CPPLOT_SURFACE_HPP_
#define _CPPLOT_SURFACE_HPP_
namespace cpplot {
    class Surface;
    /**
     * The surface_t is a shared pointer that will keep the
     * Surface object it's pointing to alive and available.
     */
    typedef boost::shared_ptr<Surface> surface_t;
    class Surface : public drawing_t_t, public boost::enable_shared_from_this<Surface> {
        public:

            virtual ~Surface(){ }; 

            dmat XData,YData,ZData,CDataIndex; ///< Datacontainers for the plot
            tcmat CData; ///< Color data container
            dvec V; ///< Data container

            std::string FaceColor;///< ColorSpec    | none | { flat}
            std::string EdgeColor;///< ColorSpec{ k} | none | flat

            std::string LineStyle;///< {-} | - - | : | -. | none
            float LineWidth; ///< Width of plotted line
            int NContour; ///< Number of contours in plot

            enum types { _2D, _3D, contourplot } type; ///< The type of plot

            /**
             * The constructor accepts as single argument a shared pointer to the axes to which the object belongs
             */
            Surface(const axes_t a)
                :   drawing_t_t(a),
                    V(),
                    FaceColor("flat"),
                    EdgeColor("b"),
                    LineStyle("-"),
                    LineWidth(0.5),
                    NContour(10),
                    type(_2D)
                {}

            void clear(); ///< Clear all data

            /**
             * Draw the surface on the axes. This is only used internally.
             */
            void draw();


            /// Color ///
            surface_t shading(const std::string c); ///< Set the shading property
            surface_t surface(const dmat& Z); ///< Plot a surface described by the argument
            surface_t surface(const dmat& Z, const dmat& C); ///< Plot a surface described by the arguments
            surface_t surface(const dmat& Z, const tcmat& C); ///< Plot a surface described by the arguments
            surface_t surface(const dvec& x, const dvec& y, const dmat& Z); ///< Plot a surface described by the arguments
            surface_t surface(const dvec& x, const dvec& y, const dmat& Z, const dmat& C); ///< Plot a surface described by the arguments
            surface_t surface(const dvec& x, const dvec& y, const dmat& Z, const tcmat& C); ///< Plot a surface described by the arguments
            surface_t surface(const dmat& X, const dmat& Y, const dmat& Z); ///< Plot a surface described by the arguments
            surface_t surface(const dmat& X, const dmat& Y, const dmat& Z, const dmat& C); ///< Plot a surface described by the arguments
            surface_t surface(const dmat& X, const dmat& Y, const dmat& Z, const tcmat& C); ///< Plot a surface described by the arguments

            /// surf
            surface_t surf(const dvec& x, const dvec& y, const dmat& Z); ///< Plot a surf surface described by the arguments

            /// create pcolor
            surface_t pcolor(dmat C); ///< Plot a pcolor surface described by the arguments
            surface_t pcolor(const tcmat& C); ///< Plot a pcolor surface described by the arguments
            surface_t pcolor(const dvec& x, const dvec& y, const dmat& C); ///< Plot a pcolor surface described by the arguments
            surface_t pcolor(const dvec& x, const dvec& y, const tcmat& C); ///< Plot a pcolor surface described by the arguments
            surface_t pcolor(const dmat& X, const dmat& Y, const dmat& C); ///< Plot a pcolor surface described by the arguments
            surface_t pcolor(const dmat& X, const dmat& Y, const tcmat& C); ///< Plot a pcolor surface described by the arguments

            /// mesh
            surface_t mesh(const dvec& x, const dvec& y, dmat Z);///< Plot a mesh described by the arguments

            /// contour
            surface_t contour(const dmat& Z); ///< Plot a surface contour described by the arguments
            surface_t contour(const dmat& Z, const int n); ///< Plot a surface contour described by the arguments
            surface_t contour(const dmat& Z, const dvec& v); ///< Plot a surface contour described by the arguments
            surface_t contour(const dvec& x, const dvec& y, const dmat& Z); ///< Plot a surface contour described by the arguments
            surface_t contour(const dvec& x, const dvec& y, const dmat& Z, const int n); ///< Plot a surface contour described by the arguments
            surface_t contour(const dvec& x, const dvec& y, const dmat& Z, const dvec& v); ///< Plot a surface contour described by the arguments

            surface_t set(const std::string p, const std::string v); ///< Set property
            surface_t set(const std::string p, const float v); ///< Set property


            void config(); ///< Configure the axes in which the surface is printed

        private:
            void draw2d(); ///< Draw a 2D surface
            void draw3d(); ///< Draw a 3D surface
            void contourc(const dvec& x, const dvec& y, const dmat& Z, const dvec& v, dmat& C); ///< Generate contour structure
            void draw_contour(); ///< Draw contour
            boost::mutex data_mutex; ///< Protects the data when reading/writing

            // contour
            /**
             * Internal contour representation
             */
            struct ContourPoint{
                double x,y;
                int xj,yi;
                int xy;
                int done;
            };

    };
}
#endif
