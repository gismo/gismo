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

#ifndef _CPPLOT_LINE_HPP_
#define _CPPLOT_LINE_HPP_

#include <valarray>

namespace cpplot {
    class Line;
    /**
     * The linee_t is a shared pointer that will keep the
     * Line object it's pointing to alive and available.
     */
    typedef boost::shared_ptr<Line> line_t;

    /**
     * The line object represent a single line drawn in an axes
     */
    class Line : public drawing_t_t, public boost::enable_shared_from_this<Line> {
        public:

            virtual ~Line(){ }; 

            bool Errorbar; ///< Plot errorbar or not

            void clear(); ///< Clear all data
            void color(const float r, const float g, const float b); ///< Set line color

            // Matlab oriented variables //
            unsigned int max_capacity;
            bool stop_at_max_;
            dvec XData,YData,ZData; ///< Data containers
            dvec YPData,YMData; ///< Error data containers (+- error)

            std::string Color; ///< Color of the line
            std::string LineStyle;/// {-} | - - | : | -. | none
            float  LineWidth; ///< Width of plotted line
            std::string Marker; ///< {none} | . | + | x | d | ^ | v | o | * | s
            float  MarkerSize; ///< Size of data marker
            std::string MarkerEdgeColor; ///< Color of marker edge
            std::string MarkerFaceColor; ///< Color of marker face
            bool visible; ///< Wether the line is plotted or not

            /**
             * The constructor accepts as single argument a shared pointer to the axes to which the object belongs
             */
            Line(const axes_t a);
            figure_t gcf(); ///< Get the figure to which the line belongs

            void draw(); ///< Draw the line on the line's axes. Used internally only.
            void config(); ///< Configure the axes to which the line belongs. Used internally only

            line_t set_capacity(unsigned int);
            line_t stop_at_max(bool = true);

            // vertex
            void  vertex(const double x, const double y); ///< Add a line segment to the line
            line_t line(const dvec& x, const dvec& y); ///< Draw a 2D-line with the supplied data (clearing any previous data in the line)
            line_t line(const dvec& x, const dvec& y, const dvec& z); ///< Draw a 3D-line with the supplied data (clearing any previous data in the line)
            // plot, semilogx, semilogy, loglog
            line_t plot(const dvec& y); ///< Draw a 2D-line with the supplied data and integer x-values indicating the index (clearing any previous data in the line)
            line_t plot(const dvec& x, const dvec& y); ///< Draw a 2D-line with the supplied data (clearing any previous data in the line)
            line_t plot(const dvec& x, const dvec& y, const dvec& z); ///< Draw a 3D-line with the supplied data (clearing any previous data in the line)
            line_t plot(std::valarray<double> x, std::valarray<double> y); ///< Draw a 2D-line with the supplied data (clearing any previous data in the line)
            line_t semilogx(const dvec& x, const dvec& y); ///< Draw a 2D-line with the supplied data (clearing any previous data in the line) where the x-axis is drawn logarithmically
            line_t semilogy(const dvec& x, const dvec& y); ///< Draw a 2D-line with the supplied data (clearing any previous data in the line) where the y-axis is drawn logarithmically
            line_t loglog(const dvec& x, const dvec& y); ///< Draw a 2D-line with the supplied data (clearing any previous data in the line) where both the axes are drawn logarithmically
            // errorbar
            void vertex(const double x, const double y, const double ep, const double em); ///< Add line segment with error to the plot
            line_t errorbar(const dvec& x, const dvec& y,dvec e); ///< Draw a 2D-line with the supplied data (clearing any previous data in the line) with error indicated in the plot
            line_t errorbar(const dvec& x, const dvec& y,dvec ep, const dvec& em); ///< ///< Draw a 2D-line with the supplied data (clearing any previous data in the line) with error indicated in the plot
            // 3D line
            void vertex(const double x, const double y, const double z); ///< Add a line segment to the line
            line_t set(const float v); ///< Set linewidth
            line_t set(const std::string v); ///< Set property
            line_t set(const std::string p, const std::string v); ///< Set property
            line_t set(const std::string p, const float v); ///< Set property

        private:
            std::pair<double, double> min_max(const dvec&, math::scale); ///< Find min and max of data
            boost::mutex data_mutex; ///< Protects the data when reading/writing
    };
}
#endif
