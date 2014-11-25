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

#ifndef _CPPLOT_PATCH_HPP_
#define _CPPLOT_PATCH_HPP_
namespace cpplot {
    class Patch;
    typedef boost::shared_ptr<Patch> patch_t;
    class Patch : public drawing_t_t, public boost::enable_shared_from_this<Patch> {
        public:

            virtual ~Patch(){ }; 

            enum types {_2D, _3D} type; ///< View type
            std::vector< std::vector<int> > faces; ///< Face data container
            dmat vertices; ///< Vertice data container
            dmat XData,YData,ZData; ///< Data containers
            tcvec CData; ///< Color data container

            std::string EdgeColor,FaceColor; ///{ ColorSpec}|none|flat|interp

            std::string LineStyle; /// {-} | - - | : | -. | none
            float LineWidth; ///< Width of plotted line

            /**
             * The constructor accepts as single argument a shared pointer to the axes to which the object belongs
             */
            Patch(const axes_t a)
                :   drawing_t_t(a),
                    type(_2D),
                    EdgeColor("k"),
                    FaceColor("r"),
                    LineStyle("-"),
                    LineWidth(1)
                {}

            void clear(); ///< Clear all data


            /**
             * Draw the patch on the axes. This is only used internally.
             */
            void draw();


            // bar
            patch_t bar(const dvec& y, float width = 0.8); ///< Plot a bar
            patch_t bar(const dvec& x, const dvec& y, const float width = 0.8); ///< Plot a bar
            /// patch
            patch_t patch(const dmat& X, const dmat& Y); ///< Plot patch
            patch_t patch(const dmat& X, const dmat& Y, const dvec& C); ///< Plot patch
            patch_t patch(const dmat& X, const dmat& Y, const tcvec& C); ///< Plot patch
            patch_t patch(const dmat& X, const dmat& Y, const dmat& Z); ///< Plot patch
            patch_t patch(const dmat& X, const dmat& Y, const dmat& Z, const dvec& C); ///< Plot patch
            patch_t patch(const dmat& X, const dmat& Y, const dmat& Z, const tcvec& C); ///< Plot patch

            patch_t set(const std::string p, const std::string v); ///< Set property
            patch_t set(const std::string p, const float v); ///< Set property

            tcvec Index2TrueColor(const dvec& IC); ///< Convert color index to a real color

            void config(); ///< Configure the axes in which the patch is printed

        private:
            void draw2d(); ///< Draw 2D patch
            void draw3d(); ///< Draw 3D patch
            boost::mutex data_mutex; ///< Protects the data when reading/writing
    };
}
#endif
