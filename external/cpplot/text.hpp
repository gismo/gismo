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

#ifndef _CPPLOT_TEXT_HPP_
#define _CPPLOT_TEXT_HPP_
namespace cpplot {
    class Text;
    /**
     * The text_t is a shared pointer that will keep the
     * Text object it's pointing to alive and available.
     */
    typedef boost::shared_ptr<Text> const text_t;

    /**
     * The text object represents a string that is printed on a plot
     */
    class Text : public drawing_t_t, public boost::enable_shared_from_this<Text> {
        public:

            virtual ~Text(){ }; 

            std::string String; ///< Contains the data that is printed to the plot
            float position[3]; ///< The three-dimensional position of the text. In a 2d-plot, only the first two are used.
            void draw(); ///< Draw the text to the plot. This is only used internally

            /**
             * The constructor accepts as single argument a shared pointer to the axes to which the object belongs
             */
            Text(const axes_t a)
                :   drawing_t_t(a)
                {}
            text_t text(double x, double y, const std::string s);
            void clear(); ///< Clear all data
            void config(){} ///< Configure the axes in which the text is printed
    };
}
#endif
