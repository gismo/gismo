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
#include "color.hpp"
#include <limits>

namespace cpplot {
    Line::Line(const axes_t a)
        :   drawing_t_t(a),
            Errorbar(false),
            max_capacity(std::numeric_limits<int>::max()),
            stop_at_max_(false),
            XData(), YData(), ZData(),
            YPData(), YMData(),
            Color("b"),
            LineStyle("-"),
            LineWidth(0.5),
            Marker("none"),
            MarkerSize(6),
            MarkerEdgeColor("k"),
            MarkerFaceColor("w"),
            visible(true)
        {}
    figure_t Line::gcf() { return ca->gcl()->gcf(); }
    void Line::clear() {
        boost::mutex::scoped_lock l(data_mutex);
        XData.clear();
        YData.clear();
        ZData.clear();
        YPData.clear();
        YMData.clear();
    }

    void Line::color(float r,float g,float b) {
        Color[0]=r;
        Color[1]=g;
        Color[2]=b;
    }

    void Line::draw() {
        boost::mutex::scoped_lock l(data_mutex);
        if(XData.size() == 0) return;

        float xx,yy;// transformed coordinates
        //float r;//marker size
        float rx,ry;
        std::vector<float> rgb = ColorSpec2RGB(Color);
        glColor3f(rgb[0],rgb[1],rgb[2]);

        glLineWidth(LineWidth);
        glPointSize(LineWidth);
        gl2psLineWidth(LineWidth);
        gl2psPointSize(LineWidth);
        // 2D //
        if(ca->type == axes_t_t::_2D) {
            if(LineStyle !="none") {// Line //
                if(LineStyle == "-") {
                    glDisable(GL_LINE_STIPPLE);
                    gl2psDisable(GL2PS_LINE_STIPPLE);
                }
                else if(LineStyle == "- -") {
                    glEnable(GL_LINE_STIPPLE);
                    glLineStipple(1, 0xF0F0);
                    gl2psEnable(GL2PS_LINE_STIPPLE);
                }
                else if(LineStyle == ":") {
                    glEnable(GL_LINE_STIPPLE);
                    glLineStipple(1, 0xCCCC);
                    gl2psEnable(GL2PS_LINE_STIPPLE);
                }
                else if(LineStyle == "-.") {
                    glEnable(GL_LINE_STIPPLE);
                    glLineStipple(1, 0x087F);
                    gl2psEnable(GL2PS_LINE_STIPPLE);
                }

                glBegin(GL_LINE_STRIP);
                for(unsigned int i = 0; i < XData.size(); ++i) {
                    //printf("i:%d %f %f\n",i,xx,yy);
                    xx = ctx(XData[i]);
                    yy = cty(YData[i]);
                    //~ std::cout << "(x,y), (xx, yy): (" << XData[i] << "," << YData[i] << "), (" << xx << "," << yy << ")" << std::endl;
                    glVertex2d(xx,yy);
                }
                glEnd();
            }

            if(Marker != "none") {// Marker //

                // Scale with window size
                //r  = MarkerSize/500.0;
                rx = MarkerSize/gcf()->window_w;
                ry = MarkerSize/gcf()->window_h;


                glDisable(GL_LINE_STIPPLE);
                gl2psDisable(GL2PS_LINE_STIPPLE);

                if(Marker == ".") {//.
                    for(unsigned int i = 0; i < XData.size(); ++i) {
                        xx = ctx(XData[i]);
                        yy = cty(YData[i]);

                        glPointSize(LineWidth);
                        glBegin(GL_POINTS);
                            glVertex2d(xx,yy);
                        glEnd();
                    }
                }
                else if(Marker == "+") {//+
                    for(unsigned int i = 0; i < XData.size(); ++i) {
                        xx = ctx(XData[i]);
                        yy = cty(YData[i]);

                        glBegin(GL_LINE_STRIP);
                            glVertex2d(xx-rx,yy);
                            glVertex2d(xx+rx,yy);
                        glEnd();

                        glBegin(GL_LINE_STRIP);
                            glVertex2d(xx,yy-ry);
                            glVertex2d(xx,yy+ry);
                        glEnd();
                    }
                }
                else if(Marker == "x") {//x
                    for(unsigned int i = 0; i < XData.size(); ++i) {
                        xx = ctx(XData[i]);
                        yy = cty(YData[i]);

                        glBegin(GL_LINE_STRIP);
                            glVertex2d(xx-rx,yy-ry);
                            glVertex2d(xx+rx,yy+ry);
                        glEnd();
                        glBegin(GL_LINE_STRIP);
                            glVertex2d(xx+rx,yy-ry);
                            glVertex2d(xx-rx,yy+ry);
                        glEnd();
                    }
                }
                else if(Marker == "d") {//d diamond
                    for(unsigned int i = 0; i < XData.size(); ++i) {
                        xx = ctx(XData[i]);
                        yy = cty(YData[i]);

                        glBegin(GL_LINE_LOOP);
                            glVertex2d(xx,   yy+ry);
                            glVertex2d(xx+rx,yy);
                            glVertex2d(xx,   yy-ry);
                            glVertex2d(xx-rx,yy);
                        glEnd();
                    }
                }
                else if(Marker == "^") {//^
                    for(unsigned int i = 0; i < XData.size(); ++i) {
                        xx = ctx(XData[i]);
                        yy = cty(YData[i]);

                        glBegin(GL_LINE_LOOP);
                            glVertex2d(xx,   yy+ry);
                            glVertex2d(xx+rx,yy-ry);
                            glVertex2d(xx-rx,yy-ry);
                        glEnd();
                    }
                }
                else if(Marker == "v") {//v
                    for(unsigned int i = 0; i < XData.size(); ++i) {
                        xx = ctx(XData[i]);
                        yy = cty(YData[i]);

                        glBegin(GL_LINE_LOOP);
                            glVertex2d(xx,   yy-ry);
                            glVertex2d(xx+rx,yy+ry);
                            glVertex2d(xx-rx,yy+ry);
                        glEnd();
                    }
                }
                else if(Marker == "o") {//o
                    for(unsigned int i = 0; i < XData.size(); ++i) {
                        xx = ctx(XData[i]);
                        yy = cty(YData[i]);

                        glBegin(GL_LINE_LOOP);
                            for(int i=0;i<20;i++) {
                                glVertex2d(xx+rx*cos(2*M_PI*(double)i/(double)(20)),
                                       yy+ry*sin(2*M_PI*(double)i/(double)(20)));
                            }
                        glEnd();
                    }
                }
                else if(Marker == "s") {//s :squire
                    for(unsigned int i = 0; i < XData.size(); ++i) {
                        xx = ctx(XData[i]);
                        yy = cty(YData[i]);

                        glBegin(GL_LINE_LOOP);
                            glVertex2d(xx-rx,yy-ry);
                            glVertex2d(xx-rx,yy+ry);
                            glVertex2d(xx+rx,yy+ry);
                            glVertex2d(xx+rx,yy-ry);
                        glEnd();
                    }
                }
                else { // if(Marker == "*") {//*
                    for(unsigned int i = 0; i < XData.size(); ++i) {
                        xx = ctx(XData[i]);
                        yy = cty(YData[i]);

                        glBegin(GL_LINE_STRIP);
                            glVertex2d(xx-rx,yy);
                            glVertex2d(xx+rx,yy);
                        glEnd();
                        glBegin(GL_LINE_STRIP);
                            glVertex2d(xx,yy-ry);
                            glVertex2d(xx,yy+ry);
                        glEnd();
                        glBegin(GL_LINE_STRIP);
                            glVertex2d(xx-rx,yy-ry);
                            glVertex2d(xx+rx,yy+ry);
                        glEnd();
                        glBegin(GL_LINE_STRIP);
                            glVertex2d(xx+rx,yy-ry);
                            glVertex2d(xx-rx,yy+ry);
                        glEnd();
                    }
                }
            }// Marker


            if(Errorbar) {// Errorbar //
                float xx,yy,yyp,yym;// transformed coordination

                glDisable(GL_LINE_STIPPLE);
                gl2psDisable(GL2PS_LINE_STIPPLE);
                //r=MarkerSize/500;

                for(unsigned int i = 0; i < XData.size(); ++i) {
                    xx  = ctx(XData[i]);
                    yy  = cty(YData[i]);
                    yyp = cty(YData[i] + YPData[i]);
                    yym = cty(YData[i] - YMData[i]);

                    glBegin(GL_LINE_STRIP);
                        glVertex2d(xx,yyp);
                        glVertex2d(xx,yym);
                    glEnd();
                    glBegin(GL_LINE_STRIP);
                        glVertex2d(xx-rx,yy);
                        glVertex2d(xx+rx,yy);
                    glEnd();
                    glBegin(GL_LINE_STRIP);
                        glVertex2d(xx-rx,yyp);
                        glVertex2d(xx+rx,yyp);
                    glEnd();
                    glBegin(GL_LINE_STRIP);
                        glVertex2d(xx-rx,yym);
                        glVertex2d(xx+rx,yym);
                    glEnd();
                }
            }//Errorbar
            //TODO:selection of error bar type
        }//2D

        // 3D //
        else {
            glBegin(GL_LINE_STRIP);
            for(unsigned int i = 0; i < XData.size(); ++i) {
                glVertex3d(ct3x(XData[i]),
                       ct3y(YData[i]),
                       ct3z(ZData[i]));
                }
            glEnd();
        }
    }

    std::pair<double, double> Line::min_max(const dvec& data, math::scale scale_ = math::linear_scale) {
        boost::mutex::scoped_lock l(data_mutex);
        std::pair<double, double> mm(std::numeric_limits<double>::max(), -std::numeric_limits<double>::max());
        if(scale_ == math::linear_scale) {//linear_scale
            for(dvec::const_iterator i = data.begin(); i != data.end(); ++i) {
                if(mm.first > *i) {
                    mm.first = *i;
                }
                if(mm.second < *i) {
                    mm.second = *i;
                }
            }
        } else {//log
            for(dvec::const_iterator i = data.begin(); i != data.end(); ++i) {
                if(mm.first > *i && *i > 0) {
                    mm.first = *i;
                }
                if(mm.second < *i) {
                    mm.second = *i;
                }
            }
        }
        return mm;
    }

    void Line::config() {
        std::pair<double, double> mm;
        mm = min_max(XData, ca->XScale);
        ca->xmin = std::min(mm.first, ca->xmin);
        ca->xmax = std::max(mm.second,ca->xmax);
        mm = min_max(YData, ca->YScale);
        ca->ymin = std::min(mm.first, ca->ymin);
        ca->ymax = std::max(mm.second,ca->ymax);
        mm = min_max(ZData, ca->ZScale);
        ca->zmin = std::min(mm.first, ca->zmin);
        ca->zmax = std::max(mm.second,ca->zmax);
    }

    line_t Line::set_capacity(unsigned int a) {
        if(a < max_capacity) {
            if(XData.size() > a) XData.erase(XData.begin(), XData.end() - a);
            if(YData.size() > a) YData.erase(YData.begin(), YData.end() - a);
            if(ZData.size() > a) ZData.erase(ZData.begin(), ZData.end() - a);
        }
        max_capacity = a;
        return shared_from_this();
    }

    line_t Line::stop_at_max(bool s) {
        stop_at_max_ = s;
        return shared_from_this();
    }


    /// vertex
    void  Line::vertex(const double x, const double y) {
        if(XData.size() == max_capacity && stop_at_max_) return;
        boost::mutex::scoped_lock l(data_mutex);
        if(ca->xmin > x) { ca->xmin = x; }
        if(ca->xmax < x) { ca->xmax = x; }
        if(ca->ymin > y) { ca->ymin = y; }
        if(ca->ymax < y) { ca->ymax = y; }
        XData.push_back(x);
        if(XData.size() > max_capacity) XData.erase(XData.begin(), XData.end()-max_capacity);
        YData.push_back(y);
        if(YData.size() > max_capacity) YData.erase(YData.begin(), YData.end()-max_capacity);
    }

    line_t Line::line(const dvec& x, const dvec& y) {
        boost::mutex::scoped_lock l(data_mutex);
        XData = x;
        YData = y;
        return shared_from_this();
    }
    line_t Line::line(const dvec& x, const dvec& y, const dvec& z) {
        boost::mutex::scoped_lock l(data_mutex);
        XData = x;
        YData = y;
        ZData = z;
        return shared_from_this();
    }
    /// plot, semilogx, semilogy, loglog
    line_t Line::plot(const dvec& y) {
        int n = y.size();
        dvec x(n);
        for(int i = 0; i < n; ++i) { x[i] = 1.0*i/(n-1); }
        return line(x,y);
    }
    line_t Line::plot(const dvec& x, const dvec& y) {
        return line(x,y);
    }
    line_t Line::plot(const dvec& x, const dvec& y, const dvec& z) {
        ca->type = axes_t_t::_3D;
        return line(x,y,z);
    }

    line_t Line::plot(std::valarray<double> x, std::valarray<double> y) {
        dvec xx,yy;
        for(unsigned int i = 0; i < x.size(); ++i) { xx.push_back(x[i]); }
        for(unsigned int i = 0; i < y.size(); ++i) { yy.push_back(y[i]); }
        return line(xx,yy);
    }
    line_t Line::semilogx(const dvec& x, const dvec& y) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->XScale = math::logarithmic_scale;
        XData = x;
        YData = y;
        return shared_from_this();
    }
    line_t Line::semilogy(const dvec& x, const dvec& y) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->YScale = math::logarithmic_scale;
        XData = x;
        YData = y;
        return shared_from_this();
    }
    line_t Line::loglog(const dvec& x, const dvec& y) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->XScale = math::logarithmic_scale;
        ca->YScale = math::logarithmic_scale;
        XData = x;
        YData = y;
        return shared_from_this();
    }
    /// errorbar
    void Line::vertex(const double x, const double y, const double ep, const double em) {//for errorbar
        if(XData.size() == max_capacity && stop_at_max_) return;
        boost::mutex::scoped_lock l(data_mutex);
        if(ca->xmin>x) { ca->xmin=x; }
        if(ca->xmax < x) { ca->xmax = x; }
        if(ca->ymin > y+ep) { ca->ymin = y+ep; }
        if(ca->ymax < y-em) { ca->ymax = y-em; }
        XData.push_back(x);
        if(XData.size() > max_capacity) XData.erase(XData.begin(), XData.end()-max_capacity);
        YData.push_back(y);
        if(YData.size() > max_capacity) YData.erase(YData.begin(), YData.end()-max_capacity);
        YPData.push_back(ep);
        if(YPData.size() > max_capacity) YPData.erase(XData.begin(), YPData.end()-max_capacity);
        YMData.push_back(em);
        if(YMData.size() > max_capacity) YMData.erase(YMData.begin(), YMData.end()-max_capacity);
    }
    line_t Line::errorbar(const dvec& x, const dvec& y,dvec e) {
        for(unsigned int i = 0; i < x.size(); ++i) { vertex(x[i],y[i],e[i],e[i]); }
        Errorbar = 1;
        return shared_from_this();
    }
    line_t Line::errorbar(const dvec& x, const dvec& y,dvec ep, const dvec& em) {
        for(unsigned int i = 0; i < x.size(); ++i) { vertex(x[i],y[i],ep[i],em[i]); }
        Errorbar = 1;
        return shared_from_this();
    }
    /// 3D line
    void Line::vertex(const double x, const double y, const double z) {
        if(XData.size() == max_capacity && stop_at_max_) return;
        boost::mutex::scoped_lock l(data_mutex);
        if(ca->xmin > x) { ca->xmin = x; }
        if(ca->xmax < x) { ca->xmax = x; }
        if(ca->ymin > y) { ca->ymin = y; }
        if(ca->ymax < y) { ca->ymax = y; }
        if(ca->zmin > z) { ca->zmin = z; }
        if(ca->zmax < z) { ca->zmax = z; }
        XData.push_back(x);
        if(XData.size() > max_capacity) XData.erase(XData.begin(), XData.end()-max_capacity);
        YData.push_back(y);
        if(YData.size() > max_capacity) YData.erase(YData.begin(), YData.end()-max_capacity);
        ZData.push_back(z);
        if(ZData.size() > max_capacity) ZData.erase(ZData.begin(), ZData.end()-max_capacity);
    }

    line_t Line::set(const float v) {
        LineWidth = v;
        MarkerSize = v;
        return shared_from_this();
    }

    line_t Line::set(const std::string p, const std::string v) {
             if(p == "COLOR") { Color = v; MarkerEdgeColor = v; }
        else if(p == "Color") { Color = v; }
        else if(p == "Marker") { Marker = v; }
        else if(p == "LineStyle") { LineStyle = v; }
        else if(p == "MarkerEdgeColor") { MarkerEdgeColor = v; }
        else if(p == "MarkerFaceColor") { MarkerFaceColor = v; }

        return shared_from_this();
    }

    line_t Line::set(const std::string p, const float v) {
             if(p == "LineWidth") { LineWidth  = v; }
        else if(p == "MarkerSize") { MarkerSize = v; }

        return shared_from_this();
    }

    line_t Line::set(const std::string v) {
             if( v == "k" ) { set("COLOR", "k"); }
        else if( v == "r" ) { set("COLOR", "r"); }
        else if( v == "b" ) { set("COLOR", "b"); }
        else if( v == "g" ) { set("COLOR", "g"); }
        else if( v == "c" ) { set("COLOR", "c"); }
        else if( v == "m" ) { set("COLOR", "m"); }
        else if( v == "y" ) { set("COLOR", "y"); }
        else if( v == "w" ) { set("COLOR", "w"); }

        else if( v == "dr" ) { set("COLOR", "dr"); }
        else if( v == "db" ) { set("COLOR", "db"); }
        else if( v == "dg" ) { set("COLOR", "dg"); }
        else if( v == "dc" ) { set("COLOR", "dc"); }
        else if( v == "dm" ) { set("COLOR", "dm"); }
        else if( v == "dy" ) { set("COLOR", "dy"); }

        else if( v == "lr" ) { set("COLOR", "lr"); }
        else if( v == "lb" ) { set("COLOR", "lb"); }
        else if( v == "lg" ) { set("COLOR", "lg"); }
        else if( v == "lc" ) { set("COLOR", "lc"); }
        else if( v == "lm" ) { set("COLOR", "lm"); }
        else if( v == "ly" ) { set("COLOR", "ly"); }

        else if( v == "ur" ) { set("COLOR", "ur"); }
        else if( v == "ub" ) { set("COLOR", "ub"); }
        else if( v == "ug" ) { set("COLOR", "ug"); }
        else if( v == "uy" ) { set("COLOR", "uy"); }
        else if( v == "uc" ) { set("COLOR", "uc"); }
        else if( v == "up" ) { set("COLOR", "up"); }
        else if( v == "uo" ) { set("COLOR", "uo"); }
        else if( v == "um" ) { set("COLOR", "um"); }
        else if( v == "ubr" ) { set("COLOR", "ubr"); }

        else if( v == "-"  ) { set("LineStyle", "-");   set("Marker", "none"); }
        else if( v == "- -") { set("LineStyle", "- -"); set("Marker", "none"); }
        else if( v == ":"  ) { set("LineStyle", ":");   set("Marker", "none"); }
        else if( v == "-." ) { set("LineStyle", "-.");  set("Marker", "none"); }

        else if( v == "." ) { set("Marker", "."); set("LineStyle", "none"); }
        else if( v == "+" ) { set("Marker", "+"); set("LineStyle", "none"); }
        else if( v == "x" ) { set("Marker", "x"); set("LineStyle", "none"); }
        else if( v == "d" ) { set("Marker", "d"); set("LineStyle", "none"); }
        else if( v == "^" ) { set("Marker", "^"); set("LineStyle", "none"); }
        else if( v == "v" ) { set("Marker", "v"); set("LineStyle", "none"); }
        else if( v == "o" ) { set("Marker", "o"); set("LineStyle", "none"); }
        else if( v == "*" ) { set("Marker", "*"); set("LineStyle", "none"); }
        else if( v == "s" ) { set("Marker", "s"); set("LineStyle", "none"); }
        else if( v == ">" ) { set("Marker", ">"); set("LineStyle", "none"); }
        else if( v == "<" ) { set("Marker", "<"); set("LineStyle", "none"); }
        else if( v == "p" ) { set("Marker", "p"); set("LineStyle", "none"); }
        else if( v == "h" ) { set("Marker", "h"); set("LineStyle", "none"); }

        return shared_from_this();
    }
}
