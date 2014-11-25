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

namespace cpplot {
    void Patch::draw() {
        if(type == _2D) { draw2d(); }
        if(type == _3D) { draw3d(); }
    }
    void Patch::clear() {
        boost::mutex::scoped_lock l(data_mutex);
        XData.clear();
        YData.clear();
        ZData.clear();
        CData.clear();
        faces.clear();
        vertices.clear();
    }

    void Patch::draw2d() {
        boost::mutex::scoped_lock l(data_mutex);

        std::vector<float> v(3);
        std::vector<int> f(3);
        float x,y;
        for(unsigned int i = 0; i < faces.size(); ++i) {
            f = faces[i];
            glBegin(GL_TRIANGLES);
                x = ctx(vertices[f[0]][0]);
                y = cty(vertices[f[0]][1]);
                glVertex2d(x,y);
                    x=ctx(vertices[f[1]][0]);
                    y=cty(vertices[f[1]][1]);
                glVertex2d(x,y);
                    x = ctx(vertices[f[2]][0]);
                    y = cty(vertices[f[2]][1]);
                glVertex2d(x,y);
            glEnd();
        }


        // XYZ Data //
        int nf, nv; //number of faces and vertex
        nf=XData.size();
        std::vector<float> rgb;

        for(int i = 0; i < nf; ++i) {
            nv = XData[i].size();

            // Edge
            if(EdgeColor!="none") {

                glLineWidth(LineWidth);
                gl2psLineWidth(LineWidth);

                rgb = ColorSpec2RGB(EdgeColor);
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glBegin(GL_LINE_LOOP);
                    for(int iv=0;iv<nv;++iv) {
                        glVertex2d( ctx(XData[i][iv]),
                                    cty(YData[i][iv]) );
                    }
                glEnd();
            }

            // Face
            if(FaceColor!="none") {
                rgb = ColorSpec2RGB(FaceColor);
                glColor3f(rgb[0], rgb[1], rgb[2]);

                if(CData.size()) {
                    rgb = CData[i];
                    glColor3d(rgb[0], rgb[1], rgb[2]);
                }

                glBegin(GL_POLYGON);
                    for(int iv = 0; iv < nv; ++iv) {
                        glVertex2d( ctx(XData[i][iv]),
                                    cty(YData[i][iv]) );
                    }
                glEnd();

            }
        }
    }


    void Patch::draw3d() {
        boost::mutex::scoped_lock l(data_mutex);
        // XYZ Data //
        unsigned int nf, nv; //number of faces and vertex
        nf = XData.size();
        std::vector<float> rgb;

        for(unsigned int i = 0; i < nf; ++i) {
            nv = XData[i].size();

            // Edge
            if(EdgeColor!="none") {
                glLineWidth(LineWidth);
                gl2psLineWidth(LineWidth);

                rgb = ColorSpec2RGB(EdgeColor);
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glBegin(GL_LINE_LOOP);
                    for(unsigned int iv = 0; iv < nv; ++iv) {
                        glVertex3d( ct3x(XData[i][iv]),
                                    ct3y(YData[i][iv]),
                                    ct3z(ZData[i][iv]) );
                    }
                glEnd();
            }

            // Face
            if(FaceColor!="none") {
                rgb = ColorSpec2RGB(FaceColor);
                glColor3f(rgb[0], rgb[1], rgb[2]);

                if(CData.size() > 0) {
                    rgb = CData[i];
                    glColor3d(rgb[0], rgb[1], rgb[2]);
                }

                glBegin(GL_POLYGON);
                    for(unsigned int iv = 0; iv < nv; ++iv) {
                        glVertex3d( ct3x(XData[i][iv]),
                                    ct3y(YData[i][iv]),
                                    ct3z(ZData[i][iv]));
                    }
                glEnd();

            }
        }
    }

    /// bar

    patch_t Patch::bar(const dvec& y, float width) {
        dvec x(y.size());
        for(unsigned int i = 0; i < y.size(); ++i) {
            x[i] = (double)(1+i);
        }
        return bar(x,y,width);
    }

    patch_t Patch::bar(const dvec& x, const dvec& y, const float width) {
        boost::mutex::scoped_lock l(data_mutex);
        type = _2D;
        XData.clear(); YData.clear(); ZData.clear();

        double wx = width*( math::max(x) - math::min(x) ) / x.size();

        dvec X(4),Y(4);
        for(unsigned int i = 0; i < x.size(); ++i) {
            X[0] = x[i] - wx/2.0; Y[0] = 0;
            X[1] = x[i] + wx/2.0; Y[1] = 0;
            X[2] = x[i] + wx/2.0; Y[2] = y[i];
            X[3] = x[i] - wx/2.0; Y[3] = y[i];
            XData.push_back(X);
            YData.push_back(Y);
        }

        return shared_from_this();
    }


    /// patch
    patch_t Patch::patch(const dmat& X, const dmat& Y) {
        boost::mutex::scoped_lock l(data_mutex);
        // Single color
        type = _2D;

        XData = X;
        YData = Y;
        ZData.clear();
        CData.clear();

        return shared_from_this();
    }
    patch_t Patch::patch(const dmat& X, const dmat& Y, const dvec& C) {
        boost::mutex::scoped_lock l(data_mutex);
        // One color per face with index color
        type = _2D;
        XData = X;
        YData = Y;
        ZData.clear();
        CData = Index2TrueColor(C);

        return shared_from_this();
    }
    patch_t Patch::patch(const dmat& X, const dmat& Y, const tcvec& C) {
        boost::mutex::scoped_lock l(data_mutex);
        // One color per face with true color
        type = _2D;
        XData = X;
        YData = Y;
        ZData.clear();
        CData = C;

        return shared_from_this();
    }
    patch_t Patch::patch(const dmat& X, const dmat& Y, const dmat& Z) {
        boost::mutex::scoped_lock l(data_mutex);
        // Single color
        ca->type = axes_t_t::_3D;
        type = _3D;
        XData = X;
        YData = Y;
        ZData = Z;
        CData.clear();

        return shared_from_this();
    }
    patch_t Patch::patch(const dmat& X, const dmat& Y, const dmat& Z, const dvec& C) {
        boost::mutex::scoped_lock l(data_mutex);
        // One color per face
        ca->type = axes_t_t::_3D;
        type = _3D;

        XData = X;
        YData = Y;
        ZData = Z;
        CData = Index2TrueColor(C);

        return shared_from_this();
    }
    patch_t Patch::patch(const dmat& X, const dmat& Y, const dmat& Z, const tcvec& C) {
        boost::mutex::scoped_lock l(data_mutex);
        // One color per face
        ca->type = axes_t_t::_3D;
        type = _3D;

        XData = X;
        YData = Y;
        ZData = Z;
        CData = C;

        return shared_from_this();
    }

    void Patch::config() {
        boost::mutex::scoped_lock l(data_mutex);
        ca->xmax = std::max( math::max(XData), ca->xmax );
        ca->xmin = std::min( math::min(XData), ca->xmin );
        ca->ymax = std::max( math::max(YData), ca->ymax );
        ca->ymin = std::min( math::min(YData), ca->ymin );
        ca->zmax = std::max( math::max(ZData), ca->zmax );
        ca->zmin = std::min( math::min(ZData), ca->zmin );
    }

    patch_t Patch::set(const std::string p, const std::string v) {
             if(p == "COLOR")    { EdgeColor = v; }
        else if(p == "LineStyle") { LineStyle = v; }
        else if(p == "EdgeColor") { EdgeColor = v; }
        else if(p == "FaceColor") { FaceColor = v; }

        return shared_from_this();
    }

    patch_t Patch::set(const std::string p, const float v) {
        if(p == "LineWidth") { LineWidth = v; }
        return shared_from_this();
    }

    tcvec Patch::Index2TrueColor(const dvec& IC) {
        if(ca->CLim[0] == ca->CLim[1]) {
            ca->CLim[0] = math::min(IC);
            ca->CLim[1] = math::max(IC);
            //~ ca->CLim[0] = std::min( math::min(IC), math::min(IC) ); ???????WHAT????
            //~ ca->CLim[1] = std::max( math::max(IC), math::max(IC) );
        }
        std::vector<float> rgb;
        tcvec tc;
        for(unsigned int j = 0; j < IC.size(); ++j) {
            rgb = ca->map2color(IC[j]);
            tc.push_back(rgb);
        }
        return tc;
    }
}
