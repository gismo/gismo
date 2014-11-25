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
#include <deque>

namespace cpplot {
    void Surface::draw() {
        //if(type==0) { display_pcolor(); }
        if(type == contourplot) {
            draw_contour();
        }
        else if(type == _3D) {
            draw3d();
        }
        else {
            draw2d();
        }
    }

    void Surface::clear() {
        boost::mutex::scoped_lock l(data_mutex);
        XData.clear();
        YData.clear();
        ZData.clear();
        CDataIndex.clear();
        CData.clear();
        V.clear();

    }

    void Surface::draw2d() {
        boost::mutex::scoped_lock l(data_mutex);
        unsigned int nxi,nxj,nyi,nyj,nzi,nzj;
        std::vector<float> rgb;
        nxi = XData.size(); if(nxi) { nxj = XData[0].size(); }
        nyi = YData.size(); if(nyi) { nyj = YData[0].size(); }
        nzi = ZData.size(); if(nzi) { nzj = ZData[0].size(); }

        //printf("%s %s:%d  (%d %d) (%d %d) (%d %d) \n", __func__, __FILE__, __LINE__,nxi,nxj,nyi,nyj,nzi,nzj);

        // (Z) // do not use
        if(nxi == 0) {
            printf("%s %s:%d\n", __func__, __FILE__, __LINE__);
            // Edge
            if(EdgeColor != "none") {
                glLineWidth(LineWidth);
                rgb = ColorSpec2RGB(EdgeColor);
                glColor3d(rgb[0], rgb[1], rgb[2]);

                for(unsigned int i = 0; i < nzi; ++i) {
                    glBegin(GL_LINE_STRIP); //TODO add more style
                        for(unsigned int j = 0; j < nzj; ++j) {
                            glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j)/(nzj-1)),
                                   cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i)/(nzi-1)));
                        }
                    glEnd();
                }
            }
            // Face
            if(FaceColor != "none") {
                for(unsigned int i = 0; i < nzi-1; ++i) {
                    for(unsigned int j = 0; j < nzj-1; ++j) {

                        rgb = ColorSpec2RGB(FaceColor);
                        if(FaceColor == "flat") { rgb = CData[i][j]; }
                        glColor3f(rgb[0], rgb[1], rgb[2]);

                        glBegin(GL_QUADS);
                            glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j  )/(float)(nzj-1)),
                                   cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i  )/(float)(nzi-1)));
                            glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j+1)/(float)(nzj-1)),
                                   cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i  )/(float)(nzi-1)));
                            glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j+1)/(float)(nzj-1)),
                                   cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i+1)/(float)(nzi-1)));
                            glVertex2d(ctx(ca->XLim[0]+(ca->XLim[1]-ca->XLim[0])*(float)(j  )/(float)(nzj-1)),
                                   cty(ca->YLim[0]+(ca->YLim[1]-ca->YLim[0])*(float)(i+1)/(float)(nzi-1)));
                        glEnd();
                    }
                }
            }//Face
        }

        // (x,y,Z) //
        else if(nxi==1) {
            //printf("%s %s:%d  %d %d \n", __func__, __FILE__, __LINE__,nxj,nyj);
            // Edge
            if(EdgeColor != "none") {
                glLineWidth(LineWidth);
                rgb=ColorSpec2RGB(EdgeColor);
                glColor3d(rgb[0], rgb[1], rgb[2]);

                for(unsigned int i = 0;i < nyj; ++i) {
                    glBegin(GL_LINE_STRIP); //TODO add more style
                        glVertex2d( ctx(XData[0][    0]),cty(YData[0][i]) );
                        glVertex2d( ctx(XData[0][nxj-1]),cty(YData[0][i]) );
                    glEnd();
                }

                for(unsigned int j = 0; j < nxj; ++j) {
                    glBegin(GL_LINE_STRIP);
                        glVertex2d( ctx(XData[0][j]),cty(YData[0][0])    );
                        glVertex2d( ctx(XData[0][j]),cty(YData[0][nyj-1]));
                    glEnd();
                }
            }
            // Face
            if(FaceColor != "none") {
                for(unsigned int i = 0; i < nyj-1; ++i) {
                    for(unsigned int j = 0; j < nxj-1; ++j) {
                        // color
                        rgb = ColorSpec2RGB(FaceColor);
                        if(FaceColor == "flat") { rgb = CData[i][j]; }
                        glColor3f(rgb[0], rgb[1], rgb[2]);

                        glBegin(GL_QUADS);
                            glVertex2d( ctx(XData[0][j]),
                                        cty(YData[0][i]));
                            glVertex2d( ctx(XData[0][j]),
                                        cty(YData[0][i+1]));
                            glVertex2d( ctx(XData[0][j+1]),
                                        cty(YData[0][i+1]));
                            glVertex2d( ctx(XData[0][j+1]),
                                        cty(YData[0][i]));
                        glEnd();
                    }
                }
            }
        }//nxi==1

        // (X,Y,C) //
        else if(nxi>1) {
            // Edge
            //printf("%s %s:%d\n", __func__, __FILE__, __LINE__);
            if(EdgeColor != "none") {
                glLineWidth(LineWidth);
                rgb = ColorSpec2RGB(EdgeColor);
                glColor3d(rgb[0], rgb[1], rgb[2]);

                for(unsigned int i = 0; i < nxi; ++i) {
                    glBegin(GL_LINE_STRIP); //TODO add more style
                        for(unsigned int j = 0; j < nxj; ++j) {
                            glVertex2d(ctx(XData[i][j]),
                                   cty(YData[i][j]));
                        }
                    glEnd();
                }
                for(unsigned int j = 0; j < nxi; ++j) {
                    glBegin(GL_LINE_STRIP);
                        for(unsigned int i = 0; i < nxj; ++i) {
                            glVertex2d(ctx(XData[i][j]),
                                   cty(YData[i][j]));
                        }
                    glEnd();
                }
            }
            // Face
            if(FaceColor != "none") {
                for(unsigned int i = 0; i < nxi-1; ++i) {
                    for(unsigned int j = 0; j < nxj-1; ++j) {
                        // color
                        rgb = ColorSpec2RGB(FaceColor);
                        if(FaceColor == "flat") { rgb = CData[i][j]; }
                        glColor3f(rgb[0], rgb[1], rgb[2]);

                        glBegin(GL_QUADS);
                            glVertex2d( ctx(XData[i  ][j]),
                                        cty(YData[i  ][j]) );
                            glVertex2d( ctx(XData[i  ][j+1]),
                                        cty(YData[i  ][j+1]));
                            glVertex2d( ctx(XData[i+1][j+1]),
                                        cty(YData[i+1][j+1]));
                            glVertex2d( ctx(XData[i+1][j]),
                                        cty(YData[i+1][j]));
                        glEnd();
                    }
                }
            }
        }
    }

    /// display 3d
    void Surface::draw3d() {
        boost::mutex::scoped_lock l(data_mutex);
        std::vector<float> rgb;
        unsigned int ny = ZData.size();
        unsigned int nx = ZData[0].size();

        if(XData.size() == 1) {
            //(x,y,Z);
            //Face
            if(FaceColor != "none") {
                for(unsigned int i = 0; i < ny-1; ++i) {
                    for(unsigned int j = 0; j < nx-1; ++j) {
                        rgb = ColorSpec2RGB(FaceColor);
                        glColor3d(rgb[0], rgb[1], rgb[2]);
                        if(FaceColor=="flat") {
                            rgb = CData[i][j];
                            glColor3f(rgb[0], rgb[1], rgb[2]);
                        }
                        glBegin(GL_TRIANGLES);
                            glVertex3d( ct3x(XData[0][j]),
                                        ct3y(YData[0][i]),
                                        ct3z(ZData[i][j]));
                            glVertex3d( ct3x(XData[0][j]),
                                        ct3y(YData[0][i+1]),
                                        ct3z(ZData[i+1][j]));
                            glVertex3d( ct3x(XData[0][j+1]),
                                        ct3y(YData[0][i+1]),
                                        ct3z(ZData[i+1][j+1]));
                        glEnd();
                        glBegin(GL_TRIANGLES);
                            glVertex3d( ct3x(XData[0][j]),
                                        ct3y(YData[0][i]),
                                        ct3z(ZData[i][j]));
                            glVertex3d( ct3x(XData[0][j+1]),
                                        ct3y(YData[0][i]),
                                        ct3z(ZData[i][j+1]));
                            glVertex3d( ct3x(XData[0][j+1]),
                                        ct3y(YData[0][i+1]),
                                        ct3z(ZData[i+1][j+1]));
                        glEnd();
                    }
                }
            }
            // Edge
            if(EdgeColor != "none") {
                glLineWidth(LineWidth);
                rgb = ColorSpec2RGB(EdgeColor);
                glColor3d(rgb[0], rgb[1], rgb[2]);

                for(unsigned int i = 0; i < ny; ++i) {
                    glBegin(GL_LINE_STRIP);
                        for(unsigned int j = 0; j < nx; ++j) {
                            glVertex3d( ct3x(XData[0][j]),
                                        ct3y(YData[0][i]),
                                        ct3z(ZData[i][j]));
                        }
                    glEnd();
                }
                for(unsigned int j = 0; j < nx; ++j) {
                    glBegin(GL_LINE_STRIP);
                        for(unsigned int i = 0; i < ny; ++i) {
                            glVertex3d( ct3x(XData[0][j]),
                                        ct3y(YData[0][i]),
                                        ct3z(ZData[i][j]));
                        }
                    glEnd();
                }
            }
        }// (x,y,Z)

        else { // XData.size() > 1
            // (X,Y,Z) //

            //Face
            if(FaceColor != "none") {
                for(unsigned int i = 0; i < ny-1; ++i) {
                    for(unsigned int j = 0; j < nx-1; ++j) {
                        rgb = ColorSpec2RGB(FaceColor);
                        glColor3d(rgb[0], rgb[1], rgb[2]);
                        if(FaceColor=="flat") {
                            rgb = CData[i][j];
                            glColor3f(rgb[0],rgb[1],rgb[2]);
                        }

                        glBegin(GL_TRIANGLES);
                            glVertex3d( ct3x(XData[i][j]),
                                        ct3y(YData[i][j]),
                                        ct3z(ZData[i][j]));
                            glVertex3d( ct3x(XData[i][j+1]),
                                        ct3y(YData[i][j+1]),
                                        ct3z(ZData[i][j+1]));
                            glVertex3d( ct3x(XData[i+1][j+1]),
                                        ct3y(YData[i+1][j+1]),
                                        ct3z(ZData[i+1][j+1]));
                        glEnd();

                        glBegin(GL_TRIANGLES);
                            glVertex3d(ct3x(XData[i][j]),
                                   ct3y(YData[i][j]),
                                   ct3z(ZData[i][j]) );
                            glVertex3d(ct3x(XData[i+1][j]),
                                   ct3y(YData[i+1][j]),
                                   ct3z(ZData[i+1][j]) );
                            glVertex3d(ct3x(XData[i+1][j+1]),
                                   ct3y(YData[i+1][j+1]),
                                   ct3z(ZData[i+1][j+1]) );
                        glEnd();
                    }
                }
            }
            // Edge
            if(EdgeColor != "none") {
                glLineWidth(LineWidth);
                rgb = ColorSpec2RGB(EdgeColor);
                glColor3d(rgb[0], rgb[1], rgb[2]);

                for(unsigned int i = 0; i < ny; ++i) {
                    glBegin(GL_LINE_STRIP);
                    for(unsigned int j = 0; j < nx; ++j) {
                    glVertex3d( ct3x(XData[i][j]),
                                ct3y(YData[i][j]),
                                ct3z(ZData[i][j]));
                    }
                    glEnd();
                }
                for(unsigned int j = 0; j < nx; ++j) {
                    glBegin(GL_LINE_STRIP);
                    for(unsigned int i = 0; i < ny; ++i) {
                    glVertex3d( ct3x(XData[i][j]),
                                ct3y(YData[i][j]),
                                ct3z(ZData[i][j]));
                    }
                    glEnd();
                }
            }
        }//(X,Y,Z)
    }



    /// display contour
    void Surface::contourc(const dvec& x, const dvec& y, const dmat& Z, const dvec& v, dmat& C) {
        //Z(i,j), x(j),y(i)
        double x0, y0, z0;
        unsigned int ny = Z.size();
        unsigned int nx = Z[0].size();
        ContourPoint c;
        std::vector<ContourPoint> vc;
        std::deque<ContourPoint> ac;
        C.resize(2);

        for(unsigned int iv = 0; iv < v.size(); ++iv) {
            z0 = v[iv];

            // find contour points
            vc.clear();
            for(unsigned int i = 0; i < ny; ++i) {
                for(unsigned int j = 0; j < nx; ++j) {
                    if( ( j < nx-1 ) && ( (Z[i][j+1] - z0)*(Z[i][j] - z0) < 0 ) ) {
                        x0 = x[j] + ( x[j+1] - x[j] )*( z0 - Z[i][j] )/( Z[i  ][j+1]-Z[i][j] );
                        c.xj = j; c.yi = i; c.xy = 0; c.done = 0;
                        c.x = x0; c.y = y[i];
                        vc.push_back(c);
                    }
                    if( ( i < ny-1 ) && ( (Z[i+1][j] - z0)*(Z[i][j] - z0) < 0 ) ) {
                        y0 = y[i] + ( y[i+1] - y[i] )*( z0 - Z[i][j] )/( Z[i+1][j  ]-Z[i][j] );
                        c.xj = j; c.yi = i; c.xy = 1; c.done = 0;
                        c.x = x[j]; c.y = y0;
                        vc.push_back(c);
                    }
                }
            }
            // sort contour points
            int is = 0;
            unsigned int m, kk; // k
            int mode, mode_next;

            //k = 0;
            mode = 0;
            while( mode < 5 ) {

                if(mode == 0) {// set initial contour point
                    ac.clear();
                    is = 0; m = 0;
                    while( !is && (m < vc.size()) ) {
                        if(!vc[m].done) { is=1; kk=m; }
                        ++m;
                    }
                    if(is) {
                        vc[kk].done = 2;
                        c = vc[kk];
                        ac.push_back(c);
                        mode_next = 1;
                    }else{
                        mode_next = 5;
                    }
                }
                else if( mode == 1 || mode == 3 ) {//search next contour point
                    is = 0;
                    m  = 0;
                    while( !is && (m < vc.size()) ) {
                        is = 0;
                        if( (!vc[m].done) || ((vc[m].done==2)&&(ac.size()>2)) ) {
                            if( (c.xy == 0) && (vc[m].xy == 0) && (vc[m].xj == c.xj  ) && (vc[m].yi == c.yi-1) ) { is = 1; }
                            if( (c.xy == 0) && (vc[m].xy == 0) && (vc[m].xj == c.xj  ) && (vc[m].yi == c.yi+1) ) { is = 2; }
                            if( (c.xy == 0) && (vc[m].xy == 1) && (vc[m].xj == c.xj  ) && (vc[m].yi == c.yi  ) ) { is = 3; }
                            if( (c.xy == 0) && (vc[m].xy == 1) && (vc[m].xj == c.xj+1) && (vc[m].yi == c.yi  ) ) { is = 4; }
                            if( (c.xy == 0) && (vc[m].xy == 1) && (vc[m].xj == c.xj  ) && (vc[m].yi == c.yi-1) ) { is = 5; }
                            if( (c.xy == 0) && (vc[m].xy == 1) && (vc[m].xj == c.xj+1) && (vc[m].yi == c.yi-1) ) { is = 6; }
                            if( (c.xy == 1) && (vc[m].xy == 1) && (vc[m].xj == c.xj+1) && (vc[m].yi == c.yi  ) ) { is = 7; }
                            if( (c.xy == 1) && (vc[m].xy == 1) && (vc[m].xj == c.xj-1) && (vc[m].yi == c.yi  ) ) { is = 8; }
                            if( (c.xy == 1) && (vc[m].xy == 0) && (vc[m].xj == c.xj  ) && (vc[m].yi == c.yi  ) ) { is = 9; }
                            if( (c.xy == 1) && (vc[m].xy == 0) && (vc[m].xj == c.xj  ) && (vc[m].yi == c.yi+1) ) { is = 10; }
                            if( (c.xy == 1) && (vc[m].xy == 0) && (vc[m].xj == c.xj-1) && (vc[m].yi == c.yi  ) ) { is = 11; }
                            if( (c.xy == 1) && (vc[m].xy == 0) && (vc[m].xj == c.xj-1) && (vc[m].yi == c.yi+1) ) { is = 12; }
                        }
                        if(is) { kk = m; }
                        m++;
                    }
                    if(is) {
                        vc[kk].done = 1;
                        c = vc[kk];
                    }
                    if(mode == 1) {
                        if(is) {
                            ac.push_back(vc[kk]);
                            mode_next = 1;
                        } else {
                            mode_next=2;
                        }
                    }
                    else if(mode == 3) {
                        if(is) {
                            ac.push_front(vc[kk]);
                            mode_next = 3;
                        } else {
                            mode_next = 4;
                        }
                    }
                }
                else if(mode == 2) { // set first terminal
                    c = ac[0];
                    mode_next = 3;
                }
                else if(mode == 4) {// move the accumlated points
                    if(ac.size()) {
                        C[0].push_back(z0);
                        C[1].push_back(ac.size());
                        for(unsigned int i = 0; i < ac.size(); ++i) {
                            C[0].push_back( ac[i].x );
                            C[1].push_back( ac[i].y );
                        }
                    }
                    mode_next = 0;
                }

                mode = mode_next;
            }

        }//iv
    }


    void Surface::draw_contour() {
        boost::mutex::scoped_lock l(data_mutex);
        //vector<float> rgb;
        //int ny=ZData.size();
        //int nx=ZData[0].size();
        dmat C;
        double xx,yy;

        contourc(XData[0],YData[0],ZData,V, C);

        glDisable(GL_LINE_STIPPLE);
        gl2psDisable(GL2PS_LINE_STIPPLE);

        //glLineWidth(cl->LineWidth);
        glLineWidth(2);
        gl2psLineWidth(2);

        //rgb=CData[i][j];
        //glColor3f(rgb[0],rgb[1],rgb[2]);
        glColor3f(0,0,0);

        //TODO: adjust line color and properties
        unsigned int k = 0, nk;
        for(unsigned int i = 0; i < C[0].size(); ++i) {
            if(k == 0) {
                nk = (unsigned int)C[1][i];
                glBegin(GL_LINE_STRIP);
            } else{
                if(k <= nk) {
                    xx = ctx( C[0][i] );
                    yy = cty( C[1][i] );
                    glVertex2d(xx,yy);
                }
            }
            ++k;
            if(k > nk) {
                k = 0;
                glEnd();
            }
        }
    }

    void Surface::config() {
        boost::mutex::scoped_lock l(data_mutex);
        // check data size
        unsigned int nzi,nzj;

        nzi = ZData.size();
        if(nzi) { nzj = ZData[0].size(); }

        unsigned int nci,ncj;
        nci = CDataIndex.size();
        if(nci) { ncj = CDataIndex[0].size(); }

        // generate x and y data
        unsigned int nx = 0,ny = 0;
        if(nzi) { ny = nzi; nx=nzj; }
        if(nci) { ny = nci; nx=ncj; }

        //printf("%s %s:%d: %d %d %d %d \n", __func__, __FILE__, __LINE__,nzi,nci,nx,ny);
        if(XData.size() == 0) {
            XData.resize(1);
            XData[0] = math::linspace<std::vector<double> >(1.0, (double)nx, nx);
        }
        if(YData.size() == 0) {
            YData.resize(1);
            YData[0] = math::linspace<std::vector<double> >(1.0, (double)ny, ny);
        }

        // config data range
        ca->xmax = std::max(math::max(XData),ca->xmax);
        ca->xmin = std::min(math::min(XData),ca->xmin);
        ca->ymax = std::max(math::max(YData),ca->ymax);
        ca->ymin = std::min(math::min(YData),ca->ymin);
        ca->zmax = std::max(math::max(ZData),ca->zmax);
        ca->zmin = std::min(math::min(ZData),ca->zmin);

        // set CLim
        //note: first called surface effects CLim
        if(ca->CLim[0] == ca->CLim[1]) {
            ca->CLim[0] = math::min(CDataIndex);
            ca->CLim[1] = math::max(CDataIndex);
        }

        // CData !!!
        if( (CData.size()==0) && (CDataIndex.size()) ) {
            std::vector<float> rgb;
            //tcmat cdata(ny,nx);
            tcmat cdata(ny);

            for(unsigned int i = 0; i < ny; ++i) {
                cdata[i].resize(nx);
                for(unsigned int j = 0; j < nx; ++j) {
                    rgb = ca->map2color(CDataIndex[i][j]);
                    cdata[i][j] = rgb;
                }
            }
            CData = cdata;
        }

        // contour plot
        if(V.size() == 0) {
            if(NContour < 1) {
                NContour = 10;
            }
            V = math::linspace<std::vector<double> >(math::min(ZData),math::max(ZData),NContour);
        }
    }

    // Color
    surface_t Surface::shading(const std::string c) {
        if(c == "faceted") {
            EdgeColor="k";
        }
        else if(c == "flat") {
            EdgeColor="none";
        }
        return shared_from_this();
    }

        /// create surface
    surface_t Surface::surface(const dmat& Z) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        ZData = Z;
        CDataIndex = Z;
        CData.clear();
        return shared_from_this();
    }
    surface_t Surface::surface(const dmat& Z, const dmat& C) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        ZData = Z;
        CDataIndex = C;
        CData.clear();

        return shared_from_this();
    }
    surface_t Surface::surface(const dmat& Z, const tcmat& C) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        ZData = Z;
        CDataIndex.clear();
        CData = C;

        return shared_from_this();
    }
    surface_t Surface::surface(const dvec& x, const dvec& y, const dmat& Z) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData = Z;
        CDataIndex = Z;
        CData.clear();
        return shared_from_this();
    }
    surface_t Surface::surface(const dvec& x, const dvec& y, const dmat& Z, const dmat& C) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData = Z;
        CDataIndex = C;
        CData.clear();
        return shared_from_this();
    }
    surface_t Surface::surface(const dvec& x, const dvec& y, const dmat& Z, const tcmat& C) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData = Z;
        CDataIndex.clear();
        CData = C;

        return shared_from_this();
    }

    surface_t Surface::surface(const dmat& X, const dmat& Y, const dmat& Z) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        XData = X;
        YData = Y;
        ZData = Z;
        CDataIndex = Z;
        CData.clear();
        return shared_from_this();
    }
    surface_t Surface::surface(const dmat& X, const dmat& Y, const dmat& Z, const dmat& C) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        XData = X;
        YData = Y;
        ZData = Z;
        CDataIndex = C;
        CData.clear();
        return shared_from_this();
    }
    surface_t Surface::surface(const dmat& X, const dmat& Y, const dmat& Z, const tcmat& C) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        XData = X;
        YData = Y;
        ZData = Z;
        CDataIndex.clear();
        CData = C;
        return shared_from_this();
    }
    // surf
    surface_t Surface::surf(const dvec& x, const dvec& y, const dmat& Z) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData = Z;
        CDataIndex = Z;
        CData.clear();
        EdgeColor = "k";
        FaceColor = "flat";
        return shared_from_this();
    }

    // create pcolor
    surface_t Surface::pcolor(dmat C) {
        boost::mutex::scoped_lock l(data_mutex);
        type = _2D;
        XData.clear();
        YData.clear();
        ZData.clear();
        CDataIndex = C;
        CData.clear();
        return shared_from_this();
    }
    surface_t Surface::pcolor(const tcmat& C) {
        boost::mutex::scoped_lock l(data_mutex);
        type = _2D;
        XData.clear();
        YData.clear();
        ZData.clear();
        CDataIndex.clear();
        CData = C;
        return shared_from_this();
    }
    surface_t Surface::pcolor(const dvec& x, const dvec& y, const dmat& C) {
        boost::mutex::scoped_lock l(data_mutex);
        type = _2D;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData.clear();
        CDataIndex = C;
        CData.clear();
        return shared_from_this();
    }
    surface_t Surface::pcolor(const dvec& x, const dvec& y, const tcmat& C) {
        boost::mutex::scoped_lock l(data_mutex);
        type = _2D;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData.clear();
        CDataIndex.clear();
        CData = C;
        return shared_from_this();
    }
    surface_t Surface::pcolor(const dmat& X, const dmat& Y, const dmat& C) {
        type = _2D;
        XData = X;
        YData = Y;
        ZData.clear();
        CDataIndex = C;
        CData.clear();
        return shared_from_this();
    }
    surface_t Surface::pcolor(const dmat& X, const dmat& Y, const tcmat& C) {
        type = _2D;
        XData = X;
        YData = Y;
        ZData.clear();
        CDataIndex.clear();
        CData = C;
        return shared_from_this();
    }

    /// mesh
    surface_t Surface::mesh(const dvec& x, const dvec& y, dmat Z) {
        boost::mutex::scoped_lock l(data_mutex);
        ca->type = axes_t_t::_3D;
        type = _3D;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData = Z;
        CDataIndex = Z;
        CData.clear();
        EdgeColor = "k";
        FaceColor = "w";
        return shared_from_this();
    }

    /// contour
    surface_t Surface::contour(const dmat& Z) {
        boost::mutex::scoped_lock l(data_mutex);
        type = contourplot;
        XData.clear();
        YData.clear();
        ZData = Z;
        NContour = 0;
        V.clear();
        return shared_from_this();
    }
    surface_t Surface::contour(const dmat& Z, const int n) {
        boost::mutex::scoped_lock l(data_mutex);
        type = contourplot;
        XData.clear();
        YData.clear();
        ZData = Z;
        NContour = n;
        V.clear();
        return shared_from_this();
    }
    surface_t Surface::contour(const dmat& Z, const dvec& v) {
        boost::mutex::scoped_lock l(data_mutex);
        type = contourplot;
        XData.clear();
        YData.clear();
        ZData = Z;
        NContour = v.size();
        V = v;
        return shared_from_this();
    }
    surface_t Surface::contour(const dvec& x, const dvec& y, const dmat& Z) {
        boost::mutex::scoped_lock l(data_mutex);
        type = contourplot;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData = Z;
        NContour = 0;
        V.clear();
        return shared_from_this();
    }
    surface_t Surface::contour(const dvec& x, const dvec& y, const dmat& Z, const int n) {
        boost::mutex::scoped_lock l(data_mutex);
        type = contourplot;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData = Z;
        NContour = n;
        V.clear();
        return shared_from_this();
    }
    surface_t Surface::contour(const dvec& x, const dvec& y, const dmat& Z, const dvec& v) {
        boost::mutex::scoped_lock l(data_mutex);
        type = contourplot;
        XData.resize(1); XData[0] = x;
        YData.resize(1); YData[0] = y;
        ZData = Z;
        NContour = v.size();
        V = v;
        return shared_from_this();
    }


    surface_t Surface::set(const std::string p, const std::string v) {
             if(p == "COLOR")    { EdgeColor=v; }
        else if(p == "LineStyle") { LineStyle=v; }
        else if(p == "EdgeColor") { EdgeColor=v; }
        else if(p == "FaceColor") { FaceColor=v; }

        return shared_from_this();
    }

    surface_t Surface::set(const std::string p, const float v) {
        if(p == "LineWidth") { LineWidth = v; }

        return shared_from_this();
    }
}
