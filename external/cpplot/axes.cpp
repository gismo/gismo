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
#include <limits>

namespace cpplot {
    axes_t_t::axes_t_t(layer_t l)
        :   layer(l),
            cta0(30),
            phi0(30),
            cta(cta0),
            phi(cta0),
            Mouse(false),
            mouse_callback(NULL),
            xmin(std::numeric_limits<double>::max()),    xmax(-std::numeric_limits<double>::max()),
            ymin(std::numeric_limits<double>::max()),    ymax(-std::numeric_limits<double>::max()),
            zmin(std::numeric_limits<double>::max()),    zmax(-std::numeric_limits<double>::max()),
            type(_2D),
            Box(true),
            GridLineStyle(":"),
            LineWidth(1),
            TickDir("in"),
            visible(true),
            XGrid(false),
            YGrid(false),
            ZGrid(false),
            selected(false),
            XAxisLocation("bottom"),
            YAxisLocation("left"),
            XLimMode(automatic),
            YLimMode(automatic),
            ZLimMode(automatic),
            XScale(math::linear_scale),
            YScale(math::linear_scale),
            ZScale(math::linear_scale),
            TickLabel(true)
    {
        position[0] = 0.13;
        position[1] = 0.11;
        position[2] = 0.775;
        position[3] = 0.815;

        viewport3d[0] = 0.0;
        viewport3d[1] = 0.0;
        viewport3d[2] = 1.0;
        viewport3d[3] = 1.0;

        camera_position[0] = 1; camera_position[1] = 1; camera_position[2] = 1;
        camera_target[0]   = 0.;camera_target[1]   = 0; camera_target[2]   = 0;
        camera_up_vector[0] = 0; camera_up_vector[1] = 0; camera_up_vector[2] = 1;

        XLim[0] = 0;    XLim[1] = 10;
        YLim[0] = 0;    YLim[1] = 10;
        ZLim[0] = 0;    ZLim[1 ]= 10;




        CLim[0] = 0; CLim[1] = 0;
        jet();
    };

    int axes_t_t::window_h() { return layer->figure->window_h; }
    int axes_t_t::window_w() { return layer->figure->window_w; }

    void axes_t_t::reset_limits() {
        xmin = std::numeric_limits<double>::max(); xmax = -std::numeric_limits<double>::max();
        ymin = std::numeric_limits<double>::max(); ymax = -std::numeric_limits<double>::max();
        zmin = std::numeric_limits<double>::max(); zmax = -std::numeric_limits<double>::max();
    }
    void axes_t_t::config() {
        reset_limits();
        for(drawings_t::iterator d = children.begin(); d != children.end(); ++d) {
            (*d)->config();
        }
        float extent = 0, extent_linear = 0.03;
        if((XLimMode == automatic) && (xmax > xmin)) {
            extent = (XScale == math::linear_scale) ? extent_linear : 0;
            XLim[0] = xmin - extent*(xmax - xmin);
            XLim[1] = xmax + extent*(xmax - xmin);
        }
        if((YLimMode == automatic) && (ymax > ymin)) {
            extent = (YScale == math::linear_scale) ? extent_linear : 0;
            YLim[0] = ymin - extent*(ymax - ymin);
            YLim[1] = ymax + extent*(ymax - ymin);
        }
        if((ZLimMode == automatic) && (zmax > zmin)) {
            extent = (ZScale == math::linear_scale) ? extent_linear : 0;
            ZLim[0] = zmin - extent*(zmax - zmin);
            ZLim[1] = zmax + extent*(zmax - zmin);
        }
        //printf("Z: %d,%f,%f\n",ZLimMode,ZLim[0],ZLim[1]);
        //if(num_child) { visible=1; }else{ visible=0; }

        XTick = make_tick(XLim[0], XLim[1]);
        YTick = make_tick(YLim[0], YLim[1]);
        ZTick = make_tick(ZLim[0], ZLim[1]);
    }

    bool axes_t_t::mouse(const int button, const int state, const int x, const int y) {
        if(!Mouse) return false;
        if ( button == GLUT_LEFT_BUTTON && state == GLUT_DOWN ) {// Left Click
            xButtonDown = x;
            yButtonDown = y;
            if(type == _3D) {//3D
                ctaButtonDown = cta;
                phiButtonDown = phi;
            }
            XMouse = x;
            YMouse = y;
            return true;
        }
        if(mouse_callback) mouse_callback(button, state, x, y);
        return false;
    }

    bool axes_t_t::motion(const int x, const int y) {
        if(type == _3D) {
            cta = (ctaButtonDown - (float) (x - xButtonDown)*1);
            phi = (phiButtonDown + (float) (y - yButtonDown)*1);
            if(phi >=  90) { phi =  90; }
            if(phi <= -90) { phi = -90; }
            if(cta >  360) { cta -= 360; }
            if(cta <    0) { cta += 360; }

            return true;
        }
        return false;
    }

    dvec axes_t_t::make_tick(const double min, const double max) {
        int i,j;
        double dg;
        double x,y;
        int z;
        x = std::abs(max - min);
        z = (int)std::log10(x);
        y = x/std::pow((double)10, (double)z);
        dg = std::pow((double)10, (double)z);
             if(y < 2) { dg *= 0.2; }
        else if(y < 5) { dg *= 0.5; }

        double min0 = min - std::fmod(min,dg);
        j = 0;

        dvec tick;
        if(max > min) { i=-2; while(max >= min0 + dg*i) { if(min <= min0 + dg*i) { tick.push_back(min0 + dg*i); ++j; } ++i; } }
        if(max < min) { i=-2; while(max <= min0 - dg*i) { if(min >= min0 - dg*i) { tick.push_back(min0 - dg*i); ++j; } ++i; } }
        return tick;
    }

    void axes_t_t::draw() {
        boost::mutex::scoped_lock l(children_mutex);
        if(children.size() > 0) {
            if(type == _2D) {
                draw2d();
            } else if(type == _3D) {
                draw3d();
            }
        }
        if(type == color_bar) {//colorbar
            draw_colorbar();
        }

        config();
        // children //
        for(drawings_t::iterator d = children.begin(); d != children.end(); ++d) {
            (*d)->draw();
        }

        if(color_bar_axes) color_bar_axes->draw();
    }

    /// display
    void axes_t_t::draw2d() {
        char ctmp[100];
        float l,b,w,h;//left,bottom,width,height
        float r = 0.01;

        l = position[0];
        b = position[1];
        w = position[2];
        h = position[3];

        // viewport figure (VpF) for drawing axes
        glViewport(0,0, (int)(window_w()), (int)(window_h()));
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );

        glDisable(GL_LINE_STIPPLE);
        gl2psDisable(GL2PS_LINE_STIPPLE);

        //float x_axis_location = 1, y_axis_location = 1;
        // if(XAxisLocation=="top"  ) { x_axis_location = 1; } else { x_axis_location = -1; }
        // if(YAxisLocation=="right") { y_axis_location = 1; } else { y_axis_location = -1; }

        int char_w = 6, char_h = 12;
        float offset = 0.01;
        int num_char = 4;

        int gridlinestyle;

        int tickdir = 1;//1:in, -1:out

        if(Box) {
            // box //
            glLineWidth(LineWidth);
            gl2psLineWidth(LineWidth);
            if(selected) {
                glLineWidth(2*LineWidth);
                gl2psLineWidth(LineWidth);
            }
            glColor3f(0,0,0);
            glBegin(GL_LINE_LOOP);
                glVertex2d(l,  b);
                glVertex2d(l+w,b);
                glVertex2d(l+w,b+h);
                glVertex2d(l,  b+h);
            glEnd();

            // mouse capture //
            if(selected) {
                sprintf(ctmp,"Mouse: (%f,%f)",XMouse,YMouse);
                ptext(l, b+h+r, ctmp);
            }

            // Grid //
            gridlinestyle = 3;
                 if(GridLineStyle == "-"  ) { gridlinestyle = 1; }
            else if(GridLineStyle == "- -") { gridlinestyle = 2; }
            else if(GridLineStyle == ":"  ) { gridlinestyle = 3; }
            else if(GridLineStyle == "-." ) { gridlinestyle = 4; }

            if(XGrid) {
                glLineWidth(LineWidth);
                gl2psLineWidth(LineWidth);
                for(unsigned int i = 0; i < XTick.size(); ++i) {

                    //cout <<"grid "<<gridlinestyle<<" "<<XTick[i]<<endl;

                    if(gridlinestyle == 1) {// -
                        glDisable(GL_LINE_STIPPLE);
                        gl2psDisable(GL2PS_LINE_STIPPLE);
                    }
                    else if(gridlinestyle == 2) {//- -
                        glEnable(GL_LINE_STIPPLE);
                        glLineStipple(1, 0xF0F0);
                        gl2psEnable(GL2PS_LINE_STIPPLE);
                    }
                    else if(gridlinestyle == 3) {//:
                        glEnable(GL_LINE_STIPPLE);
                        glLineStipple(1, 0xCCCC);
                        gl2psEnable(GL2PS_LINE_STIPPLE);
                    }
                    else if(gridlinestyle == 4) {//-.
                        glEnable(GL_LINE_STIPPLE);
                        glLineStipple(1, 0x087F);
                        gl2psEnable(GL2PS_LINE_STIPPLE);
                    }
                    glBegin(GL_LINE_STRIP);
                        glVertex2d( ctx(XTick[i]),b );
                        glVertex2d( ctx(XTick[i]),b+h );//!! TODO
                    glEnd();
                }
            }
            if(YGrid) {
                for(unsigned int i = 0; i < XTick.size(); ++i) {
                    if(gridlinestyle == 1) {// -
                        glDisable(GL_LINE_STIPPLE);
                        gl2psDisable(GL2PS_LINE_STIPPLE);
                    }
                    else if(gridlinestyle == 2) {//- -
                        glEnable(GL_LINE_STIPPLE);
                        glLineStipple(1, 0xF0F0);
                        gl2psEnable(GL2PS_LINE_STIPPLE);
                    }
                    else if(gridlinestyle == 3) {//:
                        glEnable(GL_LINE_STIPPLE);
                        glLineStipple(1, 0xCCCC);
                        gl2psEnable(GL2PS_LINE_STIPPLE);
                    }
                    else if(gridlinestyle == 4) {//-.
                        glEnable(GL_LINE_STIPPLE);
                        glLineStipple(1, 0x087F);
                        gl2psEnable(GL2PS_LINE_STIPPLE);
                    }
                    glBegin(GL_LINE_STRIP);
                        glVertex2d( l,  cty(YTick[i]) );
                        glVertex2d( l+w,cty(YTick[i]) );
                    glEnd();
                }
            }

            // Ticks //
            if(TickDir=="in") { tickdir =  1; }
            if(TickDir=="out") { tickdir = -1; }

            glDisable(GL_LINE_STIPPLE);
            gl2psDisable(GL2PS_LINE_STIPPLE);
            //TODO precise adjustment of tick location
            // x tick
            for(unsigned int i = 0; i < XTick.size(); ++i) {
                glBegin(GL_LINE_STRIP);
                    glVertex2d( ctx(XTick[i]),b );
                    glVertex2d( ctx(XTick[i]),b+tickdir*0.01 );//b-0.02*h
                glEnd();
            }
            // x tick label
            if(TickLabel) {
                for(unsigned int i = 0; i < XTick.size(); ++i) {
                    sprintf(ctmp,"%4.1f",XTick[i]);
                    //ptext( ctx(XTick[i])-0.02, b-0.025,ctmp );//b-0.05*h
                    ptext( ctx(XTick[i])-(float)num_char*char_w/window_w()/2.0,
                       b-offset-1.0*char_h/window_h(),ctmp );//b-0.05*h
                }
            }
            // y tick
            for(unsigned int i = 0; i < YTick.size(); ++i) {
                glBegin(GL_LINE_STRIP);
                    glVertex2d( l,             cty(YTick[i]) );
                    glVertex2d( l+tickdir*0.01,cty(YTick[i]) );
                glEnd();
            }
            // y tick label
            if(TickLabel) {
                for(unsigned int i = 0; i < YTick.size(); ++i) {
                    sprintf(ctmp,"%4.1f", YTick[i]);
                    //ptext( l-0.05,cty(YTick[i])-0.0,ctmp );
                    ptext( l-(float)num_char*char_w/window_w()-offset,
                       cty(YTick[i])-0.5*char_h/window_h(),ctmp );
                }
            }
        }//Box

        //Title
        num_char = Title.length();
        ptext( l+w/2.0-(float)num_char*char_w/window_w()/2.0,
           b+h+offset,
           Title );

        //XLabel
        num_char = XLabel.length();
        ptext( l+w/2.0-(float)num_char*char_w/window_w()/2.0,
           b-offset-2.0*char_h/window_h(),
           XLabel );

        //YLabel
        num_char = YLabel.length();
        ptext(l, b+h+offset, YLabel);

        // viewport Axes (VpA) for drawing lines and surfaces
        glViewport((int)(position[0]*window_w() ),
                   (int)(position[1]*window_h() ),
                   (int)(position[2]*window_w() ),
                   (int)(position[3]*window_h() ));
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );
    }

    void axes_t_t::draw3d() {
        char ctmp[100];

        // float l,b,w,h;//left,bottom,width,height
        // l = position[0];
        // b = position[1];
        // w = position[2];
        // h = position[3];

        // viewport Axes
        glViewport((int)(viewport3d[0]*window_w() ),
                   (int)(viewport3d[1]*window_h() ),
                   (int)(viewport3d[2]*window_w() ),
                   (int)(viewport3d[3]*window_h() ));

        glLoadIdentity();
        //glOrtho(-1.7, 1.7, -1.7, 1.7, -1.5, 3);
        glOrtho(-1.8, 1.8, -1.8, 1.8, -1.5, 3);

        gluLookAt(cos(cta*M_PI/180)*cos(phi*M_PI/180),
              sin(cta*M_PI/180)*cos(phi*M_PI/180),
              sin(phi*M_PI/180),
              //gluLookAt(camera_position[0],camera_position[1],camera_position[2],
              camera_target[0],  camera_target[1],  camera_target[2],
              camera_up_vector[0],camera_up_vector[1],camera_up_vector[2]);

        if(Box) {
            // tick
            float cta0;
            float r1 = 1.05;//tick width
            float r2 = 1.2;
            float r3 = 1.4;
            int signx,signy;
            cta0 = cta; cta0 = std::fmod(cta,360);
            if((  0<=cta0) && (cta0< 90)) { signx = 1;signy = 1; }
            if(( 90<=cta0) && (cta0<190)) { signx =-1;signy = 1; }
            if((180<=cta0) && (cta0<270)) { signx =-1;signy =-1; }
            if((270<=cta0) && (cta0<360)) { signx = 1;signy =-1; }

            glColor3f(0,0,0);

            // axes //
            // x
            glBegin(GL_LINE_STRIP);
                glVertex3d(-1,signy,-1);
                glVertex3d( 1,signy,-1);
            glEnd();
            // y
            glBegin(GL_LINE_STRIP);
                glVertex3d(signx,-1,-1);
                glVertex3d(signx, 1,-1);
            glEnd();
            // z
            glBegin(GL_LINE_STRIP);
                glVertex3d(signy,-signx,-1);
                glVertex3d(signy,-signx, 1);
            glEnd();

            // Tick //
            //x
            for(unsigned int i = 0; i < XTick.size(); ++i) {
                glBegin(GL_LINE_STRIP);
                    glVertex3d( ct3x(XTick[i]),signy   ,-1 );
                    glVertex3d( ct3x(XTick[i]),signy*r1,-1 );
                glEnd();
            }
            // y
            for(unsigned int i = 0; i < YTick.size(); ++i) {
                glBegin(GL_LINE_STRIP);
                    glVertex3d( signx   ,ct3y(YTick[i]),-1 );
                    glVertex3d( signx*r1,ct3y(YTick[i]),-1 );
                glEnd();
            }
            // z
            for(unsigned int i = 0; i < YTick.size(); ++i) {
                glBegin(GL_LINE_STRIP);
                    glVertex3d( signy   , -signx, ct3z(ZTick[i]) );
                    glVertex3d( signy*r1, -signx, ct3z(ZTick[i]) );
                glEnd();
            }
            // Tick Label //
            if(TickLabel) {
                //x
                for(unsigned int i = 0; i < XTick.size(); ++i) {
                    sprintf(ctmp,"%4.1f",XTick[i]);
                    ptext3c( ct3x(XTick[i]),signy*r2 ,-1,ctmp );
                }
                // y
                for(unsigned int i = 0; i < YTick.size(); ++i) {
                    sprintf(ctmp,"%4.1f",YTick[i]);
                    ptext3c( signx*r2,ct3y(YTick[i]),-1,ctmp );
                }
                // z
                for(unsigned int i = 0; i < ZTick.size(); ++i) {
                    sprintf(ctmp,"%4.1f",ZTick[i]);
                    ptext3c( signy*r2,-signx,ct3z(ZTick[i]),ctmp );
                }
            }
            // xyz Label //
            ptext3c(0,signy*r3,-1,"x");
            ptext3c(signx*r3,0,-1,"y");
            ptext3c(signy*r3,-signx,0,"z");

        }//box
    }

    /// colorbar
    axes_t axes_t_t::colorbar() {
        axes_t p(new axes_t_t(layer));
        color_bar_axes = p;
        float l,b,w,h;
        l = position[0];
        b = position[1];
        w = position[2];
        h = position[3];
        float zmin = ZLim[0];
        float zmax = ZLim[1];

        // TODO use in 3D

        p->cmap = cmap;
        p->type = color_bar;
        p->position[0] = l + w*1.01;
        p->position[1] = b;
        p->position[2] = w*0.05;
        p->position[3] = h;
        p->ZLim[0] = zmin;
        p->ZLim[1] = zmax;
        p->YLim[0] = zmin;
        p->YLim[1] = zmax;
        return color_bar_axes;
    }

    /// colorbar
    void axes_t_t::draw_colorbar() {
        char ctmp[100];
        float l,b,w,h;//left,bottom,width,height

        l = position[0];
        b = position[1];
        w = position[2];
        h = position[3];

        // viewport figure (VpF) for drawing axes
        glViewport(0,0, (int)(window_w()), (int)(window_h()));
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );

        glDisable(GL_LINE_STIPPLE);
        gl2psDisable(GL2PS_LINE_STIPPLE);

        if(Box) {
            // box
            glLineWidth(LineWidth);
            gl2psLineWidth(LineWidth);
            glColor3f(0, 0, 0);
            glBegin(GL_LINE_LOOP);
                glVertex2d(l,  b);
                glVertex2d(l+w,b);
                glVertex2d(l+w,b+h);
                glVertex2d(l,  b+h);
            glEnd();

            // z tick
            for(unsigned int i = 0; i < ZTick.size(); ++i) {
                glBegin(GL_LINE_STRIP);
                    glVertex2d( l+w,       cty(ZTick[i]) );
                    glVertex2d( l+w+0.01,  cty(ZTick[i]) );
                glEnd();
            }
            // z tick number
            for(unsigned int i = 0; i < ZTick.size(); ++i) {
                sprintf(ctmp,"%4.1f",ZTick[i]);
                ptext( l+w+0.01,cty(ZTick[i]),ctmp );
            }
        }//Box

        std::vector<float> rgb;
        unsigned int n = cmap.size();
        for(unsigned int i = 0; i < n; ++i) {
            rgb = cmap[i];
            glColor3f(rgb[0], rgb[1], rgb[2]);

            glBegin(GL_QUADS);
                glVertex2d(l  ,b+h*i/n);
                glVertex2d(l+w,b+h*i/n);
                glVertex2d(l+w,b+h*(i+1)/n);
                glVertex2d(l  ,b+h*(i+1)/n);
            glEnd();
        }
    }

    // figure coordination
    float axes_t_t::ctx(const double x) {
        double t;
        if(XScale == math::linear_scale) {//math::linear_scale
            return position[0] + position[2]*( ( x - XLim[0] ) / ( XLim[1] - XLim[0]) );
        } else {//log
            t = ( std::log10(x) - std::log10(XLim[0]) ) / ( std::log10(XLim[1]) - std::log10(XLim[0]) );
            if(x <= 0) { t = -1; }
            return position[0] + position[2]*t;
        }

    }
    float axes_t_t::cty(const double y) {
        if(YScale == math::linear_scale) {//math::linear_scale
            return position[1] + position[3]*( ( y - YLim[0] ) / ( YLim[1] - YLim[0] ) );
        } else {//log
            return position[1] + position[3]*( std::log10(y) - std::log10(YLim[0]) ) / ( std::log10(YLim[1]) - std::log10(YLim[0]) );
        }
    }


    float axes_t_t::ct3x(const double x) {
        return -1 + 2*( x - XLim[0] ) / ( XLim[1] - XLim[0] );
    }
    float axes_t_t::ct3y(const double y) {
        return -1 + 2*( y - YLim[0] ) / ( YLim[1] - YLim[0] );
    }
    float axes_t_t::ct3z(const double z) {
        return -1 + 2*( z - ZLim[0] ) / ( ZLim[1] - ZLim[0] );
    }

    /// axis
    axes_t axes_t_t::axis(const double xMin, const double xMax, const double yMin, const double yMax) {
        if(xMin != xMax) {
            XLim[0]  = xMin;
            XLim[1]  = xMax;
            XLimMode = manual;
        }
        if(yMin != yMax) {
            YLim[0]  = yMin;
            YLim[1]  = yMax;
            YLimMode = manual;
        }
        type = axes_t_t::_2D;
        return shared_from_this();
    }

    axes_t axes_t_t::axis(const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax) {
        if(xMin != xMax) {
            XLim[0]  = xMin;
            XLim[1]  = xMax;
            XLimMode = manual;
        }
        if(yMin != yMax) {
            YLim[0]  = yMin;
            YLim[1]  = yMax;
            YLimMode = manual;
        }
        if(zMin != zMax) {
            ZLim[0]  = zMin;
            ZLim[1]  = zMax;
            ZLimMode = manual;
        }
        type = axes_t_t::_3D;
        return shared_from_this();
    }

    axes_t axes_t_t::axis(const bool s)                 { Box = s; return shared_from_this(); }
    axes_t axes_t_t::axis(const std::string s)          { axis(s == "on"); return shared_from_this(); }
    axes_t axes_t_t::grid(bool s)                       { XGrid = YGrid = ZGrid = s; return shared_from_this(); }
    axes_t axes_t_t::grid(const std::string s)          { grid(s == "on"); return shared_from_this(); }
    axes_t axes_t_t::ticklabel(const bool s)            { TickLabel = s; return shared_from_this(); }
    axes_t axes_t_t::title(const std::string s)         { Title = s; return shared_from_this(); }
    axes_t axes_t_t::xlabel(const std::string s)        { XLabel = s; return shared_from_this(); }
    axes_t axes_t_t::ylabel(const std::string s)        { YLabel = s; return shared_from_this(); }
    axes_t axes_t_t::mouse_capture(const bool y)        { Mouse = y; return shared_from_this(); }

    /// ptext
    void axes_t_t::ptext(float x, float y, const std::string s) {
        // viewport figure_t_t::
        glViewport(0,0, (int)(window_w()), (int)(window_h()));
        glLoadIdentity();
        gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );

        glColor3f(0,0,0);
        glRasterPos2f( x, y );
        gl2psText(s.c_str(), "Arial", 12);
        //gl2psText(test, "Times-Roman", 24);

        for(int i = 0; i < (int)s.size(); ++i) {
            glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, s[i] );
        }
    }
    void axes_t_t::ptext3(float x, float y, float z, const std::string s) {
        glColor3f(0,0,0);
        glRasterPos3d(x,y,z);
        gl2psText(s.c_str(), "Arial", 12);
        for(int i = 0; i < (int)s.size(); ++i) {
            glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, s[i] );
        }
    }
    void axes_t_t::ptext3c(float x, float y, float z, const std::string s) {
        int char_w = 6,char_h = 12;
        int num_char = s.length();
        glColor3f(0,0,0);
        glRasterPos3d(x,y,z);
        glBitmap(0,0,0,0, -char_w*num_char/2, -char_h/2, NULL);
        gl2psText(s.c_str(), "Arial", 12);

        for(int i=0; i<(int)s.size(); ++i) {
            glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, s[i] );
        }
    }

    axes_t axes_t_t::set(const std::string p, const std::string v) {
        if(p == "TickDir") { TickDir = v; }
        return shared_from_this();
    }


    std::vector<float> axes_t_t::colormap(const std::string c, float t) {

        std::vector<float> rgb(3);
        if(t > 1) { t = 1; }
        if(t < 0) { t = 0; }

        if( c == "Gray" ) {
            rgb[0] = t;
            rgb[1] = t;
            rgb[2] = t;
            return rgb;
        }
        else if( c == "HSV" ) {
            t*=6;
            if(  0<=t && t<=1.0) { rgb[0] = 1;        rgb[1] = t;        rgb[2] = 0; }
            if(1.0<=t && t<=2.0) { rgb[0] = 1-(t-1);  rgb[1] = 1;        rgb[2] = 0; }
            if(2.0<=t && t<=3.0) { rgb[0] = 0;        rgb[1] = 1;        rgb[2] = t-2; }
            if(3.0<=t && t<=4.0) { rgb[0] = 0;        rgb[1] = 1-(t-3);  rgb[2] = 1; }
            if(4.0<=t && t<=5.0) { rgb[0] = t-4;      rgb[1] = 0;        rgb[2] = 1; }
            if(5.0<=t && t<=6.0) { rgb[0] = 1;        rgb[1] = 0;        rgb[2] = 1-(t-5); }

            return rgb;
        }
        else if( c == "Jet" ) {
            t*=8;
            if(  0<=t && t<=1.0) { rgb[0] = 0;        rgb[1] = 0;          rgb[2] = 0.5+0.5*t; }
            if(1.0<=t && t<=3.0) { rgb[0] = 0;        rgb[1] = 0.5*(t-1);  rgb[2] = 1; }
            if(3.0<=t && t<=5.0) { rgb[0] = 0.5*(t-3);rgb[1] = 1;          rgb[2] = 1-0.5*(t-3); }
            if(5.0<=t && t<=7.0) { rgb[0] = 1;        rgb[1] = 1-0.5*(t-5);rgb[2] = 0; }
            if(7.0<=t && t<=8.0) { rgb[0] = 1-0.5*(t-7);rgb[1] = 0;        rgb[2] = 0; }
            return rgb;
        }
        else if( c == "Hot" ) {
            t*=3;
            if(  0<=t && t<=1.0) { rgb[0] = t; rgb[1] = 0;   rgb[2] = 0; }
            if(1.0<=t && t<=2.0) { rgb[0] = 1; rgb[1] = t-1; rgb[2] = 0; }
            if(2.0<=t && t<=3.0) { rgb[0] = 1; rgb[1] = 1;   rgb[2] = t-2; }
            return rgb;
        }
        else if( c == "Cool" ) {
            rgb[0] = t;
            rgb[1] = 1-t;
            rgb[2] = 1;
            return rgb;
        }
        else if( c == "Spring" ) {// Magenta - Yellow
            rgb[0] = 1;
            rgb[1] = t;
            rgb[2] = 1-t;
            return rgb;
        }
        else if( c == "Summer" ) {// Green Yellow
            rgb[0] = t;
            rgb[1] = 1;
            rgb[2] = 0;
            return rgb;
        }
        else if( c == "Autumn" ) {
            rgb[0] = 1;
            rgb[1] = t;
            rgb[2] = 0;
            return rgb;
        }
        else if( c == "Winter" ) {
            rgb[0] = 0;
            rgb[1] = t;
            rgb[2] = 1-t;
            return rgb;
        }
        else { // if( c == "Bone" ) {
            rgb[0] = t;
            rgb[1] = t < 0.8 ? t : 0.8;
            rgb[2] = t;
            return rgb;
        }
    }

    void axes_t_t::colormap(const std::string c) {
        //if(is_debug1) { printf("colormap %s \n",c.c_str()); }
        const int n = 64;

        cmap.clear();
        for(int i = 0; i < n; ++i) {
            cmap.push_back(colormap(c, (float)i/(n-1)));
        }

        if(color_bar_axes) color_bar_axes->cmap = cmap;
    }

    void axes_t_t::colormap(const std::vector<std::vector<float> >& c) {
        cmap = c;
    }

    std::vector<float> axes_t_t::map2color(const double x) {
        xmin = CLim[0]; xmax = CLim[1];
        int n = cmap.size();
        float normx;
        std::vector<float> rgb(3);

        normx = (x - xmin)/(xmax - xmin);
        if(x > xmax) { normx = 1; }
        if(x < xmin) { normx = 0; }
        rgb = cmap[(int)(normx*(n-1))];
        //cout << "c: "<<(int)(normx*n) <<endl;
        //cout << "rgb: "<<rgb[0]<<" "<<endl;
        return rgb;
    }
}
