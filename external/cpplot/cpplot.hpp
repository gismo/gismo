/****************************************************************************
Copyright (c) 2011 Jonatan Olofsson All Rites Reversed
License: Gnu Public license (GPL) v3
* Author: Jonatan Olofsson (jonatan.olofsson@gmail.com)
* Version: 0.1
* Based on
Author: Yuichi Katori (yuichi.katori@gmail.com)
Project:MATPLOT++ (MATLAB-like plotting tool in C++).
Version:0.3.13
****************************************************************************/

/**
 * In this file, the main API of cpplot is held. Most of it is shorthand
 * redirects to functions defined in the main class structure of figures,
 * layers, axes and drawings (here in logical order). This file
 * may also be viewed as reference to access more advanced features
 * of cpplot which may not be available as shorthand functions.
 */

#ifndef _LIB_CPPLOT_
#define _LIB_CPPLOT_

#include "cpplot_common.hpp"

namespace cpplot {
    extern figure_t cf;

    /// Initialize glut
    inline void figure_start()
    {
        char *myargv [1];
        int myargc=1;
        myargv [0]=strdup("gsPlot");
        glut::init(myargc, myargv);
    };

    /// After making a plot, make the thread wait..
    inline void figure_end()
    { 
        //while(true) boost::this_thread::sleep(boost::posix_time::seconds(1000));
        // boost::mutex mymut;
        // boost::mutex::scoped_lock mylock(mymut);
        // mylock.lock();        
        while(true) boost::this_thread::yield();
    };
    

    /**
     * Obtain a figure with a given name. The name will, as well as
     * being a unique identifier for the figure, also be the title
     * of the window
     */
    figure_t figure(const std::string name);
    figure_t figure(); ///< Create a new, unnamed, figure
    figure_t figure(const int i); ///< Obtain, and if nescessary create, a figure with a given numerical ID



    inline figure_t gcf() { return cf ? cf : figure(0); } ///< Get a pointer to the currently active figure
    inline axes_t subplot(const int m, const int n, const int p) { return gcf()->gcl()->subplot(m,n,p); } ///< Get a new axes object, placed as a subplot in the current layer
    //void legend(string s,int N);

    // interface ///
    inline layer_t layer(const std::string name = "plot", const bool viz = true) { return gcf()->layer(name, viz); } ///< See layer.hpp
    inline axes_t gca() { return gcf()->gcl()->gca(); } ///< See layer.hpp
    inline dmat peaks(int n) { return math::peaks(n); } ///< See math.hpp

    /**
     * Obtain a vector with n linearly spaced doubles between [min, max]
     *
     * See math.hpp
     */
    inline std::vector<double> linspace(const double min, const double max, int n) { return math::linspace<std::vector<double> >(min, max, n); }

    template<typename T>
    inline void set(const std::string v) { gca()->gco<T>()->set(v); } ///< Set property of current drawing object
    inline void set(const std::string v) { gca()->gco<Line>()->set(v); } ///< Set property of current drawing object
    template<typename T>
    inline void set(const float v) { gca()->gco<T>()->set(v); } ///< Set property of current drawing object
    inline void set(const float v) { gca()->gco<Line>()->set(v); } ///< Set property of current drawing object
    template<typename T>
    inline void set(const std::string p, const std::string v) { gca()->gco<T>()->set(p,v); } ///< Set property of current drawing object
    inline void set(const std::string p, const std::string v) { gca()->gco<Line>()->set(p,v); } ///< Set property of current drawing object
    template<typename T>
    inline void set(const std::string p, float v) { gca()->gco<T>()->set(p,v); } ///< Set property of current drawing object
    inline void set(const std::string p, float v) { gca()->gco<Line>()->set(p,v); } ///< Set property of current drawing object

    // Axes
    inline void axis(const double xMin, const double xMax, const double yMin, const double yMax) {
        gca()->axis(xMin, xMax, yMin, yMax);
    } ///< See axes.hpp
    inline void axis(const double xMin, const double xMax, const double yMin, const double yMax, const double zMin, const double zMax) {
        gca()->axis(xMin, xMax, yMin, yMax, zMin, zMax);
    } ///< See axes.hpp
    inline void axis(const std::string s) { gca()->axis(s); } ///< See axes.hpp
    inline void axis(const bool s) { gca()->axis(s); } ///< See axes.hpp
    inline void grid(const std::string s) { gca()->grid(s); } ///< See axes.hpp
    inline void grid(const bool s) { gca()->grid(s); } ///< See axes.hpp
    inline void ticklabel(const bool s) { gca()->ticklabel(s); } ///< See axes.hpp
    inline void title(const std::string s) { gca()->title(s); } ///< See axes.hpp
    inline void xlabel(const std::string s) { gca()->xlabel(s); } ///< See axes.hpp
    inline void ylabel(const std::string s) { gca()->ylabel(s); } ///< See axes.hpp
    inline void mouse_capture(const bool y) { gca()->mouse_capture(y); } ///< See axes.hpp

    // Line
    inline void vertex(const double x, const double y) { gca()->gco<Line>()->vertex(x,y); } ///< See line.hpp
    inline void vertex(const double x, const double y, const double z) { gca()->gco<Line>()->vertex(x,y,z); } ///< See axes.hpp

    inline line_t plot(const dvec& y) { return gca()->add<Line>()->plot(y); } ///< See axes.hpp
    inline line_t plot(const dvec& x,const dvec& y) { return gca()->add<Line>()->plot(x,y); } ///< See axes.hpp
    inline line_t plot(const dvec& x, const dvec& y, const dvec& z) { return gca()->add<Line>()->plot(x,y,z); } ///< See axes.hpp

    inline line_t semilogx(const dvec& x, const dvec& y) { return gca()->add<Line>()->semilogx(x,y); } ///< See axes.hpp
    inline line_t semilogy(const dvec& x, const dvec& y) { return gca()->add<Line>()->semilogy(x,y); } ///< See axes.hpp
    inline line_t loglog(const dvec& x, const dvec& y)   { return gca()->add<Line>()->loglog(x,y); } ///< See axes.hpp

    inline void vertex(const double x, const double y, const double ep, const double em)
        { gca()->gco<Line>()->vertex(x,y,ep,em); } ///< See axes.hpp
    inline void errorbar(const dvec& x, const dvec& y, const dvec& e)
        { gca()->gco<Line>()->errorbar(x,y,e); } ///< See axes.hpp
    inline void errorbar(const dvec& x, const dvec& y, const dvec& ep, const dvec& em)
        { gca()->gco<Line>()->errorbar(x,y,ep, em); } ///< See axes.hpp


    // Surface, Contour
    inline surface_t surface(const dmat& Z) { return gca()->add<Surface>()->surface(Z); } ///< See surface.hpp
    inline surface_t surface(const dmat& Z, const dmat& C) { return gca()->add<Surface>()->surface(Z, C); } ///< See surface.hpp
    inline surface_t surface(const dmat& Z, const tcmat& C) { return gca()->add<Surface>()->surface(Z, C); } ///< See surface.hpp
    inline surface_t surface(const dvec& x, const dvec& y, const dmat& Z) { return gca()->add<Surface>()->surface(x,y,Z); } ///< See surface.hpp
    inline surface_t surface(const dvec& x, const dvec& y, const dmat& Z, const dmat& C) { return gca()->add<Surface>()->surface(x,y,Z,C); } ///< See surface.hpp
    inline surface_t surface(const dvec& x, const dvec& y, const dmat& Z, const tcmat& C) { return gca()->add<Surface>()->surface(x,y,Z,C); } ///< See surface.hpp
    inline surface_t surface(const dmat& X, const dmat& Y, const dmat& Z) { return gca()->add<Surface>()->surface(X,Y,Z); } ///< See surface.hpp
    inline surface_t surface(const dmat& X, const dmat& Y, const dmat& Z, const dmat& C) { return gca()->add<Surface>()->surface(X,Y,Z,C); } ///< See surface.hpp
    inline surface_t surface(const dmat& X, const dmat& Y, const dmat& Z, const tcmat& C) { return gca()->add<Surface>()->surface(X,Y,Z,C); } ///< See surface.hpp

    inline surface_t pcolor(const dmat& C) { return gca()->add<Surface>()->pcolor(C); } ///< See surface.hpp
    inline surface_t pcolor(const tcmat& C) { return gca()->add<Surface>()->pcolor(C); } ///< See surface.hpp
    inline surface_t pcolor(const dvec& x, const dvec& y, const dmat& C) { return gca()->add<Surface>()->pcolor(x,y,C); } ///< See surface.hpp
    inline surface_t pcolor(const dvec& x, const dvec& y, const tcmat& C) { return gca()->add<Surface>()->pcolor(x,y,C); } ///< See surface.hpp
    inline surface_t pcolor(const dmat& X, const dmat& Y, const dmat& C) { return gca()->add<Surface>()->pcolor(X,Y,C); } ///< See surface.hpp
    inline surface_t pcolor(const dmat& X, const dmat& Y, const tcmat& C) { return gca()->add<Surface>()->pcolor(X,Y,C); } ///< See surface.hpp

    inline surface_t contour(const dmat& Z) { return gca()->add<Surface>()->contour(Z); } ///< See surface.hpp
    inline surface_t contour(const dmat& Z,int n) { return gca()->add<Surface>()->contour(Z, n); } ///< See surface.hpp
    inline surface_t contour(const dmat& Z, const dvec& v) { return gca()->add<Surface>()->contour(Z,v); } ///< See surface.hpp
    inline surface_t contour(const dvec& x, const dvec& y, const dmat& Z) { return gca()->add<Surface>()->contour(x,y,Z); } ///< See surface.hpp
    inline surface_t contour(const dvec& x, const dvec& y, const dmat& Z, const int n) { return gca()->add<Surface>()->contour(x,y,Z,n); } ///< See surface.hpp
    inline surface_t contour(const dvec& x, const dvec& y, const dmat& Z, const dvec& v) { return gca()->add<Surface>()->contour(x,y,Z,v); } ///< See surface.hpp

    inline surface_t mesh(const dvec& x, const dvec& y, const dmat& Z) { return gca()->add<Surface>()->mesh(x,y,Z); } ///< See surface.hpp
    inline surface_t surf(const dvec& x, const dvec& y, const dmat& Z) { return gca()->add<Surface>()->surf(x,y,Z); } ///< See surface.hpp

    inline void shading(const std::string c) { gca()->gco<Surface>()->shading(c); } ///< See surface.hpp

    // Patch ///
    inline patch_t patch(const dmat& X, const dmat& Y) { return gca()->add<Patch>()->patch(X,Y); } ///< See patch.hpp
    inline patch_t patch(const dmat& X, const dmat& Y, const dvec& C) { return gca()->add<Patch>()->patch(X,Y,C); } ///< See patch.hpp
    inline patch_t patch(const dmat& X, const dmat& Y, const tcvec& C) { return gca()->add<Patch>()->patch(X,Y,C); } ///< See patch.hpp
    inline patch_t patch(const dmat& X, const dmat& Y, const dmat& Z) { return gca()->add<Patch>()->patch(X,Y,Z); } ///< See patch.hpp
    inline patch_t patch(const dmat& X, const dmat& Y, const dmat& Z, const dvec& C) { return gca()->add<Patch>()->patch(X,Y,Z,C); } ///< See patch.hpp
    inline patch_t patch(const dmat& X, const dmat& Y, const dmat& Z, const tcvec& C) { return gca()->add<Patch>()->patch(X,Y,Z,C); } ///< See patch.hpp

    inline patch_t bar(const dvec& y) { return gca()->add<Patch>()->bar(y); } ///< See patch.hpp
    inline patch_t bar(const dvec& y, const float width) { return gca()->add<Patch>()->bar(y, width); } ///< See patch.hpp
    inline patch_t bar(const dvec& x, const dvec& y) { return gca()->add<Patch>()->bar(x,y); } ///< See patch.hpp
    inline patch_t bar(const dvec& x, const dvec& y, const float width) { return gca()->add<Patch>()->bar(x,y,width); } ///< See patch.hpp

    // Text ///
    //TODO: more fonts
    inline text_t text(const double x, const double y, const std::string s) { return gca()->add<Text>()->text(x,y,s); } ///< See text.hpp

    inline axes_t colorbar() { return gca()->colorbar(); } ///< See axes.hpp
    inline void gray() { gca()->gray(); } ///< See axes.hpp
    inline void jet() { gca()->jet(); } ///< See axes.hpp
    inline void hsv() { gca()->hsv(); } ///< See axes.hpp
    inline void hot() { gca()->hot(); } ///< See axes.hpp
    inline void cool() { gca()->hot(); } ///< See axes.hpp
    inline void spring() { gca()->spring(); } ///< See axes.hpp
    inline void summer() { gca()->summer(); } ///< See axes.hpp
    inline void autumn() { gca()->autumn(); } ///< See axes.hpp
    inline void winter() { gca()->winter(); } ///< See axes.hpp


    // print
    inline void print(const std::string name = "gismo_plot.eps") { if(cf) cf->print(name); } ///< See figure.hpp

    inline void operator<<(line_t ln, const std::pair<double, double> p) { ln->vertex(p.first, p.second); } ///< Add point to line
    inline void operator<<(axes_t a, const std::pair<double, double> p) { a->gco<Line>()->vertex(p.first, p.second); } ///< Add point to current line in the axes object
    inline void operator<<(layer_t a, const std::pair<double, double> p) { a->gca()->gco<Line>()->vertex(p.first, p.second); } ///< Add point to current line in the layer object
    inline void operator<<(figure_t a, const std::pair<double, double> p) { a->gca()->gco<Line>()->vertex(p.first, p.second); } ///< Add point to current line in the figure object
}
#endif
