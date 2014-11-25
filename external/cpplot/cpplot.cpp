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

#include "cpplot.hpp"
#include <algorithm>
#include <cstdio>
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <iostream>

namespace cpplot {
    figures_t figures;
    figure_t cf;
    typedef std::map<std::string, figure_t> figuremap;
    figuremap named_figures;
    boost::mutex figures_mutex;
    boost::mutex named_figures_mutex;

    figure_t figure(const std::string name) {
        boost::mutex::scoped_lock l(named_figures_mutex);
        figure_t f = named_figures[name];
        if(f == NULL) {
            std::cout << "New figure: " << name << std::endl;
            named_figures[name] = f = figure();
            f->set_window_name(name);
        }
        return f;
    }

    int max_figure_number = 0;
    figure_t figure() {
        return figure(max_figure_number + 1);
    }

    figure_t figure(const int i) {
        if(i > max_figure_number) max_figure_number = i;
        boost::mutex::scoped_lock l(figures_mutex);
        cf = figures[i];
        if(cf == NULL) {
            figure_t p(new figure_t_t());
            cf = p;
            glut::register_figure(p);
            figures[i] = p;
            return cf;
        } else {
            return cf;
        }
    }
}
