/** @file gsAxelPluginExport.h

    @brief This file is part of the G+Smo plugin for Axel modeler.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#ifndef GSAXELPLUGINEXPORT_H
#define GSAXELPLUGINEXPORT_H

#ifdef WIN32
#ifdef GSAXELPLUGIN_EXPORT
#define GSAXELPLUGIN_EXPORT __declspec(dllexport)
#else
#define GSAXELPLUGIN_EXPORT
#endif
#else
#define GSAXELPLUGIN_EXPORT
#endif

#endif

