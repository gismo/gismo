/* (C) MyCompany */


/* Put a short description of your plugin here */

/* MyCompany-contact@mycompany.com-http://www.mycompany.com */

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

