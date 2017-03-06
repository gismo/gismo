=== THB Splines plug-ins for Rhino ===

This is the source for the plug-ins for Rhino that support truncated hierarchical B-Splines.

Please read the following before proceeding.

1. Visual Studio Versions
The source can be built for Rhino V5 (64-bit only) using Visual Studio 2010 and for Rhino V6 (currently WIP) using Visual Studio 2015. Other Visual Studio versions are not supported because building plug-ins for Rhino is only supported in these two versions.

2. Rhino V5 and Rhino V6 support
Going from Rhino V5 to V6 there have been numerous breaking changes in the C++ SDK of Rhino. To facilitate this, you must #define the RHINO_V6_READY flag when building for Rhino6. This is done automatically in the stdafx.h header file. Furthermore the RHINO_V6_READY flag is used to distinguish code for V5 and V6.

3. Rhino SDK installation
The Rhino V5 SDK can be downloaded from http://www.rhino3d.com/download/rhino-sdk/5.0/release
The Rhino V6 SDK can be downloaded from http://www.rhino3d.com/download/rhino-sdk/wip
The build expects that you have installed the Rhino SDK. For Rhino V5, the SDK must be installed in C:\Program Files (x86)\Rhino 5.0 x64 SDK. For Rhino V6, the SDK must be installed in C:\Program Files\Rhino 6.0 SDK. 

4. G+SMO build and installation
G+SMO must be built using the same Visual Studio version. Currently, the path to the install of G+SMO is hard-coded as C:\Program Files\gismo\ with the standard 'include\gismo\' and 'lib' folders. If your installation of G+SMO is not at C:\Program Files\gismo, you need to change the paths in the GismoSupport.vcxproj file.

5. Design
5a. plug-ins
The plug-ins are all designed in such a way that they do not depend on G+SMO directly, but depend on a DLL called GismoSupport.dll. This dll contains all functionality, and the plug-ins delegate their command, import and export routines to the GismoSupport DLL. For commands this means that each command has a command implementation class, that implements the ICommandImplementation pure virtual class (interface). For the import and export routines, the declaration of supported file types and the reading/writing of files is also handled by implementation classes (using interfaces IFileImportImplementation/IFileExportImplementation).

5b. Custom surface objects
The truncated hierarchical B-Spline surface is represented in Rhino as a sub-class of the CRhinoBrepObject class called CThbSurfaceObject. This object is responsible for drawing itself and enabling the colored control points. The gismo truncated hierarchical B-Spline surface object is stored in user-data (CThbSurfaceUserData). This is done to enable transformations of the surface in the CThbSurfaceUserData::Transform(const ON_Xform&) function. The user-data is also responsible for reading and writing the gismo object to the Rhino 3dm file; this is done by converting it to the XML definition and writing that to the 3dm file. (it would be nice if gismo objects could be serialized in binary form!)

5c. Control points
The control points of the surface can be displayed by enabling control points using F10 or the PointsOn command in Rhino. The custom object then uses a grips-enabler to create and enable the grips (grips == control points). The control points are managed in the CThbSurfaceGrips class, that makes individual CThbSurfaceGrip objects. Each level is drawn in a different color by selectively showing and hiding grips. When the user interacts with a grip (for instance by dragging it), the CThbSurfaceGrip::NewLocation function is called, each time the control point gets a new location. This change is propagated to the CThbSurfaceGrips::UpdateGrips(int) function, which updates the ON_Brep representation of the gismo surface object. This ON_Brep representation is used to draw the surface _while it is being edited_.

6. To truncate or not to truncate?
Currently, only gsTHBSpline<2> is supported.

7. Reduce compile time for template instantiation
Non-templated subclasses of gsTHBSpline<2> and gsTensorNurbs<2> were created in the pre-compiled header to prevent template instantiation for each compile unit. My C++ fu is not so strong that I could prevent this (CLASS_TEMPLATE_INST did not work, at least not in VS2010). This speeds up the compile time tremendously.

8. Supporting VS2010 and VS2015 with the same VCXPROJ file
To prevent having to create and support two .vcxproj files, one for each version of Visual Studio, the root CMakeLists.txt file contains a configure_file command that writes the version dependent input to a Platform.targets file. This file is included by each vcxproj file to switch between the v100 and v140 toolchain.