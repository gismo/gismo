/** @file gsReadParasolid.h

    @brief Provides declaration of gsReadParasolid functions.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsIO/gsXmlUtils.h>

#include <gsParasolid/gsPKSession.h>

typedef int PK_BSURF_t;
typedef int PK_SURF_t;
typedef int PK_PART_t;
typedef int PK_BODY_t;

namespace gismo {

namespace extensions {
    
    /// Reads a parasolid file
    /// \param fname filename (without extension)
    bool gsReadParasolid( const char * fname, internal::gsXmlTree & data  );
    
    /// Extracts a Parasolid part
    bool readPK_PART( const PK_PART_t & part, internal::gsXmlTree & data  );

    /// Extracts all geometries from a Parasolid part
    bool readPK_PART_geoms( const PK_PART_t & part, internal::gsXmlTree & data  );

    /// Extracts a Parasolid body
    bool readPK_BODY( const PK_BODY_t & body, internal::gsXmlTree & data  );

    /// Extracts a surface from Parasolid
    bool readPK_SURF( const PK_SURF_t & pkbs, internal::gsXmlTree & data  );

    /// Extracts a b-surface from Parasolid
    bool readPK_BSURF( const PK_BSURF_t & pkbs, internal::gsXmlTree & data  );

}//extensions

}//gismo

