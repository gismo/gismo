/* gsWriter.cpp ---
 *
 *
 */

#include "gsWriter.h"

#include <dtkCoreSupport/dtkAbstractData.h>
#include <dtkCoreSupport/dtkAbstractDataFactory.h>

#include<axlCore/axlAbstractDataConverter.h>
#include<axlCore/axlMesh.h>
#include<axlCore/axlPoint.h>
#include<axlCore/axlPointConverter.h>


// /////////////////////////////////////////////////////////////////
// gsWriterPrivate
// /////////////////////////////////////////////////////////////////

class gsWriterPrivate
{
public:
    QFile *filename;

};


// /////////////////////////////////////////////////////////////////
// gsWriter
// /////////////////////////////////////////////////////////////////

gsWriter::gsWriter(void) : dtkAbstractDataWriter(), d(new gsWriterPrivate)
{
    this->setObjectName(this->description());
    d->filename = new QFile();

}

gsWriter::~gsWriter(void)
{

    if(d->filename){
        delete d->filename;
        d->filename = NULL;
    }

    delete d;
    d = NULL;


}



//! Return the identifier "gsWriter".
/*!
 *
 */
QString gsWriter::identifier(void) const
{
    return "gsWriter";
}


//! Return a description of the writer.
/*!
 *
 */
QString gsWriter::description(void) const
{
    return "gsWriter";
}


//! Return the appropriate extension of files (.off) written by this writer.
/*!
 *
 */
QStringList gsWriter::handled(void) const
{
    return QStringList() << ".off";
}



//! Register the writer in the factory.
/*!
 *
 */
bool gsWriter::registered(void)
{
    return dtkAbstractDataFactory::instance()->registerDataWriterType("gsWriter", QStringList(), creategsWriter);
}


//! Check if the writer can write the file.
/*! file must have extension .off
 *
 */
bool gsWriter::canWrite(const QString& file){
    return (file.endsWith(".off"));
}

//! Write the .off file.
/*!
 *
 */
bool gsWriter::write(const QString& file)
{

    return true;
}



dtkAbstractDataWriter *creategsWriter(void)
{
    return new gsWriter;
}

