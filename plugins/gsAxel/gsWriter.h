/* gsWriter.h ---
 *
 */

#ifndef AXLOFFWRITER_H
#define AXLOFFWRITER_H

#include "gsAxelPluginExport.h"

#include <dtkCore/dtkAbstractDataWriter.h>
#include <axlCore/axlCoreExport.h>


class dtkAbstractData;
class gsWriterPrivate;


class GSAXELPLUGIN_EXPORT gsWriter : public dtkAbstractDataWriter
{
    Q_OBJECT

public :
    gsWriter(void);
    ~gsWriter(void);

public:
    QString identifier(void) const;
    QString description(void) const;
    QStringList handled(void) const;

    static bool registered(void);

public:
    bool canWrite(const QString& file);
    bool write(const QString& file);


private:
    gsWriterPrivate *d;

};

dtkAbstractDataWriter *creategsWriter(void);

#endif // AXLOFFWRITER_H
