#ifndef GSGEOMETRYCONVERTER_H
#define GSGEOMETRYCONVERTER_H

#include <axlCore/axlAbstractDataConverter.h>

class gsBasisConverterPrivate;


// public axlAbstactData
class gsBasisConverter : public axlAbstractDataConverter
{
    Q_OBJECT

public:
     gsBasisConverter(void);
    ~gsBasisConverter(void);

    QString  description (void) const;
    QStringList fromTypes(void) const;
    QString       toType (void) const;

public:
    static bool registered(void);

public slots:
    axlMesh *toMesh(void);

public:
    void setData(dtkAbstractData *data);

private:
    gsBasisConverterPrivate *d;
};

dtkAbstractDataConverter *createGoGeometryConverter(void);

#endif
