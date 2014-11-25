#ifndef GSGEOMETRYCONVERTER_H
#define GSGEOMETRYCONVERTER_H

#include <axlCore/axlAbstractDataConverter.h>

class gsGeometryConverterPrivate;

class gsGeometryConverter : public axlAbstractDataConverter
{
    Q_OBJECT

public:
     gsGeometryConverter(void);
    ~gsGeometryConverter(void);

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
    gsGeometryConverterPrivate *d;
};

dtkAbstractDataConverter *createGoGeometryConverter(void);

#endif
