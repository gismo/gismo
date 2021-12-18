
#ifndef GOSURFACEBSPLINEREADER_H
#define GOSURFACEBSPLINEREADER_H

#include <axlCore/axlAbstractDataReader.h>

#include "gsAxelPluginExport.h"

class dtkAbstractData;

class GSAXELPLUGIN_EXPORT gsNativeAxelReader : public axlAbstractDataReader
{
    Q_OBJECT

public :
      gsNativeAxelReader(void);
     ~gsNativeAxelReader(void);

public:
    QString identifier(void) const;
    QString description(void) const;
    QStringList handled(void) const;

    static bool registered(void);

public:
    bool accept(const QDomNode& node);
    bool reject(const QDomNode& node);
    
    axlAbstractData *read(const QDomNode& node);
};

dtkAbstractDataReader *creategsNativeAxelReader(void);

#endif //GOSURFACEBSPLINEREADER_H
