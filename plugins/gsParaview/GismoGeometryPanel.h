
#ifndef _GismoGeometryPanel_h
#define _GismoGeometryPanel_h

#include "pqLoadedFormObjectPanel.h"
#include "pqObjectPanelInterface.h"

class GismoGeometryPanel : public pqLoadedFormObjectPanel
{
  Q_OBJECT
public:
  GismoGeometryPanel(pqProxy* proxy, QWidget* p);
  ~GismoGeometryPanel();
};

#endif

