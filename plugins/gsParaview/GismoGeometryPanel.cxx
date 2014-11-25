
#include "GismoGeometryPanel.h"

#include <QLayout>
#include <QLabel>

#include "vtkSMProxy.h"
#include "pqProxy.h"

GismoGeometryPanel::GismoGeometryPanel(pqProxy* pxy, QWidget* p)
  : pqLoadedFormObjectPanel("GismoGeometryPanel.ui", pxy, p)
{
  this->layout()->addWidget(new QLabel("This is GISMO plugin", this));
}

GismoGeometryPanel::~GismoGeometryPanel()
{
}


