
#pragma once

#include <QActionGroup>
/// This example illustrates adding a toolbar to ParaView to create a sphere and
/// a cylinder source.
class GismoToolbar : public QActionGroup
{
  Q_OBJECT
public:
  GismoToolbar(QObject* p);

public slots:
  /// Callback for each action triggerred.
  void onAction(QAction* a);
};

