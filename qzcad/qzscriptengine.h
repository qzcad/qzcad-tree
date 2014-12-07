#ifndef QZSCRIPTENGINE_H
#define QZSCRIPTENGINE_H

#include <QScriptEngine>

class QZScriptEngine : public QScriptEngine
{
    Q_OBJECT
public:
    explicit QZScriptEngine(QObject *parent = 0);

signals:

public slots:

};

#endif // QZSCRIPTENGINE_H
