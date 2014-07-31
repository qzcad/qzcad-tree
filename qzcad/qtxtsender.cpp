#include "qtxtsender.h"
#include <QApplication>

QTxtSender::QTxtSender(QObject *parent) :
    QObject(parent)
{
}

void QTxtSender::printMessage(const char *msg, int size)
{
    emit messageChanged(QByteArray(msg, size)); // Вызов сигнала
    qApp->processEvents(); // Обработка событий приложением
}
