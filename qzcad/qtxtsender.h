/**
  * @author Сергей Чопоров
  * @date 31/07/08
  * @version 1.0.0
  **/
#ifndef QTXTSENDER_H
#define QTXTSENDER_H

#include <QObject>
/**
 * @brief Класс для отправки сигналов с текстовыми сообщениями
 * @see QObject
 */
class QTxtSender : public QObject
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор по умолчанию
     * @param parent Указатель на родительский объект
     */
    explicit QTxtSender(QObject *parent = 0);

signals:
    /**
     * @brief Сигнал - изменено сообщение
     * @param text Текст сообщения
     */
    void messageChanged(const QString & text);
protected:
    /**
     * @brief Напечатать сообщение
     * @param msg Указатель на сообщение
     * @param size Длина сообщения (размер массива msg)
     */
    void printMessage(const char *msg, int size);

};

#endif // QTXTSENDER_H
