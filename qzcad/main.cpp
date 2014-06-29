/**
  * Система автоматизированного проектирования qzcad
  * Основной модуль
  * Цель разарботки: автоматизация проктирования конструкций путем объединения контурных примитивов и R-функций
  * Результат работы системы: геометрическая + дискретная модели
  * @author Чопоров Сергей
  * @date 30/12/2013
  * @version 1.0
  */
#include "mainwindow.h"
#include <QApplication>
/**
 * @brief main Точка входа в программу (основная функция)
 * @param argc Количество параметров командной строки
 * @param argv Массив параметров командной строки
 * @return Статус завершения программы
 */
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    
    return a.exec();
}
