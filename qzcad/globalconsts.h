/**
  * Модуль для хранения глобальных констант управления GUI (устойчивых к изменеиям данным)
  * @author Чопоров Сергей
  * @date 30/12/2013
  * @version 1.0
  */
#ifndef GLOBALCONSTS_H
#define GLOBALCONSTS_H
// Коды категорий
// от 100 до 199 - операции
#define _OPERATION_MIN_ID_ 100
#define _OPERATION_MAX_ID_ 199
#define _UNION_ 100 // объединение
#define _INTERSECTION_ 101 // пересечение
#define _DIFFERENCE_ 102 // разность
#define _NEGATION_ 103 // отрицание
// от 200 до 999 - примитивы
#define _PRIMITIVE_MIN_ID_ 200
#define _PRIMITIVE_MAX_ID_ 999
#define _PATH_ 200 // свободный контур
// от 1000 - модели пользователя
#define _MODEL_ 1000
#endif // GLOBALCONSTS_H
