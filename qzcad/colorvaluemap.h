/**
  * @author Сергей Чопоров
  * @date 13/04/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef COLORVALUEMAP_H
#define COLORVALUEMAP_H

#include <QColor>
/**
 * @brief Класс ColorValueMap служит для получения цвета, соответсвующего точке.
 * Схема получения цвета точки определятеся цветовой картой
 */
class ColorValueMap
{
public:
    /**
      * @enum Перечень доступных цветовых карт
      */
    enum ColorMapName{WINTER = 0, SPRING = 1, SUMMER = 2, AUTUMN = 3, GREYSCALE = 4, HOT = 5, COOL = 6};
    /**
     * @brief Конструктор
     * @param min Левая граница диапазона значений
     * @param max Правая граница диапазона значений
     * @param colorMapName Цветовая карта
     */
    ColorValueMap(double min = -1.0, double max = 1.0, ColorMapName colorMapName = WINTER);
    /**
     * @brief Получить минимально допустимое значение
     * @return Минимально допустимое значение
     */
    double min() const;
    /**
     * @brief Установить минимально допустимое значение
     * @param min Новое минимально допустимое значение
     */
    void setMin(double min);
    /**
     * @brief Получить максимально допустимое значение
     * @return Максимально допустимое значение
     */
    double max() const;
    /**
     * @brief Установить максимально допустимое значение
     * @param max Новое максимально допустимое значение
     */
    void setMax(double max);
    /**
     * @brief Получить цветову схему (карту)
     * @return Значение индикатора цветовой схемы (карты)
     */
    ColorMapName colorMap() const;
    /**
     * @brief Задать цветовую схему (карту)
     * @param colorMap
     */
    void setColorMap(const ColorMapName &colorMap);
    /**
     * @brief Получить цвет, соответствующий значению
     * @param value Значение, для которого необходимо найти цвет
     * @return Цвет, соответствующий значению
     */
    QColor color(const double &value, int colors = 8);

protected:
    /**
     * @brief Отобразить значение из интеравала [min_; max_] на интервал [a; b]
     * @param value Значение для отображения
     * @param a Левая граница интервала
     * @param b Правая граница интервала
     * @param colors Количество цветов (значение по умолчанию: 8)
     * @return Точка из интервала [a; b]
     */
    double scaleOn(const double &value, const double &a, const double &b, int colors = 8);
    /**
     * @brief Отобразить значение из интеравала [min_; max_] на интервал [a; b]
     * @param value Значение для отображения
     * @param min Левое дополнительное ограничение: если value < min, то функция вернет a
     * @param max Правое дополнительное ограничение: если value > max, то функция вернет b
     * @param a Левая граница интервала
     * @param b Правая граница интервала
     * @param colors Количество цветов (значение по умолчанию: 8)
     * @return Точка из интервала [a; b]
     */
    double scaleOn(const double &value, const double &min, const double &max, const double &a, const double &b, int colors = 8);

private:
    double min_; //!< Минимально допустимое значение
    double max_; //!< Максимально допустимое значение
    ColorMapName colorMap_; //!< Индикатор цветовой схемы
};

#endif // COLORVALUEMAP_H
