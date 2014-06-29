/**
  * @author Сергей Чопоров
  * @date 13/04/2014
  * @version 1.0.1
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
    enum ColorMapName{WINTER = 0, SPRING = 1, SUMMER = 2, AUTUMN = 3, GREYSCALE = 4};
    /**
     * @brief ColorValueMap Конструктор
     * @param min Левая граница диапазона значений
     * @param max Правая граница диапазона значений
     * @param colorMapName Цветовая карта
     */
    ColorValueMap(double min = -1.0, double max = 1.0, ColorMapName colorMapName = WINTER);

    double min() const;
    void setMin(double min);

    double max() const;
    void setMax(double max);

    ColorMapName colorMap() const;
    void setColorMap(const ColorMapName &colorMap);
    /**
     * @brief color Получить цвет, соответствующий значению
     * @param value Значение, для которого необходимо найти цвет
     * @return Цвет, соответствующий значению
     */
    QColor color(const double &value);

protected:
    /**
     * @brief scaleOn Отобразить значение из интеравала [min_; max_] на интервал [a; b]
     * @param value Значение для отображения
     * @param a Левая граница интервала
     * @param b Правая граница интервала
     * @return Точка из интервала [a; b]
     */
    double scaleOn(const double &value, const double &a, const double &b);

private:
    double min_;
    double max_;
    ColorMapName colorMap_;
};

#endif // COLORVALUEMAP_H
