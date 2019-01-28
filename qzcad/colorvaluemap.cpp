#include "colorvaluemap.h"
#include <math.h>

ColorValueMap::ColorValueMap(double min, double max, ColorMapName colorMapName)
{
    min_ = min;
    max_ = max;
    if (min_ == max_) min_ -= 1.0; // защита от деления на ноль
    colorMap_ = colorMapName;
}

ColorValueMap::ColorMapName ColorValueMap::colorMap() const
{
    return colorMap_;
}

void ColorValueMap::setColorMap(const ColorMapName &colorMap)
{
    colorMap_ = colorMap;
}

QColor ColorValueMap::color(const double &value, int colors)
{
    QColor result;
    switch (colorMap_) {
    case GREYSCALE:
        result.setRgbF(scaleOn(value, 0.1, 0.9, colors),
                       scaleOn(value, 0.1, 0.9, colors),
                       scaleOn(value, 0.1, 0.9, colors));
        break;
    case SPRING:
        result.setRgbF(1.0,
                       scaleOn(value, 0.0, 1.0, colors),
                       scaleOn(value, 1.0, 0.0, colors));
        break;
    case SUMMER:
        result.setRgbF(scaleOn(value, 0.0, 1.0, colors),
                       scaleOn(value, 0.5, 1.0, colors),
                       0.4);
        break;
    case AUTUMN:
        result.setRgbF(1.0,
                       scaleOn(value, 0.0, 1.0, colors),
                       0.0);
        break;
    case WINTER:
        result.setRgbF(0.0,
                       scaleOn(value, 0.0, 1.0, colors),
                       scaleOn(value, 1.0, 0.5, colors));
        break;
    case HOT:
        result.setRgbF(scaleOn(value, min_, min_ + (max_ - min_) / 3.0, 0.0, 1.0, colors / 2),
                       scaleOn(value, min_ + (max_ - min_) / 3.0, min_ + 2.0 * (max_ - min_) / 3.0, 0.0, 1.0, colors / 2),
                       scaleOn(value, min_ + 2.0 * (max_ - min_) / 3.0, max_, 0.0, 1.0, colors / 2));
        break;
    case COOL:
        result.setRgbF(scaleOn(value, 0.0, 1.0, colors),
                       scaleOn(value, 1.0, 0.0, colors),
                       1.0);
        break;
    default:
        break;
    }
    return result;
}

double ColorValueMap::scaleOn(const double &value, const double &a, const double &b, int colors)
{
    double h = (b - a) / static_cast<double>(colors);
    int j = static_cast<int>(round(((b - a) * (value - min_) / (max_ - min_)) / h));
    return a + static_cast<double>(j) * h;
}

double ColorValueMap::scaleOn(const double &value, const double &min, const double &max, const double &a, const double &b, int colors)
{
    if (value < min) return a;
    if (value > max) return b;
    double h = (b - a) / static_cast<double>(colors);
    int j = static_cast<int>(round(((b - a) * (value - min) / (max - min)) / h));
    return a + static_cast<double>(j) * h;
}

double ColorValueMap::max() const
{
    return max_;
}

void ColorValueMap::setMax(double max)
{
    max_ = max;
    if (min_ == max_) min_ -= 1.0; // защита от деления на ноль
}

double ColorValueMap::min() const
{
    return min_;
}

void ColorValueMap::setMin(double min)
{
    min_ = min;
    if (min_ == max_) max_ += 1.0; // защита от деления на ноль
}
