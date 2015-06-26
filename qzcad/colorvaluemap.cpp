#include "colorvaluemap.h"

ColorValueMap::ColorValueMap(double min, double max, ColorMapName colorMapName)
{
    min_ = min;
    max_ = max;
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

QColor ColorValueMap::color(const double &value)
{
    QColor result;
    switch (colorMap_) {
    case GREYSCALE:
        result.setRgbF(scaleOn(value, 0.0, 1.0), scaleOn(value, 0.0, 1.0), scaleOn(value, 0.0, 1.0));
        break;
    case SPRING:
        result.setRgbF(1.0, scaleOn(value, 0.0, 1.0), scaleOn(value, 1.0, 0.0));
        break;
    case SUMMER:
        result.setRgbF(scaleOn(value, 0.0, 1.0), scaleOn(value, 0.5, 1.0), 0.4);
        break;
    case AUTUMN:
        result.setRgbF(1.0, scaleOn(value, 0.0, 1.0), 0.0);
        break;
    case WINTER:
        result.setRgbF(0.0, scaleOn(value, 0.0, 1.0), scaleOn(value, 1.0, 0.5));
        break;
    default:
        break;
    }
    return result;
}

double ColorValueMap::scaleOn(const double &value, const double &a, const double &b)
{
    int n = 15;
    double h = (b - a) / (double)(n - 1);
    int j = (int)round(((b - a) * (max_ - value) / (max_ - min_)) / h);
    return a + (double)j * h;
//    return a + (b - a) * (max_ - value) / (max_ - min_);
}

double ColorValueMap::max() const
{
    return max_;
}

void ColorValueMap::setMax(double max)
{
    max_ = max;
}

double ColorValueMap::min() const
{
    return min_;
}

void ColorValueMap::setMin(double min)
{
    min_ = min;
}
