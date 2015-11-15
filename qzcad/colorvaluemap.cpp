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

QColor ColorValueMap::color(const double &value, int colors)
{
    QColor result;
    switch (colorMap_) {
    case GREYSCALE:
        result.setRgbF(scaleOn(value, 0.0, 1.0, colors),
                       scaleOn(value, 0.0, 1.0, colors),
                       scaleOn(value, 0.0, 1.0, colors));
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
    default:
        break;
    }
    return result;
}

double ColorValueMap::scaleOn(const double &value, const double &a, const double &b, int colors)
{
    double h = (b - a) / (double)colors;
    int j = (int)round(((b - a) * (max_ - value) / (max_ - min_)) / h);
    return a + (double)j * h;
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
