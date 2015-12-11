#ifndef NAMEDDOUBLEVECTOR_H
#define NAMEDDOUBLEVECTOR_H

#include <vector>
#include <string>

namespace msh {
/**
 * @brief Класс для хранения именованных векторов данных
 */
class NamedDoubleVector
{
public:
    /**
     * @brief Конструктор
     * @param name Имя
     * @param vector Значения
     */
    NamedDoubleVector(const std::string &name, const std::vector<double> &vector);
    // Getters and setters
    std::vector<double> vector() const;
    void setVector(const std::vector<double> &vector);

    std::string name() const;
    void setName(const std::string &name);

private:
    std::vector<double> vector_;
    std::string name_;
};

}

#endif // NAMEDDOUBLEVECTOR_H
