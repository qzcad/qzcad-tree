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
    NamedDoubleVector();
    /**
     * @brief Конструктор
     * @param name Имя
     * @param vector Значения
     */
    NamedDoubleVector(const std::string &name, const std::vector<double> &vector);
    /**
     * @brief Конструктор копирования
     * @param ndv Экземпляр объекта для копирования
     */
    NamedDoubleVector(const NamedDoubleVector &ndv);

    /// @brief Методы для получения и установки значений полей класса
    /// @{
    std::vector<double> vector() const;
    void setVector(const std::vector<double> &vector);

    std::string name() const;
    void setName(const std::string &name);
    /// @}

    /**
     * @brief Локальное переопределение типа счетчика в векторе
     * @see std::vector::size_type
     */
    typedef std::vector<double>::size_type size_type;
    /**
     * @brief Получить размер вектора чисел
     * @return Размер вектора (количество элементов)
     */
    size_type size() const;
    /**
     * @brief Оператор доступа по индексу к элементам ветора
     * @param i Номер элемента
     * @return Значение элемента с номером i
     */
    double &operator[](const size_type &i);
    /**
     * @brief Константный оператор доступа по индексу к элементам ветора
     * @param i Номер элемента
     * @return Значение элемента с номером i
     */
    const double &operator[](const size_type &i) const;
    /**
     * @brief Найти значение минимального элемента вектора
     * @return Значение минимального элемента вектора
     */
    double min() const;
    /**
     * @brief Найти индекс и значение минимального элемента вектора
     * @param index Индекс минимального элемента вектора
     * @return Значение минимального элемента вектора
     */
    double min(size_type &index) const;
    /**
     * @brief Найти значение максимального элемента вектора
     * @return Значение максимального элемента вектора
     */
    double max() const;
    /**
     * @brief max Найти индекс и значение максимального элемента вектора
     * @param index Индекс максимального элемента вектора
     * @return Значение максимального элемента вектора
     */
    double max(size_type &index) const;

private:
    std::vector<double> vector_;
    std::string name_;
};

}

#endif // NAMEDDOUBLEVECTOR_H
