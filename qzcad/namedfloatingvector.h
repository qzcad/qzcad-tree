/**
  * @author Сергей Чопоров
  * @date 12/09/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef NAMEDFLOATINGVECTOR_H
#define NAMEDFLOATINGVECTOR_H

#include <QString>
#include <vector>
/**
 * @brief Класс для хранения вектора действительных чисел с именем
 */
class NamedFloatingVector
{
public:
    /**
     * @brief Конструктор по умолчанию
     */
    NamedFloatingVector();
    /**
     * @brief Конструктор
     * @param name Имя вектор
     * @param data Вектор
     */
    NamedFloatingVector(const QString &name, const std::vector<double> &data);
    /**
     * @brief Конструктор копирования
     * @param nfv Экземпляр объекта для копирования
     */
    NamedFloatingVector(const NamedFloatingVector &nfv);
    /**
     * @brief Локальное переопределение типа счетчика в векторе
     * @see std::vector::size_type
     */
    typedef std::vector<double>::size_type size_type;
    /**
     * @brief Получить размер вектора чисел
     * @return Размер вектора (количество элементов)
     */
    size_type size();
    /**
     * @brief data Значение вектора с заданным номеров
     * @param i Номер значения
     * @return Значение с номером i
     */
    double data(const size_type &i);
    /**
     * @brief Оператор доступа по индексу к элементам ветора
     * @param i Номер элемента
     * @return Значение элемента с номером i
     */
    double &operator[](const size_type &i);
    /**
     * @brief Найти значение минимального элемента вектора
     * @return Значение минимального элемента вектора
     */
    double min() const;
    /**
     * @brief Найти значение максимального элемента вектора
     * @return Значение максимального элемента вектора
     */
    double max() const;
    /// @brief Методы для получения и установки значений полей класса
    /// @{
    QString name() const;
    void setName(const QString &name);

    std::vector<double> data() const;
    void setData(const std::vector<double> &data);
    /// @}
private:
    QString name_; //!< Имя вектора
    std::vector<double> data_; //!< Вектор действительный чисел
};

#endif // NAMEDFLOATINGVECTOR_H
