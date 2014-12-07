/**
  * @author Сергей Чопоров
  * @date 23/09/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef DOUBLEVECTOR_H
#define DOUBLEVECTOR_H

#include "mtx_types.h"

namespace mtx {
/**
 * @brief Класс для хранения массивов (векторов) действительных чисел (двойная точность)
 */
class DoubleVector
{
public:
    /**
     * @brief Конструктор по умолчанию. Массив нулевой длины
     */
    DoubleVector();
    /**
     * @brief Конструктор для формирования массива заданной длины (размера)
     * @param size Размер массива
     */
    DoubleVector(size_type size);
    /**
     * @brief Конструктор для формирования массива заданной длины (размера) с инициализацией элементов заданным значением
     * @param size Размер массива
     * @param initValue Начальное значение для инициализации элементов массива
     */
    DoubleVector(size_type size, const_reference initValue);
    /**
     * @brief Конструктор копирования
     * @param dv Экземпляр объекта для копирования
     */
    DoubleVector(const DoubleVector &dv);
    ~DoubleVector();
    /**
     * @brief Метод для получения указателя на первый элемент массива
     * @return Указатель на первый элемент массива
     */
    pointer begin();
    /**
     * @brief Метод для получения константного указателя на первый элемент массива
     * @return Константный указатель на первый элемент массива
     */
    const_pointer begin() const;
    /**
     * @brief Метод для получения указателя на гипотетический конец массива (этот указатель не моэет быть разыменован)
     * @return Указатель на гипотетический конец массива
     */
    pointer end();
    /**
     * @brief Метод для получения константного указателя на гипотетический конец массива (этот указатель не моэет быть разыменован)
     * @return Константный указатель на гипотетический конец массива
     */
    const_pointer end() const;
    /**
     * @brief Метод для получения размра массива
     * @return Размер массива
     */
    size_type size() const;
    /**
     * @brief Метод для получения копии значения элемента массива
     * @param i Индекс элемента массива
     * @return Кпоию значения элемента массива м индексом i
     */
    double data(size_type i) const;
    /**
     * @brief Присвоить всем элементам массива заданное
     * @param initValue
     */
    void set(const_reference initValue);
    /**
     * @brief Изменить размер массива
     * @param size Новый размер массива
     */
    void resize(size_type size);
    /**
     * @brief Метод для поиска значения минимального элемента
     * @return Значение минимального элемента
     */
    double min() const;
    /**
     * @brief Метод для поиска значения максимального элемента
     * @return Значение максимального элемента
     */
    double max() const;
    /**
     * @brief Вычислить эвклидову норму вектора
     * @return Эвклидова норма вектора
     */
    double norm_2() const;
    /**
     * @brief Оператор для доступа по индексу
     * @param i Индекс элемента
     * @return Ссылка на элемент с заданным индексом
     */
    reference operator ()(size_type i);
    /**
     * @brief Оператор для доступа к элементу по индексу
     * @param i Индекс элемента
     * @return Ссылка на элемент с заданным индексом
     */
    reference operator [](size_type i);
    /**
     * @brief Константный оператор для доступа к элементу по индексу
     * @param i Индекс элемента
     * @return Коснтантная ссылка на элемент с заданным индексом
     */
    const_reference operator [](size_type i) const;
    /**
     * @brief Оператор присвоения
     * @param dv Экземпляр объекта для присвоения
     * @return Ссылка на копию
     */
    DoubleVector &operator =(const DoubleVector &dv);
    /**
     * @brief Оператор для нахожения скалярного (внутреннего) произведения с заданным вектором
     * @param dv Экземляр вектора, с котором необходимо найти произведение
     * @return Значение векторного произведения
     */
    double operator *(const DoubleVector &dv) const;
    /**
     * @brief Метод для умножения всех элементов массива на число
     * @param d Число на которое умножается массив
     */
    void scale(const double &d);
    /**
     * @brief Оператор для нахождения произведения числа на вектор
     * @param d Число, которое необходимо умножить на вектор
     * @param vec Вектор для умножения
     * @return Вектор - результат поэлементного умножения на число
     */
    friend DoubleVector operator *(const double &d, const DoubleVector &vec);
    /**
     * @brief Оператор присвоения с прибавлением
     * @param vec Вектор, который необходимо присвоить с прибавлением
     * @return Ссылка на результат
     */
    DoubleVector &operator +=(const DoubleVector &vec);
    /**
     * @brief Оператор присвоения с вычитание
     * @param vec Вектор, который необходимо вычесть
     * @return Ссылка на результат
     */
    DoubleVector &operator -=(const DoubleVector &vec);
    /**
     * @brief Метод для печати элементов массива на консоль (стандартный вывод)
     * @param separator Разделитель между элементами
     */
    void print(char separator = '\t') const;
protected:
    /**
     * @brief Метод для выделения памяти под массив без инициализации элементов
     */
    void alloc();
    /**
     * @brief Метод для выделения памяти под массив с инизиализацией элементов
     * @param initValue Исходное знаяение для инициализации элементов массива
     */
    void alloc(const_reference initValue);
    /**
     * @brief Очистить массив
     */
    void clear();
private:
    pointer data_; //!< Массив данных
    size_type size_; //!< Размер (длина) массива
};
}

#endif // DOUBLEVECTOR_H
