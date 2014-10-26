/**
  * @author Сергей Чопоров
  * @date 23/09/2014
  * @version 1.0.0
  * */
#ifndef DOUBLEMATRIX_H
#define DOUBLEMATRIX_H

#include "doublevector.h"

namespace mtx {
/**
 * @brief Класс для хранения матриц действительных чисел (двойная точность)
 */
class DoubleMatrix
{
public:
    /**
     * @brief Конструктор по умолчанию создает матрицу 0x0
     */
    DoubleMatrix();
    /**
     * @brief Конструктор для создания квадратной матрицы
     * @param size Размер матрицы
     */
    DoubleMatrix(size_type size);
    /**
     * @brief Конструктор для создания квадратной матрицы с инициализацией элементов
     * @param size Размер матрицы
     * @param initValue Начальное значение
     */
    DoubleMatrix(size_type size, const_reference initValue);
    /**
     * @brief Конструктор для создания прямоугольной матрицы
     * @param rows Количество строк
     * @param cols Количество столбцов
     */
    DoubleMatrix(size_type rows, size_type cols);
    /**
     * @brief Конструктор для создания прямоугольной матрицы с инициализацией элементов
     * @param rows Количество строк
     * @param cols Количество столбцов
     * @param initValue Начальное значение
     */
    DoubleMatrix(size_type rows, size_type cols, const_reference initValue);
    /**
     * @brief Конструктор копирования
     * @param dm Экземпляр объекта для копирования
     */
    DoubleMatrix(const DoubleMatrix &dm);
    ~DoubleMatrix();
    /**
     * @brief Метод для получения количества строк в матрице
     * @return Количество строк в матрице
     */
    size_type rowCount() const;
    /**
     * @brief Метод для получения количества столбцов в матрице
     * @return Количество столбцов в матрице
     */
    size_type colCount() const;
    /**
     * @brief Метод для присвоения всем элементам матрицы заданного значения
     * @param initValue Значение, которое необходимо присвоить элементам матрицы
     */
    void set(const_reference initValue);
    /**
     * @brief Метод для присвоения всем элементам заданной строки матрицы указанного значения
     * @param i Номер строки
     * @param value Значение присваемое строке матрицы
     */
    void setRow(size_type i, const_reference value);
    /**
     * @brief Метод для присвоения всем элементам заданного солбца матрицы указанного значения
     * @param j Номер столбца
     * @param value Значение присваемое столбцу матрицы
     */
    void setCol(size_type j, const_reference value);
    /**
     * @brief Метод для получения значения элемента, заданного индексами
     * @param i Номер строки
     * @param j Номер столбца
     * @return Значение элемента, заданного индексами
     */
    double data(size_type i, size_type j) const;
    /**
     * @brief Метод для изменения размеров матрицы
     * @param rows Количество строк
     * @param cols Количество столбцов
     */
    void resize(size_type rows, size_type cols);
    /**
     * @brief Оператор доступа по индексу
     * @param i Номер строки
     * @return Ссылка на элементы строки
     */
    DoubleVector &operator [](size_type i);
    /**
     * @brief Оператор доступа по индексу (для константных ссылок-параметров)
     * @param i Номер строки
     * @return Константная ссылка на элементы строки
     */
    const DoubleVector &operator [](size_type i) const;
    /**
     * @brief Оператор доступа по паре индексов
     * @param i Номер строки
     * @param j Номер столбца
     * @return Ссылка на элемент с индексами i, j
     */
    reference operator ()(size_type i, size_type j);
    /**
     * @brief Оператор доступа по паре индексов (для константных ссылок-параметров)
     * @param i Номер строки
     * @param j Номер столбца
     * @return Константная ссылка на элемент с индексами i, j
     */
    const_reference operator ()(size_type i, size_type j) const;
    /**
     * @brief Оператор присвоения
     * @param dm Экземпляр объекта для присвоения
     * @return Ссылка на копию
     */
    DoubleMatrix &operator =(const DoubleMatrix &dm);
    /**
     * @brief Оператор присвоения с прибавлением
     * @param dm Матрица, которою нужно прибавить
     * @return Результат суммирования
     */
    DoubleMatrix &operator +=(const DoubleMatrix &dm);
    /**
     * @brief Оператор для умножения матрицы на число
     * @param d Число
     * @param dm Матрица
     * @return Результат поэлементного умножения
     */
    friend DoubleMatrix operator *(const double &d, const DoubleMatrix &dm);
    /**
     * @brief print Метод для печати матрицы на консоль (стандартный вывод)
     * @param separator Разделитель столбцов
     */
    void print(char separator = '\t') const;
    /**
     * @brief Метод для получения траспонированной копии матрицы
     * @return Транспонированная копия матрицы
     */
    DoubleMatrix transpose();
    /**
     * @brief Оператор умножения матриц
     * @param a Первый операнд произведения
     * @param b Второй операнд произведения
     * @return Матрица c = a * b
     */
    friend DoubleMatrix operator * (const DoubleMatrix &a, const DoubleMatrix &b);
    /**
     * @brief Оператор умножения матрицы на вектор
     * @param a Матрица
     * @param b Вектор
     * @return Вектор c = a * b
     */
    friend DoubleVector operator * (const DoubleMatrix &a, const DoubleVector &b);
    /**
     * @brief Метод для вычисления определеителя в частном случае матриц 2x2
     * @return Значение определителя
     */
    double det2x2() const;
    /**
     * @brief Метод для вычисления определеителя в частном случае матриц 3x3
     * @return Значение определителя
     */
    double det3x3() const;
    /**
     * @brief Метод для получения обратной матрицы в частном случае 2x2
     * @return Обратная матрица
     */
    DoubleMatrix inverted2x2() const;
    /**
     * @brief Метод для получения обратной матрицы в частном случае 3x3
     * @return Обратная матрица
     */
    DoubleMatrix inverted3x3() const;
    /**
     * @brief Метод для получения обратной матрицы методом Жордана-Гаусса
     * @return Обратная матрица
     */
    DoubleMatrix inverted() const;
protected:
    /**
     * @brief Метод выделения памяти для хранения матрицы
     * @param rows Количество строк
     * @param cols Количество столбцов
     */
    void alloc(size_type rows, size_type cols);
    /**
     * @brief Очистить память
     */
    void clear();
private:
    DoubleVector *data_; //!< Массив векторов для хранения строк матрицы
    size_type rowCount_; //!< Количество строк в матрице
};
}

#endif // DOUBLEMATRIX_H
