/**
  * @author Сергей Чопоров
  * @date 24/09/2014
  * @version 1.0.0
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  * */
#ifndef MAPPEDDOUBLEMATRIX_H
#define MAPPEDDOUBLEMATRIX_H
#include <map>
#include "mtx_types.h"
#include "doublevector.h"
namespace mtx {
/**
 * @brief Сжатый массив на основе отображения ключ-значение
 */
typedef std::map<size_type, double> MappedDoubleVector;
/**
 * @brief Класс для работы с разреженными квадтратными матрицами
 * Матрица храница построчно. Строки сжимаются путем использования отображения ключ-значение
 * @see std::map
 */
class MappedDoubleMatrix
{
public:
    /**
     * @brief Конструктор по умолчанию
     */
    MappedDoubleMatrix();
    /**
     * @brief Конструктор для создания матрицы заданного размера
     * @param size Размер матрицы
     */
    MappedDoubleMatrix(size_type size);
    /**
     * @brief Конструктор копирования
     * @param mdm Экземпляр объекта для копирования
     */
    MappedDoubleMatrix(const MappedDoubleMatrix &mdm);
    ~MappedDoubleMatrix();
    /**
     * @brief Метод для получения размера матрицы
     * @return Размер матрицы
     */
    size_type size() const;
    /**
     * @brief Метод для получения количества не нулевых (хранящихся) элементов в строке
     * @param i Номер строки
     * @return Количество элементов, под которые выделена память в строке
     */
    size_type size(size_type i) const;
    /**
     * @brief Метод для получения значения элемента, заданного индексами
     * @param i Номер строки
     * @param j Номер столбца
     * @return Значение элемента, заданного индексами
     */
    double data(size_type i, size_type j) const;
    /**
     * @brief Метод для изменения размера матрицы
     * @param size Новый размер матрицы
     */
    void resize(size_type size);
    /**
     * @brief Оператор для доступа по индексу
     * @param i Индекс строки
     * @return Ссылку на вектор данных строки
     * Использование квадратных скобок для доступа по индексу вернет ссылку на значение, если оно было определено ранее.
     * Если значение с указанным индексом не было определено ранее, то будет возвращена ссылка новый элемент.
     * @code{.cpp}
     * MappedDoubleVector mpd(5);
     * mpd[0][1] = 1.3; // создает новый объект, которому присваевается значение 1.3
     * @endcode
     */
    MappedDoubleVector &operator [](size_type i);
    /**
     * @brief Оператор для доступа по паре индексов
     * @param i Номер строки
     * @param j Номер стобца
     * @return Ссылку на значение с заданными индексами
     * Если значение с указанным индексом не было определено ранее, то будет возвращена ссылка новый элемент.
     * @code{.cpp}
     * MappedDoubleVector mpd(5);
     * mpd(0, 1) = 1.3; // создает новый объект, которому присваевается значение 1.3
     * @endcode
     */
    reference operator () (size_type i, size_type j);
    /**
     * @brief Оператор присвоения
     * @param mdm Экземпляр объекта для копирования
     * @return Ссылку на копию
     */
    MappedDoubleMatrix &operator =(const MappedDoubleMatrix &mdm);
    /**
     * @brief print Метод для печати матрицы на консоль (стандартный вывод)
     * @param separator Разделитель столбцов
     */
    void print(char separator = '\t') const;
    /**
     * @brief Оператор умножения матрицы на столбец
     * @param mdm Экземпляр матрицы
     * @param dv Экземпляр вектора
     * @return Вектор mul = mdm * dv
     * Умножение "быстрое": перебераются только ненулевые элементы строки
     */
    friend DoubleVector operator *(const MappedDoubleMatrix &mdm, const DoubleVector dv);
    /**
     * @brief Метод сопряженных градиентов для решения СЛАУ, определенного матрицей и свободным вектором-столбцом
     * @param B Вектор-столбец СЛАУ
     * @param epsilon Точность меотда
     * @param niter Максимальное количество итераций
     * @param printMessages Флаг пеачти сообщений
     * @param messageStep Щаг печати сообщений (через сколько итеаций выводить сообщение
     * @return Решение СЛАУ
     */
    DoubleVector conjugateGradient(const DoubleVector &B, double epsilon = 1.0E-6, unsigned int niter = 100000, bool printMessages = true, unsigned int messageStep = 25) const;
    /**
     * @brief Метод для получения результата умножения текущей матрицы на вектор
     * @param dv Вектор, на который умножается матрица
     * @return Результат умножения матрицы на вектор
     */
    DoubleVector product(const DoubleVector &dv) const;
    /**
     * @brief Метод для получения результата умножения текущей матрицы на вектор (процедурный стиль)
     * @param dv Вектор, на который умножается матрица
     * @param res Ссылка для записи результата умножения матрицы на вектор
     */
    void product(const DoubleVector &dv, DoubleVector &res) const;
    /**
     * @brief Метод для обнуления указанной строки матрицы
     * @param i Номер строки
     */
    void zeroRow(size_type i);
    /**
     * @brief Метод для обнуления указанного столбца
     * @param j Номер столбца
     */
    void zeroCol(size_type j);
    /**
     * @brief Метод симметричного обнуления
     * @param i строка(столбец) для обнуления
     */
    void zeroSym(size_type i);
    /**
     * @brief Метод для очистки данных заданной строки
     * @param i Номер строки для очистки
     */
    void clear(size_type i);
    /**
     * @brief Метод для получения константного итератора строки матрицы, указывающего на первый элемент
     * @param i Номер строки
     * @return Константный итератор
     */
    MappedDoubleVector::const_iterator begin(size_type i) const;
    /**
     * @brief end Указатель на "гипотетический" последний элемент (адрес) строки матрицы (итератор)
     * @param i Номер строки
     * @return Константный итератор
     */
    MappedDoubleVector::const_iterator end(size_type i) const;
    /**
     * @brief Решение СЛАУ методом Холецкого (трбует для работы n * (n - 1) / 2 паямти)
     * @param B Вектор столбец
     * @return Вектор корней СЛАУ
     */
    DoubleVector cholesky(DoubleVector &B);
protected:
    /**
     * @brief Метод для выделения памяти под строки матрицы
     * @param size Размер матрицы
     */
    void alloc(size_type size);
    /**
     * @brief Метод для очистки использованной памяти
     */
    void clear();
    /**
     * @brief Метод для нахождения невязки для решения СЛАУ, определенного матрицей
     * @param X Корни СЛАУ
     * @param B Вектор-столбец СЛАУ
     * @param R Невязка - выходной параметр
     */
    void residual(const DoubleVector &X, const DoubleVector &B, DoubleVector &R) const;
private:
    MappedDoubleVector *data_; //!< Массив векторов для хранения строк матрицы
    size_type size_; //!< Размер матрицы
};
}


#endif // MAPPEDDOUBLEMATRIX_H
