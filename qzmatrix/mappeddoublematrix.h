/**
  * @author Сергей Чопоров
  * @date 24/09/2014
  * @version 1.0.0
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
private:
    MappedDoubleVector *data_; //!< Массив векторов для хранения строк матрицы
    size_type size_; //!< Размер матрицы
};
}


#endif // MAPPEDDOUBLEMATRIX_H
