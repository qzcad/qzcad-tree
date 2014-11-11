/**
  * @author Сергей Чопоров
  * @date 12/10/2014
  * @version 1.0.0
  * */
#ifndef ROWDOUBLEMATRIX_H
#define ROWDOUBLEMATRIX_H

#include "mtx_types.h"
#include "doublevector.h"
#include "mappeddoublematrix.h"

namespace mtx {
/**
 * @brief Класс построчного сжатого хранения матриц
 */
class RowDoubleMatrix
{
public:
    /**
     * @brief Конструктор по умолчанию. Создает пустую матрицу.
     */
    RowDoubleMatrix();
    /**
     * @brief Конструктор создает копию матрицы с построчной очисткой
     * @param M Экземпляр матрицы для копирования
     */
    RowDoubleMatrix(MappedDoubleMatrix &M);
    ~RowDoubleMatrix();
    /**
     * @brief Метод создает копию матрицы с построчной очисткой
     * @param M Экземпляр матрицы для копирования
     */
    void assign(MappedDoubleMatrix &M);
    /**
     * @brief Метод для получения результата умножения текущей матрицы на вектор (процедурный стиль)
     * @param dv Вектор, на который умножается матрица
     * @param res Ссылка для записи результата умножения матрицы на вектор
     */
    void product(const DoubleVector &dv, DoubleVector &res) const;
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
protected:
    /**
     * @brief Метод для очистки памяти
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
    /**
     * @brief Структура для хранения пар номер_столбца-значение
     */
    typedef struct
    {
        size_type column; //!< номер столбца
        double value; //!< значение матрицы
    } ColumnValuePair;
    size_type size_; //!< размер матрицы
    size_type *row_size_; //!< массив размеров строки
    ColumnValuePair **data_; //!< матрица данных (ленточного типа)
};

}

#endif // ROWDOUBLEMATRIX_H
