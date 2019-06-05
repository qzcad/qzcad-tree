/**
  * @author Сергей Чопоров
  * @date 02/01/2015
  * @version 1.0.0
  * @copyright Copyright 2015 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef CONSOLEPROGRESS_H
#define CONSOLEPROGRESS_H

#include <string>
#include <iostream>

using namespace std;

/**
 * @brief Класс ConsoleProgress показывает индикацию прогресса на основе текстового представления путем вывода в указанный поток.
 * За основу взята реализация из проекта boost (http://www.boost.org/doc/libs/1_57_0/boost/progress.hpp)
 * Зависимости сокращены до стандартной библиотеки c++
 */
class ConsoleProgress
{
public:
    /**
     * @brief Конструктор
     * @param expectedCount Максимальное (ожидаемое) количество итераций
     * @param stream Ссылка на поток для вывода прогресса
     * @param s1 Строка, которая будет выведена перед цифрами индикации
     * @param s2 Строка, которая будет выведена перед отметками индикации
     * @param s3 Строка, которая будет выведена перед заполнением прогресса
     * @param spacer Символ-заполнитель прогресса
     */
    explicit ConsoleProgress(unsigned long expectedCount,
                             ostream &stream = cout,
                             const string &s1 = "\n",
                             const string &s2 = "",
                             const string &s3 = "",
                             char spacer = '*');
    /**
     * @brief Перезапустить индикатор прогресса
     * @param expectedCount Максимальное (ожидаемое) количество итераций
     * Индикатор прогресса сбрасывается в исходное (начальное) состояние и заново выводится в поток.
     */
    void restart(unsigned long expectedCount);
    /**
     * @brief Увеличивает количество выполненных итераций на заданный шаг. Если необходимо, то выводит новый символ в строку прогресса.
     * @param increment Шаг для прироста итераций
     * @return Актуальное количество итераций
     */
    unsigned long operator+=(unsigned long increment);
    /**
     * @brief Инкрементация количества итераций. Если необходимо, то выводит новый символ в строку прогресса.
     * @return  Актуальное количество итераций
     */
    unsigned long operator++();
    /**
     * @brief Получить количество выполненных итераций
     * @return Количество выполненных итераций
     */
    unsigned long count() const;
    /**
     * @brief Получить максимальное количество итераций
     * @return Максимальное количество итераций
     */
    unsigned long expectedCount() const;
    /**
     * @brief Проверить выполнено ли ожидаемое количество итераций
     * @return true, если выполнено ожидаемое число итераций, false - иначе
     */
    bool isExpectedCount() const;
    /**
     * @brief Увеличивает количество выполненных итераций на заданный шаг. Если необходимо, то выводит новый символ в строку прогресса.
     * @param increment Шаг для прироста итераций (по умолчанию 1)
     * @return Актуальное количество итераций
     */
    unsigned long inc(unsigned long increment=1);
protected:
    void displayTic();
private:
    ostream &stream_; //!< Поток для вывода
    const string s1_; //!< Строка, которая будет выведена перед цифрами
    const string s2_; //!< Строка, которая будет выведена перед отметками шкалы
    const string s3_; //!< Строка, которая будет выведена перед заполнением прогресса
    char spacer_; //!< Заполнитель индикации прогресса
    unsigned long expectedCount_; //!< Ожидаемое количество итераций во время прогресса
    unsigned long count_; //!< Пройденное количество итераций
    unsigned long nextTicCount_; //!< Количество итераций до следующего тика индикатора
    unsigned int tic_; //!< Количество тиков индикатора
};

#endif // CONSOLEPROGRESS_H
