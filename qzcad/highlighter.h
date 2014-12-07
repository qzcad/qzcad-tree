/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef HIGHLIGHTER_H
#define HIGHLIGHTER_H

#include <QSyntaxHighlighter>
/**
 * @brief Класс Highlighter для подсветки синтаксиса в редакторе кода модели
 */
class Highlighter : public QSyntaxHighlighter
{
    Q_OBJECT
public:
    explicit Highlighter(QTextDocument *parent = 0);

protected:
    void highlightBlock(const QString &text);
    
private:
    /**
     * @brief Структура HighlightingRule для определения правила (шаблона и формата) подсветки синтаксиса
     */
    struct HighlightingRule
    {
        QRegExp pattern; //!< Регулярное выражения для определения шаблона
        QTextCharFormat format; //!< Формат, соответствующий шаблону
    };
    QVector<HighlightingRule> highlightingRules; //!< Массив правил подсветки синтаксиса

    QRegExp commentStartExpression; //!< Начало комментария
    QRegExp commentEndExpression; //!< Конец комментария

    QTextCharFormat keywordFormat; //!< Формат ключевых слов
    QTextCharFormat classFormat; //!< Формат классов
    QTextCharFormat singleLineCommentFormat; //!< Формат однострочного комментария
    QTextCharFormat multiLineCommentFormat; //!< Формат многострочного комментария
    QTextCharFormat quotationFormat; //!< Формат строк (выражений, заключенных в кавычки)
    QTextCharFormat functionFormat; //!< Формат функций
    
};

#endif // HIGHLIGHTER_H
