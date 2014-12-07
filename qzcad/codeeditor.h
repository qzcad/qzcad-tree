/**
  * @author Сергей Чопоров
  * @date 7/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef CODEEDITOR_H
#define CODEEDITOR_H

#include <QPlainTextEdit>
#include <QObject>
/**
 * @brief Класс редактора кода
 */
class CodeEditor : public QPlainTextEdit
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор
     * @param parent Указатель на родительский виджет
     */
    explicit CodeEditor(QWidget *parent = 0);
    void lineNumberAreaPaintEvent(QPaintEvent *event);
    int lineNumberAreaWidth();
    int tabSymbols() const;

    bool autoindent() const;


public slots:
    void setTabSymbols(int tabSymbols);
    void setAutoindent(bool autoindent);

protected:
    void resizeEvent(QResizeEvent *e);
    void keyReleaseEvent(QKeyEvent *e);
private slots:
    void updateLineNumberAreaWidth(int newBlockCount);
    void highlightCurrentLine();
    void updateLineNumberArea(const QRect &rect, int dy);
    void setTabStopSymbols(int symbols);
private:
    QWidget *lineNumberArea_;
    QPair<int, int> countCache_;
    bool autoindent_; //!< Включатель авто отступа
    int tabSymbols_; //!< Количество символов - ширина табуляции
};
/**
 * @brief Класс площади для вывода номеров строк
 */
class LineNumberArea : public QWidget
{
public:
    LineNumberArea(CodeEditor *editor) : QWidget(editor) {
        codeEditor_ = editor;
    }

    QSize sizeHint() const {
        return QSize(codeEditor_->lineNumberAreaWidth(), 0);
    }

protected:
    void paintEvent(QPaintEvent *event) {
        codeEditor_->lineNumberAreaPaintEvent(event);
    }

private:
    CodeEditor *codeEditor_;
};

#endif // CODEEDITOR_H
