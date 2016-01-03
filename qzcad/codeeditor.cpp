#include <QtGui>
#include "codeeditor.h"

CodeEditor::CodeEditor(QWidget *parent) :
    QPlainTextEdit(parent)
{
    countCache_.first = -1;
    countCache_.second = -1;
    lineNumberArea_ = new LineNumberArea(this);
    tabSymbols_ = 3; // ширина по умолчани - 3 символа

    connect(this, SIGNAL(blockCountChanged(int)), this, SLOT(updateLineNumberAreaWidth(int)));
    connect(this, SIGNAL(updateRequest(QRect,int)), this, SLOT(updateLineNumberArea(QRect,int)));
    connect(this, SIGNAL(cursorPositionChanged()), this, SLOT(highlightCurrentLine()));

    updateLineNumberAreaWidth(0);
    highlightCurrentLine();
    QFont font("Monospace");
#ifdef Q_OS_WIN
    font.setStyleHint(QFont::TypeWriter);
#endif
    setFont(font);
    setTabStopWidth(fontMetrics().width(QLatin1Char(' ')) * tabSymbols_);
    autoindent_ = true;
}

void CodeEditor::lineNumberAreaPaintEvent(QPaintEvent *event)
{
    QPainter painter(lineNumberArea_);
    painter.fillRect(event->rect(), Qt::lightGray);


    QTextBlock block = firstVisibleBlock();
    int blockNumber = block.blockNumber();
    int top = (int) blockBoundingGeometry(block).translated(contentOffset()).top();
    int bottom = top + (int) blockBoundingRect(block).height();

    while (block.isValid() && top <= event->rect().bottom())
    {
        if (block.isVisible() && bottom >= event->rect().top())
        {
            QString number = QString::number(blockNumber + 1);
            painter.setPen(Qt::black);
            painter.drawText(0, top, lineNumberArea_->width(), fontMetrics().height(),
                             Qt::AlignRight, number);
        }

        block = block.next();
        top = bottom;
        bottom = top + (int) blockBoundingRect(block).height();
        ++blockNumber;
    }
}

int CodeEditor::lineNumberAreaWidth()
{
    int digits = 1;
    int max = qMax(1, blockCount());
    while (max >= 10)
    {
        max /= 10;
        ++digits;
    }

    int space = 3 + fontMetrics().width(QLatin1Char('9')) * digits;

    return space;
}

void CodeEditor::resizeEvent(QResizeEvent *e)
{
    QPlainTextEdit::resizeEvent(e);

    QRect cr = contentsRect();
    lineNumberArea_->setGeometry(QRect(cr.left(), cr.top(), lineNumberAreaWidth(), cr.height()));
}

void CodeEditor::keyReleaseEvent(QKeyEvent *e)
{
    if (autoindent_ && (e->key() == Qt::Key_Return || e->key() == Qt::Key_Enter))
    {
        QStringList lines = toPlainText().split('\n');
        int indexCurrentLine = textCursor().blockNumber();

        if (indexCurrentLine == -1)
            indexCurrentLine = lines.size()-1;

        int numberOfTabulations = lines[indexCurrentLine-1].count('\t');

        if(lines[indexCurrentLine-1].contains('{'))
            numberOfTabulations++;
        if(lines[indexCurrentLine-1].contains('}'))
            --numberOfTabulations;

        for (register int i = 0; i < numberOfTabulations; i++)
            insertPlainText(QString('\t'));
    }
}

void CodeEditor::updateLineNumberAreaWidth(int newBlockCount)
{
    setViewportMargins(lineNumberAreaWidth(), 0, 0, 0);
    Q_UNUSED(newBlockCount)
}

void CodeEditor::highlightCurrentLine()
{
    QList<QTextEdit::ExtraSelection> extraSelections;

    if (!isReadOnly())
    {
        QTextEdit::ExtraSelection selection;
        QColor lineColor = QColor(Qt::yellow).lighter(160);
        selection.format.setBackground(lineColor);
        selection.format.setProperty(QTextFormat::FullWidthSelection, true);
        selection.cursor = textCursor();
        selection.cursor.clearSelection();
        extraSelections.append(selection);
    }

    setExtraSelections(extraSelections);
}

void CodeEditor::updateLineNumberArea(const QRect &rect, int dy)
{
    if (dy)
    {
        lineNumberArea_->scroll(0, dy);
    } else if (countCache_.first != blockCount()
               || countCache_.second != textCursor().block().lineCount())
    {
        lineNumberArea_->update(0, rect.y(), lineNumberArea_->width(), rect.height());
        countCache_.first = blockCount();
        countCache_.second = textCursor().block().lineCount();
    }

    if (rect.contains(viewport()->rect()))
        updateLineNumberAreaWidth(0);
}

void CodeEditor::setTabStopSymbols(int symbols)
{
    setTabStopWidth(fontMetrics().width(QLatin1Char(' ')) * symbols);
}
bool CodeEditor::autoindent() const
{
    return autoindent_;
}

void CodeEditor::setAutoindent(bool autoindent)
{
    autoindent_ = autoindent;
}

int CodeEditor::tabSymbols() const
{
    return tabSymbols_;
}

void CodeEditor::setTabSymbols(int tabSymbols)
{
    tabSymbols_ = tabSymbols;
    setTabStopWidth(fontMetrics().width(QLatin1Char(' ')) * tabSymbols_);
}

