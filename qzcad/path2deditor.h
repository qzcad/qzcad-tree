#ifndef PATH2DEDITOR_H
#define PATH2DEDITOR_H

#include <QDialog>

namespace Ui {
class Path2DEditor;
}

class Path2DEditor : public QDialog
{
    Q_OBJECT
    
public:
    explicit Path2DEditor(QWidget *parent = 0);
    ~Path2DEditor();
    
private slots:
    void on_addNodeToolButton_clicked();

private:
    Ui::Path2DEditor *ui;
};

#endif // PATH2DEDITOR_H
