/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef RECTMESHDIALOG_H
#define RECTMESHDIALOG_H

#include <QDialog>

namespace Ui {
class RectMeshDialog;
}

class RectMeshDialog : public QDialog
{
    Q_OBJECT

public:
    explicit RectMeshDialog(QWidget *parent = 0, int dimension = 2);
    ~RectMeshDialog();
    double xMin();
    double yMin();
    double zMin();
    double rectWidth();
    double rectHeight();
    double rectDepth();
    int xCount();
    int yCount();
    int zCount();
private:
    Ui::RectMeshDialog *ui;
};

#endif // RECTMESHDIALOG_H
