/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef POLYGONALMODELDIALOG_H
#define POLYGONALMODELDIALOG_H

#include <QDialog>
#include <QDomDocument>
#include <QTreeWidgetItem>

#include "highlighter.h"

#include "point2d.h"

namespace Ui {
class PolygonalModelDialog;
}

class PolygonalModelDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit PolygonalModelDialog(QWidget *parent = 0);
    ~PolygonalModelDialog();
    /**
     * @brief quadsCount Количество четырехугольных областей
     * @return Количество четырехугольных областей в моделе
     */
    int quadsCount() const;
    /**
     * @brief quad Четрыхугольная область с заданным номером
     * @param num Номер области
     * @param xiCount Количество узлов по первому направлению
     * @param etaCount Количество узлов по второму напрвлению
     * @return Массив координат узлов
     */
    std::vector<msh::Point2D> quad(int num, int &xiCount, int &etaCount) const;
    /**
     * @brief trianglesCount Количество треугольных областей
     * @return Количество треугольных областей в моделе
     */
    int trianglesCount() const;
    /**
     * @brief triangle Треугольная область с заданным номером
     * @param num Номер области
     * @param nodesCount Количество узлов на сторону области
     * @return Массив координат узлов
     */
    std::vector<msh::Point2D> triangle(int num, int &nodesCount) const;

public slots:
    void addQuad(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, int xiCount, int etaCount);
    void addTriangle(double x0, double y0, double x1, double y1, double x2, double y2, int nodesCount);

private slots:
    void on_newModelButton_clicked();

    void on_openModelButton_clicked();

    void on_saveModelButton_clicked();

    void on_treeWidget_itemDoubleClicked(QTreeWidgetItem *item, int column);

    void on_addQuadButton_clicked();

    void on_addTriangleButton_clicked();

    void on_buildModelButton_clicked();

    void on_newParametricModelButton_clicked();

    void on_openParametricModelButton_clicked();

    void on_saveParametricModelButton_clicked();

private:
    bool readFile(const QString &fileName);
    void parsePolygons(const QDomElement &element);
    void parseTriangle(const QDomElement &triangle, QTreeWidgetItem *parent);
    void parseQuad(const QDomElement &quad, QTreeWidgetItem *parent);
    void parseVertex(const QDomElement &vertex, QTreeWidgetItem *parent);
    bool writeFile(const QString &fileName);
    void clearAll();
    bool writeParametricModel(const QString &fileName);
    void saveAsParametricModel();
    void saveParametricModel();
    
private:
    Ui::PolygonalModelDialog *ui;
    Highlighter *highlighter;
    QString scriptFileName;
};

#endif // POLYGONALMODELDIALOG_H
