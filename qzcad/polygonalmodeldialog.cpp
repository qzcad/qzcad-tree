#include "polygonalmodeldialog.h"
#include "ui_polygonalmodeldialog.h"
#include <QFile>
#include <QFileDialog>
#include <QTextStream>
#include <iostream>
#include "qzscriptengine.h"
#include <QMessageBox>
#include "structuredisomesh2ddialog.h"
#include "pointeditordialog.h"
#include "baryquadsdialog.h"

PolygonalModelDialog::PolygonalModelDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PolygonalModelDialog)
{
    ui->setupUi(this);
    highlighter = new Highlighter(ui->textEdit->document());
    scriptFileName = "";
    setWindowFlags(Qt::Window | Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint);
    showMaximized();
}

PolygonalModelDialog::~PolygonalModelDialog()
{
    delete ui;
}

int PolygonalModelDialog::quadsCount() const
{
    int count = 0;
    for (int i = 0; i < ui->treeWidget->topLevelItemCount(); i++)
    {
        QTreeWidgetItem *item = ui->treeWidget->topLevelItem(i);
        if (item->text(0) == "quad")
        {
            count++;
        }
    }
    return count;
}

std::vector<msh::Point2D> PolygonalModelDialog::quad(int num, int &xiCount, int &etaCount) const
{
    int count = 0;
    std::vector<msh::Point2D> nodes(4);
    for (int i = 0; i < ui->treeWidget->topLevelItemCount(); i++)
    {
        QTreeWidgetItem *item = ui->treeWidget->topLevelItem(i);
        if (item->text(0) == "quad")
        {
            if (count == num)
            {

                // извлечение координат четырехугольника
                for (int j = 0; j < 4; j++)
                {
                    QStringList coords = item->child(j)->text(1).split(";");
                    double x = coords.at(0).toDouble();
                    double y = coords.at(1).toDouble();
                    nodes[j].set(x, y);
                }
                QStringList nodesList = item->text(1).split(";");
                xiCount = nodesList.at(0).toInt();
                etaCount = nodesList.at(1).toInt();
                break;
            }
            count++;
        }
    }
    return nodes;
}

int PolygonalModelDialog::trianglesCount() const
{
    int count = 0;
    for (int i = 0; i < ui->treeWidget->topLevelItemCount(); i++)
    {
        QTreeWidgetItem *item = ui->treeWidget->topLevelItem(i);
        if (item->text(0) == "triangle")
        {
            count++;
        }
    }
    return count;
}

std::vector<msh::Point2D> PolygonalModelDialog::triangle(int num, int &nodesCount) const
{
    int count = 0;
    std::vector<msh::Point2D> nodes(3);
    for (int i = 0; i < ui->treeWidget->topLevelItemCount(); i++)
    {
        QTreeWidgetItem *item = ui->treeWidget->topLevelItem(i);
        if (item->text(0) == "triangle")
        {
            if (count == num)
            {

                // извлечение координат четырехугольника
                for (int j = 0; j < 3; j++)
                {
                    QStringList coords = item->child(j)->text(1).split(";");
                    double x = coords.at(0).toDouble();
                    double y = coords.at(1).toDouble();
                    nodes[j].set(x, y);
                }
                nodesCount = item->text(1).toInt();
                break;
            }
            count++;
        }
    }
    return nodes;
}

void PolygonalModelDialog::addQuad(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, int xiCount, int etaCount)
{
    QTreeWidgetItem *item = new QTreeWidgetItem(ui->treeWidget->invisibleRootItem());
    double x[] = {x0, x1, x2, x3};
    double y[] = {y0, y1, y2, y3};
    item->setText(0, "quad");
    item->setText(1, QString("%1;%2").arg(xiCount).arg(etaCount));
    for (int i = 0; i < 4; i++)
    {
        QTreeWidgetItem *child = new QTreeWidgetItem(item);
        child->setText(0, "vertex");
        child->setText(1, QString::number(x[i], 'g', 10) + ";" + QString::number(y[i], 'g', 10));
    }
}

void PolygonalModelDialog::addTriangle(double x0, double y0, double x1, double y1, double x2, double y2, int nodesCount)
{
    QTreeWidgetItem *item = new QTreeWidgetItem(ui->treeWidget->invisibleRootItem());
    double x[] = {x0, x1, x2};
    double y[] = {y0, y1, y2};
    item->setText(0, "triangle");
    item->setText(1, QString("%1").arg(nodesCount));
    for (int i = 0; i < 3; i++)
    {
        QTreeWidgetItem *child = new QTreeWidgetItem(item);
        child->setText(0, "vertex");
        child->setText(1, QString::number(x[i], 'g', 10) + ";" + QString::number(y[i], 'g', 10));
    }
}

bool PolygonalModelDialog::readFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text))
    {
        std::cerr << "Ошибка: Невозможно прочитать файл " << qPrintable(fileName) << ": " << qPrintable(file.errorString()) << std::endl;
        return false;
    }
    QString errorStr;
    int errorLine;
    int errorColumn;
    QDomDocument doc;
    if (!doc.setContent(&file, false, &errorStr, &errorLine, &errorColumn))
    {
        std::cerr << "Ошибка: При разборе файла в сточке " << errorLine << ", столбце " << errorColumn << ": " << qPrintable(errorStr) << std::endl;
        return false;
    }
    QDomElement root = doc.documentElement();
    if (root.tagName() != "polygons")
    {
        std::cerr << "Ошибка: Файл не является полигональной моделью" << std::endl;
        return false;
    }
    parsePolygons(root);
    return true;
}

void PolygonalModelDialog::parsePolygons(const QDomElement &element)
{
    QDomNode child = element.firstChild();
    while (!child.isNull())
    {
        if (child.toElement().tagName() == "quad")
        {
            parseQuad (child.toElement(), ui->treeWidget->invisibleRootItem());
        }
        if (child.toElement().tagName() == "triangle")
        {
            parseTriangle(child.toElement(), ui->treeWidget->invisibleRootItem());
        }
        child = child.nextSibling();
    }
}

void PolygonalModelDialog::parseTriangle(const QDomElement &triangle, QTreeWidgetItem *parent)
{
    QTreeWidgetItem *item = new QTreeWidgetItem(parent);
    item->setText(0, "triangle");
    item->setText(1, triangle.attribute("nodes"));
    QDomNode vertex = triangle.firstChild();
    int count = 0;
    while (!vertex.isNull())
    {
        if (vertex.toElement().tagName() == "vertex")
        {
            parseVertex(vertex.toElement(), item);
        }
        vertex = vertex.nextSibling();
        count++;
    }
    if (count != 3)
    {
        std::cerr << "Ошибка: в треугольнике количество вершин не равно 3" << std::endl;
    }
}

void PolygonalModelDialog::parseQuad(const QDomElement &quad, QTreeWidgetItem *parent)
{
    QTreeWidgetItem *item = new QTreeWidgetItem(parent);
    item->setText(0, "quad");
    item->setText(1, quad.attribute("nodes1") + ";" + quad.attribute("nodes2"));
    QDomNode vertex = quad.firstChild();
    int count = 0;
    while (!vertex.isNull())
    {
        if (vertex.toElement().tagName() == "vertex")
        {
            parseVertex(vertex.toElement(), item);
        }
        vertex = vertex.nextSibling();
        count++;
    }
    if (count != 4)
    {
        std::cerr << "Ошибка: в четырехугольнике количество вершин не равно 4" << std::endl;
    }
}

void PolygonalModelDialog::parseVertex(const QDomElement &vertex, QTreeWidgetItem *parent)
{
    QTreeWidgetItem *item = new QTreeWidgetItem(parent);
    item->setText(0, "vertex");
    item->setText(1, vertex.attribute("x") + ";" + vertex.attribute("y"));
}

bool PolygonalModelDialog::writeFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text))
    {
        std::cerr << "Ошибка: Невозможно открыть файл для записи " << qPrintable(fileName) << ": " << qPrintable(file.errorString()) << std::endl;
        return false;
    }
    QDomDocument doc;
    QDomElement polygons = doc.createElement("polygons");
    for (int i = 0; i < ui->treeWidget->topLevelItemCount(); i++)
    {
        QTreeWidgetItem *item = ui->treeWidget->topLevelItem(i);
        QDomElement polygon = doc.createElement(item->text(0));
        // свойства многоугольника
        if (item->text(0) == "triangle")
        {
            polygon.setAttribute("nodes", item->text(1));
        }
        if (item->text(0) == "quad")
        {
            QStringList nodesList = item->text(1).split(";");
            polygon.setAttribute("nodes1", nodesList.at(0));
            polygon.setAttribute("nodes2", nodesList.at(1));
        }
        // вершины
        for (int j = 0; j < item->childCount(); j++)
        {
            QTreeWidgetItem *child = item->child(j);
            QDomElement vertex = doc.createElement(child->text(0));
            QStringList coords = child->text(1).split(";");
            vertex.setAttribute("x", coords.at(0));
            vertex.setAttribute("y", coords.at(1));
            polygon.appendChild(vertex);
        }
        polygons.appendChild(polygon);
    }
    doc.appendChild(polygons);
    QTextStream out(&file);
    out << doc.toString();
    return true;
}

void PolygonalModelDialog::clearAll()
{
    scriptFileName.clear();
    ui->textEdit->clear();
    ui->treeWidget->clear();
}

bool PolygonalModelDialog::writeParametricModel(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text))
    {
        std::cerr << "Ошибка: Невозможно открыть файл для записи " << qPrintable(fileName) << ": " << qPrintable(file.errorString()) << std::endl;
        return false;
    }
    QTextStream outStream(&file);
    outStream << ui->textEdit->toPlainText();
    ui->textEdit->document()->setModified(false);
    return true;
}

void PolygonalModelDialog::saveAsParametricModel()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Сохранить параметрическую модель"), "", tr("js-файлы (*.js)"));
    if(!fileName.isEmpty())
        if (writeParametricModel(fileName)) scriptFileName = fileName;
}

void PolygonalModelDialog::saveParametricModel()
{
    if (scriptFileName.isEmpty())
        saveAsParametricModel();
    else
        writeParametricModel(scriptFileName);
}

void PolygonalModelDialog::on_newModelButton_clicked()
{
    ui->treeWidget->clear();
}

void PolygonalModelDialog::on_openModelButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Открыть полигональную модель"), "", tr("xml-файлы (*.xml);;Любой файл (*)"));
    if(!fileName.isEmpty())
    {
        ui->treeWidget->clear();
        readFile(fileName);
    }
}

void PolygonalModelDialog::on_saveModelButton_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Сохранить полигональную модель"), "", tr("xml-файлы (*.xml)"));
    if(!fileName.isEmpty()) writeFile(fileName);
}

void PolygonalModelDialog::on_treeWidget_itemDoubleClicked(QTreeWidgetItem *item, int column)
{
    if (item->text(0) == "quad")
    {
        QStringList nodesList = item->text(1).split(";");
        int xiCount = nodesList.at(0).toInt();
        int etaCount = nodesList.at(1).toInt();

        QStringList coords0 = item->child(0)->text(1).split(";");
        double x0 = coords0.at(0).toDouble();
        double y0 = coords0.at(1).toDouble();

        QStringList coords1 = item->child(1)->text(1).split(";");
        double x1 = coords1.at(0).toDouble();
        double y1 = coords1.at(1).toDouble();

        QStringList coords2 = item->child(2)->text(1).split(";");
        double x2 = coords2.at(0).toDouble();
        double y2 = coords2.at(1).toDouble();

        QStringList coords3 = item->child(3)->text(1).split(";");
        double x3 = coords3.at(0).toDouble();
        double y3 = coords3.at(1).toDouble();

        StructuredIsoMesh2DDialog dialog(this, x0, y0, x1, y1, x2, y2, x3, y3, xiCount, etaCount);
        dialog.exec();
        if (dialog.result() == QDialog::Accepted)
        {
            x0 = dialog.x0();
            y0 = dialog.y0();
            x1 = dialog.x1();
            y1 = dialog.y1();
            x2 = dialog.x2();
            y2 = dialog.y2();
            x3 = dialog.x3();
            y3 = dialog.y3();
            xiCount = dialog.xiCount();
            etaCount = dialog.etaCount();
            item->setText(1, QString("%1;%2").arg(xiCount).arg(etaCount));
            item->child(0)->setText(1, QString::number(x0, 'g', 10) + ";" + QString::number(y0, 'g', 10));
            item->child(1)->setText(1, QString::number(x1, 'g', 10) + ";" + QString::number(y1, 'g', 10));
            item->child(2)->setText(1, QString::number(x2, 'g', 10) + ";" + QString::number(y2, 'g', 10));
            item->child(3)->setText(1, QString::number(x3, 'g', 10) + ";" + QString::number(y3, 'g', 10));
        }
    }
    else if (item->text(0) == "triangle")
    {
        int nodesCount = item->text(1).toInt();

        QStringList coords0 = item->child(0)->text(1).split(";");
        double x0 = coords0.at(0).toDouble();
        double y0 = coords0.at(1).toDouble();

        QStringList coords1 = item->child(1)->text(1).split(";");
        double x1 = coords1.at(0).toDouble();
        double y1 = coords1.at(1).toDouble();

        QStringList coords2 = item->child(2)->text(1).split(";");
        double x2 = coords2.at(0).toDouble();
        double y2 = coords2.at(1).toDouble();

        BaryQuadsDialog dialog(this, x0, y0, x1, y1, x2, y2, nodesCount);
        dialog.exec();
        if (dialog.result() == QDialog::Accepted)
        {
            x0 = dialog.x0();
            y0 = dialog.y0();
            x1 = dialog.x1();
            y1 = dialog.y1();
            x2 = dialog.x2();
            y2 = dialog.y2();
            nodesCount = dialog.nodesCount();
            item->setText(1, QString("%1").arg(nodesCount));
            item->child(0)->setText(1, QString::number(x0, 'g', 10) + ";" + QString::number(y0, 'g', 10));
            item->child(1)->setText(1, QString::number(x1, 'g', 10) + ";" + QString::number(y1, 'g', 10));
            item->child(2)->setText(1, QString::number(x2, 'g', 10) + ";" + QString::number(y2, 'g', 10));
        }
    }
    else if (item->text(0) == "vertex")
    {
        QStringList coords = item->text(1).split(";");
        double x = coords.at(0).toDouble();
        double y = coords.at(1).toDouble();
        PointEditorDialog dialog(this, x, y);
        dialog.exec();
        if (dialog.result() == QDialog::Accepted)
        {
            x = dialog.xValue();
            y = dialog.yValue();
            item->setText(1, QString::number(x, 'g', 10) + ";" + QString::number(y, 'g', 10));
        }
    }
    Q_UNUSED(column)
}

void PolygonalModelDialog::on_addQuadButton_clicked()
{
    StructuredIsoMesh2DDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        QTreeWidgetItem *item = new QTreeWidgetItem(ui->treeWidget->invisibleRootItem());
        item->setText(0, "quad");
        item->setText(1, QString("%1;%2").arg(dialog.xiCount()).arg(dialog.etaCount()));
        for (int i = 0; i < 4; i++)
        {
            QTreeWidgetItem *child = new QTreeWidgetItem(item);
            child->setText(0, "vertex");
            child->setText(1, QString::number(dialog.x(i), 'g', 10) + ";" + QString::number(dialog.y(i), 'g', 10));
        }
    }
}

void PolygonalModelDialog::on_addTriangleButton_clicked()
{
    BaryQuadsDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        QTreeWidgetItem *item = new QTreeWidgetItem(ui->treeWidget->invisibleRootItem());
        item->setText(0, "triangle");
        item->setText(1, QString("%1").arg(dialog.nodesCount()));
        for (int i = 0; i < 3; i++)
        {
            QTreeWidgetItem *child = new QTreeWidgetItem(item);
            child->setText(0, "vertex");
            child->setText(1, QString::number(dialog.x(i), 'g', 10) + ";" + QString::number(dialog.y(i), 'g', 10));
        }
    }
}

void PolygonalModelDialog::on_buildModelButton_clicked()
{
    QZScriptEngine interpreter;
    QScriptValue qsModel = interpreter.newQObject(this);

    interpreter.globalObject().setProperty("Model", qsModel);

    ui->treeWidget->clear();
    interpreter.evaluate(ui->textEdit->toPlainText());

    if (interpreter.hasUncaughtException())
    {
        QMessageBox::warning(this, tr("Ошибка разбора параметрической модели"), interpreter.uncaughtException().toString());
    }
}

void PolygonalModelDialog::on_newParametricModelButton_clicked()
{
    clearAll();
}

void PolygonalModelDialog::on_openParametricModelButton_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Открыть параметрическую модель"), "", tr("js-файлы (*.js);;Любой файл (*)"));
    if(!fileName.isEmpty())
    {
        QFile file(fileName);
        if (!file.open(QFile::ReadOnly | QFile::Text))
        {
            std::cerr << "Ошибка: Невозможно прочитать файл " << qPrintable(fileName) << ": " << qPrintable(file.errorString()) << std::endl;
            return;
        }
        clearAll();

        QTextStream textStream(&file);
        ui->textEdit->setText(textStream.readAll());
        scriptFileName = fileName;
        ui->textEdit->document()->setModified(false);
    }

}

void PolygonalModelDialog::on_saveParametricModelButton_clicked()
{
    saveParametricModel();
}
