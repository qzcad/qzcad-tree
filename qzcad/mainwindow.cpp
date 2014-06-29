#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QMessageBox>
#include <QTextCodec>
#include <QFileDialog>
#include <QTextStream>

#include "globalconsts.h"

#include "rectmeshdialog.h"
#include "structuredisomesh2ddialog.h"
#include "path2deditor.h"
#include "polygonalmodeldialog.h"
#include "baryquadsdialog.h"

#include "quadrilateralmesh2d.h"
#include "quadrilateralunion2d.h"
#include "hexahedralmesh3d.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    // устанавливаем UTF-8 по умолчанию
    QTextCodec::setCodecForLocale(QTextCodec::codecForName("UTF-8"));
#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
    QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8"));
    QTextCodec::setCodecForTr(QTextCodec::codecForName("UTF-8"));
#endif
    showMaximized();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionUnion_triggered()
{
    addElement(_UNION_, "Операция", "Объединение", ":/icons/qzicons/union.png");
}

void MainWindow::addObjectTreeItem(QTreeWidgetItem * current,
                                   const int &id,
                                   const QString &category,
                                   const QString &type,
                                   const QString &iconName,
                                   QString parameters)
{
    current->setText(0, type); // тип
    current->setText(1, category); // категория
    current->setText(2, parameters); // сводка параметров
    current->setText(3, QString::number(id)); // идентификатор
    current->setIcon(0, QPixmap(iconName));
    ui->objectTree->setCurrentItem(current);
}

void MainWindow::addElement(const int &id,
                            const QString &category,
                            const QString &type,
                            const QString &iconName,
                            QString parameters)
{
    QTreeWidgetItem *current = ui->objectTree->currentItem();
    if (!current)
    {
        // добавляется первый элемент
        addObjectTreeItem(new QTreeWidgetItem(ui->objectTree), id, category, type, iconName, parameters);
    }
    else
    {
        // вложеные элементы могут быть только в операциях
        int element_id = current->text(3).toInt(); // идентификатор текущего элемента
        if (_OPERATION_MIN_ID_ <= element_id && element_id <= _OPERATION_MAX_ID_)
        {
            addObjectTreeItem(new QTreeWidgetItem(current), id, category, type, iconName, parameters);
        }
    }
}

void MainWindow::on_actionIntersection_triggered()
{
    addElement(_INTERSECTION_, "Операция", "Пересечение", ":/icons/qzicons/intersection.png");
}

void MainWindow::on_actionDifference_triggered()
{
    addElement(_DIFFERENCE_, "Операция", "Разность", ":/icons/qzicons/difference.png");
}

void MainWindow::on_actionPath_triggered()
{
    Path2DEditor dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        addElement(_PATH_, "Примитив", "Свободный контур", ":/icons/qzicons/path.png");
    }
}

void MainWindow::on_actionModel_triggered()
{
    addElement(_MODEL_, "Аналитическая модель", "модель", ":/icons/qzicons/model.png");
}

void MainWindow::on_actionDeleteObject_triggered()
{
    QTreeWidgetItem *current = ui->objectTree->currentItem();
    if (current) delete current;
}

void MainWindow::on_actionStructQuads_triggered()
{
    RectMeshDialog dialog(this, 2);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::QuadrilateralMesh2D * mesh = new msh::QuadrilateralMesh2D(dialog.xCount(), dialog.yCount(),
                                                                       dialog.xMin(), dialog.yMin(),
                                                                       dialog.rectWidth(), dialog.rectHeight());
        msh::MeshPointer meshPtr(mesh);
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionStructIsoQuads_triggered()
{
    StructuredIsoMesh2DDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::QuadrilateralMesh2D * mesh = new msh::QuadrilateralMesh2D(dialog.xiCount(),
                                                                       dialog.etaCount(),
                                                                       msh::Point2D(dialog.x0(), dialog.y0()),
                                                                       msh::Point2D(dialog.x1(), dialog.y1()),
                                                                       msh::Point2D(dialog.x2(), dialog.y2()),
                                                                       msh::Point2D(dialog.x3(), dialog.y3()));
        //        QuadrilateralMesh2D * mesh = new QuadrilateralMesh2D(dialog.xiCount(),
        //                                                             Point2D(dialog.x0(), dialog.y0()),
        //                                                             Point2D(dialog.x1(), dialog.y1()),
        //                                                             Point2D(dialog.x2(), dialog.y2()));
        msh::MeshPointer meshPtr(mesh);
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionBaryQuads_triggered()
{
    BaryQuadsDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::QuadrilateralMesh2D * mesh = new msh::QuadrilateralMesh2D(dialog.nodesCount(),
                                                                       msh::Point2D(dialog.x0(), dialog.y0()),
                                                                       msh::Point2D(dialog.x1(), dialog.y1()),
                                                                       msh::Point2D(dialog.x2(), dialog.y2()));
        msh::MeshPointer meshPtr(mesh);
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionPolygonalModel_triggered()
{
    PolygonalModelDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::QuadrilateralUnion2D * mesh = new msh::QuadrilateralUnion2D();

        int quadsCount = dialog.quadsCount();
        for (int i = 0; i < quadsCount; i++)
        {
            int xiCount = 0;
            int etaCount = 0;
            std::vector<msh::Point2D> quad = dialog.quad(i, xiCount, etaCount);
            msh::QuadrilateralMesh2D qmesh (xiCount, etaCount, quad[0], quad[1], quad[2], quad[3]);
            mesh->addMesh(&qmesh);
        }

        int trianglesCount = dialog.trianglesCount();
        for (int i = 0; i < trianglesCount; i++)
        {
            int nodesCount = 0;
            std::vector<msh::Point2D> tri = dialog.triangle(i, nodesCount);
            msh::QuadrilateralMesh2D qmesh (nodesCount, tri[0], tri[1], tri[2]);
            mesh->addMesh(&qmesh);
        }
//        msh::HexahedralMesh3D *hmesh = new msh::HexahedralMesh3D(mesh, 0.0, -1990.0, 45.0, 350, true);
        msh::MeshPointer meshPtr(mesh);
//        msh::MeshPointer meshPtr(hmesh);
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionStructHex_triggered()
{
    RectMeshDialog dialog(this, 3);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::HexahedralMesh3D * mesh = new msh::HexahedralMesh3D(dialog.xCount(), dialog.yCount(), dialog.zCount(),
                                                                 dialog.xMin(), dialog.yMin(), dialog.zMin(),
                                                                 dialog.rectWidth(), dialog.rectHeight(), dialog.rectDepth());
        msh::MeshPointer meshPtr(mesh);
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionSaveMesh_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getMesh();
    if (mesh)
    {
        QString fileName = QFileDialog::getSaveFileName(this, "qMesher: Сохранить дискретную модель", "", tr("Текстовые файлы (*.txt);;Любой файл (*)"));
        QFile file(fileName);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text))
            return;

        QTextStream out(&file);
        int dim = mesh->dimesion();
        out << dim << '\n';
        out << mesh->nodesCount() << '\n';
        for (msh::UInteger i = 0; i < mesh->nodesCount(); i++)
        {
            msh::PointPointer point = mesh->node(i);
            out << point->x() << " " << point->y() << " ";
            if (dim == 3)
                out << point->z() << " ";
            out << mesh->nodeType(i) << '\n';
        }
        out << mesh->elementsCount() << '\n';
        for (msh::UInteger i = 0; i < mesh->elementsCount(); i++)
        {
            msh::ElementPointer element = mesh->element(i);
            for (int j = 0; j < element->verticesCount(); j++)
            {
                out << element->vertexNode(j) << " ";
            }
            out << '\n';
        }
    }
}
