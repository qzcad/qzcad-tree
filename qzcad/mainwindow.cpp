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
#include "rotationbodymeshdialog.h"
#include "elasticfemdialog.h"

#include "quadrilateralmesh2d.h"
#include "quadrilateralunion2d.h"
#include "hexahedralmesh3d.h"
#include "exportmeshdialog.h"
#include "qstdredirector.h"

#include "hexahedralfem.h"
#include "qtscriptfemcondition3d.h"
#include "qtscriptforcecondition3d.h"

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
    stdRedirector = new QStdRedirector<>(std::cout, this);
    connect(stdRedirector, SIGNAL(messageChanged(QString)), this, SLOT(onConsoleMessage(QString)));
    std::cout << "Система успешно запущена и готова к использованию..." << std::endl;
}

MainWindow::~MainWindow()
{
    msh::MeshPointer meshPtr = ui->pictureControl->releaseMesh();
    if (meshPtr != NULL) delete meshPtr; // чистим старые данные
    delete stdRedirector;
    delete ui;
}

void MainWindow::onConsoleMessage(QString message)
{
    ui->terminalText->insertPlainText(message);
    QTextCursor c =  ui->terminalText->textCursor();
    c.movePosition(QTextCursor::End);
    ui->terminalText->setTextCursor(c);
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

void MainWindow::clearMesh(MeshPointer mesh)
{
    if (mesh)
    {
        delete mesh;
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
        msh::MeshPointer meshPtr = ui->pictureControl->releaseMesh();
        clearMesh(meshPtr);
        msh::QuadrilateralMesh2D * mesh = new msh::QuadrilateralMesh2D(dialog.xCount(), dialog.yCount(),
                                                                       dialog.xMin(), dialog.yMin(),
                                                                       dialog.rectWidth(), dialog.rectHeight());
        meshPtr = mesh;
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionStructIsoQuads_triggered()
{
    StructuredIsoMesh2DDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::MeshPointer meshPtr = ui->pictureControl->releaseMesh();
        clearMesh(meshPtr);
        msh::QuadrilateralMesh2D * mesh = new msh::QuadrilateralMesh2D(dialog.xiCount(),
                                                                       dialog.etaCount(),
                                                                       msh::Point2D(dialog.x0(), dialog.y0()),
                                                                       msh::Point2D(dialog.x1(), dialog.y1()),
                                                                       msh::Point2D(dialog.x2(), dialog.y2()),
                                                                       msh::Point2D(dialog.x3(), dialog.y3()));
        meshPtr = mesh;
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionBaryQuads_triggered()
{
    BaryQuadsDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::MeshPointer meshPtr = ui->pictureControl->releaseMesh();
        clearMesh(meshPtr);
        msh::QuadrilateralMesh2D * mesh = new msh::QuadrilateralMesh2D(dialog.nodesCount(),
                                                                       msh::Point2D(dialog.x0(), dialog.y0()),
                                                                       msh::Point2D(dialog.x1(), dialog.y1()),
                                                                       msh::Point2D(dialog.x2(), dialog.y2()));
        meshPtr = mesh;
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionPolygonalModel_triggered()
{
    PolygonalModelDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::MeshPointer meshPtr = ui->pictureControl->releaseMesh();
        clearMesh(meshPtr);
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
        meshPtr = mesh;
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionStructHex_triggered()
{
    RectMeshDialog dialog(this, 3);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::MeshPointer meshPtr = ui->pictureControl->releaseMesh();
        clearMesh(meshPtr);
        msh::HexahedralMesh3D * mesh = new msh::HexahedralMesh3D(dialog.xCount(), dialog.yCount(), dialog.zCount(),
                                                                 dialog.xMin(), dialog.yMin(), dialog.zMin(),
                                                                 dialog.rectWidth(), dialog.rectHeight(), dialog.rectDepth());
        meshPtr = mesh;
        ui->pictureControl->setMesh(meshPtr);
    }
}

void MainWindow::on_actionSaveMesh_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->releaseMesh();
    if (mesh)
    {
        QString fileName = QFileDialog::getSaveFileName(this, "qMesher: Сохранить дискретную модель", "", tr("Текстовые файлы (*.txt);;Любой файл (*)"));
        if (fileName.isEmpty())
            return;
        QFile file(fileName);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate | QIODevice::Text))
            return;
        ExportMeshDialog dialog(this);
        dialog.exec();
        if (dialog.result() == QDialog::Accepted)
        {
            QTextStream out(&file);
            int dim = mesh->dimesion();
            int elementNodes = mesh->element(0)->verticesCount(); // определяем количество узлов в элементе
            out << dim << ' ' << elementNodes << '\n';
            out << mesh->nodesCount() << ' ' << dialog.isNodeValue() << '\n';
            for (msh::UInteger i = 0; i < mesh->nodesCount(); i++)
            {
                msh::PointPointer point = mesh->node(i);
                out << point->x() << " " << point->y() << " ";
                if (dim == 3)
                    out << point->z() << " ";
                out << mesh->nodeType(i);
                if (dialog.isNodeValue()) out << " " << mesh->nodeValue(i);
                out << '\n';
            }
            out << mesh->elementsCount() << ' ' << dialog.isElementValue() << '\n';
            for (msh::UInteger i = 0; i < mesh->elementsCount(); i++)
            {
                msh::ElementPointer element = mesh->element(i);
                for (int j = 0; j < element->verticesCount(); j++)
                {
                    out << element->vertexNode(j) << " ";
                }
                if (dialog.isElementValue()) out << mesh->elementValue(i);
                out << '\n';
            }
        }

        ui->pictureControl->setMesh(mesh);
    }
}

void MainWindow::on_actionRotationBodyMesh_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getMesh();
    if (mesh)
    {
        msh::QuadrilateralMesh2D *qmesh = dynamic_cast<msh::QuadrilateralMesh2D *>(mesh);
        if (qmesh)
        {
            // сетка четырехугольников
            RotationBodyMeshDialog dialog(this);
            dialog.exec();
            if (dialog.result() == QDialog::Accepted)
            {
                ui->pictureControl->resetMesh();
                msh::HexahedralMesh3D *hmesh;
                int axe = dialog.axe();
                double radius = dialog.radius();
                int layersCount = dialog.layersCount();
                bool isClosedBody = dialog.isClosedBody();
                double angle = dialog.rotationAngle();
                if (isClosedBody)
                    hmesh = new msh::HexahedralMesh3D(qmesh, (axe == 0) ? 0.0 : radius, (axe == 0) ? radius : 0.0, layersCount, (axe == 0) ? true : false );
                else
                    hmesh = new msh::HexahedralMesh3D(qmesh, (axe == 0) ? 0.0 : radius, (axe == 0) ? radius : 0.0, angle, layersCount, (axe == 0) ? true : false );
                msh::MeshPointer mesh3d(hmesh);
                ui->pictureControl->setMesh(mesh3d);
                // !!! удаляем старые данные
                clearMesh(qmesh);

            }
        }
        else
        {
            QString message = tr("Построение модели тела вращения для данного типа сетки не предусмотрено");
            QMessageBox::warning(this, message, message);
        }
    }
}

void MainWindow::on_actionFlipVertically_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getMesh();
    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            ui->pictureControl->resetMesh();
            mesh2d->flipVertically();
            ui->pictureControl->setMesh(mesh2d);
        }
        else
        {
            QString message = tr("Отражение по вертикали для данного типа сетки не предусмотрено");
            QMessageBox::warning(this, message, message);
        }
    }
}

void MainWindow::on_actionFlipHorizontally_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getMesh();
    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            ui->pictureControl->resetMesh();
            mesh2d->flipHorizontally();
            ui->pictureControl->setMesh(mesh2d);
        }
        else
        {
            QString message = tr("Отражение по горизонтали для данного типа сетки не предусмотрено");
            QMessageBox::warning(this, message, message);
        }
    }
}

void MainWindow::on_actionMirrorVertically_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getMesh();
    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            ui->pictureControl->resetMesh();
            mesh2d->mirrorVertically();
            ui->pictureControl->setMesh(mesh2d);
        }
        else
        {
            QString message = tr("Зеракльное отражение по вертикали для данного типа сетки не предусмотрено");
            QMessageBox::warning(this, message, message);
        }
    }
}

void MainWindow::on_actionMirrorrHorizontally_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getMesh();
    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            ui->pictureControl->resetMesh();
            mesh2d->mirrorHorizontally();
            ui->pictureControl->setMesh(mesh2d);
        }
        else
        {
            QString message = tr("Зеркальное отражение по горизонтали для данного типа сетки не предусмотрено");
            QMessageBox::warning(this, message, message);
        }
    }
}

void MainWindow::on_actionArea_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getMesh();
    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            QString message = tr("Площадь модели равна %1").arg(mesh2d->area());
            QMessageBox::information(this, message, message);
        }
        else
        {
            msh::Mesh3D *mesh3d = dynamic_cast<msh::Mesh3D*>(mesh);
            if (mesh3d)
            {
                QString message = tr("Площадь поверхности модели равна %1").arg(mesh3d->surfaceArea());
                QMessageBox::information(this, message, message);
            }
            else
            {
                QString message = tr("Вычисление площади для данного типа сетки не предусмотрено");
                QMessageBox::warning(this, message, message);
            }
        }
    }
}

void MainWindow::on_actionLoadMesh_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->releaseMesh();

    clearMesh(mesh);

    QString fileName = QFileDialog::getOpenFileName(this, "qMesher: Загрузить дискретную модель", "", tr("Текстовые файлы (*.txt);;Любой файл (*)"));
    if (fileName.isEmpty())
        return;
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    std::cout << "Загррузка дискретной модели из файла: " << fileName.toAscii().constData() << std::endl;
    QTextStream in(&file);
    int dim;
    int elementNodes;
    msh::UInteger nodesCount;
    int isNodeValue;
    msh::UInteger elementsCount;
    int isElementValue;
    in >> dim;
    in >> elementNodes;
    in >> nodesCount;
    in >> isNodeValue;
    std::cout << "Размерность: " << dim << std::endl << "Количество узлов: " << nodesCount << std::endl;
    if (dim == 2 && elementNodes == 4) // четырехугольники
    {
        msh::QuadrilateralMesh2D *qMesh = new msh::QuadrilateralMesh2D();
        for (msh::UInteger i = 0; i < nodesCount; i++)
        {
            msh::Floating x, y, val;
            msh::Point2D point;
            int nodeType;
            in >> x;
            in >> y;
            in >> nodeType;
            point.set(x, y);
            qMesh->pushNode(point, static_cast<msh::NodeType>(nodeType));
            if (isNodeValue)
            {
                in >> val;
                qMesh->pushNodeValue(val);
            }
        }
        in >> elementsCount;
        in >> isElementValue;
        std::cout << "Количнство элементов: " << elementsCount << std::endl;
        for (msh::UInteger i = 0; i < elementsCount; i++)
        {
            msh::UInteger p[elementNodes];
            msh::Floating val;
            for (int j = 0; j < elementNodes; j++)
                in >> p[j];
            qMesh->addElement(p[0], p[1], p[2], p[3]);
            if (isElementValue)
            {
                in >> val;
                qMesh->pushElementValue(val);
            }
        }
        qMesh->updateDomain();
        ui->pictureControl->setMesh(qMesh);
    }
    if (dim == 3 && elementNodes == 8) // шестигранники
    {
        msh::HexahedralMesh3D *hMesh = new msh::HexahedralMesh3D();
        for (msh::UInteger i = 0; i < nodesCount; i++)
        {
            msh::Floating x, y, z, val;
            msh::Point3D point;
            int nodeType;
            in >> x;
            in >> y;
            in >> z;
            in >> nodeType;
            point.set(x, y, z);
            hMesh->pushNode(point, static_cast<msh::NodeType>(nodeType));
            if (isNodeValue == 1)
            {
                in >> val;
                hMesh->pushNodeValue(val);
            }
        }
        in >> elementsCount;
        in >> isElementValue;
        std::cout << "Количнство элементов: " << elementsCount << std::endl;
        for (msh::UInteger i = 0; i < elementsCount; i++)
        {
            msh::UInteger p[elementNodes];
            msh::Floating val;
            for (int j = 0; j < elementNodes; j++)
                in >> p[j];
            hMesh->addElement(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
            if (isElementValue == 1)
            {
                in >> val;
                hMesh->pushElementValue(val);
            }
        }
        hMesh->updateDomain();
        ui->pictureControl->setMesh(hMesh);
    }
}

void MainWindow::on_actionElasticFem_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->releaseMesh();
    ElasticFemDialog dialog(this);

    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            dialog.set2dMode();
            dialog.exec();
            // to do
        }
        else
        {
            msh::HexahedralMesh3D *hexahedralMesh = dynamic_cast<msh::HexahedralMesh3D*>(mesh);
            if (hexahedralMesh)
            {
                dialog.exec();
                if (dialog.result() == QDialog::Accepted)
                {
                    MechanicalParameters3D *params;
                    if (dialog.isIsotropy())
                        params = new MechanicalParameters3D(dialog.e(), dialog.nu());
                    else if (dialog.isOrthotropyENuG())
                        params = new MechanicalParameters3D(dialog.e(), dialog.nu(), dialog.g());
                    else
                        params = new MechanicalParameters3D(dialog.e1(), dialog.e2(), dialog.e3(),
                                                            dialog.nu12(), dialog.nu13(), dialog.nu23(),
                                                            dialog.g12(), dialog.g13(), dialog.g23());
                    std::vector<FEMCondition3DPointer> boundaryConditions;
                    for (int i = 0; i < dialog.boundaryCount(); i++)
                    {
                        boundaryConditions.push_back(new QtScriptFemCondition3D(dialog.boundaryCondition(i),
                                                                                dialog.boundaryIsU(i),
                                                                                dialog.boundaryU(i).toDouble(),
                                                                                dialog.boundaryIsV(i),
                                                                                dialog.boundaryV(i).toDouble(),
                                                                                dialog.boundaryIsW(i),
                                                                                dialog.boundaryW(i).toDouble()));
                    }
                    std::vector<FEMCondition3DPointer> forces;
                    for (int i = 0; i < dialog.forcesCount(); i++)
                    {
                        forces.push_back(new QtScriptForceCondition3D(dialog.forceCondition(i),
                                                                                  dialog.forceU(i),
                                                                                  dialog.forceV(i),
                                                                                  dialog.forceW(i)));
                    }
                    HexahedralFEM fem(hexahedralMesh, *params, forces, boundaryConditions);

                    // очистка памяти
                    delete params;
                    for (int i = 0; i < dialog.boundaryCount(); i++)
                        delete boundaryConditions[i];
                    for (int i = 0; i < dialog.forcesCount(); i++)
                        delete forces[i];

                    ui->pictureControl->setMesh(hexahedralMesh);
                    ui->pictureControl->pushNodeValuesVector(NamedFloatingVector("u", fem.u()));
                    ui->pictureControl->pushNodeValuesVector(NamedFloatingVector("v", fem.v()));
                    ui->pictureControl->pushNodeValuesVector(NamedFloatingVector("w", fem.w()));
                    ui->pictureControl->pushElementValuesVector(NamedFloatingVector("Sigma u", fem.sigmaX()));
                    ui->pictureControl->pushElementValuesVector(NamedFloatingVector("Sigma v", fem.sigmaY()));
                    ui->pictureControl->pushElementValuesVector(NamedFloatingVector("Sigma w", fem.sigmaZ()));
                    ui->pictureControl->pushElementValuesVector(NamedFloatingVector("Tau uv", fem.tauXY()));
                    ui->pictureControl->pushElementValuesVector(NamedFloatingVector("Tau vw", fem.tauYZ()));
                    ui->pictureControl->pushElementValuesVector(NamedFloatingVector("Tau wu", fem.tauZX()));
                }
            }
        }
    }

}
