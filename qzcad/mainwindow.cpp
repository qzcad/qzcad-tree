#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <float.h>
#include <iomanip>

#include <QMessageBox>
#include <QTextCodec>
#include <QFileDialog>
#include <QTextStream>
#include <QInputDialog>
#include <QTime>
#include <QFileInfo>
#include <QTime>
#include <QColorDialog>

#include "globalconsts.h"

#include "rectmeshdialog.h"
#include "structuredisomesh2ddialog.h"
#include "path2deditor.h"
#include "polygonalmodeldialog.h"
#include "baryquadsdialog.h"
#include "rotationbodymeshdialog.h"
#include "elasticfemdialog.h"

#include "segmentmesh2d.h"
#include "quadrilateralmesh2d.h"
#include "quadrilateralunion2d.h"
#include "quadrilateralmesh3d.h"
#include "trianglemesh2d.h"
#include "trianglemesh3d.h"
#include "hexahedralmesh3d.h"

#include "exportmeshdialog.h"
#include "qstdredirector.h"

#include "hexahedralfem.h"
#include "qtscriptfemcondition3d.h"
#include "qtscriptforcecondition3d.h"
#include "qzscriptengine.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // по умолчанию дерево с операциями скрыто
    ui->dockWidgetStruct->hide();

    // устанавливаем UTF-8 по умолчанию
    QTextCodec::setCodecForLocale(QTextCodec::codecForName("UTF-8"));
#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
    QTextCodec::setCodecForCStrings(QTextCodec::codecForName("UTF-8"));
    QTextCodec::setCodecForTr(QTextCodec::codecForName("UTF-8"));
#endif

#ifdef Q_OS_WIN
    ui->actionDoubleBufferGL->setChecked(false); // disable double buffer in Win32
    ui->terminalText->setFont(QFont("Courier"));
#else
    ui->terminalText->setFont(QFont("Monospace"));
#endif

    showMaximized();

    highlighter = new Highlighter(ui->codeEditor->document());

    stdRedirector = new QStdRedirector<>(std::cout, this);
    connect(stdRedirector, SIGNAL(messageChanged(QString)), this, SLOT(onConsoleMessage(QString)));

    std::cout << std::scientific << std::setprecision(9);

    std::cout << QTime::currentTime().toString("HH:mm:ss").toStdString() << ": система успешно запущена и готова к использованию..." << std::endl;
}

MainWindow::~MainWindow()
{
    msh::MeshPointer meshPtr = ui->pictureControl->getGlMeshPicture()->releaseMesh();
    if (meshPtr != NULL) delete meshPtr; // чистим старые данные
    delete highlighter;
    delete stdRedirector;
    delete ui;
}

void MainWindow::closeEvent(QCloseEvent *event)
{
    if (maybeSaveScript())
        event->accept();
    else
        event->ignore();
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

void MainWindow::setCurrentScriptName(const QString &fileName)
{
    scriptFileName = fileName;
    ui->codeEditor->document()->setModified(false);
    setWindowModified(false);
    setWindowFilePath(scriptFileName);
}

bool MainWindow::writeScript(const QString &fileName)
{
    QFile file(fileName);
    if(!file.open(QFile::WriteOnly | QFile::Text))
    {
        QString fileWriteErrorBoxTitle = tr("Ошибка");
        QString fileWriteErrorBoxMessage = tr("Невозможно записать файл %1:\n%2").arg(fileName).arg(file.errorString());
        QMessageBox::warning(this, fileWriteErrorBoxTitle, fileWriteErrorBoxMessage);
        return false;
    }
    QTextStream out(&file);

    out << ui->codeEditor->toPlainText();
    setCurrentScriptName(fileName);

    return true;
}

bool MainWindow::readScript(const QString &fileName)
{
    QFile file(fileName);
    if(!file.open(QFile::ReadOnly | QFile::Text))
    {
        QString readFileErrorBoxTitle = tr("Ошибка");
        QString readFileErrorBoxMessage = tr("Невозможно прочитать файл %1:\n%2").arg(fileName).arg(file.errorString());
        QMessageBox::warning(this, readFileErrorBoxTitle, readFileErrorBoxMessage);
        return false;
    }
    QTextStream in(&file);
    ui->codeEditor->setPlainText(in.readAll());
    setCurrentScriptName(fileName);

    return true;
}

bool MainWindow::saveScriptAs()
{
    QString saveAsDialogTitle = tr("Сохранить как");
    QString dir = (scriptFileName.isEmpty()) ? "" : QFileInfo(scriptFileName).absolutePath();
    QString saveAsDialogVariants = tr("Qt Script файлы (*.js);;Текстовые файлы (*.txt);;Любой файл (*)");
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    saveAsDialogTitle,
                                                    dir,
                                                    saveAsDialogVariants);

    if(fileName.isEmpty())
        return false;

    return writeScript(fileName);
}

bool MainWindow::saveScript()
{
    if (scriptFileName.isEmpty())
        return saveScriptAs();

    return writeScript(scriptFileName);
}

bool MainWindow::maybeSaveScript()
{
    if(ui->codeEditor->document()->isModified())
    {
        QMessageBox::StandardButton ret;
        QString maybeSaveDialogTitile = tr("Сохранить изменения?");
        QString maybeSaveDialogMessage =  tr("Код модели модифицирован. Сохранить изменения?");
        ret = QMessageBox::question(this,
                                    maybeSaveDialogTitile,
                                    maybeSaveDialogMessage,
                                    QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
        if (ret == QMessageBox::Save)
            return saveScript();
        else if (ret == QMessageBox::Cancel)
            return false;
    }
    return true;
}

bool MainWindow::openScript()
{
    if(maybeSaveScript())
    {
        QString openFileDialogTitle = tr("Открыть файл");
        QString dir = (scriptFileName.isEmpty()) ? "" : QFileInfo(scriptFileName).absolutePath();
        QString openFileDialogVariants = tr("Qt Script файлы (*.js);;Текстовые файлы (*.txt);;Любой файл (*)");
        QString fileName = QFileDialog::getOpenFileName(this,
                                                        openFileDialogTitle,
                                                        dir,
                                                        openFileDialogVariants);
        if(!fileName.isEmpty())
            return readScript(fileName);
    }
    return false;
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
        msh::MeshPointer meshPtr = ui->pictureControl->getGlMeshPicture()->releaseMesh();
        clearMesh(meshPtr);
        msh::QuadrilateralMesh2D * mesh = new msh::QuadrilateralMesh2D();
        mesh->rectangleDomain(dialog.xCount(), dialog.yCount(), dialog.xMin(), dialog.yMin(), dialog.rectWidth(), dialog.rectHeight());
        meshPtr = mesh;
        ui->pictureControl->getGlMeshPicture()->setMesh(meshPtr);
    }
}

void MainWindow::on_actionStructIsoQuads_triggered()
{
    StructuredIsoMesh2DDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::MeshPointer meshPtr = ui->pictureControl->getGlMeshPicture()->releaseMesh();
        clearMesh(meshPtr);
        msh::QuadrilateralMesh2D * mesh = new msh::QuadrilateralMesh2D();
        mesh->quadDomain(dialog.xiCount(),
                         dialog.etaCount(),
                         msh::Point2D(dialog.x0(), dialog.y0()),
                         msh::Point2D(dialog.x1(), dialog.y1()),
                         msh::Point2D(dialog.x2(), dialog.y2()),
                         msh::Point2D(dialog.x3(), dialog.y3()));
        meshPtr = mesh;
        ui->pictureControl->getGlMeshPicture()->setMesh(meshPtr);
    }
}

void MainWindow::on_actionBaryQuads_triggered()
{
    BaryQuadsDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::MeshPointer meshPtr = ui->pictureControl->getGlMeshPicture()->releaseMesh();
        clearMesh(meshPtr);
        msh::QuadrilateralMesh2D * mesh = new msh::QuadrilateralMesh2D();
        mesh->triangleDomain(dialog.nodesCount(),
                             msh::Point2D(dialog.x0(), dialog.y0()),
                             msh::Point2D(dialog.x1(), dialog.y1()),
                             msh::Point2D(dialog.x2(), dialog.y2()));
        meshPtr = mesh;
        ui->pictureControl->getGlMeshPicture()->setMesh(meshPtr);
    }
}

void MainWindow::on_actionPolygonalModel_triggered()
{
    PolygonalModelDialog dialog(this);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::MeshPointer meshPtr = ui->pictureControl->getGlMeshPicture()->releaseMesh();
        clearMesh(meshPtr);
        msh::QuadrilateralUnion2D * mesh = new msh::QuadrilateralUnion2D();
        int quadsCount = dialog.quadsCount();
        for (int i = 0; i < quadsCount; i++)
        {
            int xiCount = 0;
            int etaCount = 0;
            std::vector<msh::Point2D> quad = dialog.quad(i, xiCount, etaCount);
            msh::QuadrilateralMesh2D qmesh;
            qmesh.quadDomain(xiCount, etaCount, quad[0], quad[1], quad[2], quad[3]);
            mesh->addMesh(&qmesh);
        }

        int trianglesCount = dialog.trianglesCount();
        for (int i = 0; i < trianglesCount; i++)
        {
            int nodesCount = 0;
            std::vector<msh::Point2D> tri = dialog.triangle(i, nodesCount);
            msh::QuadrilateralMesh2D qmesh;
            qmesh.triangleDomain(nodesCount, tri[0], tri[1], tri[2]);
            mesh->addMesh(&qmesh);
        }
        meshPtr = mesh;
        ui->pictureControl->getGlMeshPicture()->setMesh(meshPtr);
    }
}

void MainWindow::on_actionStructHex_triggered()
{
    RectMeshDialog dialog(this, 3);
    dialog.exec();
    if (dialog.result() == QDialog::Accepted)
    {
        msh::MeshPointer meshPtr = ui->pictureControl->getGlMeshPicture()->releaseMesh();
        clearMesh(meshPtr);
        msh::HexahedralMesh3D * mesh = new msh::HexahedralMesh3D();
        mesh->prismDomain(dialog.xCount(), dialog.yCount(), dialog.zCount(),
                          dialog.xMin(), dialog.yMin(), dialog.zMin(),
                          dialog.rectWidth(), dialog.rectHeight(), dialog.rectDepth());
        meshPtr = mesh;
        ui->pictureControl->getGlMeshPicture()->setMesh(meshPtr);
    }
}

void MainWindow::on_actionSaveMesh_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->releaseMesh();
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
            out.setRealNumberPrecision(DBL_DIG);
            int dim = mesh->dimesion();
            int elementNodes = mesh->element(0)->verticesCount(); // определяем количество узлов в элементе
            out << dim << ' ' << elementNodes;
            out << ' ' << mesh->element(0)->facesCount(); // количество граней для устранения неоднозначности (3d: у четырехугольника и тетраэдра одинаковое количество узлов)
            out << '\n';
            out << mesh->nodesCount() << '\n';
            for (msh::UInteger i = 0; i < mesh->nodesCount(); i++)
            {
                msh::PointPointer point = mesh->node(i);
                out << point->x() << " " << point->y() << " ";
                if (dim == 3)
                    out << point->z() << " ";
                out << mesh->nodeType(i);
                out << '\n';
            }
            out << mesh->elementsCount() << ' ' << dialog.isLayers() << '\n';
            for (msh::UInteger i = 0; i < mesh->elementsCount(); i++)
            {
                msh::ElementPointer element = mesh->element(i);
                for (int j = 0; j < element->verticesCount(); j++)
                {
                    out << element->vertexNode(j) << " ";
                }
                if (dialog.isLayers()) out << mesh->layer(i);
                out << '\n';
            }
        }

        ui->pictureControl->getGlMeshPicture()->setMesh(mesh);
    }
}

void MainWindow::on_actionRotationBodyMesh_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
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
                ui->pictureControl->getGlMeshPicture()->resetMesh();
                msh::HexahedralMesh3D *hmesh;
                int axe = dialog.axe();
                double radius = dialog.radius();
                int layersCount = dialog.layersCount();
                bool isClosedBody = dialog.isClosedBody();
                double angle = dialog.rotationAngle();
                hmesh = new msh::HexahedralMesh3D();
                if (isClosedBody)
                    hmesh->rotateBaseMesh(qmesh, (axe == 0) ? 0.0 : radius, (axe == 0) ? radius : 0.0, layersCount, (axe == 0) ? true : false);
                else
                    hmesh->rotateBaseMesh(qmesh, (axe == 0) ? 0.0 : radius, (axe == 0) ? radius : 0.0, angle, layersCount, (axe == 0) ? true : false );
                msh::MeshPointer mesh3d(hmesh);
                ui->pictureControl->getGlMeshPicture()->setMesh(mesh3d);
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
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            ui->pictureControl->getGlMeshPicture()->resetMesh();
            mesh2d->flipVertically();
            ui->pictureControl->getGlMeshPicture()->setMesh(mesh2d);
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
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            ui->pictureControl->getGlMeshPicture()->resetMesh();
            mesh2d->flipHorizontally();
            ui->pictureControl->getGlMeshPicture()->setMesh(mesh2d);
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
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            ui->pictureControl->getGlMeshPicture()->resetMesh();
            mesh2d->mirrorVertically();
            ui->pictureControl->getGlMeshPicture()->setMesh(mesh2d);
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
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh)
    {
        msh::Mesh2D *mesh2d = dynamic_cast<msh::Mesh2D*>(mesh);
        if (mesh2d)
        {
            ui->pictureControl->getGlMeshPicture()->resetMesh();
            mesh2d->mirrorHorizontally();
            ui->pictureControl->getGlMeshPicture()->setMesh(mesh2d);
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
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
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
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->releaseMesh();

    clearMesh(mesh);

    QString fileName = QFileDialog::getOpenFileName(this, "qMesher: Загрузить дискретную модель", "", tr("Текстовые файлы (*.txt);;Любой файл (*)"));
    if (fileName.isEmpty())
        return;
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return;
    std::cout << std::endl;
    std::cout << "Загррузка дискретной модели из файла: " << fileName.toStdString() << std::endl;
    QTextStream in(&file);
    int dim;
    int elementNodes;
    msh::UInteger nodesCount;
    msh::UInteger elementsCount;
    msh::UInteger facesCount;
    int isLayers;
    in >> dim;
    in >> elementNodes;
    in >> facesCount;
    in >> nodesCount;
    std::cout << "Размерность: " << dim << std::endl << "Количество узлов: " << nodesCount << std::endl;
    if (dim == 2)
    {
        switch (elementNodes) {
        case 2:
            mesh = new msh::SegmentMesh2D();
            break;
        case 3:
            mesh = new msh::TriangleMesh2D();
            break;
        case 4:
            mesh = new msh::QuadrilateralMesh2D();
            break;
        default:
            std::cout << "Error:Mesh type is not recognized." << std::endl;
            return;
        }
    }
    else if (dim == 3)
    {
        switch (elementNodes) {
        case 3:
            mesh = new msh::TriangleMesh3D();
            std::cout << "3d t";
            break;
        case 4:
            if (facesCount == 1)
                mesh = new msh::QuadrilateralMesh3D();
//            else
//                mesh = new msh::TetrahedralMeshD(); // TODO
            break;
        case 8:
            mesh = new msh::HexahedralMesh3D();
            break;
        default:
            std::cout << "Error:Mesh type is not recognized." << std::endl;
            return;
        }
    }
    for (msh::UInteger i = 0; i < nodesCount; i++)
    {
        double x = 0.0, y = 0.0, z = 0.0;
        int nodeType;
        Point3D p;
        in >> x;
        in >> y;
        if (dim == 3)
            in >> z;
        in >> nodeType;
        p.set(x, y, z);
        mesh->pushNode(&p, static_cast<msh::NodeType>(nodeType));
    }
    in >> elementsCount;
    in >> isLayers;
    std::cout << "Количнство элементов: " << elementsCount << std::endl;
    for (msh::UInteger i = 0; i < elementsCount; i++)
    {
        std::vector<msh::UInteger> nodes_ref(elementNodes);
        int l = 0;
        for (int j = 0; j < elementNodes; j++)
            in >> nodes_ref[j];
        mesh->addElement(nodes_ref);
        if (isLayers)
        {
            in >> l;
            mesh->pushLayer(l);
        }
    }
    mesh->updateDomain();
    ui->pictureControl->getGlMeshPicture()->setMesh(mesh);
}

void MainWindow::on_actionElasticFem_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->releaseMesh();
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
                    std::vector<ForceCondition3DPointer> forces;
                    for (int i = 0; i < dialog.forcesCount(); i++)
                    {
                        forces.push_back(new QtScriptForceCondition3D(dialog.forceCondition(i),
                                                                                  dialog.forceU(i),
                                                                                  dialog.forceV(i),
                                                                                  dialog.forceW(i), static_cast<ForceType>(dialog.forceTypeIndex(i))));
                    }
                    QTime beginTime = QTime::currentTime();
                    std::cout << beginTime.toString("HH:mm:ss").toStdString() << ": запущен расчет методом конечных элементов" << std::endl;
                    HexahedralFEM fem(hexahedralMesh, *params, forces, boundaryConditions);
                    QTime endTime = QTime::currentTime();
                    std::cout << endTime.toString("HH:mm:ss").toStdString() << ": завершен расчет методом конечных элементов; продолжительность рачета: " << beginTime.secsTo(endTime) << " секунд." << std::endl;
                    // очистка памяти
                    delete params;
                    for (int i = 0; i < dialog.boundaryCount(); i++)
                        delete boundaryConditions[i];
                    for (int i = 0; i < dialog.forcesCount(); i++)
                        delete forces[i];

                    hexahedralMesh->addDataVector("u", fem.u());
                    hexahedralMesh->addDataVector("v", fem.v());
                    hexahedralMesh->addDataVector("w", fem.w());
                    hexahedralMesh->addDataVector("Sigma u", fem.sigmaX());
                    hexahedralMesh->addDataVector("Sigma v", fem.sigmaY());
                    hexahedralMesh->addDataVector("Sigma w", fem.sigmaZ());
                    hexahedralMesh->addDataVector("Tau uv", fem.tauXY());
                    hexahedralMesh->addDataVector("Tau vw", fem.tauYZ());
                    hexahedralMesh->addDataVector("Tau wu", fem.tauZX());
                    hexahedralMesh->addDataVector("sigma", fem.sigma());
                    ui->pictureControl->getGlMeshPicture()->setMesh(hexahedralMesh);
                }
            }
        }
    }

}

void MainWindow::on_actionLoadNodeValue_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh != NULL)
    {
        QString fileName = QFileDialog::getOpenFileName(this, "Загрузить значение в узле", "", tr("Текстовые файлы (*.txt);;Любой файл (*)"));
        if (fileName.isEmpty())
            return;
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return;
        bool ok;
        QString name = QInputDialog::getText(this, tr("QInputDialog::getText()"),
                                             tr("Название вектора значений:"), QLineEdit::Normal,
                                             QFileInfo(fileName).baseName(), &ok);
        if (ok && !name.isEmpty())
        {
            QTextStream in(&file);
            msh::UInteger count = 0;
            std::vector<double> v;
            in >> count;
            if (mesh->nodesCount() != count)
            {
                QString msg = tr("Ошибка: Количество элементов вектора не равно количству узлов сетки.");
                QMessageBox::warning(this, msg, msg);
                return;
            }
            for (msh::UInteger i = 0; i < count; i++)
            {
                double val = 0.0;
                in >> val;
                v.push_back(val);
            }
            mesh->addDataVector(name.toStdString(), v);
        }
    }
}

void MainWindow::on_actionLoadElementValue_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh != NULL)
    {
        QString fileName = QFileDialog::getOpenFileName(this, "Загрузить значение на элементе", "", tr("Текстовые файлы (*.txt);;Любой файл (*)"));
        if (fileName.isEmpty())
            return;
        QFile file(fileName);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return;
        bool ok;
        QString name = QInputDialog::getText(this, tr("QInputDialog::getText()"),
                                             tr("Название вектора значений:"), QLineEdit::Normal,
                                             QFileInfo(fileName).baseName(), &ok);
        if (ok && !name.isEmpty())
        {
            QTextStream in(&file);
            msh::UInteger count = 0;
            std::vector<double> v;
            in >> count;
            if (mesh->elementsCount() != count)
            {
                QString msg = tr("Ошибка: Количество элементов вектора не равно количству элементов сетки.");
                QMessageBox::warning(this, msg, msg);
                return;
            }
            for (msh::UInteger i = 0; i < count; i++)
            {
                double val = 0.0;
                in >> val;
                v.push_back(val);
            }
            mesh->addDataVector(name.toStdString(), v);
        }
    }
}

void MainWindow::on_actionExtremeValuesStatistica_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh != NULL)
    {
        ui->tabWidget->setCurrentIndex(2); // switch to terminal's tab
        msh::NamedDoubleVector ndv = ui->pictureControl->getGlMeshPicture()->dataVector();
        std::cout << std::endl;
        if (mesh->elementsCount() == ndv.size())
        {
            double val;
            msh::NamedDoubleVector::size_type index = 0;
            std::cout << ndv.name() << ": " << std::endl;
            val = ndv.min(index);
            std::cout << "Минимальное значение: " << val << " в элементе с номером: " << index << std::endl;
            val = ndv.max(index);
            std::cout << "Максимальное значение: " << val << " в элементе с номером: " << index << std::endl;
        }
        else if (mesh->nodesCount() == ndv.size())
        {
            double val;
            msh::NamedDoubleVector::size_type index = 0;
            msh::PointPointer node;
            std::cout << ndv.name() << ": " << std::endl;
            val = ndv.max(index);
            node = mesh->node(index);
            std::cout << "Минимальное значение: " << val << " в узле с номером: " << index << std::endl;
            std::cout << "( " << node->x() << "; " << node->y() << "; " << node->z() << " )" << std::endl;
            val = ndv.max(index);
            node = mesh->node(index);
            std::cout << "Максимальное значение: " << val << " в узле с номером: " << index << std::endl;
            std::cout << "( " << node->x() << "; " << node->y() << "; " << node->z() << " )" << std::endl;
        }
    }
}

void MainWindow::on_actionNewScript_triggered()
{
    if (maybeSaveScript())
    {
        ui->codeEditor->clear();
        setCurrentScriptName("");
    }
}

void MainWindow::on_actionOpenScript_triggered()
{
    openScript();
    ui->tabWidget->setCurrentIndex(1); // switch to script's tab
}

void MainWindow::on_actionSaveScript_triggered()
{
    saveScript();
}

void MainWindow::on_actionSaveAsScript_triggered()
{
    saveScriptAs();
}

void MainWindow::on_actionRunScript_triggered()
{
    QZScriptEngine engine(this);
    msh::Mesh *mesh = ui->pictureControl->getGlMeshPicture()->releaseMesh();
    QTime time;
    ui->tabWidget->setCurrentIndex(2); // switch to terminal's tab
    std::cout << std::endl;
    std::cout << "QZScriptEngine started..." << std::endl;
    time.start();
    engine.setMesh(mesh);
    QScriptValue result = engine.evaluate(ui->codeEditor->toPlainText());
    if (engine.hasUncaughtException()) {
        int line = engine.uncaughtExceptionLineNumber();
        std::cout << "Uncaught exception at line " << line << ": " << result.toString().toStdString() << std::endl;
    }
    std::cout << "QZScriptEngine finished in " << time.elapsed() << " ms." << std::endl;
    if (engine.mesh() != NULL)
    {
//        clearMesh(ui->pictureControl->getGlMeshPicture()->releaseMesh());
        ui->pictureControl->getGlMeshPicture()->setMesh(engine.mesh());

//        for (unsigned i = 0; i < engine.getNodeValuesSize(); i++)
//            ui->pictureControl->getGlMeshPicture()->pushNodeValuesVector(engine.getNodeValues(i));
//        for (unsigned i = 0; i < engine.getElementValuesSize(); i++)
//            ui->pictureControl->getGlMeshPicture()->pushElementValuesVector(engine.getElementValues(i));
        ui->tabWidget->setCurrentIndex(0); // switch to picture's tab
    }
}

void MainWindow::on_actionJacobianMetric_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh != NULL)
    {
        ui->tabWidget->setCurrentIndex(2); // switch to terminal's tab
        std::cout << std::endl;
        std::cout << "Вчисление значений якобиана... ";
        double min = 0.0, max = 0.0;
        if (dynamic_cast<msh::TriangleMesh2D*>(mesh))
        {
            msh::TriangleMesh2D *triangles = dynamic_cast<msh::TriangleMesh2D*>(mesh);

            if (triangles->elementsCount() > 0)
            {
                std::vector<double> j(triangles->elementsCount());
                for (msh::UInteger i = 0; i < triangles->elementsCount(); i++)
                {
                    j[i] = triangles->jacobian(i);
                    if (i == 0)
                    {
                        min = max = j[i];
                    }
                    else
                    {
                        if (min > j[i])
                            min = j[i];
                        if (max < j[i])
                            max = j[i];
                    }
                }
                mesh->addDataVector("Jacobian", j);
            }
        }
        std::cout << "Выполнено: " << min << " <= J <= " << max << std::endl;
    }
}

void MainWindow::on_actionLengthAspect_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh != NULL)
    {
        ui->tabWidget->setCurrentIndex(2); // switch to terminal's tab
        std::cout << std::endl;
        std::cout << "Вычисление соотношений длин сторон элементов... ";
        double min = 0.0, max = 0.0;
        if (dynamic_cast<msh::TriangleMesh2D*>(mesh))
        {
            msh::TriangleMesh2D *triangles = dynamic_cast<msh::TriangleMesh2D*>(mesh);
            if (triangles->elementsCount() > 0)
            {
                std::vector<double> j(triangles->elementsCount());
                for (msh::UInteger i = 0; i < triangles->elementsCount(); i++)
                {
                    j[i] = triangles->lengthAspect(i);
                    if (i == 0)
                    {
                        min = max = j[i];
                    }
                    else
                    {
                        if (min > j[i])
                            min = j[i];
                        if (max < j[i])
                            max = j[i];
                    }
                }
                mesh->addDataVector("length ratio", j);
            }
        }
        std::cout << "Выполнено: " << min << "<= aspect(l) <= " << max << std::endl;
    }
}

void MainWindow::on_actionMinAngleMetric_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh != NULL)
    {
        ui->tabWidget->setCurrentIndex(2); // switch to terminal's tab
        std::cout << std::endl;
        std::cout << "Вычисление значений минимальных углов... ";
        double min = 0.0, max = 0.0;
        if (dynamic_cast<msh::TriangleMesh2D*>(mesh))
        {
            msh::TriangleMesh2D *triangles = dynamic_cast<msh::TriangleMesh2D*>(mesh);
            if (triangles->elementsCount() > 0)
            {
                std::vector<double> j(triangles->elementsCount());
                for (msh::UInteger i = 0; i < triangles->elementsCount(); i++)
                {
                    j[i] = triangles->minAngle(i) * 180.0 / M_PI;
                    if (i == 0)
                    {
                        min = max = j[i];
                    }
                    else
                    {
                        if (min > j[i])
                            min = j[i];
                        if (max < j[i])
                            max = j[i];
                    }
                }
                mesh->addDataVector("min angle, gradus", j);
            }
        }
        if (dynamic_cast<msh::TriangleMesh3D*>(mesh))
        {
            msh::TriangleMesh3D *triangles = dynamic_cast<msh::TriangleMesh3D*>(mesh);
            if (triangles->elementsCount() > 0)
            {
                std::vector<double> j(triangles->elementsCount());
                for (msh::UInteger i = 0; i < triangles->elementsCount(); i++)
                {
                    j[i] = triangles->minAngle(i) * 180.0 / M_PI;
                    if (i == 0)
                    {
                        min = max = j[i];
                    }
                    else
                    {
                        if (min > j[i])
                            min = j[i];
                        if (max < j[i])
                            max = j[i];
                    }
                }
                mesh->addDataVector("min angle, gradus", j);
            }
        }
        std::cout << "Выполнено: " << min << " <= min(alpha) <= " << max << std::endl;
    }
}

void MainWindow::on_actionAngleAspect_triggered()
{
    msh::MeshPointer mesh = ui->pictureControl->getGlMeshPicture()->getMesh();
    if (mesh != NULL)
    {
        ui->tabWidget->setCurrentIndex(2); // switch to terminal's tab
        std::cout << std::endl;
        std::cout << "Вычисление соотношений углов... ";
        double min = 0.0, max = 0.0;
        if (dynamic_cast<msh::TriangleMesh2D*>(mesh))
        {
            msh::TriangleMesh2D *triangles = dynamic_cast<msh::TriangleMesh2D*>(mesh);
            if (triangles->elementsCount() > 0)
            {
                std::vector<double> j(triangles->elementsCount());
                for (msh::UInteger i = 0; i < triangles->elementsCount(); i++)
                {
                    j[i] = triangles->angleAspect(i);
                    if (i == 0)
                    {
                        min = max = j[i];
                    }
                    else
                    {
                        if (min > j[i])
                            min = j[i];
                        if (max < j[i])
                            max = j[i];
                    }
                }
                mesh->addDataVector("angle ratio", j);
            }
        }
        std::cout << "Выполнено: " << min << "<= aspect(alpha) <= " << max << std::endl;
    }
}

void MainWindow::on_actionSegments_triggered()
{
    QString command = "var mesh = new Segments2D(xCount: Integer, yCount: Integer, origin: new Point2D(x, y), width: Floating, height: Floating, function (x, y) { return value; } [, points: Array]);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionTrianglesScript_triggered()
{
    QString command = "var mesh = new Triangles2D(xCount: Integer, yCount: Integer, origin: new Point2D(x, y), width: Floating, height: Floating, function (x, y) { return value; } [, points: Array]);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionQuadsScript_triggered()
{
    QString command = "var mesh = new Quads2D(xCount: Integer, yCount: Integer, origin: new Point2D(x, y), width: Floating, height: Floating, function (x, y) { return value; } [, points: Array]);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionDelaunayTriScript_triggered()
{
    QString command = "var mesh = new delaunay(xCount: Integer, yCount: Integer, origin: new Point2D(x, y), width: Floating, height: Floating, function (x, y) { return value; } [, points: Array]);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionRuppertTriScript_triggered()
{
    QString command = "var mesh = new ruppert(xCount: Integer, yCount: Integer, origin: new Point2D(x, y), width: Floating, height: Floating, function (x, y) { return value; } [, points: Array]);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionRectangularTriScript_triggered()
{
    QString command = "var mesh = new Triangles2D(xCount: Integer, yCount: Integer, origin: new Point2D(x, y), width: Floating, height: Floating);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionRectangularQuadsScript_triggered()
{
    QString command = "var mesh = new Quads2D(xCount: Integer, yCount: Integer, origin: new Point2D(x, y), width: Floating, height: Floating);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionSumScript_triggered()
{
    QString command = "sum(a, b, c, ...)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionConjunctionScript_triggered()
{
    QString command = "con(a, b, c, ...)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionDisjunctionScript_triggered()
{
    QString command = "dis(a, b, c, ...)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionDifferenceScript_triggered()
{
    QString command = "diff(a, b, c, ...)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionCylinderTriScript_triggered()
{
    QString command = "var mesh = CylinderTriangles(rCount: Integer, lCount: Integer, radius: Floating, length: Floating[, function (x, y, z) { return value; }]);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionCylinderQuadsScript_triggered()
{
    QString command = "var mesh = CylinderQuads(rCount: Integer, lCount: Integer, radius: Floating, length: Floating[, function (x, y, z) { return value; }]);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionConeTriScript_triggered()
{
    QString command = "var mesh = ConeTriangles(rCount: Integer, lCount: Integer, radiusBottom: Floating, radiusTop: Floating, length: Floating[, function (x, y, z) { return value; }]);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionConeQuadsScript_triggered()
{
    QString command = "var mesh = ConeQuads(rCount: Integer, lCount: Integer, radiusBottom: Floating, radiusTop: Floating, length: Floating[, function (x, y, z) { return value; }]);";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionCircleScript_triggered()
{
    QString command = "circle(x, y, r)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionEllipseScript_triggered()
{
    QString command = "ellipse(x, y, a, b)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionBandScript_triggered()
{
    QString command = "band(x, w)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionLineScript_triggered()
{
    QString command = "line(x, y, x1, y1, x2, y2)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionRectangleScript_triggered()
{
    QString command = "rectangle(x, y, w, h)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionRoundRectangleScript_triggered()
{
    QString command = "rectangle(x, y, w, h, r)";
    ui->codeEditor->insertPlainText(command);
}

void MainWindow::on_actionChangePictureBackground_triggered()
{
    QColor color = QColorDialog::getColor(ui->pictureControl->getGlMeshPicture()->backgroundColor());
    if (color.isValid())
        ui->pictureControl->getGlMeshPicture()->setBackgroundColor(color);
}

void MainWindow::on_actionChangeMeshColor_triggered()
{
    QColor color = QColorDialog::getColor(ui->pictureControl->getGlMeshPicture()->meshColor());
    if (color.isValid())
        ui->pictureControl->getGlMeshPicture()->setMeshColor(color);
}

void MainWindow::on_actionChangeElementColor_triggered()
{
    QColor color = QColorDialog::getColor(ui->pictureControl->getGlMeshPicture()->elementColor());
    if (color.isValid())
        ui->pictureControl->getGlMeshPicture()->setElementColor(color);
}
