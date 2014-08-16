#include "glmeshpicture.h"

GLMeshPicture::GLMeshPicture(QWidget *parent) :
    QGLWidget(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer), parent)
{
    mesh_ = NULL;
    showMesh_ = showCube_ = true;
    backgroundColor_ = QColor (Qt::white);
    elementColor_ = QColor (Qt::green);
    meshColor_ = QColor (Qt::black);
    textColor_ = QColor (Qt::blue);
    mouseMode_ = ROTATION;
    visualizationMode_ = USER_COLOR;
    //    colorMapName_ = ColorValueMap::WINTER;
    setDefault();
    isMousePressed = false;
    showColorBar_ = false;
}

GLMeshPicture::~GLMeshPicture()
{
}

QSize GLMeshPicture::minimumSizeHint() const
{
    return QSize (50, 50);
}

QSize GLMeshPicture::sizeHint() const
{
    return QSize (500, 500);
}

msh::MeshPointer GLMeshPicture::getMesh()
{
    return mesh_;
}

void GLMeshPicture::setMesh(msh::MeshPointer mesh)
{
    setDefault();
    mesh_ = mesh;
    //    xMin_ = mesh->xMin ();
    //    xMax_ = mesh->xMax ();
    //    yMin_ = mesh->yMin ();
    //    yMax_ = mesh->yMax ();
    //    zMin_ = mesh->zMin ();
    //    zMax_ = mesh->zMax ();
    double dx = mesh->xMax () - mesh->xMin ();
    double dy = mesh->yMax () - mesh->yMin ();
    double dz = mesh->zMax () - mesh->zMin ();
    xMin_ = -1.0; xMax_ = 1.0;
    yMin_ = -1.0 * dy / dx;
    yMax_ = 1.0 * dy / dx;
    zMin_ = -1.0 * dz / dx;
    zMax_ = 1.0 * dz / dx;

    //    emit xMinChanged(xMin_);
    //    emit xMaxChanged(xMax_);
    //    emit yMinChanged(yMin_);
    //    emit yMaxChanged(yMax_);
    //    emit zMinChanged(zMin_);
    //    emit zMaxChanged(zMax_);
    resetProjectionMatrix ();
    setVisualisationMode(visualizationMode_); // необходимо обновить цветовую карту
    updateGL ();
}

void GLMeshPicture::resetMesh()
{
    //    if (mesh_)mesh_.reset();
    mesh_ = NULL;
    setDefault ();
    updateGL ();
}

msh::MeshPointer GLMeshPicture::releaseMesh()
{
    msh::MeshPointer meshPointer = mesh_;
    resetMesh();
    return meshPointer;
}

void GLMeshPicture::setDefault()
{
    xRot_ = yRot_ = zRot_ = 0;
    zoom_ = 1.0;
    xTranslate_ = yTranslate_ = zTranslate_ = 0.0;
    translateStep_ = 0.01;
    xMin_ = yMin_ = zMin_ = -1.0;
    xMax_ = yMax_ = zMax_ = 1.0;
    emit xRotationChanged(xRot_);
    emit yRotationChanged(yRot_);
    emit zRotationChanged(zRot_);
    emit zoomChanged(zoom_);
    emit translateStepChanged(translateStep_);
    resetProjectionMatrix ();
}

int GLMeshPicture::normalizedViewAngle(const int &angle)
{
    int a = angle;
    // отрицательный угол
    while (a < 0)
        a += 360;
    // поворот на 360 или больше градусов
    while (a >= 360)
        a -= 360;

    return a;
}

void GLMeshPicture::drawRegionBorders()
{
    // куб видимости
    glBegin (GL_LINE_STRIP);
    glColor3f (0.0f, 0.0f, 0.0f);
    glVertex3d (xMin_, yMin_, zMin_);
    glColor3f (1.0f, 0.0f, 0.0f);
    glVertex3d (xMax_, yMin_, zMin_);
    glColor3f (1.0f, 1.0f, 0.0f);
    glVertex3d (xMax_, yMax_, zMin_);
    glColor3f (0.0f, 1.0f, 0.0f);
    glVertex3d (xMin_, yMax_, zMin_);
    glColor3f (0.0f, 0.0f, 0.0f);
    glVertex3d (xMin_, yMin_, zMin_);
    glEnd ();
    glBegin (GL_LINE_STRIP);
    glColor3f (0.0f, 0.0f, 1.0f);
    glVertex3d (xMin_, yMin_, zMax_);
    glColor3f (1.0f, 0.0f, 1.0f);
    glVertex3d (xMax_, yMin_, zMax_);
    glColor3f (1.0f, 1.0f, 1.0f);
    glVertex3d (xMax_, yMax_, zMax_);
    glColor3f (0.0f, 1.0f, 1.0f);
    glVertex3d (xMin_, yMax_, zMax_);
    glColor3f (0.0f, 0.0f, 1.0f);
    glVertex3d (xMin_, yMin_, zMax_);
    glEnd ();
    glBegin (GL_LINES);
    glColor3f (0.0f, 0.0f, 0.0f);
    glVertex3d (xMin_, yMin_, zMin_);
    glColor3f (0.0f, 0.0f, 1.0f);
    glVertex3d (xMin_, yMin_, zMax_);
    glColor3f (0.0f, 1.0f, 0.0f);
    glVertex3d (xMin_, yMax_, zMin_);
    glColor3f (0.0f, 1.0f, 1.0f);
    glVertex3d (xMin_, yMax_, zMax_);
    glEnd ();
    glBegin (GL_LINES);
    glColor3f (1.0f, 0.0f, 0.0f);
    glVertex3d (xMax_, yMin_, zMin_);
    glColor3f (1.0f, 0.0f, 1.0f);
    glVertex3d (xMax_, yMin_, zMax_);
    glColor3f (1.0f, 1.0f, 0.0f);
    glVertex3d (xMax_, yMax_, zMin_);
    glColor3f (1.0f, 1.0f, 1.0f);
    glVertex3d (xMax_, yMax_, zMax_);
    glEnd ();
    // Отрисовка подписей в углах куба видимости
    //setFont (font().setPointSize(18));
    qglColor (textColor_);
    if (mesh_)
    {
        renderText (xMin_, yMin_, zMin_,
                    QString::number(mesh_->xMin()) +
                    "; " +
                    QString::number(mesh_->yMin()) +
                    "; " +
                    QString::number(mesh_->zMin()));
        renderText (xMax_, yMin_, zMin_,
                    QString::number(mesh_->xMax()) +
                    "; " +
                    QString::number(mesh_->yMin()) +
                    "; " +
                    QString::number(mesh_->zMin()));
        renderText (xMin_, yMax_, zMin_,
                    QString::number(mesh_->xMin()) +
                    "; " +
                    QString::number(mesh_->yMax()) +
                    "; " +
                    QString::number(mesh_->zMin()));
        renderText (xMin_, yMin_, zMax_,
                    QString::number(mesh_->xMin()) +
                    "; " +
                    QString::number(mesh_->yMin()) +
                    "; " +
                    QString::number(mesh_->zMax()));
        renderText (xMax_, yMax_, zMin_,
                    QString::number(mesh_->xMax()) +
                    "; " +
                    QString::number(mesh_->yMax()) +
                    "; " +
                    QString::number(mesh_->zMin()));
        renderText (xMax_, yMin_, zMax_,
                    QString::number(mesh_->xMax()) +
                    "; " +
                    QString::number(mesh_->yMin()) +
                    "; " +
                    QString::number(mesh_->zMax()));
        renderText (xMin_, yMax_, zMax_,
                    QString::number(mesh_->xMin()) +
                    "; " +
                    QString::number(mesh_->yMax()) +
                    "; " +
                    QString::number(mesh_->zMax()));
        renderText (xMax_, yMax_, zMax_,
                    QString::number(mesh_->xMax()) +
                    "; " +
                    QString::number(mesh_->yMax()) +
                    "; " +
                    QString::number(mesh_->zMax()));
    }else
    {
        renderText (xMin_, yMin_, zMin_, QString::number(xMin_) + "; " + QString::number(yMin_) + "; " + QString::number(zMin_));
        renderText (xMax_, yMin_, zMin_, QString::number(xMax_) + "; " + QString::number(yMin_) + "; " + QString::number(zMin_));
        renderText (xMin_, yMax_, zMin_, QString::number(xMin_) + "; " + QString::number(yMax_) + "; " + QString::number(zMin_));
        renderText (xMin_, yMin_, zMax_, QString::number(xMin_) + "; " + QString::number(yMin_) + "; " + QString::number(zMax_));
        renderText (xMax_, yMax_, zMin_, QString::number(xMax_) + "; " + QString::number(yMax_) + "; " + QString::number(zMin_));
        renderText (xMax_, yMin_, zMax_, QString::number(xMax_) + "; " + QString::number(yMin_) + "; " + QString::number(zMax_));
        renderText (xMin_, yMax_, zMax_, QString::number(xMin_) + "; " + QString::number(yMax_) + "; " + QString::number(zMax_));
        renderText (xMax_, yMax_, zMax_, QString::number(xMax_) + "; " + QString::number(yMax_) + "; " + QString::number(zMax_));
    }
}

void GLMeshPicture::resetProjectionMatrix()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
#ifdef QT_OPENGL_ES_1
    glOrthof(xMin_, xMax_,
             yMin_, yMax_,
             zMin_, zMax_);
#else
    double minVal = qMin(xMin_, qMin(yMin_, zMin_));
    double maxVal = qMax(xMax_, qMax(yMax_, zMax_));
    glOrtho(minVal, maxVal,
            minVal, maxVal,
            minVal, maxVal);
#endif
    glMatrixMode(GL_MODELVIEW);
}

void GLMeshPicture::pointToGLVertex(const msh::PointPointer &point) const
{
    double x = xMin_ + (xMax_ - xMin_) * (point->x() - mesh_->xMin()) / (mesh_->xMax() - mesh_->xMin());
    double y = yMin_ + (yMax_ - yMin_) * (point->y() - mesh_->yMin()) / (mesh_->yMax() - mesh_->yMin());
    double z = zMin_ + (zMax_ - zMin_) * (point->z() - mesh_->zMin()) / (mesh_->zMax() - mesh_->zMin());
    glVertex3d(x, y, z);
}

void GLMeshPicture::drawColorBar()
{
    double min = map_.min();
    double max = map_.max();
    double minVal = -1.0;
    double maxVal = 1.0;
    const double width = 0.04 * (maxVal - minVal);
    const double length = 0.2 * (maxVal - minVal);
    const double top = -1.0 + length;
    const double right = minVal + width;
    // изменение режима проекции для отрисовки контрольной полосы
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
#ifdef QT_OPENGL_ES_1
    glOrthof(minVal, maxVal,
             minVal, maxVal,
             minVal, maxVal);
#else
    glOrtho(minVal, maxVal,
            minVal, maxVal,
            minVal, maxVal);
#endif
    // отрисовка котрольной полосы
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    qglColor (textColor_);
    renderText (right, top - length, maxVal, QString::number(min));
    renderText (right, top, maxVal, QString::number(max));
    glDisable(GL_LIGHTING);
    glBegin(GL_POLYGON);
    qglColor(map_.color(min));
    glVertex3d(right - width, top - length, maxVal);
    glVertex3d(right, top - length, maxVal);
    qglColor(map_.color(max));
    glVertex3d(right, top, maxVal);
    glVertex3d(right - width, top, maxVal);
    glEnd();
    glEnable(GL_LIGHTING);
    // востановление исходного режима проекции
    resetProjectionMatrix();
}

void GLMeshPicture::setXRotation(int angle)
{
    int normalized = normalizedViewAngle(angle);
    if (xRot_ != normalized)
    {
        xRot_ = normalized;
        emit xRotationChanged(xRot_);
        updateGL();
    }
}

void GLMeshPicture::setYRotation(int angle)
{
    int normalized = normalizedViewAngle(angle);
    if (yRot_ != normalized)
    {
        yRot_ = normalized;
        emit yRotationChanged(yRot_);
        updateGL();
    }
}

void GLMeshPicture::setZRotation(int angle)
{
    int normalized = normalizedViewAngle(angle);
    if (zRot_ != normalized)
    {
        zRot_ = normalized;
        emit zRotationChanged(zRot_);
        updateGL();
    }
}

void GLMeshPicture::setZoom(double zoom)
{
    if (zoom_ != zoom)
    {
        zoom_ = zoom;
        emit zoomChanged(zoom_);
        updateGL();
    }
}

void GLMeshPicture::setXTranslation(double pos)
{
    if (xTranslate_ != pos)
    {
        xTranslate_ = pos;
        emit xTranslationChanged(xTranslate_);
        updateGL();
    }
}

void GLMeshPicture::setYTranslation(double pos)
{
    if (yTranslate_ != pos)
    {
        yTranslate_ = pos;
        emit yTranslationChanged(yTranslate_);
        updateGL();
    }
}

void GLMeshPicture::setZTranslation(double pos)
{
    if (zTranslate_ != pos)
    {
        zTranslate_ = pos;
        emit zTranslationChanged(zTranslate_);
        updateGL();
    }
}

void GLMeshPicture::setTranslateStep(double step)
{
    if (translateStep_ != step)
    {
        translateStep_ = step;
        emit translateStepChanged(translateStep_);
    }
}

void GLMeshPicture::setShowMesh(bool show)
{
    showMesh_ = show;
    updateGL();
}

void GLMeshPicture::setShowCube(bool show)
{
    showCube_ = show;
    updateGL();
}

void GLMeshPicture::setShowColorBar(bool show)
{
    showColorBar_ = show;
    updateGL();
}

void GLMeshPicture::setBackgroundColor(QColor color)
{
    backgroundColor_ = color;
    updateGL();
}

void GLMeshPicture::setTextColor(QColor color)
{
    textColor_ = color;
    updateGL();
}

void GLMeshPicture::setMouseMode(int mode)
{
    switch(mode)
    {
    case 0:
        mouseMode_ = ROTATION;
        break;
    case 1:
        mouseMode_ = TRANSLATION;
        break;
    case 2:
        mouseMode_ = SELECTION;
        break;
    default:
        mouseMode_ = ROTATION;
    }
}

void GLMeshPicture::activateDoubleBufferGL(bool activate)
{
    if (activate)
        setFormat(QGLFormat(QGL::DoubleBuffer | QGL::DepthBuffer));
    else
        setFormat(QGLFormat(QGL::SingleBuffer | QGL::DepthBuffer));
}

void GLMeshPicture::setVisualisationMode(int mode)
{
    visualizationMode_ = static_cast<VisualizationMode>(mode);
    if(mesh_)
    {
        if (visualizationMode_ == ELEMENT_VALUE)
        {
            msh::Floating minEV = mesh_->elementValue(0);
            msh::Floating maxEV = mesh_->elementValue(0);
            for (msh::UInteger i = 1; i < mesh_->elementsCount(); i++)
            {
                if (minEV > mesh_->elementValue(i)) minEV = mesh_->elementValue(i);
                if (maxEV < mesh_->elementValue(i)) maxEV = mesh_->elementValue(i);
            }
            map_.setMin(minEV);
            map_.setMax(maxEV);
        }
        else if (visualizationMode_ == NODE_VALUE)
        {
            msh::Floating minEV = mesh_->nodeValue(0);
            msh::Floating maxEV = mesh_->nodeValue(0);
            for (msh::UInteger i = 1; i < mesh_->nodesCount(); i++)
            {
                if (minEV > mesh_->nodeValue(i)) minEV = mesh_->nodeValue(i);
                if (maxEV < mesh_->nodeValue(i)) maxEV = mesh_->nodeValue(i);
            }
            map_.setMin(minEV);
            map_.setMax(maxEV);
        }
    }
    updateGL();
}

void GLMeshPicture::setColormapName(int mapID)
{
    map_.setColorMap(static_cast<ColorValueMap::ColorMapName>(mapID));
    updateGL();
}

void GLMeshPicture::initializeGL()
{
    qglClearColor(backgroundColor_);
    // включаем освещение
    glEnable(GL_LIGHTING);
    // вкдючаем нулевой источник света
    glEnable(GL_LIGHT0);
    //    glEnable(GL_LIGHT1);

    // включаем пересчет нормалей при масштабировании
#if defined(GL_RESCALE_NORMAL)
    glEnable(GL_RESCALE_NORMAL);
#else
    glEnable(GL_NORMALIZE);
#endif

    glEnable ( GL_COLOR_MATERIAL );
    glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;

    //glFrontFace(GL_CCW);
    glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glHint( GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST );

    glEnable(GL_DEPTH_TEST);

    glShadeModel(GL_SMOOTH);

#if defined(GL_MULTISAMPLE)
    glEnable(GL_MULTISAMPLE);
#endif
}

void GLMeshPicture::paintGL()
{
    // координаты источника света
    float light_position0[] = {0.0f, 0.0f, -10.0f * (float)qMax(xMax_, qMax(yMax_, zMax_)), 0.0f};
    // Для трехмерных объектов оставляем камеру неподвижной
    if (mesh_ && mesh_->dimesion() == 3)
        glLightfv (GL_LIGHT0, GL_POSITION, light_position0);
    // сброс буферов
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
    // загрузка единичной матрицы
    glLoadIdentity();
    // зум
    glScaled(zoom_, zoom_, zoom_);
    // вращение
    glRotatef(xRot_, 1.0, 0.0, 0.0);
    glRotatef(yRot_, 0.0, 1.0, 0.0);
    glRotatef(zRot_, 0.0, 0.0, 1.0);
    // движение
    glTranslated(xTranslate_, yTranslate_, zTranslate_);
    // для двумерных объетов камера двигается одновременно со сценой
    if (mesh_ && mesh_->dimesion() == 2)
        glLightfv (GL_LIGHT0, GL_POSITION, light_position0);
    // границы области
    if (showCube_) drawRegionBorders();
    if (mesh_)
    {
        for (msh::UInteger i = 0; i < mesh_->elementsCount(); i++)
        {
            if (mesh_->dimesion() == 2 || mesh_->isBorderElement(i))
            {
                for (int p = 0; p < mesh_->element(i)->facesCount(); p++)
                {
                    msh::UIntegerVector face = mesh_->element(i)->face(p);
                    msh::PointPointer a = mesh_->node(face[1]);
                    msh::PointPointer b = mesh_->node(face[0]);
                    msh::PointPointer c = mesh_->node(face[2]);

                    if (showMesh_ && !isMousePressed)
                    {
                        glEnable(GL_POLYGON_OFFSET_FILL);
                        glPolygonOffset(1.0, 1.0);
                    }
                    if (!isMousePressed)
                    {
                        // вычисление нормали к грани
                        msh::Floating nx = (b->y() - a->y()) * (c->z() - a->z()) - (b->z() - a->z()) * (c->y() - a->y());
                        msh::Floating ny = (b->z() - a->z()) * (c->x() - a->x()) - (b->x() - a->x()) * (c->z() - a->z());
                        msh::Floating nz = (b->x() - a->x()) * (c->y() - a->y()) - (b->y() - a->y()) * (c->x() - a->x());
                        msh::Floating nn = sqrt(nx * nx + ny * ny + nz * nz);
                        // нормализация
                        nx = nx / nn;
                        ny = ny / nn;
                        nz = nz / nn;
                        if (visualizationMode_ == ELEMENT_VALUE)
                            qglColor(map_.color(mesh_->elementValue(i)));
                        else if (visualizationMode_ == USER_COLOR)
                            qglColor(elementColor_);
                        // нормаль к многоугольнику
                        if(mesh_->dimesion() == 2)
                            glNormal3d(nx, ny, nz);
                        else
                            glNormal3d(-nx, -ny, -nz);

                        glBegin(GL_POLYGON);
                        for (msh::UInteger j = 0; j < face.size(); j++)
                        {
                            if (visualizationMode_ == NODE_VALUE)
                                qglColor(map_.color(mesh_->nodeValue(face[j])));
                            pointToGLVertex(mesh_->node(face[j]));
                        }
                        glEnd();
                    }
                    if(showMesh_ || isMousePressed)
                    {
                        glDisable(GL_POLYGON_OFFSET_FILL);
                        qglColor(meshColor_);
                        glBegin(GL_LINE_LOOP);
                        for (msh::UInteger j = 0; j < face.size(); j++)
                        {
                            pointToGLVertex(mesh_->node(face[j]));
                        }
                        glEnd();
                    }
                }
            }
        }
    }
    glPopMatrix();

    if (showColorBar_) drawColorBar();

    glFlush();
    swapBuffers();
}

void GLMeshPicture::resizeGL(int width, int height)
{
    int side = qMin(width, height);
    glViewport((width - side) / 2, (height - side) / 2, side, side);
    resetProjectionMatrix();
}

void GLMeshPicture::mousePressEvent(QMouseEvent *event)
{
    mousePressPos_ = event->pos();
    isMousePressed = true;
}

void GLMeshPicture::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - mousePressPos_.x();
    int dy = event->y() - mousePressPos_.y();
    if (mouseMode_ == ROTATION)
    {
        if (event->buttons() & Qt::LeftButton)
        {
            setXRotation(xRot_ + dy);
            setYRotation(yRot_ + dx);
        }else if (event->buttons() & Qt::RightButton)
        {
            setXRotation(xRot_ + dy);
            setZRotation(zRot_ + dx);
        }
    }else if (mouseMode_ == TRANSLATION)
    {
        if (event->buttons() & Qt::LeftButton)
        {
            setXTranslation(xTranslate_ + translateStep_ * dx);
            setYTranslation(yTranslate_ - translateStep_ * dy);
        }else if (event->buttons() & Qt::RightButton)
        {
            setXTranslation(xTranslate_ + translateStep_ * dx);
            setZTranslation(zTranslate_ + translateStep_ * dy);
        }
    }
    mousePressPos_ = event->pos();
}

void GLMeshPicture::wheelEvent(QWheelEvent *event)
{
    double numDegrees = -event->delta() / 8.0;
    double numSteps = numDegrees / 15.0;
    setZoom(zoom_ * pow(1.25, numSteps));
}

void GLMeshPicture::mouseReleaseEvent(QMouseEvent *event)
{
    isMousePressed = false;
    updateGL();
}
