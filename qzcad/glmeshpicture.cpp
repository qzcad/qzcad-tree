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
    isMousePressed_ = false;
    showColorBar_ = false;
    isLighting_ = true;
    isUseVector_ = false;
    setDefault();
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
    xMin_ = mesh->xMin ();
    xMax_ = mesh->xMax ();
    yMin_ = mesh->yMin ();
    yMax_ = mesh->yMax ();
    zMin_ = mesh->zMin ();
    zMax_ = mesh->zMax ();
    centerX_ = (xMax_ + xMin_) / 2.0;
    centerY_ = (yMax_ + yMin_) / 2.0;
    centerZ_ = (zMax_ + zMin_) / 2.0;

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

void GLMeshPicture::pushNodeValuesVector(const NamedFloatingVector &vector)
{
    nodeValues_.push_back(vector);
    if (nodeValues_.size() == 1)
    {
        setVisualisationMode(visualizationMode_); // необходимо обновить цветовую карту
    }
}

void GLMeshPicture::clearNodeValues()
{
    nodeValues_.clear();
    valueIndex_ = 0;
}

void GLMeshPicture::pushElementValuesVector(const NamedFloatingVector &vector)
{
    elementValues_.push_back(vector);
    if (elementValues_.size() == 1)
    {
        setVisualisationMode(visualizationMode_); // необходимо обновить цветовую карту
    }
}

void GLMeshPicture::clearElementValues()
{
    elementValues_.clear();
    valueIndex_ = 0;
}

void GLMeshPicture::setDefault()
{
    xRot_ = yRot_ = zRot_ = 0;
    zoom_ = 1.0;
    xTranslate_ = yTranslate_ = zTranslate_ = 0.0;
    translateStep_ = 0.01;
    xMin_ = yMin_ = zMin_ = -1.0;
    xMax_ = yMax_ = zMax_ = 1.0;
    centerX_ = centerY_ = centerZ_ = 0.0;
    vectorScale_ = 1.0;
    clearNodeValues();
    clearElementValues();
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
    const double xmin = xMin_ - centerX_;
    const double ymin = yMin_ - centerY_;
    const double zmin = zMin_ - centerZ_;
    const double xmax = xMax_ - centerX_;
    const double ymax = yMax_ - centerY_;
    const double zmax = zMax_ - centerZ_;
    glBegin (GL_LINE_STRIP);
    glColor3f (0.0f, 0.0f, 0.0f);
    glVertex3d (xmin, ymin, zmin);
    glColor3f (1.0f, 0.0f, 0.0f);
    glVertex3d (xmax, ymin, zmin);
    glColor3f (1.0f, 1.0f, 0.0f);
    glVertex3d (xmax, ymax, zmin);
    glColor3f (0.0f, 1.0f, 0.0f);
    glVertex3d (xmin, ymax, zmin);
    glColor3f (0.0f, 0.0f, 0.0f);
    glVertex3d (xmin, ymin, zmin);
    glEnd ();
    glBegin (GL_LINE_STRIP);
    glColor3f (0.0f, 0.0f, 1.0f);
    glVertex3d (xmin, ymin, zmax);
    glColor3f (1.0f, 0.0f, 1.0f);
    glVertex3d (xmax, ymin, zmax);
    glColor3f (1.0f, 1.0f, 1.0f);
    glVertex3d (xmax, ymax, zmax);
    glColor3f (0.0f, 1.0f, 1.0f);
    glVertex3d (xmin, ymax, zmax);
    glColor3f (0.0f, 0.0f, 1.0f);
    glVertex3d (xmin, ymin, zmax);
    glEnd ();
    glBegin (GL_LINES);
    glColor3f (0.0f, 0.0f, 0.0f);
    glVertex3d (xmin, ymin, zmin);
    glColor3f (0.0f, 0.0f, 1.0f);
    glVertex3d (xmin, ymin, zmax);
    glColor3f (0.0f, 1.0f, 0.0f);
    glVertex3d (xmin, ymax, zmin);
    glColor3f (0.0f, 1.0f, 1.0f);
    glVertex3d (xmin, ymax, zmax);
    glEnd ();
    glBegin (GL_LINES);
    glColor3f (1.0f, 0.0f, 0.0f);
    glVertex3d (xmax, ymin, zmin);
    glColor3f (1.0f, 0.0f, 1.0f);
    glVertex3d (xmax, ymin, zmax);
    glColor3f (1.0f, 1.0f, 0.0f);
    glVertex3d (xmax, ymax, zmin);
    glColor3f (1.0f, 1.0f, 1.0f);
    glVertex3d (xmax, ymax, zmax);
    glEnd ();
    // Отрисовка подписей в углах куба видимости
    //setFont (font().setPointSize(18));
    qglColor (textColor_);
    if (mesh_)
    {
        renderText (xmin, ymin, zmin,
                    QString::number(mesh_->xMin()) +
                    "; " +
                    QString::number(mesh_->yMin()) +
                    "; " +
                    QString::number(mesh_->zMin()));
        renderText (xmax, ymin, zmin,
                    QString::number(mesh_->xMax()) +
                    "; " +
                    QString::number(mesh_->yMin()) +
                    "; " +
                    QString::number(mesh_->zMin()));
        renderText (xmin, ymax, zmin,
                    QString::number(mesh_->xMin()) +
                    "; " +
                    QString::number(mesh_->yMax()) +
                    "; " +
                    QString::number(mesh_->zMin()));
        renderText (xmin, ymin, zmax,
                    QString::number(mesh_->xMin()) +
                    "; " +
                    QString::number(mesh_->yMin()) +
                    "; " +
                    QString::number(mesh_->zMax()));
        renderText (xmax, ymax, zmin,
                    QString::number(mesh_->xMax()) +
                    "; " +
                    QString::number(mesh_->yMax()) +
                    "; " +
                    QString::number(mesh_->zMin()));
        renderText (xmax, ymin, zmax,
                    QString::number(mesh_->xMax()) +
                    "; " +
                    QString::number(mesh_->yMin()) +
                    "; " +
                    QString::number(mesh_->zMax()));
        renderText (xmin, ymax, zmax,
                    QString::number(mesh_->xMin()) +
                    "; " +
                    QString::number(mesh_->yMax()) +
                    "; " +
                    QString::number(mesh_->zMax()));
        renderText (xmax, ymax, zmax,
                    QString::number(mesh_->xMax()) +
                    "; " +
                    QString::number(mesh_->yMax()) +
                    "; " +
                    QString::number(mesh_->zMax()));
    }
}

void GLMeshPicture::resetProjectionMatrix()
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double dx = fabs(xMax_ - xMin_);
    double dy = fabs(yMax_ - yMin_);
    double dz = fabs(zMax_ - zMin_);
    double max = qMax(dx, qMax(dy, dz));
#ifdef QT_OPENGL_ES_1
    glOrthof(-max, max,
             -max, max,
             -max, max);
#else
    glOrtho(-max, max,
            -max, max,
            -max, max);
#endif
    glMatrixMode(GL_MODELVIEW);
}

void GLMeshPicture::pointToGLVertex(const msh::PointPointer &point, double dx, double dy, double dz) const
{
    glVertex3d(point->x() + dx - centerX_, point->y() + dy - centerY_, point->z() + dz - centerZ_);
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
    if (visualizationMode_ == NODE_VALUE && valueIndex_ < nodeValues_.size())
    {
        renderText (right, top + 0.07, maxVal, nodeValues_[valueIndex_].name());
    }
    else if (visualizationMode_ == ELEMENT_VALUE && valueIndex_ < elementValues_.size())
    {
        renderText (right, top + 0.07, maxVal, elementValues_[valueIndex_].name());
    }
    if (isLighting_) glDisable(GL_LIGHTING);
    glBegin(GL_POLYGON);
    qglColor(map_.color(min));
    glVertex3d(right - width, top - length, maxVal);
    glVertex3d(right, top - length, maxVal);
    qglColor(map_.color(max));
    glVertex3d(right, top, maxVal);
    glVertex3d(right - width, top, maxVal);
    glEnd();
    if (isLighting_) glEnable(GL_LIGHTING);
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
        valueIndex_ = 0;
        if (visualizationMode_ == ELEMENT_VALUE  && valueIndex_ < elementValues_.size())
        {
            map_.setMin(elementValues_[valueIndex_].min());
            map_.setMax(elementValues_[valueIndex_].max());
        }
        else if (visualizationMode_ == NODE_VALUE && valueIndex_ < nodeValues_.size())
        {
            map_.setMin(nodeValues_[valueIndex_].min());
            map_.setMax(nodeValues_[valueIndex_].max());
        }
    }
    updateGL();
}

void GLMeshPicture::setColormapName(int mapID)
{
    map_.setColorMap(static_cast<ColorValueMap::ColorMapName>(mapID));
    updateGL();
}

void GLMeshPicture::setIsLighting(bool isLighting)
{
    isLighting_ = isLighting;
    if (isLighting_)
        glEnable(GL_LIGHTING);
    else
        glDisable(GL_LIGHTING);
    updateGL();
}

void GLMeshPicture::nextValueIndex()
{
    valueIndex_ = valueIndex_ + 1;
    if ((visualizationMode_ == NODE_VALUE && valueIndex_ >= nodeValues_.size()) || (visualizationMode_ == ELEMENT_VALUE && valueIndex_ >= elementValues_.size()))
        valueIndex_= 0;
    if (visualizationMode_ == NODE_VALUE)
    {
        map_.setMin(nodeValues_[valueIndex_].min());
        map_.setMax(nodeValues_[valueIndex_].max());
    }
    else if (visualizationMode_ == ELEMENT_VALUE)
    {
        map_.setMin(elementValues_[valueIndex_].min());
        map_.setMax(elementValues_[valueIndex_].max());
    }
    updateGL();
}

void GLMeshPicture::setIsUseVector(bool isUseVector)
{
    isUseVector_ = isUseVector;
    updateGL();
}

void GLMeshPicture::setVectorScale(double vectorScale)
{
    vectorScale_ = vectorScale;
    updateGL();
}

void GLMeshPicture::initializeGL()
{
    qglClearColor(backgroundColor_);
    // включаем освещение
    if (isLighting_) glEnable(GL_LIGHTING);
    // вкдючаем нулевой источник света
    if (isLighting_) glEnable(GL_LIGHT0);
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


                    if (!isMousePressed_)
                    {
                        if (showMesh_)
                        {
                            glEnable(GL_POLYGON_OFFSET_FILL);
                            glPolygonOffset(1.0, 1.0);
                        }
                        // вычисление нормали к грани
                        double nx = (b->y() - a->y()) * (c->z() - a->z()) - (b->z() - a->z()) * (c->y() - a->y());
                        double ny = (b->z() - a->z()) * (c->x() - a->x()) - (b->x() - a->x()) * (c->z() - a->z());
                        double nz = (b->x() - a->x()) * (c->y() - a->y()) - (b->y() - a->y()) * (c->x() - a->x());
                        double nn = sqrt(nx * nx + ny * ny + nz * nz);
                        // нормализация
                        nx = nx / nn;
                        ny = ny / nn;
                        nz = nz / nn;
                        if (visualizationMode_ == ELEMENT_VALUE  && valueIndex_ < elementValues_.size())
                            qglColor(map_.color(elementValues_[valueIndex_][i]));
                        else if (visualizationMode_ == USER_COLOR || visualizationMode_ == ELEMENT_VALUE)
                            qglColor(elementColor_);
                        // нормаль к многоугольнику
                        if(mesh_->dimesion() == 2)
                            glNormal3d(nx, ny, nz);
                        else
                            glNormal3d(-nx, -ny, -nz);

                        glBegin(GL_POLYGON);
                        for (msh::UInteger j = 0; j < face.size(); j++)
                        {
                            if (visualizationMode_ == NODE_VALUE && valueIndex_ < nodeValues_.size())
                                qglColor( map_.color( nodeValues_[valueIndex_][face[j]] ) );
                            else if (visualizationMode_ == NODE_VALUE)
                                qglColor(elementColor_); // если индекс вне диапазона, то цветом пользователя
                            if (isUseVector_ && nodeValues_.size() >= 3)
                                pointToGLVertex(mesh_->node(face[j]),
                                                vectorScale_ * nodeValues_[0][face[j]],
                                        vectorScale_ * nodeValues_[1][face[j]],
                                        vectorScale_ * nodeValues_[2][face[j]]);
                            else
                                pointToGLVertex(mesh_->node(face[j]));
                        }
                        glEnd();
                        if(showMesh_)
                        {
                            glDisable(GL_POLYGON_OFFSET_FILL);
                            qglColor(meshColor_);
                            glBegin(GL_LINE_LOOP);
                            for (msh::UInteger j = 0; j < face.size(); j++)
                            {
                                if (isUseVector_ && nodeValues_.size() >= 3)
                                    pointToGLVertex(mesh_->node(face[j]),
                                                    vectorScale_ * nodeValues_[0][face[j]],
                                            vectorScale_ * nodeValues_[1][face[j]],
                                            vectorScale_ * nodeValues_[2][face[j]]);
                                else
                                    pointToGLVertex(mesh_->node(face[j]));
                            }
                            glEnd();
                        }
                    }
                    else
                    {
                        qglColor(meshColor_);
                        glBegin(GL_POINTS);
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
    isMousePressed_ = true;
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
    isMousePressed_ = false;
    updateGL();
}
