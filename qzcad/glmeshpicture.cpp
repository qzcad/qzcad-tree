#include "glmeshpicture.h"
#include "segmentmesh2d.h"
#include "point3d.h"

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
    visualizationMode_ = COLOR;
    isMousePressed_ = false;
    showColorBar_ = false;
    isLighting_ = true;
    isUseVector_ = false;
    isShowInitialFrames = false;

#ifdef Q_OS_WIN
    setFont(QFont("Courier New", 14, QFont::Bold));
#else
    setFont(QFont("Monospace", 14, QFont::Bold));
#endif
    setDefault();
#ifdef Q_OS_WIN
    setFormat(QGLFormat(QGL::SingleBuffer | QGL::DepthBuffer));
#endif
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
//    updateGL ();
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
    centerX_ = centerY_ = centerZ_ = 0.0;
    vectorScale_ = 1.0;
    emit xRotationChanged(xRot_);
    emit yRotationChanged(yRot_);
    emit zRotationChanged(zRot_);
    emit zoomChanged(zoom_);
    emit translateStepChanged(translateStep_);
//    resetProjectionMatrix (); // на некоторых системах приводит к краху при вызове из конструткора
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
        glDisable(GL_DEPTH_TEST);
        if (mesh_->dimesion() == 2)
        {
            renderText (xmin, ymin, 0, "(" + QString::number(mesh_->xMin()) + "; " + QString::number(mesh_->yMin()) + ")" );
            renderText (xmax, ymin, 0, "(" + QString::number(mesh_->xMax()) + "; " + QString::number(mesh_->yMin()) + ")" );
            renderText (xmax, ymax, 0, "(" + QString::number(mesh_->xMax()) + "; " + QString::number(mesh_->yMax()) + ")" );
            renderText (xmin, ymax, 0, "(" + QString::number(mesh_->xMin()) + "; " + QString::number(mesh_->yMax()) + ")" );
        }
        else
        {
            renderText (xmin, ymin, zmin, "(" + QString::number(mesh_->xMin()) + "; " + QString::number(mesh_->yMin()) + "; " + QString::number(mesh_->zMin()) + ")");
            renderText (xmax, ymin, zmin, "(" + QString::number(mesh_->xMax()) + "; " + QString::number(mesh_->yMin()) + "; " + QString::number(mesh_->zMin()) + ")");
            renderText (xmin, ymax, zmin, "(" + QString::number(mesh_->xMin()) + "; " + QString::number(mesh_->yMax()) + "; " + QString::number(mesh_->zMin()) + ")");
            renderText (xmin, ymin, zmax, "(" + QString::number(mesh_->xMin()) + "; " + QString::number(mesh_->yMin()) + "; " + QString::number(mesh_->zMax()) + ")");
            renderText (xmax, ymax, zmin, "(" + QString::number(mesh_->xMax()) + "; " + QString::number(mesh_->yMax()) + "; " + QString::number(mesh_->zMin()) + ")");
            renderText (xmax, ymin, zmax, "(" + QString::number(mesh_->xMax()) + "; " + QString::number(mesh_->yMin()) + "; " + QString::number(mesh_->zMax()) + ")");
            renderText (xmin, ymax, zmax, "(" + QString::number(mesh_->xMin()) + "; " + QString::number(mesh_->yMax()) + "; " + QString::number(mesh_->zMax()) + ")");
            renderText (xmax, ymax, zmax, "(" + QString::number(mesh_->xMax()) + "; " + QString::number(mesh_->yMax()) + "; " + QString::number(mesh_->zMax()) + ")");
        }
        glEnable(GL_DEPTH_TEST);
    }
}

void GLMeshPicture::resetProjectionMatrix()
{
    double dx = fabs(xMax_ - xMin_);
    double dy = fabs(yMax_ - yMin_);
    double dz = fabs(zMax_ - zMin_);
    double max = qMax(dx, qMax(dy, dz));
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
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

msh::Point3D GLMeshPicture::pointToScenePoint(const msh::PointPointer &point, double dx, double dy, double dz) const
{
    return msh::Point3D(point->x() + dx - centerX_, point->y() + dy - centerY_, point->z() + dz - centerZ_);
}

void GLMeshPicture::drawColorBar()
{
    double min = map_.min();
    double max = map_.max();
    double minVal = -1.0;
    double maxVal = 1.0;
    const double width = 0.04 * (maxVal - minVal);
    const double length = 0.2 * (maxVal - minVal);
    const double top = -1.0 + length + 0.1;
    const double right = 0.01 + minVal + width;
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
    if (isLighting_) glDisable(GL_LIGHTING);

    int nc = 9;
    double bar_step = length / (double)nc;
    double value_step = (max - min) / (double)nc;
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    for (int i = 0; i < nc; i++)
    {
        qglColor(map_.color(min + (double)(i+1) * value_step));
        glBegin(GL_POLYGON);
        glVertex3d(right - width, top - length + (double)i * bar_step, maxVal);
        glVertex3d(right, top - length + (double)i * bar_step, maxVal);
        glVertex3d(right, top - length + (double)(i+1) * bar_step, maxVal);
        glVertex3d(right - width, top - length + (double)(i+1) * bar_step, maxVal);
        glEnd();
    }
    glDisable(GL_POLYGON_OFFSET_FILL);
    qglColor(meshColor_);
    for (int i = 0; i < nc; i++)
    {
        glBegin(GL_LINE_LOOP);
        glVertex3d(right - width, top - length + (double)i * bar_step, maxVal);
        glVertex3d(right, top - length + (double)i * bar_step, maxVal);
        glVertex3d(right, top - length + (double)(i+1) * bar_step, maxVal);
        glVertex3d(right - width, top - length + (double)(i+1) * bar_step, maxVal);
        glEnd();
    }
    qglColor (textColor_);
    glDisable(GL_DEPTH_TEST);
    for (int i = 0; i < nc; i++)
        renderText (right+0.01, top - length + (double)i * bar_step - 0.4 * bar_step, maxVal, QString::number(min + (double)(i) * value_step));
    renderText (right+0.01, top - 0.4 * bar_step, maxVal, QString::number(max));
    if (mesh_ && valueIndex_ < mesh_->dataVectorsCount())
        renderText (right + 0.01, top + 0.6 * bar_step, maxVal, QString::fromStdString(mesh_->data(valueIndex_).name()));
    glEnable(GL_DEPTH_TEST);
    if (isLighting_) glEnable(GL_LIGHTING);
    // востановление исходного режима проекции
    resetProjectionMatrix();
}

void GLMeshPicture::drawAxesDirection()
{
    const double border = 1.0;
    const double spacing = 0.3;
    const double axesLength = 0.2;
    const double c = border - spacing - axesLength;
    const double e = border - spacing;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
#ifdef QT_OPENGL_ES_1
    glOrthof(-border, border,
             -border, border,
             -border, border);
#else
    glOrtho(-border, border,
            -border, border,
            -border, border);
#endif
    // отрисовка котрольной полосы
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glRotatef(xRot_, 1.0, 0.0, 0.0);
    glRotatef(yRot_, 0.0, 1.0, 0.0);
    glRotatef(zRot_, 0.0, 0.0, 1.0);

    // x
    glColor3ub(255, 0, 0);
    glBegin(GL_LINES);
        glVertex3d(c, c, c);
        glVertex3d(e, c, c);
    glEnd();

    // y
    qglColor(Qt::green);
    glBegin(GL_LINES);
        glVertex3d(c, c, c);
        glVertex3d(c, e, c);
    glEnd();

    // z
    qglColor(Qt::blue);
    glBegin(GL_LINES);
        glVertex3d(c, c, c);
        glVertex3d(c, c, e);
    glEnd();
    // подписи к осям
    qglColor (textColor_);
    glDisable(GL_DEPTH_TEST);
    renderText (e, c, c, "x");
    renderText (c, e, c, "y");
    renderText (c, c, e, "z");
    glEnable(GL_DEPTH_TEST);
    // востановление исходного режима проекции
    resetProjectionMatrix();
}

void GLMeshPicture::drawFace(const msh::UIntegerVector &face, GLenum mode, GLfloat width, GLfloat size, bool useNodeColors, msh::NodeType filter)
{
    glLineWidth(width);
    glPointSize(size);
    msh::UInteger face_size = face.size();
    if (useNodeColors && (face_size == 3 || face_size == 4) && visualizationMode_ == VALUE && valueIndex_ < mesh_->dataVectorsCount() && mesh_->data(valueIndex_).size() == mesh_->nodesCount())
    {
        bool usable = true;
        msh::Point3D p[face_size];
        double v[face_size];
        for (msh::UInteger j = 0; j < face_size; j++)
        {
            if (!(filter == msh::UNDEFINED || filter == mesh_->nodeType(face[j])))
            {
                usable = false;
                break;
            }

            if (isUseVector_ && mesh_->dataVectorsCount() >= 2 && mesh_->dimesion() == 2 && mesh_->data(0).size() == mesh_->nodesCount() && mesh_->data(1).size() == mesh_->nodesCount())
                p[j] = pointToScenePoint(mesh_->node(face[j]),
                                      vectorScale_ * mesh_->data(0)[face[j]],
                                      vectorScale_ * mesh_->data(1)[face[j]]);
            else if (isUseVector_ && mesh_->dataVectorsCount() >= 3 && mesh_->dimesion() == 3 && mesh_->data(0).size() == mesh_->nodesCount() && mesh_->data(1).size() == mesh_->nodesCount() && mesh_->data(2).size() == mesh_->nodesCount())
                p[j] = pointToScenePoint(mesh_->node(face[j]),
                                      vectorScale_ * mesh_->data(0)[face[j]],
                                      vectorScale_ * mesh_->data(1)[face[j]],
                                      vectorScale_ * mesh_->data(2)[face[j]]);
            else
                p[j] = pointToScenePoint(mesh_->node(face[j]));

            v[j] = mesh_->data(valueIndex_)[face[j]];
        }
        if (!usable)
            return;
        drawTriangle(p[0], p[1], p[2], v[0], v[1], v[2], 1);
        if (face_size == 4)
            drawTriangle(p[0], p[2], p[3], v[0], v[2], v[3], 1);
    }
    else
    {
        glBegin(mode);
        for (msh::UInteger j = 0; j < face.size(); j++)
        {
            if (filter == msh::UNDEFINED || filter == mesh_->nodeType(face[j]))
            {
                if (useNodeColors && visualizationMode_ == VALUE && valueIndex_ < mesh_->dataVectorsCount() && mesh_->data(valueIndex_).size() == mesh_->nodesCount())
                    qglColor( map_.color( mesh_->data(valueIndex_)[face[j]] ) );
                if (isUseVector_ && mesh_->dataVectorsCount() >= 2 && mesh_->dimesion() == 2 && mesh_->data(0).size() == mesh_->nodesCount() && mesh_->data(1).size() == mesh_->nodesCount())
                    pointToGLVertex(mesh_->node(face[j]),
                                    vectorScale_ * mesh_->data(0)[face[j]],
                            vectorScale_ * mesh_->data(1)[face[j]]);
                else if (isUseVector_ && mesh_->dataVectorsCount() >= 3 && mesh_->dimesion() == 3 && mesh_->data(0).size() == mesh_->nodesCount() && mesh_->data(1).size() == mesh_->nodesCount() && mesh_->data(2).size() == mesh_->nodesCount())
                    pointToGLVertex(mesh_->node(face[j]),
                                    vectorScale_ * mesh_->data(0)[face[j]],
                            vectorScale_ * mesh_->data(1)[face[j]],
                            vectorScale_ * mesh_->data(2)[face[j]]);
                else
                    pointToGLVertex(mesh_->node(face[j]));
            }
        }
        glEnd();
    }
    glLineWidth(1.0);
    glPointSize(1.0);
}

void GLMeshPicture::drawTriangle(const msh::Point3D &p0, const msh::Point3D &p1, const msh::Point3D &p2, double v0, double v1, double v2, int level)
{
    QColor c0 = map_.color( v0 );
    QColor c1 = map_.color( v1 );
    QColor c2 = map_.color( v2 );
    if ((c0 == c1 && c1 == c2) || level >= 6)
    {
        glBegin(GL_TRIANGLES);
        qglColor(c0);
        glVertex3d(p0.x(), p0.y(), p0.z());
        qglColor(c1);
        glVertex3d(p1.x(), p1.y(), p1.z());
        qglColor(c2);
        glVertex3d(p2.x(), p2.y(), p2.z());
        glEnd();
    }
    else
    {
        msh::Point3D p01 = 0.5 * (p0 + p1);
        double v01 = 0.5 * (v0 + v1);
        msh::Point3D p12 = 0.5 * (p1 + p2);
        double v12 = 0.5 * (v1 + v2);
        msh::Point3D p20 = 0.5 * (p2 + p0);
        double v20 = 0.5 * (v2 + v0);
        drawTriangle(p0, p01, p20, v0, v01, v20, level + 1);
        drawTriangle(p1, p12, p01, v1, v12, v01, level + 1);
        drawTriangle(p2, p20, p12, v2, v20, v12, level + 1);
        drawTriangle(p01, p12, p20, v01, v12, v20, level + 1);
    }
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
    initializeGL();
    updateGL();
}

void GLMeshPicture::setMeshColor(QColor color)
{
    meshColor_ = color;
    updateGL();
}

void GLMeshPicture::setElementColor(QColor color)
{
    elementColor_ = color;
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
        if (visualizationMode_ == VALUE && mesh_->dataVectorsCount() > 0)
        {
            map_.setMin(mesh_->data(valueIndex_).min());
            map_.setMax(mesh_->data(valueIndex_).max());
        }
    }
    updateGL();
}

void GLMeshPicture::setColormap(int mapID)
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
    if (mesh_)
    {
        if (mesh_->dataVectorsCount() > 0)
        {
            ++valueIndex_;
            if (valueIndex_ >= mesh_->dataVectorsCount())
                valueIndex_ = 0;

            map_.setMin(mesh_->data(valueIndex_).min());
            map_.setMax(mesh_->data(valueIndex_).max());
            updateGL();
        }
    }
}

void GLMeshPicture::prevValueIndex()
{
    if (mesh_)
    {
        if (mesh_->dataVectorsCount() > 0)
        {
            if (valueIndex_ == 0)
                valueIndex_ = mesh_->dataVectorsCount() - 1;
            else
                --valueIndex_;

            map_.setMin(mesh_->data(valueIndex_).min());
            map_.setMax(mesh_->data(valueIndex_).max());
            updateGL();
        }
    }
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

void GLMeshPicture::setShowInitialFrames(bool isShow)
{
    isShowInitialFrames = isShow;
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
            msh::ElementPointer element = mesh_->element(i);
            if (dynamic_cast<const msh::Segment *>(element) != NULL)
            {
                glNormal3d(0.0, 0.0, -1.0);
                msh::UIntegerVector face = element->face(0);
                if (visualizationMode_ == VALUE  && valueIndex_ < mesh_->dataVectorsCount() && mesh_->data(valueIndex_).size() == mesh_->elementsCount())
                    qglColor(map_.color(mesh_->data(valueIndex_)[i], 256));
                else
                    qglColor(elementColor_);
                drawFace(face, GL_LINES, 2.0, 1.0);
                drawFace(face, GL_POINTS, 1.0, 6.0);
            }
            else if (mesh_->dimesion() == 2 || mesh_->isBorderElement(element))
            {
                for (int p = 0; p < element->facesCount(); p++)
                {
                    msh::UIntegerVector face = element->face(p);
                    if (mesh_->dimesion() == 3 && !mesh_->isBorderFace(face))
                    {
                        continue; // для трехмерных объектов рисуются только наружные грани
                    }

                    if (!isMousePressed_)
                    {
                        if (showMesh_)
                        {
                            glEnable(GL_POLYGON_OFFSET_FILL);
                            glPolygonOffset(1.0, 1.0);
                        }
                        if (visualizationMode_ == VALUE  && valueIndex_ < mesh_->dataVectorsCount() && mesh_->data(valueIndex_).size() == mesh_->elementsCount())
                            qglColor(map_.color(mesh_->data(valueIndex_)[i], 1024));
                        else
                            qglColor(elementColor_);
                        // нормаль к многоугольнику
                        msh::Point3D n = mesh_->normal(face);
                        glNormal3d(n.x(), n.y(), n.z());;
                        drawFace(face, GL_POLYGON);
                        if(showMesh_)
                        {
                            glDisable(GL_POLYGON_OFFSET_FILL);
                            qglColor(meshColor_);
                            drawFace(face, GL_LINE_LOOP, 1.0, 1.0, false);
                            drawFace(face, GL_POINTS, 1.0, 6.0, false, msh::CHARACTER);
                        }
                        if (isUseVector_ && isShowInitialFrames)
                        {
                            // отрисовка исходного каркаса
                            glColor3f(0.5, 0.5, 0.5);
                            glBegin(GL_LINE_LOOP);
                            for (msh::UInteger j = 0; j < face.size(); j++)
                            {
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

    drawAxesDirection();

    if (showColorBar_) drawColorBar();

    glPopMatrix();

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
    Q_UNUSED(event);
}

msh::NamedDoubleVector GLMeshPicture::dataVector() const
{
    msh::NamedDoubleVector data;
    if(mesh_ && valueIndex_ < mesh_->dataVectorsCount()) data = mesh_->data(valueIndex_);
    return data;
}

QColor GLMeshPicture::backgroundColor() const
{
    return backgroundColor_;
}

QColor GLMeshPicture::meshColor() const
{
    return meshColor_;
}

QColor GLMeshPicture::elementColor() const
{
    return elementColor_;
}

