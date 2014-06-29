/**
  * @author Сергей Чопоров
  * @version 1.0.3
  * @date 30/12/2013
  */
#ifndef GLMESHPICTURE_H
#define GLMESHPICTURE_H

#include <QtOpenGL>
#include <QWidget>

#include "meshpointer.h"
#include "colorvaluemap.h"

class GLMeshPicture : public QGLWidget
{
    Q_OBJECT
public:
    /**
     * @brief Конструктор
     * @param parent Указатель на родительский виджет
     */
    explicit GLMeshPicture (QWidget *parent = 0);
    /// Деструктор
    ~GLMeshPicture ();
    /**
     * @brief Минимальный размер виджета
     * @return 50x50
     */
    virtual QSize minimumSizeHint () const;
    /**
     * @brief Рекомендованный по умолчанию размер виджета
     * @return 500x500
     */
    virtual QSize sizeHint () const;
    /**
     * @brief Получить указатель на сетку
     * @return "Умный" указатель на актуальную сетку
     */
    msh::MeshPointer getMesh();
    /**
     * @brief Установить сетку для рисования
     * @param mesh Указатель на сетку
     */
    void setMesh (msh::MeshPointer mesh);
    /**
     * @brief Спросить указатель на сетку в NULL
     */
    void resetMesh();

signals:
    /**
     * @brief Сигнал об изменении угла поворота вокруг оси Ox
     * @param angle Угол поворота вокруг оси Ox
     */
    void xRotationChanged (int angle);
    /**
     * @brief Сигнал об изменении угла поворота вокруг оси Oy
     * @param angle Угол поворота вокруг оси Oy
     */
    void yRotationChanged (int angle);
    /**
     * @brief Сигнал об изменении угла поворота вокруг оси Oz
     * @param angle Угол поворота вокруг оси Oz
     */
    void zRotationChanged (int angle);
    /**
     * @brief Сигнал об изменении коэффициента приближения
     * @param zoom Значение коэффициента приближения
     */
    void zoomChanged (double zoom);
    /**
     * @brief Сигнал об изменении смещения вдоль оси Ox
     * @param step Значение смещения вдоль оси Ox
     */
    void xTranslationChanged (double step);
    /**
     * @brief Сигнал об изменении смещения вдоль оси Oy
     * @param step Значение смещения вдоль оси Oy
     */
    void yTranslationChanged (double step);
    /**
     * @brief Сигнал об изменении смещения вдоль оси Oz
     * @param step Значение смещения вдоль оси Oz
     */
    void zTranslationChanged (double step);
    /**
     * @brief Сигнал об изменении шага смещения
     * @param step Значение шага смещения
     */
    void translateStepChanged (double step);
//    void xMinChanged (double val);
//    void xMaxChanged (double val);
//    void yMinChanged (double val);
//    void yMaxChanged (double val);
//    void zMinChanged (double val);
//    void zMaxChanged (double val);

public slots:
    /**
     * @brief Установить угол вращения вокруг оси Ox
     * @param angle Значение угла вращения
     */
    void setXRotation (int angle);
    /**
     * @brief Установить угол вращения вокруг оси Oy
     * @param angle Значение угла вращения
     */
    void setYRotation (int angle);
    /**
     * @brief Установить угол вращения вокруг оси Oz
     * @param angle Значение угла вращения
     */
    void setZRotation (int angle);
    /**
     * @brief Установить коэффициент приближения
     * @param zoom Значение коэффициента приближения
     */
    void setZoom (double zoom);
    /**
     * @brief Установить смещение вдоль оси Ox
     * @param pos Координата точки смщения
     */
    void setXTranslation (double pos);
    /**
     * @brief Установить смещение вдоль оси Oy
     * @param pos Координата точки смщения
     */
    void setYTranslation (double pos);
    /**
     * @brief Установить смещение вдоль оси Oz
     * @param pos Координата точки смщения
     */
    void setZTranslation (double pos);
    /**
     * @brief Установить значение шага смещения
     * @param step Шаг смещения
     */
    void setTranslateStep (double step);
    /**
     * @brief Установить значение флага видимости сетки
     * @param show Значение флага видимости
     */
    void setShowMesh (bool show);
    /**
     * @brief Установить значение флага видимости границ области
     * @param show Значение флага
     */
    void setShowCube (bool show);
    /**
     * @brief Установить цвет фона
     * @param color Цвет фона
     */
    void setBackgroundColor (QColor color);
    /**
     * @brief Установить цвет текста
     * @param color Цвет текста
     */
    void setTextColor(QColor color);
    /**
     * @brief Установить режим обработки мыши
     * @param mode Код режима обработки мыши
     * @see MouseMode
     */
    void setMouseMode(int mode);
    /**
     * @brief Актививровать двойной буффер в OpenGL
     * @param activate Флаг активации двойного буффера. Если ativate = true, то включается режим двойного буффера OpenGL, иначе используется один буффер
     */
    void activateDoubleBufferGL(bool activate);
    /**
     * @brief Установить режим цветовой визуализации
     * @param mode Режим визуализации
     * @see VisaulizationMode
     */
    void setVisualisationMode(int mode);
    /**
     * @brief Установить номер схемы цветовой визуализации
     * @param mapID Номер схемы
     */
    void setColormapName(int mapID);

protected:
    /**
     * @brief Инициализация состояния openGL
     */
    virtual void initializeGL ();
    /**
     * @brief Рисование изображения
     */
    virtual void paintGL ();
    /**
     * @brief Изменить размеры изображения
     * @param width Ширина
     * @param height Высота
     */
    virtual void resizeGL ( int width, int height );
    /**
     * @brief Обработчик события нажатия кнопки мыши
     * @param event Указатель на событие
     */
    virtual void mousePressEvent(QMouseEvent *event);
    /**
     * @brief Обработчик события движения указателя мыши
     * @param event Указатель на событие
     */
    virtual void mouseMoveEvent(QMouseEvent *event);
    /**
     * @brief Обработчик события прокрутки
     * @param event Указатель на событие
     */
    virtual void wheelEvent(QWheelEvent *event);
    /**
     * @brief Обработчик события отпускания кнопки мыши
     * @param event Указатель на событие
     */
    virtual void mouseReleaseEvent(QMouseEvent * event);

private:
    msh::MeshPointer mesh_; //!< Указатель на сетку
    int xRot_; //!< Угол вращения вокруг оси Ox
    int yRot_; //!< Угол вращения вокруг оси Oy
    int zRot_; //!< Угол вращения вокруг оси Oz
    double zoom_; //!< Коэффициент масштабирования изображения
    double xTranslate_; //!< Смещение вдоль оси Ox
    double yTranslate_; //!< Смещение вдоль оси Oy
    double zTranslate_; //!< Смещение вдоль оси Oz
    double translateStep_; //!< Шаг смещения
    double xMin_; //!< Минимальное значение ординаты области видимости
    double xMax_; //!< Максимальное значение ординаты области видимости
    double yMin_; //!< Минимальное значение абсциссы области видимости
    double yMax_; //!< Максимальное значение абсциссы области видимости
    double zMin_; //!< Минимальное значение аппликаты области видимости
    double zMax_; //!< Максимальное значение аппликаты области видимости
    bool showMesh_; //!< Показать или нет сетку элементов
    bool showCube_; //!< Показать или нет куб видимости
    QColor backgroundColor_; //!< Цвет фона
    QColor textColor_; //!< Цвет текста
    QColor elementColor_; //!< Цвет элемента
    QColor meshColor_; //!< Цвет сетки
    QPoint mousePressPos_; //!< Точка нажатия кнопки мыши
    /**
     * @brief Возможные режимы обработки движения мыши
     */
    enum MouseMode { ROTATION = 0, TRANSLATION = 1, SELECTION = 2 };
    MouseMode mouseMode_; //!< Режим обработки движения мыши
    ColorValueMap map_; //!< Карта отображеня значения из диапазона на цветовую схему
    /**
     * @brief Список режимов визуализации:
     * @value USER_COLOR цветовая визуализация на основе цветов, указанных пользователем
     * @value NODE_VALUE цветовая вмзуализация на основе значений в узлах сетки и цветовых схем
     * @value ELEMENT_VALUE цветовая вмзуализация на основе значений на элементах и цветовых схем
     */
    enum VisualizationMode { USER_COLOR = 0, NODE_VALUE = 1, ELEMENT_VALUE = 2 };
    VisualizationMode visualizationMode_; //!< Индикатор включения цветовой визуализации значений, определенных на элементе
//    ColorValueMap::ColorMapName colorMapName_; //!< Используемая цветовая схема
    bool isMousePressed;
private:
    /**
     * @brief Установить настройки изображения по умолчанию
     */
    void setDefault ();
    /**
     * @brief Нормализация угла обзора
     * @param angle Значение угла длянормализации
     * @return Эквивалентное значение угла в диапазоне от 0 до 360 градусов
     */
    int normalizedViewAngle (const int &angle);
    /**
     * @brief Нарисовать границы области
     */
    void drawRegionBorders ();
    /**
     * @brief Сбросить видовую матрицу
     */
    void resetProjectionMatrix ();
    /**
     * @brief Преобразование точки сетки в точку openGL
     * @param point Указатель на точку сетки
     */
    void pointToGLVertex(const msh::PointPointer &point) const;
};

#endif // GLMESHPICTURE_H
