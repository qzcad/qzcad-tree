/**
  * @author Сергей Чопоров
  * @date 30/12/2013
  * @version 1.0.5
  * @copyright Copyright 2013 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef GLMESHPICTURE_H
#define GLMESHPICTURE_H

#include <QtOpenGL>
#include <QWidget>

#include "meshpointer.h"
#include "colorvaluemap.h"
#include "namedfloatingvector.h"

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
     * @return Указатель на актуальную сетку
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
    /**
     * @brief Освободить указатель на сетку (перестать использовать и вернуть его значение, память не очищается)
     * @return Значение указателя на сетку
     */
    msh::MeshPointer releaseMesh();
    /**
     * @brief Добавить вектор значений, определенных в узле
     * @param vector Вектор значений, определенных в узле
     */
    void pushNodeValuesVector(const NamedFloatingVector &vector);
    /**
     * @brief Очистить вектора значений, определенных в узле
     */
    void clearNodeValues();
    /**
     * @brief Добавить вектор значений, определенных на элементе
     * @param vector Вектор значений определенных на элементе
     */
    void pushElementValuesVector(const NamedFloatingVector &vector);
    /**
     * @brief Очистить вектор значений, определенных на элементе
     */
    void clearElementValues();
    /**
     * @brief Получить текущий (выбранный) вектор значений
     * @return Вектор значений, который выбрал пользователь для визуализации
     * Если нет векторов в соответствующих массивах, будет возвращен пустой вектор.
     */
    NamedFloatingVector currentValuesVector();

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
     * @brief Установить значение флага видимости контрольной цветной полосы
     * @param show Значение флага
     */
    void setShowColorBar(bool show);
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
    /**
     * @brief Включить/выключить освещение
     * @param isLighting true - освещение включено; false - освещение выключено
     */
    void setIsLighting(bool isLighting);
    /**
     * @brief Переключить на следующий индекс последовательности значений для визуализации
     */
    void nextValueIndex();
    /**
     * @brief Переключить на предыдущий индекс последовательности значений для визуализации
     */
    void prevValueIndex();
    /**
     * @brief Включить/выключить режим использования вектора (перемещений на основе узловых значений) при отрисовке конструкции
     * @param isUseVector Флаг включения/выключения режима
     */
    void setIsUseVector(bool isUseVector);
    /**
     * @brief Установить значения множителя вектора (перемещений на основе узловых значений) при отрисовке
     * @param vectorScale Значение множителя
     */
    void setVectorScale(double vectorScale);
    /**
     * @brief Включить/выключить отображение исходного каркаса
     * @param isShow true для включения
     */
    void setShowInitialFrames(bool isShow);

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
    bool showColorBar_; //!< Показать или нет контрольную цветную полосу
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
    bool isMousePressed_; //!< Флаг нажатия кнопки мышки
    bool isLighting_; //!< Флаг включения/выключения освещения
    double centerX_; //!< Абсцисса центра сцены
    double centerY_; //!< Ордината центра сцены
    double centerZ_; //!< Аппликата центра сцены
    std::vector<NamedFloatingVector> nodeValues_; //!< Массив векторов значений, определенных в узле
    std::vector<NamedFloatingVector> elementValues_; //!< Массив векторов значений, определенных в узле
    msh::UInteger valueIndex_; //!< Индекс текущего вектора для визуализации
    bool isUseVector_; //!< Применить вектор при визуализации
    double vectorScale_; //!< Множитель вектора
    bool isShowInitialFrames; //!< Фглаг включения/выключения отображения каркаса исходного состояния в режиме использования вектора

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
    void pointToGLVertex(const msh::PointPointer &point, double dx = 0.0, double dy = 0.0, double dz = 0.0) const;
    /**
     * @brief Нарисовать контрольную цветную полосу
     */
    void drawColorBar();
    /**
     * @brief Нарисовать направления осей координат
     */
    void drawAxesDirection();
    /**
     * @brief Нарисовать грань
     * @param face Ссылка на грань
     * @param mode Режим OpenGL для отрисовки
     * @param width Толщина линии (по умолчанию 1.0)
     * @param size Толщина точки (по умолчанию 1.0)
     * @param filter Фильтер типа узлов (если UNDEFINED, то используются все узлы, иначе только указанного типа)
     */
    void drawFace(const msh::UIntegerVector &face, GLenum mode, GLfloat width = 1.0, GLfloat size = 1.0, msh::NodeType filter = msh::UNDEFINED);
};

#endif // GLMESHPICTURE_H
