/**
  * @author Сергей Чопоров
  * @date 02/07/2014
  * @version 1.0.1
  */
#ifndef ROTATIONBODYMESHDIALOG_H
#define ROTATIONBODYMESHDIALOG_H

#include <QDialog>

namespace Ui {
class RotationBodyMeshDialog;
}
/**
 * @brief Диалог для создания дискретных моделей трехмерных тел путем вращения дискретной модели плоского профиля
 */
class RotationBodyMeshDialog : public QDialog
{
    Q_OBJECT

public:
    explicit RotationBodyMeshDialog(QWidget *parent = 0);
    ~RotationBodyMeshDialog();
    /**
     * @brief Индекс выбраной оси координат
     * @return 0 - Ox; 1 - Oy
     */
    int axe();
    /**
     * @brief Радиус тела вращения
     * @return Радиус тела вращения
     */
    double radius();
    /**
     * @brief Количество слоев элементов
     * @return Количество слоев элементов в теле вращения
     */
    int layersCount();
    /**
     * @brief Является ли тело вращения замкнутым
     * @return true - если тело замкнуто; false - иначе
     */
    bool isClosedBody();
    /**
     * @brief Угол поворота профиля (для моделей незамкнутых тел)
     * @return Угол поворота профиля
     */
    double rotationAngle();

private:
    Ui::RotationBodyMeshDialog *ui;
};

#endif // ROTATIONBODYMESHDIALOG_H
