/**
  * Заголовочный файл класса главного окна системы
  * @author Чопоров Сергей
  * @date 30/12/2013
  * @version 1.0
  */
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QTreeWidget>
#include "mesh.h"

namespace Ui {
class MainWindow;
}
/**
 * @brief Класс главного окна
 */
class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    /**
     * @brief Конструктор
     * @param parent Указатель на родительский виджет (объект)
     */
    explicit MainWindow(QWidget *parent = 0);
    /**
      * @brief Деструктор
      */
    ~MainWindow();
    
private slots:
    void on_actionUnion_triggered();

    void on_actionIntersection_triggered();

    void on_actionDifference_triggered();

    void on_actionPath_triggered();

    void on_actionModel_triggered();

    void on_actionDeleteObject_triggered();

    void on_actionStructQuads_triggered();

    void on_actionStructIsoQuads_triggered();

    void on_actionBaryQuads_triggered();

    void on_actionPolygonalModel_triggered();

    void on_actionStructHex_triggered();

    void on_actionSaveMesh_triggered();

    void on_actionRotationBodyMesh_triggered();

    void on_actionFlipVertically_triggered();

    void on_actionFlipHorizontally_triggered();

    void on_actionMirrorVertically_triggered();

    void on_actionMirrorrHorizontally_triggered();

private:
    Ui::MainWindow *ui; //!< Контейнер элементов графического интерфейса
private:
    /**
     * @brief addObjectTreeItem Добавить элемент в дерево на форме
     * @param current Указатель на объект-контейнер
     * @param id Идентификатор элемента дерева
     * @param category Категория
     * @param type Тип
     * @param iconName Имя иконки в файле ресурсов
     * @param parameters Сводка параметров
     */
    void addObjectTreeItem(QTreeWidgetItem * current, const int &id, const QString& category, const QString& type, const QString &iconName, QString parameters = "");
    /**
     * @brief addElement Добавить элемент модели (произодится проверка возможности вставки)
     * @param id Идентификатор
     * @param category Категория
     * @param type Тип
     * @param iconName Имя иконки в файле ресурсов
     * @param parameters Параметры
     * @see addObjectTreeItem;
     */
    void addElement(const int &id, const QString &category, const QString &type, const QString &iconName, QString parameters = "");
};

#endif // MAINWINDOW_H
