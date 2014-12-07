/**
  * @author Сергей Чопоров
  * @date 30/12/2013
  * @version 1.0.1
  * @copyright Copyright 2013 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QTreeWidget>
#include "meshpointer.h"
#include "qstdredirector.h"
#include "highlighter.h"

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

    virtual void closeEvent(QCloseEvent *event);
    
private slots:
    void onConsoleMessage(QString message);
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

    void on_actionArea_triggered();

    void on_actionLoadMesh_triggered();

    void on_actionElasticFem_triggered();

    void on_actionLoadNodeValue_triggered();

    void on_actionLoadElementValue_triggered();

    void on_actionExtremeValuesStatistica_triggered();

    void on_actionNewScript_triggered();

    void on_actionOpenScript_triggered();

    void on_actionSaveScript_triggered();

    void on_actionSaveAsScript_triggered();

    void on_actionRunScript_triggered();

private:
    Ui::MainWindow *ui; //!< Контейнер элементов графического интерфейса
    QStdRedirector<> *stdRedirector; //!< Перехватчик сообщений стандартного потока
    Highlighter *highlighter;
    QString scriptFileName;

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
    void clearMesh(msh::MeshPointer mesh);
    /**
     * @brief Установить имя текущего файла со скриптом модели
     * @param fileName Имя файла
     */
    void setCurrentScriptName(const QString &fileName);
    /**
     * @brief Метод для записи скрипта модели в текстовый файл
     * @param fileName Имя файла для записи
     * @return true, если запись успешно выполнена
     */
    bool writeScript(const QString &fileName);
    /**
     * @brief Метод для чтения скрипта модели из текстового файла
     * @param fileName Имя файла для считывания
     * @return true, если файл успешно считан
     */
    bool readScript(const QString &fileName);
    /**
     * @brief Метод для сохранения файла по пути, указаному в диалоге выбора файла
     * @return true, если сохранение файла успешно выполнено
     */
    bool saveScriptAs();
    /**
     * @brief Метод для сохранения файла, который при необходимости предлагает диалог выбора пути
     * @return true, если сохранение файла успешно выполнено
     */
    bool saveScript();
    /**
     * @brief Метод при необходимости запрашивает подтвержедение сохранения изменение в скрипте модели
     * @return true -- Ok, false -- Cancel
     */
    bool maybeSaveScript();
    bool openScript();
};

#endif // MAINWINDOW_H
