/**
  * @author Сергей Чопоров
  * @date 07/12/2014
  * @version 1.0.1
  * @copyright Copyright 2014 Sergey Choporov. All rights reserved.
  * This project is released under the GNU Lesser General Public License.
  */
#ifndef EXPORTMESHDIALOG_H
#define EXPORTMESHDIALOG_H

#include <QDialog>

namespace Ui {
class ExportMeshDialog;
}
/**
 * @brief Диалог параметров экпрота-импорта дискретной модели
 * @see QDialog
 */
class ExportMeshDialog : public QDialog
{
    Q_OBJECT

public:
    /**
     * @brief Конструктор
     * @param parent Родительское окно
     */
    explicit ExportMeshDialog(QWidget *parent = 0);
    ~ExportMeshDialog();
    /**
     * @brief Сохранять информацию о слоях
     * @return true, если сохранять
     */
    bool isLayers();
private:
    Ui::ExportMeshDialog *ui; //!< Указатель на контейнер элементов формы
};

#endif // EXPORTMESHDIALOG_H
