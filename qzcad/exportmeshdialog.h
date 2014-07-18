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
     * @brief Флаг учета значений в узле
     * @return true - начение в узле необходимо учесть
     */
    bool isNodeValue();
    /**
     * @brief Флаг учета значений на элементе
     * @return true - начение на элементе необходимо учесть
     */
    bool isElementValue();
private:
    Ui::ExportMeshDialog *ui; //!< Указатель на контейнер элементов формы
};

#endif // EXPORTMESHDIALOG_H
