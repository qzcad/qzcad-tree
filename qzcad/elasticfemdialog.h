#ifndef ELASTICFEMDIALOG_H
#define ELASTICFEMDIALOG_H

#include <QDialog>

namespace Ui {
class ElasticFemDialog;
}

class ElasticFemDialog : public QDialog
{
    Q_OBJECT

public:
    explicit ElasticFemDialog(QWidget *parent = 0);
    ~ElasticFemDialog();

private:
    Ui::ElasticFemDialog *ui;
};

#endif // ELASTICFEMDIALOG_H
