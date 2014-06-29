#include "path2deditor.h"
#include "ui_path2deditor.h"

Path2DEditor::Path2DEditor(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Path2DEditor)
{
    ui->setupUi(this);
    ui->pathTree->addAction(ui->actionAddNode);
    QTreeWidgetItemIterator it(ui->pathTree, QTreeWidgetItemIterator::All);
    ui->pathTree->setCurrentItem((*it));
}

Path2DEditor::~Path2DEditor()
{
    delete ui;
}

void Path2DEditor::on_addNodeToolButton_clicked()
{

}
