#include "ScreenshotInterface.h"
#include "Interface/CommonInterface.h"

ScreenshotInterface::ScreenshotInterface(QWidget *parent)
    : ActionInterface("screenshot", "Screenshot", "view", "Take screenshots of maps", "", parent)
{}

void ScreenshotInterface::display(const Vector3 &camPos)
{

}

void ScreenshotInterface::replay(nlohmann::json action)
{

}

QLayout *ScreenshotInterface::createGUI()
{
    return new QGridLayout;
}

void ScreenshotInterface::takeScreenshots()
{
    bool multiFiles = true;

    if (multiFiles) {
        QStringList fileNames = QFileDialog::getOpenFileNames(
            nullptr,            // Parent widget
            "Select one or more files",  // Dialog title
            QDir::homePath(),   // Initial directory
            "Images (*.png *.xpm *.jpg);;Text files (*.txt);;All files (*.*)" // File types
            );

        // Print selected file paths
        for (const QString &fileName : fileNames) {
            qDebug() << "Selected file:" << fileName;
        }
    } else {
        QFileDialog dialog;
        dialog.setFileMode(QFileDialog::Directory);
        dialog.setOption(QFileDialog::DontUseNativeDialog, true);
        dialog.setDirectory(QDir::homePath());
        dialog.setWindowTitle("Select Directories");

        // Find QListView and QTreeView inside QFileDialog
        QListView *listView = dialog.findChild<QListView*>("listView");
        if (listView) listView->setSelectionMode(QAbstractItemView::MultiSelection);
        QTreeView *treeView = dialog.findChild<QTreeView*>();
        if (treeView) treeView->setSelectionMode(QAbstractItemView::MultiSelection);

        QStringList directoryPaths;
        if (dialog.exec()) {
            directoryPaths = dialog.selectedFiles();
        }

        // Print selected directory paths
        for (const QString &path : directoryPaths) {
            qDebug() << "Selected directory:" << path;
        }
    }
}
