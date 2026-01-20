#ifndef SPOILER_H
#define SPOILER_H

#include <QObject>
#include <QWidget>
#include <QGridLayout>
#include <QToolButton>
#include <QFrame>
#include <QParallelAnimationGroup>
#include <QScrollArea>

// Taken from https://stackoverflow.com/a/37119983/9863298
class Spoiler : public QWidget {
    Q_OBJECT
private:
    QGridLayout mainLayout;
    QToolButton toggleButton;
    QFrame headerLine;
    QParallelAnimationGroup toggleAnimation;
    QScrollArea contentArea;
    int animationDuration{300};

public Q_SLOTS:
    Spoiler* toggle(const bool opening);
    Spoiler* open();
    Spoiler* close();
public:
    explicit Spoiler(const QString & title = "", const int animationDuration = 300, QWidget *parent = 0);
    Spoiler* setContentLayout(QLayout & contentLayout);

};

#endif // SPOILER_H
