#include <QApplication>
#include "glwidget.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    GLWidget w;
    w.resize(600,600);
    w.show();
    
    return a.exec();
}
