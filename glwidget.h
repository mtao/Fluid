#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <GL/gl.h>
#include <GL/glu.h>
#include <QtGui>
#include <QtOpenGL/QGLWidget>

#include "simulation.h"

class GLWidget : public QGLWidget
{
    Q_OBJECT
public:
    GLWidget(QWidget *parent = NULL);

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    void timeOut();

protected slots:
    void timeOutSlot();

private:
    QTimer *m_timer;

    Simulation sim;
    QPoint lastPos;

    float width;
    float height;
    
};

#endif // GLWIDGET_H
