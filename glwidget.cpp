#include "glwidget.h"

GLWidget::GLWidget(QWidget *parent) :
    QGLWidget(parent)
{
    m_timer = new QTimer(this);
    connect(m_timer, SIGNAL(timeout()), this, SLOT(timeOutSlot()));
    m_timer->start(50);
}


void GLWidget::initializeGL() {

    //glEnable(GL_TEXTURE_2D);                        // Enable Texture Mapping

    //glShadeModel(GL_SMOOTH);                        // Enables Smooth Shading

    glClearColor(1,1,1,1);           // Black Background
/*
    //depth buffer
    glClearDepth(1.0f);                             // Depth Buffer Setup
    glEnable(GL_DEPTH_TEST);                        // Enables Depth Testing
    glDepthFunc(GL_LEQUAL);                         // The Type Of Depth Test To Do

    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);          // Really Nice Perspective Calculations

    GLfloat LightAmbient[4] = {0.5f, 0.5f, 0.5f, 1.0f};
    GLfloat LightDiffuse[4] = {1.0f, 1.0f, 1.0f, 1.0f};
    GLfloat LightPosition[4] = {0.0f, 0.0f, 5.0f, 1.0f};

    glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);             // Setup The Ambient Light
    glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);             // Setup The Diffuse Light
    glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);            // Position The Light
    glEnable(GL_LIGHT1);                            // Enable Light One
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    */

    sim.initialization(50, 1);

}

void GLWidget::resizeGL(int w, int h) {
    if (h==0)										// Prevent A Divide By Zero By
    {
        h=1;										// Making Height Equal One
    }

    width = w;
    height = h;
    glViewport(0, 0, w, h);

    // projection mode
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity(); // reset

    // Calculate The Aspect Ratio Of The Window
    //gluPerspective(45,(GLfloat)w/(GLfloat)h,0.1f,100);
    // first term is FOV
    // second term is aspect ratio
    // third term is near plane
    // fourth term is far plane

    //glOrtho (0, w, h, 0, 0, 1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}


void draw_circle2d(const Vec2f& centre, float rad, int segs)
{
   glBegin(GL_POLYGON);
   for(int i=0;i<segs;i++){
      float cosine=rad*cos(i*2*M_PI/(float)(segs));
      float sine=rad* sin(i*2*M_PI/(float)(segs));
      Vec2f temp = Vec2f(cosine,sine) + centre;
      glVertex2f(temp[0], temp[1]);
   }
   glEnd();
}

void draw_grid2d(const Vec2f& origin, float dx, int nx, int ny) {
   float width = nx*dx;
   float height = ny*dx;

   glBegin(GL_LINES);
   for(int i = 0; i <= nx; i++) {
      Vec2f a(i*dx, 0);
      Vec2f b(i*dx, height);
      Vec2f temp = origin + a;
      glVertex2f(temp[0], temp[1]);
      temp = origin + b;
      glVertex2f(temp[0], temp[1]);
   }
   for(int j = 0; j <= ny; ++j) {
      Vec2f a(0,j*dx);
      Vec2f b(width,j*dx);
      Vec2f temp = origin + a;
      glVertex2f(temp[0], temp[1]);
      temp = origin + b;
      glVertex2f(temp[0], temp[1]);
   }
   glEnd();
}

void draw_arrow2d(const Vec2f& start, const Vec2f& end, float arrow_head_len)
{
   Vec2f direction = end - start;

   Vec2f dir_norm = direction;

   //TODO Possibly automatically scale arrowhead length based on vector magnitude
   if(dir_norm.norm() < 1e-14)
      return;

   dir_norm.normalize();
   Vec2f perp(dir_norm[1],-dir_norm[0]);

   Vec2f tip_left = end + arrow_head_len/(float)sqrt(2.0)*(-dir_norm + perp);
   Vec2f tip_right = end + arrow_head_len/(float)sqrt(2.0)*(-dir_norm - perp);

   glBegin(GL_LINES);
   glVertex2f(start[0], start[1]);
   glVertex2f(end[0], end[1]);
   glVertex2f(end[0], end[1]);
   glVertex2f(tip_left[0], tip_left[1]);
   glVertex2f(end[0], end[1]);
   glVertex2f(tip_right[0], tip_right[1]);
   glEnd();

}

void GLWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT);

    glLoadIdentity();
    glScalef(2, 2, 0);
    glTranslatef(-0.5, -0.5, 0);

    //grid
    glColor3f(0,0,0);
    glLineWidth(1);
    draw_grid2d(Vec2f(0,0), sim.dx, sim.n, sim.n);

    //particles
    glColor3f(0,0,1);
    for(int i=0; i < sim.particles.size(); i++){
        draw_circle2d(sim.particles[i], sim.particle_radius, 20);
    }

    //velocity
    glColor3f(1,0,0);
    for(int j = 0;j < sim.n; ++j) for(int i = 0; i < sim.n; ++i) {
       Vec2f pos((i+0.5)*sim.dx,(j+0.5)*sim.dx);
       draw_arrow2d(pos, pos + 0.01f*sim.get_velocity(pos), 0.01*sim.dx);
    }

}

void GLWidget::mousePressEvent(QMouseEvent *event) {
    float dx = event->pos().x();
    float dy = event->pos().y();
    lastPos = event->pos();
}
void GLWidget::mouseMoveEvent(QMouseEvent *event) {

}
void GLWidget::mouseReleaseEvent(QMouseEvent *event) {
    sim.computeVelocity(lastPos.x()/width, 1 - lastPos.y()/height,
                        event->pos().x()/width, 1 - event->pos().y()/height);
}

void GLWidget::keyPressEvent( QKeyEvent *e ){
    switch( e->key() )
      {
      /*case Qt::Key_Up:
        xrot += 0.5;

        break;

      case Qt::Key_Down:
        xrot -= 0.5;

        break;
*/
      default:
        GLWidget::keyPressEvent( e );
    }
}

void GLWidget::timeOutSlot(){
    timeOut();
}

void GLWidget::timeOut(){

    sim.step(0.002);

    //qDebug("rawr");

    updateGL();
}
