#include<stdio.h>
#include<stdlib.h>
#include<GL/glut.h>
#include<math.h>

#define pi 3.1416
#define converter (pi/180.0)

// Camera position and orientation
double eye_x = 0, eye_y = 5, eye_z = 7;          // Camera position point coordinates
double center_x = 0, center_y = 0, center_z = 0; // Look-at point coordinates
double up_x = 0, up_y = 1, up_z = 0;             // Up vector coordinates

double rotationAngle = 0.0;
double fan_rotationAngle = 0.0;
double fan_speed = 0.5;
double blade_radius = 2.0;

void drawAxes(){
    glBegin(GL_LINES);{
        glColor3f(1.0,1.0,1.0); //white
        
        // x axis
        glVertex3f(-50.0, 0.0, 0.0);
        glVertex3f(50.0, 0.0, 0.0);

        // y axis
        glVertex3f(0.0, -50.0, 0.0);
        glVertex3f(0.0, 50.0, 0.0);

        // z axis
        glVertex3f(0.0, 0.0, -50.0);
        glVertex3f(0.0, 0.0, 50.0);

    }glEnd();
}

void drawFanHolder(){
    glBegin(GL_QUADS);{
        glColor3f(1.0, 0.0, 1.0); // Purple
        
        glVertex3f(0.1, 2.1, 1.3);
        glVertex3f(0.1, 1.9, 1.3);
        glVertex3f(-0.1, 1.9, 1.3);
        glVertex3f(-0.1, 2.1, 1.3);

        glVertex3f(0.1, 2.1, -1.3);
        glVertex3f(0.1, 1.9, -1.3);
        glVertex3f(-0.1, 1.9, -1.3);
        glVertex3f(-0.1, 2.1, -1.3);

        glVertex3f(0.1, 2.1, 1.3);
        glVertex3f(0.1, 2.1, -1.3);
        glVertex3f(-0.1, 2.1, -1.3);
        glVertex3f(-0.1, 2.1, 1.3);

        glVertex3f(0.1, 2.1, 1.3);
        glVertex3f(0.1, 2.1, -1.3);
        glVertex3f(0.1, 1.9, -1.3);
        glVertex3f(0.1, 1.9, 1.3);

        glVertex3f(0.1, 1.9, 1.3);
        glVertex3f(0.1, 1.9, -1.3);
        glVertex3f(-0.1, 1.9, -1.3);
        glVertex3f(-0.1, 1.9, 1.3);

        glVertex3f(-0.1, 2.1, 1.3);
        glVertex3f(-0.1, 2.1, -1.3);
        glVertex3f(-0.1, 1.9, -1.3);
        glVertex3f(-0.1, 1.9, 1.3);

    }glEnd();

}

void drawFanBlades(){
    glBegin(GL_TRIANGLE_FAN);{
        glColor3f(0.0, 1.0, 1.0); // Cyan
        for(int i=0; i < 3; i++){
            double theta = i * 120.0 + fan_rotationAngle;
            glVertex3f(0.0, 2.0, 1.31);
            glVertex3f(blade_radius * sin(theta*converter), 2.0 + blade_radius * cos(theta*converter), 1.31);
            glVertex3f(blade_radius * sin((theta + 10.0)*converter), 2.0 + blade_radius * cos((theta + 10.0)*converter), 1.31);
        }
    }glEnd();
}


void drawCube(){
    glBegin(GL_QUADS);{
        // Front Face
        glColor3f(0.0, 1.0, 0.0); // Green
        glVertex3f(1.0, 3.0, 1.0);
        glVertex3f(1.0, -3.0, 1.0);
        glVertex3f(-1.0, -3.0, 1.0);
        glVertex3f(-1.0, 3.0, 1.0);

        // Right Face
        glColor3f(0.0, 0.0, 1.0); // Blue
        glVertex3f(1.0, 3.0, 1.0);
        glVertex3f(1.0, 3.0, -1.0);
        glVertex3f(1.0, -3.0, -1.0);
        glVertex3f(1.0, -3.0, 1.0);

        // Back Face
        glColor3f(0.0, 1.0, 0.0); // Green
        glVertex3f(1.0, 3.0, -1.0);
        glVertex3f(-1.0, 3.0, -1.0);
        glVertex3f(-1.0, -3.0, -1.0);
        glVertex3f(1.0, -3.0, -1.0);

        // Left Face
        glColor3f(0.0, 0.0, 1.0); // Blue
        glVertex3f(-1.0, 3.0, -1.0);
        glVertex3f(-1.0, 3.0, 1.0);
        glVertex3f(-1.0, -3.0, 1.0);
        glVertex3f(-1.0, -3.0, -1.0);


    }glEnd();

    drawFanHolder();
    drawFanBlades();
}

void alnumkey(unsigned char key, int x, int y){
    switch(key){
        case 'a':{
            printf("Rotate Left\n");
            rotationAngle -= 5;
            break;
        }
        case 'd':{
            printf("Rotate Right\n");
            rotationAngle += 5;
            break;
        }
        case 'w':{
            printf("Increase Fan Speed\n");
            fan_speed += 0.1;
            break;
        }
        case 's':{
            printf("Decrease Fan Speed\n");
            if(fan_speed >= 0.1)
                fan_speed -= 0.1;
            else
                fan_speed = 0.0; 
            break;
        }
        default:{
            printf("Unknown key pressed\n");
            break;
        }
    }
}

void display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    // Although, eye and center are points, up is a vector. 
    // So, look_at_vector = {center_x - eye_x, center_y - eye_y, center_z - eye_z}
    // up_vector = {up_x, up_y, up_z}
    gluLookAt(
        eye_x,eye_y,eye_z, // Camera position
        center_x,center_y,center_z, // Look-at point
        up_x,up_y,up_z // Up vector
    );
    
    drawAxes();
    glPushMatrix();
    glRotatef(rotationAngle, 0.0, 1.0, 0.0); // Rotate around Y-axis
    drawCube();
    glPopMatrix();


    glutSwapBuffers();
}

void reshape(int width, int height)
{
    double aspect = (double)width / height;
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(60.0f, aspect, 0.1f, 100.0f);
}

void timer(int value){
    glutPostRedisplay();

    if(fan_rotationAngle >= 360.0)
        fan_rotationAngle = 0.0; // Reset the angle
    else
        fan_rotationAngle += fan_speed; // Rotate the fan blades
    
    glutTimerFunc(16, timer, 0); // ~60Hz
}

int main(int argc, char** argv){
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Windmill");

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(alnumkey);
    glutTimerFunc(0, timer, 0);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glEnable(GL_DEPTH_TEST);

    glutMainLoop();

    return 0;
}