#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<GL/glut.h>


// Camera position and orientation
double eye_x = 0, eye_y = 0, eye_z = 4;          // Camera position coordinates
double center_x = 0, center_y = 0, center_z = 0; // Look-at point coordinates
double up_x = 0, up_y = 1, up_z = 0;             // Up vector coordinates

void drawCube(){
    glBegin(GL_QUADS);{
        // Front Face
        glColor3f(1.0, 0.0, 0.0); // Red
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(1.0, -1.0, 1.0);
        glVertex3f(-1.0, -1.0, 1.0);
        glVertex3f(-1.0, 1.0, 1.0);

        // Right Face
        glColor3f(0.0, 0.0, 1.0); // Blue
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, -1.0);
        glVertex3f(1.0, -1.0, -1.0);
        glVertex3f(1.0, -1.0, 1.0);

        // Back Face
        glColor3f(1.0, 1.0, 0.0); // Yellow
        glColor3f(1.0, 1.0, -1.0);
        glColor3f(-1.0, 1.0, -1.0);
        glColor3f(-1.0, -1.0, -1.0);
        glColor3f(1.0, -1.0, -1.0);

        // Left Face
        glColor3f(1.0, 0.0, 1.0);
        glVertex3f(-1.0, 1.0, -1.0);
        glVertex3f(-1.0, 1.0, 1.0);
        glVertex3f(-1.0, -1.0, 1.0);
        glVertex3f(-1.0, -1.0, -1.0);

        // Top Face
        glColor3f(0.0, 1.0, 0.0); // Green
        glVertex3f(-1.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, -1.0);
        glVertex3f(-1.0, 1.0, -1.0);

        // Bottom Face
        glColor3f(1.0, 0.5, 0.0); // Orange
        glVertex3f(1.0, -1.0, 1.0);
        glVertex3f(1.0, -1.0, -1.0);
        glVertex3f(-1.0, -1.0, -1.0);
        glVertex3f(-1.0, -1.0, 1.0);
    }glEnd();
}

void SpecialKeyListener(int key, int x, int y){
    double shift_value = 0.1;

    switch(key){
        case GLUT_KEY_UP:{
            eye_z -= shift_value;
            center_z -= shift_value;
            break;
        }
        case GLUT_KEY_DOWN:{
            eye_z += shift_value;
            center_z += shift_value;
            break;
        }
        case GLUT_KEY_LEFT:{
            eye_x -= shift_value;
            center_x -= shift_value;
            break;
        }
        case GLUT_KEY_RIGHT:{
            eye_x += shift_value;
            center_x += shift_value;
            break;
        }

        case GLUT_KEY_PAGE_UP:{
            eye_y += shift_value;
            center_y += shift_value;
            break;
        }

        case GLUT_KEY_PAGE_DOWN:{
            eye_y -= shift_value;
            center_y -= shift_value;
            break;
        }

        default:{
            printf("Unknown key pressed\n");
            break;
        }

    }

    glutPostRedisplay();
}


void display(){
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    gluLookAt(
        eye_x,eye_y,eye_z, // Camera position
        center_x,center_y,center_z, // Look-at point
        up_x,up_y,up_z // Up vector
    );

    drawCube();

    glutSwapBuffers();
}

void reshape(int width, int height)
{
    double aspect = (double)width / height;

    // Set viewport to cover entire window
    glViewport(0, 0, width, height);

    // Set up perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // 45-degree field of view, aspect ratio, near and far clipping planes
    gluPerspective(60.0f, aspect, 0.1f, 100.0f);
}


int main(int argc, char** argv){
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(1000,1000);
    glutInitWindowPosition(1300,500);

    glutCreateWindow("Bouncing Ball");

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutSpecialFunc(SpecialKeyListener);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glEnable(GL_DEPTH_TEST);

    glutMainLoop();
    return 0;
    
}