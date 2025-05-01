#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<GL/glut.h>

#define pi 3.1416

enum rotation_3d {yaw, pitch_look_at, pitch_up, roll};

// Camera position and orientation
double eye_x = 4, eye_y = 4, eye_z = 4;          // Camera position point coordinates
double center_x = 0, center_y = 0, center_z = 0; // Look-at point coordinates
double up_x = 0, up_y = 1, up_z = 0;             // Up vector coordinates

typedef struct{
    double x, y, z;
}vector_3d;

vector_3d add(vector_3d a, vector_3d b, vector_3d c){
    vector_3d result = {
        a.x + b.x + c.x,
        a.y + b.y + c.y,
        a.z + b.z + c.z
    };
    return result;
}

vector_3d multiply(vector_3d a, double scalar){
    vector_3d result = {
        a.x * scalar,
        a.y * scalar,
        a.z * scalar
    };
    return result;
}

double dot(vector_3d a, vector_3d b){
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

vector_3d cross(vector_3d a, vector_3d b){
    vector_3d result = {
        a.y * b.z - b.y * a.z,
        b.x * a.z - a.x * b.z,
        a.x * b.y - b.x * a.y
    };

    return result;
}

vector_3d normalize(vector_3d a){
    double magnitude = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    vector_3d result = {
        a.x / magnitude,
        a.y / magnitude,
        a.z / magnitude
    };
    return result;
}

vector_3d rotate(vector_3d vector, vector_3d axis, double angle){
    double radian_angle = angle * pi /180.0;
    /*
      Rodrigues' rotation formula
      ----------------------------
      v_rot = v cos(theta) + (k x v) sin(theta) + k (k . v) (1- cos(theta)) 

      where, k is the unit vector along the axis of rotation
             v is the vector to be rotated
             theta is the angle of rotation in radians
             v_rot is the rotated vector
    */

    vector_3d k = normalize(axis);
    vector_3d first_term = multiply(vector, cos(radian_angle));
    vector_3d second_term = multiply(cross(k, vector), sin(radian_angle));
    vector_3d third_term = multiply(k, dot(k, vector) * (1 - cos(radian_angle)));

    vector_3d v_rot = add(first_term, second_term, third_term);
    
    return v_rot;
}

vector_3d getUpVector(){
    vector_3d vector = {up_x, up_y, up_z};
    return vector;
}
vector_3d getLookAtVector(){
    vector_3d vector = {
        center_x - eye_x,
        center_y - eye_y,
        center_z - eye_z
    };
    return vector;
}
vector_3d getRightVector(){
    vector_3d up_vector = getUpVector();
    vector_3d look_at_vector = getLookAtVector();
    vector_3d right_vector = cross(look_at_vector, up_vector);
    return right_vector;
}
vector_3d getRotatedVector(enum rotation_3d type, double angle){
    vector_3d vector;
    vector_3d axis;
    switch(type){
        case yaw:{
            vector = getLookAtVector();
            axis = getUpVector();
            break;
        }
        case pitch_look_at:{
            vector = getLookAtVector();
            axis = getRightVector();
            break;
        }
        case pitch_up:{
            vector = getUpVector();
            axis = getRightVector();
            break;
        }
        case roll:{
            vector = getUpVector();
            axis = getLookAtVector();
            break;
        }
        default:{
            printf("Unknown rotation type\n");
            exit(1);
        }
    }
    vector_3d rotated_vector = rotate(vector, axis, angle);
    return rotated_vector;
}

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
        glVertex3f(1.0, 1.0, -1.0);
        glVertex3f(-1.0, 1.0, -1.0);
        glVertex3f(-1.0, -1.0, -1.0);
        glVertex3f(1.0, -1.0, -1.0);

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

void AlphaNumericKeyListener(unsigned char key, int x, int y){
    double shift_value = 0.1;
    double shift_angle = 1.0;
    switch(key){
        case '1':{
            // look left / yaw left 
            printf("Look Left\n");
            vector_3d rotated_look_at_vector = getRotatedVector(yaw, shift_angle);
            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;

            vector_3d right_vector = getRightVector();
            vector_3d up_vector = cross(right_vector, rotated_look_at_vector);
            up_vector = normalize(up_vector);
            
            up_x = up_vector.x;
            up_y = up_vector.y;
            up_z = up_vector.z;

            break;
        }
        case '2':{
            // look right / yaw right 
            printf("Look Right\n");
            vector_3d rotated_look_at_vector = getRotatedVector(yaw, -shift_angle);
            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;

            vector_3d right_vector = getRightVector();
            vector_3d up_vector = cross(right_vector, rotated_look_at_vector);
            up_vector = normalize(up_vector);
            
            up_x = up_vector.x;
            up_y = up_vector.y;
            up_z = up_vector.z;

            break;
        }
        case '3':{
            // look up / pitch up
            printf("Look Up\n");
            vector_3d rotated_look_at_vector = getRotatedVector(pitch_look_at, shift_angle);
            vector_3d rotated_up_vector = getRotatedVector(pitch_up, shift_angle);
            rotated_up_vector = normalize(rotated_up_vector);

            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;

            up_x = rotated_up_vector.x;
            up_y = rotated_up_vector.y;
            up_z = rotated_up_vector.z;
            
            break;
        }
        case '4':{
            // look down / pitch down
            printf("Look Down\n");
            vector_3d rotated_look_at_vector = getRotatedVector(pitch_look_at, -shift_angle);
            vector_3d rotated_up_vector = getRotatedVector(pitch_up, -shift_angle);
            rotated_up_vector = normalize(rotated_up_vector);

            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;

            up_x = rotated_up_vector.x;
            up_y = rotated_up_vector.y;
            up_z = rotated_up_vector.z;
            
            break;
        }
        case '5':{
            // tilt clockwise / roll clockwise
            printf("Tilt Clockwise\n");
            vector_3d rotated_up_vector = getRotatedVector(roll, -shift_angle);
            rotated_up_vector = normalize(rotated_up_vector);
            up_x = rotated_up_vector.x;
            up_y = rotated_up_vector.y;
            up_z = rotated_up_vector.z;
            break;
        }
        case '6':{
            // tilt counter clockwise / roll counter clockwise
            printf("Tilt Counter Clockwise\n");
            vector_3d rotated_up_vector = getRotatedVector(roll, shift_angle);
            rotated_up_vector = normalize(rotated_up_vector);
            up_x = rotated_up_vector.x;
            up_y = rotated_up_vector.y;
            up_z = rotated_up_vector.z;
            break;
        }
        case 'w':{
            // move upward without changing the reference point
            printf("Move Upward Without Changing the Reference Point\n");
            eye_y += shift_value; 
            break;
        }
        case 's':{
            // move downward without changing the reference point
            printf("Move Downward Without Changing the Reference Point\n");
            eye_y -= shift_value; 
            break;
        }
        default:{
            printf("Unknown key pressed\n");
            break;
        }
    }

    glutPostRedisplay();
}

void SpecialKeyListener(int key, int x, int y){
    double shift_value = 0.1;

    switch(key){
        case GLUT_KEY_UP:{
            // move forward
            printf("Move Forward\n");
            eye_z -= shift_value;
            center_z -= shift_value;
            break;
        }
        case GLUT_KEY_DOWN:{
            // move backward
            printf("Move Backward\n");
            eye_z += shift_value;
            center_z += shift_value;
            break;
        }
        case GLUT_KEY_LEFT:{
            // move left
            printf("Move Left\n");
            eye_x -= shift_value;
            center_x -= shift_value;
            break;
        }
        case GLUT_KEY_RIGHT:{
            // move right
            printf("Move Right\n");
            eye_x += shift_value;
            center_x += shift_value;
            break;
        }

        case GLUT_KEY_PAGE_UP:{
            // move up
            printf("Move Up\n");
            eye_y += shift_value;
            center_y += shift_value;
            break;
        }

        case GLUT_KEY_PAGE_DOWN:{
            // move down
            printf("Move Down\n");
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

    // Although, eye and center are points, up is a vector. 
    // So, look_at_vector = {center_x - eye_x, center_y - eye_y, center_z - eye_z}
    // up_vector = {up_x, up_y, up_z}
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
    glutKeyboardFunc(AlphaNumericKeyListener);
    glutSpecialFunc(SpecialKeyListener);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glEnable(GL_DEPTH_TEST);

    glutMainLoop();
    return 0;
    
}