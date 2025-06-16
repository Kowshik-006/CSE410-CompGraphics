#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<GL/glut.h>
#define pi 3.1416

enum rotation_3d {yaw, pitch, roll};

int simulate = 0;
int showVelocityArrow = 0;

// Camera position and orientation
double eye_x = 4, eye_y = 4, eye_z = 4;          // Camera position point coordinates
double center_x = 0, center_y = 0, center_z = 0; // Look-at point coordinates
double up_x = 0, up_y = 1, up_z = 0;             // Up vector coordinates

typedef struct{
    double x, y, z;
}vector_3d;

typedef struct{
    double radius;
    vector_3d position;
    vector_3d velocity;
    double rotationAngle;
    vector_3d rotationAxis;
}Ball;

Ball ball = {
             0.1,
             {0.0, 0.9, 0.0}, // position 
             {0.3, 0.1, 0.4}, // velocity
             0.0,
             {0.0, 0.0, 1.0}  // rotationAxis (around z-axis)
            };

double gravity = -0.98; // The down direction is negative Y-axis
double restitution = 0.75;
double deltaTime = 0.016; // 60 FPS

// Generate a random number between -1 and 1
double getRandom() {
    return ((double)rand() / RAND_MAX) * 2.0 - 1.0;
}

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

vector_3d rotate(vector_3d v, vector_3d k, double angle){
    /*
        Rodrigues' rotation formula
        ----------------------------
        v_rot = v cos(theta) + (k x v) sin(theta) + k (k . v) (1- cos(theta)) 
        
        where, 
            k     = the unit vector along the axis of rotation
            v     = the vector to be rotated
            theta = the angle of rotation in radians
            v_rot = the rotated vector
    */
   
    double theta = angle * pi /180.0;
    
    vector_3d first_term = multiply(v, cos(theta));
    vector_3d second_term = multiply(cross(k, v), sin(theta));
    vector_3d third_term = multiply(k, dot(k, v) * (1 - cos(theta)));

    vector_3d v_rot = add(first_term, second_term, third_term);
    
    return v_rot;
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
    vector_3d up_vector = {up_x, up_y, up_z};
    vector_3d look_at_vector = getLookAtVector();
    vector_3d right_vector = cross(look_at_vector, up_vector);
    return right_vector;
}

vector_3d getUpVector(){
    // The up vector created with the up coordinates might not be orthogonal with the look_at_vector
    // So, we need to calculate the up vector using the right vector and look_at_vector
    vector_3d look_at_vector = getLookAtVector();
    vector_3d right_vector = getRightVector();
    vector_3d up_vector = cross(right_vector, look_at_vector);
    return up_vector;
}

vector_3d getRotatedVector(enum rotation_3d type, double angle){
    vector_3d vector;
    vector_3d axis;
    switch(type){
        case yaw:{
            vector = getLookAtVector();
            axis = normalize(getUpVector());
            break;
        }
        case pitch:{
            vector = getLookAtVector();
            axis = normalize(getRightVector());
            break;
        }
        case roll:{
            vector = getUpVector();
            axis = normalize(getLookAtVector());
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
void drawBottomFace(){
    int is_white = 1;
    for(double i = -1; i <= 0.9; i += 0.1){
        for(double j = -1; j <= 0.9; j += 0.1){
            if(is_white){
                // white
                glColor3f(1.0, 1.0, 1.0);
            }
            else{
                // black
                glColor3f(0.0, 0.0, 0.0);
            }
            glBegin(GL_QUADS);{
                glVertex3f(i, -1, j);
                glVertex3f(i + 0.1, -1, j);
                glVertex3f(i + 0.1, -1, j + 0.1);
                glVertex3f(i, -1, j + 0.1);
            }glEnd();
            is_white = !is_white;
        }
        is_white = !is_white;
    }
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

    }glEnd();

    // Bottom Face
    drawBottomFace();
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

            break;
        }
        case '2':{
            // look right / yaw right 
            printf("Look Right\n");
            vector_3d rotated_look_at_vector = getRotatedVector(yaw, -shift_angle);
            
            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;

            break;
        }
        case '3':{
            // look up / pitch up
            printf("Look Up\n");
            vector_3d rotated_look_at_vector = getRotatedVector(pitch, shift_angle);
            
            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;
            
            // Updating the up vector since look_at_vector has changed
            vector_3d up_vector = getUpVector();
            up_vector = normalize(up_vector);

            up_x = up_vector.x;
            up_y = up_vector.y;
            up_z = up_vector.z;
            
            break;
        }
        case '4':{
            // look down / pitch down
            printf("Look Down\n");
            vector_3d rotated_look_at_vector = getRotatedVector(pitch, -shift_angle);
            
            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;
            
            // Updating the up vector since look_at_vector has changed
            vector_3d up_vector = getUpVector();
            up_vector = normalize(up_vector);

            up_x = up_vector.x;
            up_y = up_vector.y;
            up_z = up_vector.z;
            
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
            vector_3d up_vector = normalize(getUpVector());
            eye_x += up_vector.x * shift_value;
            eye_y += up_vector.y * shift_value;
            eye_z += up_vector.z * shift_value;

            // Updating the up vector coordinates since look_at_vector has changed
            up_vector = normalize(getUpVector());
            up_x = up_vector.x;
            up_y = up_vector.y;
            up_z = up_vector.z;
            break;
        }
        case 's':{
            // move downward without changing the reference point
            printf("Move Downward Without Changing the Reference Point\n");
            vector_3d up_vector = normalize(getUpVector());
            eye_x -= up_vector.x * shift_value;
            eye_y -= up_vector.y * shift_value;
            eye_z -= up_vector.z * shift_value;

            // Updating the up vector coordinates since look_at_vector has changed
            up_vector = normalize(getUpVector());
            up_x = up_vector.x;
            up_y = up_vector.y;
            up_z = up_vector.z;
            break;
        }
        case ' ':{ 
            if(simulate){
                simulate = 0;
                printf("Simulation Paused\n");
            }
            else{
                simulate = 1;
                printf("Simulation Started\n");
            } 
            break;
        }
        case 'r':{
            if(simulate){
                return;
            }
            printf("The position and velocity of the ball have been reset\n");
            ball.position = (vector_3d){getRandom(), getRandom(), getRandom()};
            ball.velocity = (vector_3d){getRandom(), fabs(getRandom()), getRandom()};
            ball.rotationAngle = 0.0;
            ball.rotationAxis = (vector_3d){0.0, 0.0, 1.0};
            break;
        }
        case '+':{
            if(simulate){
                return;
            }
            if(ball.velocity.y > 0.0){
                ball.velocity.y += 0.02;
            }
            else{
                ball.velocity.y -= 0.02;
            }
            printf("Ball speed increased\n");
            break;
        }
        case '-':{
            if(simulate){
                return;
            }
            if(fabs(ball.velocity.y) < 0.01){
                printf("Cannot decrease the speed further\n");
                return;
            }
            if(ball.velocity.y > 0.0){
                ball.velocity.y -= 0.02;
                printf("Ball speed decreased\n");
            }
            else{
                ball.velocity.y += 0.02;
                printf("Ball velocity decreased\n");
            }
            break;
    }
        case 'v':{ 
            if(showVelocityArrow){
                showVelocityArrow = 0;
                printf("Velocity Arrow Turned Off\n");
            }
            else{
                showVelocityArrow = 1;
                printf("Velocity Arrow Turned On\n");
            }
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
            vector_3d look_at_vector = normalize(getLookAtVector());

            eye_x += look_at_vector.x * shift_value;
            eye_y += look_at_vector.y * shift_value;
            eye_z += look_at_vector.z * shift_value;

            center_x += look_at_vector.x * shift_value;
            center_y += look_at_vector.y * shift_value;
            center_z += look_at_vector.z * shift_value;
            break;
        }
        case GLUT_KEY_DOWN:{
            // move backward
            printf("Move Backward\n");
            vector_3d look_at_vector = normalize(getLookAtVector());

            eye_x -= look_at_vector.x * shift_value;
            eye_y -= look_at_vector.y * shift_value;
            eye_z -= look_at_vector.z * shift_value;

            center_x -= look_at_vector.x * shift_value;
            center_y -= look_at_vector.y * shift_value;
            center_z -= look_at_vector.z * shift_value;
            break;
        }
        case GLUT_KEY_LEFT:{
            // move left
            printf("Move Left\n");
            vector_3d right_vector = normalize(getRightVector());

            eye_x -= right_vector.x * shift_value;
            eye_y -= right_vector.y * shift_value;
            eye_z -= right_vector.z * shift_value;

            center_x -= right_vector.x * shift_value;
            center_y -= right_vector.y * shift_value;
            center_z -= right_vector.z * shift_value;
            break;
        }
        case GLUT_KEY_RIGHT:{
            // move right
            printf("Move Right\n");
            vector_3d right_vector = normalize(getRightVector());

            eye_x += right_vector.x * shift_value;
            eye_y += right_vector.y * shift_value;
            eye_z += right_vector.z * shift_value;

            center_x += right_vector.x * shift_value;
            center_y += right_vector.y * shift_value;
            center_z += right_vector.z * shift_value;
            break;
        }

        case GLUT_KEY_PAGE_UP:{
            // move up
            printf("Move Up\n");
            vector_3d up_vector = normalize(getUpVector());

            eye_x += up_vector.x * shift_value;
            eye_y += up_vector.y * shift_value;
            eye_z += up_vector.z * shift_value;

            center_x += up_vector.x * shift_value;
            center_y += up_vector.y * shift_value;
            center_z += up_vector.z * shift_value;
            break;
        }

        case GLUT_KEY_PAGE_DOWN:{
            // move down
            printf("Move Down\n");
            vector_3d up_vector = normalize(getUpVector());

            eye_x -= up_vector.x * shift_value;
            eye_y -= up_vector.y * shift_value;
            eye_z -= up_vector.z * shift_value;

            center_x -= up_vector.x * shift_value;
            center_y -= up_vector.y * shift_value;
            center_z -= up_vector.z * shift_value;
            break;
        }

        default:{
            printf("Unknown key pressed\n");
            break;
        }

    }

    glutPostRedisplay();
}

void drawBall(double radius, int slices, int stacks) {
    int brown = 0;
    
    for (int i = 0; i < stacks; ++i) {
        // Alternate starting color for each stack
        brown = !brown;
        
        glBegin(GL_QUAD_STRIP);
        for (int j = 0; j <= slices; ++j) {
            if (brown) glColor3f(0.63, 0.32, 0.18);  // Brown
            else glColor3f(1.0, 1.0, 1.0);        // White
            
            double theta = 2.0 * pi * j / slices;  // Longitude
            double phi = pi * i / stacks;          // Latitude
            
            double x = radius * cos(theta) * sin(phi);
            double y = radius * sin(theta) * sin(phi);
            double z = radius * cos(phi);
            
            glVertex3f(x, y, z);
            
            if (brown) glColor3f(1.0, 1.0, 1.0);  // White
            else glColor3f(0.63, 0.32, 0.18);      // Brown
            
            double phi_next = (i + 1) * pi / stacks;
            x = radius * cos(theta) * sin(phi_next);
            y = radius * sin(theta) * sin(phi_next);
            z = radius * cos(phi_next);
            
            glVertex3f(x, y, z);
            
            brown = !brown;
        }
        glEnd();
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
    
    drawCube();

    glPushMatrix();
    glTranslatef(ball.position.x, ball.position.y, ball.position.z);
    glRotatef(ball.rotationAngle, ball.rotationAxis.x, ball.rotationAxis.y, ball.rotationAxis.z);
    // glColor3f(0.65, 0.5, 0.39); // Wooden ball
    // glutSolidSphere(ball.radius, 32, 32);
    drawBall(ball.radius, 32, 32); // Patterned sphere
    glPopMatrix();

    if (showVelocityArrow) {
        double scale = 0.8;
        vector_3d velocity_arrow = {
            ball.velocity.x * scale,
            ball.velocity.y * scale,
            ball.velocity.z * scale
        };
        glLineWidth(4.0);
        glColor3f(1.0f, 0.0f, 0.0f);  // Red line for velocity
        glBegin(GL_LINES);
        glVertex3f(ball.position.x, ball.position.y, ball.position.z);
        glVertex3f(ball.position.x + velocity_arrow.x, ball.position.y + velocity_arrow.y, ball.position.z + velocity_arrow.z);
        glEnd();
    }

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

void updateBall() {
    if(!simulate){
        return;
    }
    // Integrate
    ball.velocity.y += gravity * deltaTime;
    ball.position.x += ball.velocity.x * deltaTime;
    ball.position.y += ball.velocity.y * deltaTime;
    ball.position.z += ball.velocity.z * deltaTime;

    double minBound = -1.0 + ball.radius;
    double maxBound =  1.0 - ball.radius;

    // Floor bounce
    if (ball.position.y < minBound) {
        ball.position.y = minBound;
        ball.velocity.y = -ball.velocity.y * restitution;
        
        // Ball comes to rest if vertical velocity is too small
        if (fabs(ball.velocity.y) < 0.01){
            ball.velocity.y = 0.0;
        }
    }

    // Ceiling
    if (ball.position.y > maxBound) {
        ball.position.y = maxBound;
        ball.velocity.y = -ball.velocity.y * restitution;
    }

    // Side walls (X)
    if (ball.position.x < minBound) {
        ball.position.x = minBound;
        ball.velocity.x = -ball.velocity.x * restitution;
    }
    if (ball.position.x > maxBound) {
        ball.position.x = maxBound;
        ball.velocity.x = -ball.velocity.x * restitution;
    }

    // Front/back walls (Z)
    if (ball.position.z < minBound) {
        ball.position.z = minBound;
        ball.velocity.z = -ball.velocity.z * restitution;
    }
    if (ball.position.z > maxBound) {
        ball.position.z = maxBound;
        ball.velocity.z = -ball.velocity.z * restitution;
    }

    // Compute horizontal velocity (XZ plane)
    double dx = ball.velocity.x * deltaTime;
    double dz = ball.velocity.z * deltaTime;
    double distance = sqrt(dx * dx + dz * dz);

    if (distance > 0.0001f) {
        // Axis = (-Z, 0, -X) to rotate perpendicular to velocity and Y-up
        ball.rotationAxis.x = -ball.velocity.z;
        ball.rotationAxis.y = 0.0;
        ball.rotationAxis.z = -ball.velocity.x;

        // Normalize axis
        ball.rotationAxis = normalize(ball.rotationAxis);

        // Update angle (convert to degrees)
        double angleDelta = (distance / ball.radius) * (180.0 / pi);
        ball.rotationAngle += angleDelta;
        if (ball.rotationAngle > 360.0){
            ball.rotationAngle -= 360.0;
        }
    }
}

void timer(int) {
    updateBall();
    glutPostRedisplay();
    glutTimerFunc(16, timer, 0); // ~60Hz
}


int main(int argc, char** argv){
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(1000,1000);
    glutInitWindowPosition(1300,500);

    glutCreateWindow("Bouncing Ball");

    // Updating the up coordinates to correspond to the actual up vector
    vector_3d up_vector = normalize(getUpVector());
    up_x = up_vector.x;
    up_y = up_vector.y;
    up_z = up_vector.z;

    srand(555); // Seed for random number generation

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(0, timer, 0);
    glutKeyboardFunc(AlphaNumericKeyListener);
    glutSpecialFunc(SpecialKeyListener);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glEnable(GL_DEPTH_TEST);

    glutMainLoop();
    return 0;
    
}