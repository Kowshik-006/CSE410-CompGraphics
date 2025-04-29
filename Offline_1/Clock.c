#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>
#include<time.h> // Include time.h for struct tm
#include<GL/glut.h>

#define pi 3.1416

enum time_unit {hour,minute,second};

// Clock parameters
double clock_radius = 0.5;

double hour_hand_depth = 0.01;
double minute_hand_depth = 0.02;
double second_hand_depth = 0.03;

double hour_hand_length = 0.2;
double minute_hand_length = 0.3;
double second_hand_length = 0.4;

// Half length of the box
double box_size = 0.02;

struct time {
    double hour;
    double minute;
    double second;
};

struct point {
    double x;
    double y;
};

struct time getTime(){
    struct timeval tv;
    gettimeofday(&tv,NULL);
    time_t rawTime = tv.tv_sec;
    struct tm * timeInfo = localtime(&rawTime);

    struct time currentTime;
    currentTime.hour = timeInfo->tm_hour%12;
    currentTime.minute = timeInfo->tm_min;
    currentTime.second = timeInfo->tm_sec + (tv.tv_usec/1000000.0);

    return currentTime;
}

void drawCircle(double radius,int segments)
{
    struct point points[segments+1];
    glColor3f(0.5,0.5,0.5); // medium gray
    //generate points
    for(int i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(2*pi*((double)i/(double)segments));
        points[i].y=radius*sin(2*pi*((double)i/(double)segments));
    }
    //draw segments using generated points
    for(int i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex2f(points[i].x,points[i].y);
			glVertex2f(points[i+1].x,points[i+1].y);
        }
        glEnd();
    }
}

void drawHand(enum time_unit unit, double time){
    double hand_x,hand_y, box_x, box_y;
    double angle;
    double hand_length;
    double depth;

    switch(unit){
        case hour:
            angle = 2*pi*(time/12);
            hand_length = hour_hand_length;
            depth = hour_hand_depth;
            glColor3f(1.0,1.0,1.0); // white
            glLineWidth(4.0);
            break;
        case minute:
            angle = 2*pi*(time/60);
            hand_length = minute_hand_length;
            depth = minute_hand_depth;
            glColor3f(1.0,1.0,1.0); // white
            glLineWidth(3.0);
            break;
        case second:
            angle = 2*pi*(time/60);
            hand_length = second_hand_length;
            depth = second_hand_depth;
            glColor3f(1.0, 0.0, 0.0); // red
            glLineWidth(2.0);
            break;
    }

    hand_x = hand_length*sin(angle);
    hand_y = hand_length*cos(angle);

    box_x = clock_radius*sin(angle);
    box_y = clock_radius*cos(angle);

    glBegin(GL_LINES);
    {
        glVertex3f(0,0,depth);
        glVertex3f(hand_x,hand_y,depth);
    }
    glEnd();

    glBegin(GL_QUADS);
    glVertex3f(box_x - box_size, box_y + box_size, depth);
    glVertex3f(box_x + box_size, box_y + box_size, depth);
    glVertex3f(box_x + box_size, box_y - box_size, depth);
    glVertex3f(box_x - box_size, box_y - box_size, depth);
    glEnd();
}

void drawClockFace(){
    drawCircle(clock_radius,500);

    for(int i = 0; i < 60; i++){
        double angle = 2 * pi * (i/60.0);
        double x = clock_radius * sin(angle);
        double y = clock_radius * cos(angle);

        glBegin(GL_LINES);{
            glColor3f(0.5,0.5,0.5);
            glVertex3f(x*0.95,y*0.95,0.0);
            if(i % 5 == 0){
                glVertex3f(x*0.85,y*0.85,0.0);
            }
            else{
                glVertex3f(x*0.9,y*0.9,0.0);
            }
        }glEnd();
    }
}

void display()
{
    // Clear the color and depth buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Reset transformations
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0); // Adjust as needed
    glMatrixMode(GL_MODELVIEW);

    drawClockFace();
    struct time currentTime = getTime();
    // Draw the hour hand
    drawHand(hour, currentTime.hour + (currentTime.minute/60.0));
    // Draw the minute hand
    drawHand(minute, currentTime.minute + (currentTime.second/60.0));
    // Draw the second hand 
    drawHand(second, currentTime.second);

    // Display the result using double buffering
    glutSwapBuffers();
}

void timerFunction(int value)
{
    // Request a redisplay
    glutPostRedisplay();

    // Register the timer again
    glutTimerFunc(16, timerFunction, 0);
}


int main(int argc,char** argv){
    // Initialize GLUT
    glutInit(&argc, argv);

    // Set initial window size and position
    glutInitWindowSize(800, 800);
    glutInitWindowPosition(1500, 500);

    // Set display mode: RGB color, double buffering, and depth buffering
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    // Create a window with the title
    glutCreateWindow("Analog Clock");

    // Register callback functions
    glutDisplayFunc(display);            // Called when screen needs to be redrawn

    // Use timer function instead of idle function for better control of animation speed
    glutTimerFunc(0, timerFunction, 0);

    // Set the background (clear) color to dark blue
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    // Enable depth testing for proper 3D rendering
    glEnable(GL_DEPTH_TEST);

    glutMainLoop();

    return 0;
}