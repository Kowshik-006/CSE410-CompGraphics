#include "2005006_classes.hpp"
#include<GL/glut.h>
#define deg2rad M_PI/180.0

using namespace std;

vector<Object> objects;
vector<Light> lights;

string scene_file = "scene.txt";
int recursion_levels;
int image_width;
int image_height;

double eye_x = 4, eye_y = 4, eye_z = 4;          // Camera position point coordinates
double center_x = 0, center_y = 0, center_z = 0; // Look-at point coordinates
double up_x = 0, up_y = 1, up_z = 0;             // Up vector coordinates

enum rotation_type {yaw,pitch,roll};

void Sphere :: draw(){
    // lalala
}
void Triangle :: draw(){
    // lalala
}
void Floor :: draw(){
    // lalala
}

Vector getLookAtVector(){
    return Vector(
        center_x - eye_x,
        center_y - eye_y,
        center_z - eye_z
    );
}
Vector getRightVector(){
    return getLookAtVector() ^ Vector(up_x, up_y, up_z);
}

Vector getUpVector(){
    // The up vector created with the up coordinates might not be orthogonal with the look_at_vector
    // So, we need to calculate the up vector using the right vector and look_at_vector
    return getRightVector() ^ getLookAtVector();
}

Vector rotate(Vector v, double angle, Vector& axis){
    angle = angle * deg2rad;
    Vector first_term = v * cos(angle);
    Vector second_term = axis * (axis * v) * (1 - cos(angle));
    Vector third_term = (axis ^ v) * sin(angle);
    return first_term + second_term + third_term;
}

Vector getRotatedVector(enum rotation_type type, double angle){
    Vector vector;
    Vector axis;
    switch(type){
        case yaw:{
            vector = getLookAtVector();
            axis = getUpVector().normalize();
            break;
        }
        case pitch:{
            vector = getLookAtVector();
            axis = getRightVector().normalize();
            break;
        }
        case roll:{
            vector = getUpVector();
            axis = getLookAtVector().normalize();
            break;
        }
        default:{
            printf("Unknown rotation type\n");
            exit(1);
        }
    }
    Vector rotated_vector = rotate(vector, angle, axis);
    return rotated_vector;
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

void AlphaNumericKeyListener(unsigned char key, int x, int y){
    double shift_value = 0.1;
    double shift_angle = 1.0;
    switch(key){
        case '1':{
            // look left / yaw left 
            printf("Look Left\n");
            Vector rotated_look_at_vector = getRotatedVector(yaw, shift_angle);
            
            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;

            break;
        }
        case '2':{
            // look right / yaw right 
            printf("Look Right\n");
            Vector rotated_look_at_vector = getRotatedVector(yaw, -shift_angle);
            
            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;

            break;
        }
        case '3':{
            // look up / pitch up
            printf("Look Up\n");
            Vector rotated_look_at_vector = getRotatedVector(pitch, shift_angle);
            
            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;
            
            // Updating the up vector since look_at_vector has changed
            Vector up_vector = getUpVector().normalize();

            up_x = up_vector.x;
            up_y = up_vector.y;
            up_z = up_vector.z;
            
            break;
        }
        case '4':{
            // look down / pitch down
            printf("Look Down\n");
            Vector rotated_look_at_vector = getRotatedVector(pitch, -shift_angle);
            
            center_x = eye_x + rotated_look_at_vector.x;
            center_y = eye_y + rotated_look_at_vector.y;
            center_z = eye_z + rotated_look_at_vector.z;
            
            // Updating the up vector since look_at_vector has changed
            Vector up_vector = getUpVector().normalize();

            up_x = up_vector.x;
            up_y = up_vector.y;
            up_z = up_vector.z;
            
            break;
        }
        case '5':{
            // tilt clockwise / roll clockwise
            printf("Tilt Clockwise\n");
            Vector rotated_up_vector = getRotatedVector(roll, -shift_angle).normalize();

            up_x = rotated_up_vector.x;
            up_y = rotated_up_vector.y;
            up_z = rotated_up_vector.z;
            break;
        }
        case '6':{
            // tilt counter clockwise / roll counter clockwise
            printf("Tilt Counter Clockwise\n");
            Vector rotated_up_vector = getRotatedVector(roll, shift_angle).normalize();

            up_x = rotated_up_vector.x;
            up_y = rotated_up_vector.y;
            up_z = rotated_up_vector.z;
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
            Vector look_at_vector = getLookAtVector().normalize();

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
            Vector look_at_vector = getLookAtVector().normalize();

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
            Vector right_vector = getRightVector().normalize();

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
            Vector right_vector = getRightVector().normalize();

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
            Vector up_vector = getUpVector().normalize();

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
            Vector up_vector = getUpVector().normalize();

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

void load_data(){
    ifstream input_file(scene_file);
    if (!input_file.is_open()) {
        cerr << "Error opening file: " << scene_file << endl;
        return;
    }

    input_file >> recursion_levels >> image_width;
    image_height = image_width;

    int num_objects, num_point_lights, num_spot_lights;
    input_file >> num_objects;

    for(int i=0; i<num_objects; i++){
        string object_type;
        input_file >> object_type;

        Color color;
        CoEfficients coeff;
        int shine;
        Object obj;

        if(object_type == "sphere"){
            Vector center;
            double radius;
            
            input_file >> center.x >> center.y >> center.z
                       >> radius;
                       
            obj = Sphere(center, radius);
        }
        else if(object_type == "triangle"){
            Vector a, b, c;
            input_file >> a.x >> a.y >> a.z
                       >> b.x >> b.y >> b.z
                       >> c.x >> c.y >> c.z;
            obj = Triangle(a,b,c);
        }

        else if(object_type == "general"){
            double A,B,C,D,E,F,G,H,I,J;
            Vector ref_point;
            double height, width, length;
            input_file >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J
                       >> ref_point.x >> ref_point.y >> ref_point.z >> height >> width >> length; 
        }
        else{
            cerr << "Invalid  object type : " << object_type << endl;
            return;
        }

        input_file >> color.r >> color.g >> color.b
                   >> coeff.ambient >> coeff.diffuse >> coeff.specular >> coeff.reflection
                   >> shine;
        
        obj.setColor(color);
        obj.setCoefficients(coeff);
        obj.setShine(shine);

        objects.push_back(obj);
    }

    input_file >> num_point_lights;

    for(int i=0; i<num_point_lights; i++){
        Vector position;
        Color color;
        input_file >> position.x >> position.y >> position.z
                   >> color.r >> color.g >> color.b;
        
        lights.push_back(PointLight(position, color));
    }

    input_file >> num_spot_lights;

    for(int i=0; i<num_spot_lights; i++){
        Vector position, direction;
        Color color;
        double cutoff_angle;

        input_file >> position.x >> position.y >> position.z
                   >> color.r >> color.g >> color.b
                   >> direction.x >> direction.y >> direction.z
                   >> cutoff_angle;
                   
        lights.push_back(SpotLight(position, color, direction, cutoff_angle));
    }

    input_file.close();
    cout << "Successfully loaded data from " << scene_file << endl;
}

void init_openGL(int* argc, char** argv){
    glutInit(argc, argv);

    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(image_width,image_height);
    glutInitWindowPosition(750,100);

    glutCreateWindow("Ray Tracing");

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(AlphaNumericKeyListener);
    glutSpecialFunc(SpecialKeyListener);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glEnable(GL_DEPTH_TEST);

    glutMainLoop();
}

int main(int argc, char** argv){
    load_data();
    init_openGL(&argc,argv);
    return 0;
}