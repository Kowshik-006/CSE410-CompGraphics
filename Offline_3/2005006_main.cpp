#include "2005006_classes.hpp"

using namespace std;

vector<Object*> objects;
vector<Light*> lights;

string scene_file = "scene.txt";
int recursion_levels;
int image_width;
int image_height;

int window_width = 100;
int window_height = 100;
double view_angle = 60.0; // Field of view angle in degrees

int image_no = 11;


double eye_x = 90, eye_y = -125, eye_z = 90;        // Camera position point coordinates
double center_x = 45, center_y = -25, center_z = 45; // Look-at point coordinates
double up_x = 0, up_y = 0, up_z = 1;             // Up vector coordinates

bool Object :: valid(const Ray& ray,double t){
    if(t < 0){
        return false;
    }
    Vector p = ray.start + ray.direction * t;
    return (((reference_point.x <= p.x && p.x <= reference_point.x + length) || length == 0) &&
            ((reference_point.y <= p.y && p.y <= reference_point.y + width) || width == 0) &&
            ((reference_point.z <= p.z && p.z <= reference_point.z + height) || height == 0));
}
double Object :: intersect(const Ray& ray, const Color& color, int level){
    Vector R_o = ray.start;
    Vector R_d = ray.direction;

    double a = A * R_d.x * R_d.x + B * R_d.y * R_d.y + C * R_d.z * R_d.z 
              + D * R_d.x * R_d.y + E * R_d.y * R_d.z + F * R_d.z * R_d.x;

    double b = 2 * (A * R_o.x * R_d.x + B * R_o.y * R_d.y + C * R_o.z * R_d.z)
              + D * (R_o.x * R_d.y + R_o.y * R_d.x)
              + E * (R_o.y * R_d.z + R_o.z * R_d.y)
              + F * (R_o.z * R_d.x + R_o.x * R_d.z)
              + G * R_d.x + H * R_d.y + I * R_d.z;
    double c = A * R_o.x * R_o.x + B * R_o.y * R_o.y + C * R_o.z * R_o.z
              + D * R_o.x * R_o.y + E * R_o.y * R_o.z + F * R_o.z * R_o.x
              + G * R_o.x + H * R_o.y + I * R_o.z + J;

    double discriminant = b * b - 4 * a * c;
    if(discriminant < 0){
        return -1; // No intersection
    }
    double t1 = (-b - sqrt(discriminant)) / (2 * a);
    double t2 = (-b + sqrt(discriminant)) / (2 * a);
    
    return valid(ray, t1) ? (valid(ray, t2) ? min(t1, t2) : t1 ) : (valid(ray, t2) ? t2 : -1);
}

void Sphere :: draw(){
    int stacks = 32;
    int slices = 32;
    glPushMatrix();
    glTranslated(reference_point.x,reference_point.y,reference_point.z);
    for (int i = 0; i < stacks; i++) {
        glBegin(GL_QUAD_STRIP);
        glColor3d(color.r,color.g,color.b);
        for (int j = 0; j <= slices; j++) {
            double theta = 2.0 * M_PI * j / slices;  // Longitude
            double phi = M_PI * i / stacks;          // Latitude
            
            double x = length * cos(theta) * sin(phi);
            double y = length * sin(theta) * sin(phi);
            double z = length * cos(phi);
            
            glVertex3d(x, y, z);
            
            double phi_next = (i + 1) * M_PI / stacks;
            x = length * cos(theta) * sin(phi_next);
            y = length * sin(theta) * sin(phi_next);
            z = length * cos(phi_next);
            
            glVertex3d(x, y, z);
        }
        glEnd();
    }
    glPopMatrix();
}

double Sphere :: intersect(const Ray& ray, const Color& color, int level){
    // reference_point = center of the sphere
    // length = radius of the sphere
    /*
      P(t) = Ray_start + t*Ray_direction
      Sphere equation: (P - reference_point) * (P - reference_point) - length * length = 0
      t^2 * (Ray_direction . Ray_direction) +
      2 * t * (Ray_start - reference_point) . Ray_direction +
      (Ray_start - reference_point) * (Ray_start - reference_point) - length * length = 0
    */
    Vector oc = ray.start - reference_point;
    double a = ray.direction * ray.direction;
    double b = 2.0 * (oc * ray.direction);
    double c = (oc * oc) - (length * length);
    double discriminant = (b * b) - (4.0 * a * c);
    
    if(discriminant < 0){
        return -1; // No intersection
    }
    
    double t1 = (-b - sqrt(discriminant)) / (2.0 * a);
    double t2 = (-b + sqrt(discriminant)) / (2.0 * a);
    
    if(t1 < 0 && t2 < 0){
        return -1; 
    }
    if(t1 < 0){
        return t2; 
    }
    if(t2 < 0){
        return t1; 
    }
    
    return min(t1, t2); 
}
void Triangle :: draw(){
    glBegin(GL_TRIANGLES);{
        glColor3d(color.r,color.g,color.b);
        glVertex3d(a.x,a.y,a.z);
        glVertex3d(b.x,b.y,b.z);
        glVertex3d(c.x,c.y,c.z);
    }glEnd();
}

double Triangle :: intersect(const Ray& ray, const Color& color, int level){
    Matrix A = Matrix(ray.direction*(-1), b-a, c-a);
    Matrix t_m = Matrix(ray.start-a, b-a, c-a);
    Matrix beta_m = Matrix(ray.direction*(-1), ray.start-a, c-a);
    Matrix gamma_m = Matrix(ray.direction*(-1), b-a, ray.start-a);

    double t = t_m.determinant() / A.determinant();
    double beta = beta_m.determinant() / A.determinant();
    double gamma = gamma_m.determinant() / A.determinant();

    if( beta >= 0 && gamma >= 0 && (beta + gamma) <= 1){
        return t;
    }
    return -1; // No intersection
}
void Floor :: draw(){
    // cout<< "Drawing Floor" << endl;
    int horizontal_tiles = (-2)*reference_point.x / length;
    int vertical_tiles = (-2)*reference_point.y / length;
    // cout << "Horizontal Tiles: " << horizontal_tiles << ", Vertical Tiles: " << vertical_tiles << endl;
    bool black_tile = true;
    for(int i = 0; i < vertical_tiles; i++){
        for(int j = 0; j < horizontal_tiles; j++){
            double bottom_left_x = reference_point.x + j * length;
            double bottom_left_y = reference_point.y + i * length;
            // double bottom_left_z = reference_point.y + i * length;
            glBegin(GL_QUADS);{
                if(black_tile){
                    glColor3d(0.0, 0.0, 0.0);
                }
                else{
                    glColor3d(1.0, 1.0, 1.0);
                }
                black_tile = !black_tile;
                // xy plane
                glVertex3d(bottom_left_x, bottom_left_y, reference_point.z);
                glVertex3d(bottom_left_x + length, bottom_left_y, reference_point.z);
                glVertex3d(bottom_left_x + length, bottom_left_y + length, reference_point.z);
                glVertex3d(bottom_left_x, bottom_left_y + length, reference_point.z);
                
                // xz plane
                // glVertex3d(bottom_left_x, reference_point.y, bottom_left_z);
                // glVertex3d(bottom_left_x + length, reference_point.y, bottom_left_z );
                // glVertex3d(bottom_left_x + length, reference_point.y, bottom_left_z + length);
                // glVertex3d(bottom_left_x, reference_point.y, bottom_left_z + length);
            }glEnd();
        }
        if(horizontal_tiles % 2 == 0){
            black_tile = !black_tile; 
        }
    }

}

void Floor :: setColor(Vector p){
    int row = (p.y - reference_point.y) / length;
    int col = (p.x - reference_point.x) / length;
    if((row + col) % 2 == 0){
        color = Color(0, 0, 0); // Black tile
    }
    else{
        color = Color(1, 1, 1); // White tile
    }
}

double Floor :: intersect(const Ray& ray, const Color& color, int level){
    Vector normal = Vector(0, 0, 1); 
    double numerator = normal * (reference_point - ray.start);
    double denominator = normal * ray.direction;
    if(denominator == 0){
        return -1; // Ray is parallel to the floor
    }
    double t = numerator / denominator;
    if(t < 0){
        return -1; // Intersection is behind the ray start  
    }
    this->setColor(ray.start + ray.direction * t);
    return t; 
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

    for(int i=0; i<objects.size(); i++) {
        if(i == 6 || i == 7){
            continue;
        }
        objects[i]->draw();
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

    // 60-degree field of view, aspect ratio, near and far clipping planes
    gluPerspective(60.0, aspect, 0.1, 1001);
}
void getPosition(){
    cout << "Eye: ("<< eye_x << ", " << eye_y << ", " << eye_z << ")\n";
    cout << "Center: ("<< center_x << ", " << center_y << ", " << center_z << ")\n";
    cout << "Up: ("<< up_x << ", " << up_y << ", " << up_z << ")\n";
}

string getColor(Color color){
    if(color.r == 1 && color.g == 0 && color.b == 0){
        return "Red";
    }
    else if(color.r == 0 && color.g == 1 && color.b == 0){
        return "Green";
    }
    else if(color.r == 0 && color.g == 0 && color.b == 1){
        return "Blue";
    }
    else if(color.r == 1 && color.g == 1 && color.b == 0){
        return "Yellow";
    }
    return "Unknown Color";

}
void capture(){
    bitmap_image image(image_width, image_height);
    Vector l = getLookAtVector().normalize();
    Vector u = getUpVector().normalize();
    Vector r = getRightVector().normalize();

    double plane_distance = (window_height / 2) / tan(view_angle * deg2rad / 2);

    Vector top_left = Vector(
        eye_x + l.x * plane_distance - r.x * (window_width / 2) + u.x * (window_height / 2),
        eye_y + l.y * plane_distance - r.y * (window_width / 2) + u.y * (window_height / 2),
        eye_z + l.z * plane_distance - r.z * (window_width / 2) + u.z * (window_height / 2)
    );

    double du = (double)window_width / image_width;
    double dv = (double)window_height / image_height;

    top_left = top_left + r * (du/2) - u * (dv/2);

    int nearest;
    double t, t_min;

    for(int i=0; i<image_width; i++){
        for(int j=0; j<image_height; j++){
            Vector current_pixel = top_left + r * (i * du) - u * (j * dv);
            Ray ray(Vector(eye_x, eye_y, eye_z), current_pixel - Vector(eye_x, eye_y, eye_z));
            Color color = Color(0, 0, 0);
            // maximum double number
            double t_min = INFINITY;
            Object* o_min = nullptr; 
            for (int i = 0; i < objects.size(); i++) {
                Object* obj = objects[i];
                t = obj->intersect(ray, color, 0);
                if (t < t_min && t >= 0){
                    // cout<< "Found object " << i << endl;
                    // cout<<"Color: (" << obj->color.r << ", " << obj->color.g << ", " << obj->color.b << ")\n";
                    t_min = t;
                    o_min = obj;
                }
            }
            if(o_min != nullptr){
                // cout<<"Found object\n";
                // cout<<"Color: (" << o_min->color.r << ", " << o_min->color.g << ", " << o_min->color.b << ")\n";
                color = o_min->color;
                // cout<<"Color: "<< getColor(color) << endl;
            }
            image.set_pixel(i, j, color.r*255, color.g*255, color.b*255);
        }
    }

    image.save_image("Output_" + to_string(image_no) + ".bmp");
    image_no++;
}

void AlphaNumericKeyListener(unsigned char key, int x, int y){
    double shift_value = 0.1;
    double shift_angle = 1.0;
    switch(key){
        case '0':{
            capture();
            cout<<"Captured scene\n";
            break;
        }
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
        case 'p':{
            // print position
            getPosition();
            break;
        }
        case 'x':{
            cout << "Terminating program\n";
            for(Object* obj : objects){
                delete obj; // Free memory
            }
            for(Light* light : lights){
                delete light; // Free memory
            }
            exit(0);
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
    double shift_value = 0.5;

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
        Object *obj;

        if(object_type == "sphere"){
            Vector center;
            double radius;
            
            input_file >> center.x >> center.y >> center.z
                       >> radius;
                       
            obj = new Sphere(center, radius);
        }
        else if(object_type == "triangle"){
            Vector a, b, c;
            input_file >> a.x >> a.y >> a.z
                       >> b.x >> b.y >> b.z
                       >> c.x >> c.y >> c.z;
            obj = new Triangle(a,b,c);
        }

        else if(object_type == "general"){
            double A,B,C,D,E,F,G,H,I,J;
            Vector ref_point;
            double height, width, length;
            input_file >> A >> B >> C >> D >> E >> F >> G >> H >> I >> J
                       >> ref_point.x >> ref_point.y >> ref_point.z >> length >> width >> height; 
            obj = new Object(ref_point, length, width, height);
            obj->setGeneralCoefficients(A, B, C, D, E, F, G, H, I, J);
        }
        else{
            cerr << "Invalid  object type : " << object_type << endl;
            return;
        }

        input_file >> color.r >> color.g >> color.b
                   >> coeff.ambient >> coeff.diffuse >> coeff.specular >> coeff.reflection
                   >> shine;
        
        obj->setColor(color);
        obj->setCoefficients(coeff);
        obj->setShine(shine);

        objects.push_back(obj);
    }

    Floor* floor = new Floor(1000,20);
    floor->setCoefficients(CoEfficients(0.4, 0.2, 0.2, 0.2));
    floor->setShine(10);
    objects.push_back(floor);

    input_file >> num_point_lights;

    for(int i=0; i<num_point_lights; i++){
        Vector position;
        Color color;
        input_file >> position.x >> position.y >> position.z
                   >> color.r >> color.g >> color.b;
        
        lights.push_back(new PointLight(position, color));
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
                   
        lights.push_back(new SpotLight(position, color, direction, cutoff_angle));
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