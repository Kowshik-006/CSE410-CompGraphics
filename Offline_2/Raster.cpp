#include<bits/stdc++.h>
#include "bitmap_image.hpp"

#define deg2rad M_PI/180.0

#define epsilon 1e-7f

using namespace std;

random_device rd;
mt19937 gen(rd());

int getRandomInt(int min, int max) {
    std::uniform_int_distribution<int> distribution(min, max);
    return distribution(gen);
}

class Point{
    public:
    double x,y,z;
    Point(){
        x = 0;
        y = 0;
        z = 0;
    }
    Point(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    // this . other
    double operator*(const Point& other) const {
        return (x * other.x) + (y * other.y) + (z * other.z);
    }
    // this x other
    Point operator^(const Point& other) const {
        return Point(
            (y * other.z) - (other.y * z),
            (other.x * z) - (x * other.z),
            (x * other.y) - (other.x * y)
        );
    }
    // this * scalar
    Point operator*(double scalar) const {
        return Point(x * scalar, y * scalar, z * scalar);
    }
    // this + other
    Point operator+(const Point& other) const {
        return Point(x + other.x, y + other.y, z + other.z);
    }

    Point operator-(const Point& other) const {
        return Point(x - other.x, y - other.y, z - other.z);
    }

    void normalize(){
        double length = sqrt((x * x) + (y * y) + (z * z));
        if(length == 0){
            return;
        }
        x /= length;
        y /= length;
        z /= length;
    }
};

// Vector and Point are represented in the same way
typedef Point Vector;

class Color{
    public:
    int r, g, b;
    Color(int r, int g, int b){
        this->r = r;
        this->g = g;
        this->b = b;
    }
    Color(){
        // Default color is black
        r = 0;
        g = 0;
        b = 0;
    }
};

Color palette[] = {
    Color(255, 0, 0),    // Red
    Color(0, 255, 0),    // Green
    Color(0, 0, 255),    // Blue
    Color(255, 255, 0),  // Yellow
    Color(255, 0, 255),  // Magenta
    Color(255, 128, 0),  // Orange
    Color(128, 0, 255),  // Violet
    Color(255, 255, 255), // White
    Color(0, 255, 255),  // Cyan
    Color(255, 192, 203), // Pink
    Color(153, 102, 51), // Brown
    Color(102, 0, 153)   // Purple
};

int palette_size = sizeof(palette) / sizeof(palette[0]);
class Triangle{
    public:
    Point p1, p2, p3;
    Color color;
    Triangle(Point p1, Point p2, Point p3){
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
        int index = getRandomInt(0, palette_size - 1);
        this->color = palette[index];
    }
};

vector<Triangle> triangles;

class Matrix{
    public:
    int rows;
    int cols;
    bool push_matrix = false;
    vector<vector<double>> m;
    Matrix(int rows, int cols){
        this->rows = rows;
        this->cols = cols;
        m.resize(rows, vector<double>(cols, 0.0));
    }
    Matrix(vector<vector<double>> m){
        this->rows = m.size();
        this->cols = m[0].size();
        m.resize(rows, vector<double>(cols, 0.0));
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                this->m[i][j] = m[i][j];
            }
        }
    }
    double get(int i, int j){
        return m[i][j];
    }
    void set(int i, int j, double value){
        m[i][j] = value;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw invalid_argument("Matrix dimensions do not match for multiplication");
        }
        Matrix result = Matrix(rows, other.cols);
        for (int i = 0; i < result.rows; i++) {
            for (int j = 0; j < result.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    result.m[i][j] += m[i][k] * other.m[k][j];
                }
            }
        }
        return result;
    } 
};

Matrix pointToMatrix(Point& p){
    Matrix m(4,1);
    m.set(0,0,p.x);
    m.set(1,0,p.y);
    m.set(2,0,p.z);
    m.set(3,0,1.0);
    return m;
}

// Global variables
Point eye, look;
Vector up;
double fovY, aspectRatio, near, far;
int screenWidth, screenHeight;
double screenLeft, screenBottom;
double zFront, zRear;
double dx, dy, left_x, right_x, top_y, bottom_y;
double **zBuffer;
Color **frameBuffer;

Matrix getIdentityMatrix(){
    Matrix identity(4,4);
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            if(i == j){
                identity.set(i,j,1);
            }
            else{
                identity.set(i,j,0);
            }
        }
    }
    return identity;
}

void basic_setup(ifstream& scene_file){
    scene_file  >> eye.x >> eye.y >> eye.z
                >> look.x >> look.y >> look.z
                >> up.x >> up.y >> up.z
                >> fovY >> aspectRatio >> near >> far;

}

Matrix getTranslationMatrix(double tx, double ty, double tz){
    Matrix translation(4,4);
    translation.set(0,0,1);
    translation.set(1,1,1);
    translation.set(2,2,1);
    translation.set(3,3,1);
    translation.set(0,3,tx);
    translation.set(1,3,ty);
    translation.set(2,3,tz);
    return translation;
}

Matrix getScalingMatrix(double sx, double sy, double sz){
    Matrix scaling(4,4);
    scaling.set(0,0,sx);
    scaling.set(1,1,sy);
    scaling.set(2,2,sz);
    scaling.set(3,3,1);
    return scaling;
}

Vector rotateUsingRodrigues(Vector p, double angle, Vector& axis){
    // point * point -> dot product
    // point ^ point -> cross product
    // point * scalar -> scalar multiplication
    // point + point -> vector addition
    angle = angle * deg2rad;
    Vector first_term = p * cos(angle);
    Vector second_term = axis * (axis * p) * (1 - cos(angle));
    Vector third_term = (axis ^ p) * sin(angle);
    return first_term + second_term + third_term;
}

Matrix getRotationMatrix(double angle, Vector& axis){
    Vector c1 = rotateUsingRodrigues(Vector(1,0,0), angle, axis);
    Vector c2 = rotateUsingRodrigues(Vector(0,1,0), angle, axis);
    Vector c3 = rotateUsingRodrigues(Vector(0,0,1), angle, axis);

    Matrix rotation(4,4);
    rotation.set(0,0,c1.x);
    rotation.set(0,1,c2.x);
    rotation.set(0,2,c3.x);
    rotation.set(1,0,c1.y);
    rotation.set(1,1,c2.y);
    rotation.set(1,2,c3.y);
    rotation.set(2,0,c1.z);
    rotation.set(2,1,c2.z);
    rotation.set(2,2,c3.z);
    rotation.set(3,3,1);
    
    return rotation;
}


void perform_transformations(ifstream& scene_file, ofstream& stage1_file){
    stack<Matrix> matrix_stack;
    matrix_stack.push(getIdentityMatrix());
    stage1_file << fixed << setprecision(7);
    while (true){
        string command;
        getline(scene_file, command);
        stringstream ss(command);
        ss >> command;
        if(command == "triangle"){
            // cout<< "Reading triangle" << endl;
            for(int i=0;i<3;i++){
                double x,y,z;
                scene_file >> x >> y >> z;
                Point p = Point(x,y,z);
                Matrix transformed_p = matrix_stack.top() * pointToMatrix(p);
                p.x = transformed_p.get(0,0);
                p.y = transformed_p.get(1,0);
                p.z = transformed_p.get(2,0);
                stage1_file << p.x << " " << p.y << " " << p.z << endl;
            }
            stage1_file << endl;
        }
        else if (command == "translate"){
            double tx, ty, tz;
            scene_file >> tx >> ty >> tz;

            Matrix translation_matrix = getTranslationMatrix(tx, ty, tz);
            matrix_stack.push(matrix_stack.top() * translation_matrix);
        }
        else if (command == "scale"){
            double sx, sy, sz;
            scene_file >> sx >> sy >> sz;

            Matrix scaling_matrix = getScalingMatrix(sx, sy, sz);
            matrix_stack.push(matrix_stack.top() * scaling_matrix);
        }
        else if (command == "rotate"){
            double angle, ax, ay, az;
            scene_file >> angle >> ax >> ay >> az;

            Vector axis = Vector(ax, ay, az);
            axis.normalize();

            Matrix rotation_matrix = getRotationMatrix(angle, axis);
            matrix_stack.push(matrix_stack.top() * rotation_matrix);
        }
        else if (command == "push"){
            if(!matrix_stack.empty()){
                // Multiplying the top matrix with identity matrix changes the object while keeping the same values
                Matrix top_matrix = matrix_stack.top() * getIdentityMatrix();
                top_matrix.push_matrix = true;
                matrix_stack.push(top_matrix);
            }
            else{
                cerr << "Error: Stack underflow. No matrix to push." << endl;
            }
        }
        else if (command == "pop"){
            if(matrix_stack.empty()){
                cerr << "Error: Stack underflow. No matrix to pop." << endl;
                break;
            }
            while(!matrix_stack.empty()){
                if(matrix_stack.top().push_matrix){
                    matrix_stack.pop();
                    break;
                }
                else{
                    matrix_stack.pop();
                }
            }
        }
        else if(command == "end"){
            break;
        }

    }
    
}

Matrix getViewTransformationMatrix(){
    Vector l = look - eye;
    l.normalize();
    Vector r = l ^ up;
    r.normalize();
    Vector u = r ^ l;

    u.normalize();
    Matrix T(4,4);
    for(int i=0; i<4; i++){
        T.set(i,i,1);
    }
    T.set(0,3,-eye.x);
    T.set(1,3,-eye.y);
    T.set(2,3,-eye.z);

    Matrix R(4,4);
    for(int i=0; i<3; i++){
        Vector v = i == 0 ? r : (i == 1 ? u : l*(-1));
        for(int j=0; j<3; j++){
            double value = j == 0 ? v.x : (j == 1 ? v.y : v.z);
            R.set(i,j,value);
        }
    }
    R.set(3,3,1);

    Matrix V = R * T;
    return V;
}

Matrix getProjectionTransformationMatrix(){
    double fovX = fovY * aspectRatio;
    double t = near * tan(fovY * deg2rad / 2);
    double r = near * tan(fovX * deg2rad / 2);

    Matrix P(4,4);
    P.set(0,0,near/r);
    P.set(1,1,near/t);
    P.set(2,2,-(far + near)/(far - near));
    P.set(2,3,-(2 * far * near)/(far - near));
    P.set(3,2,-1);

    return P;
}

void transformPoints(Matrix& transformMatrix, ifstream& input_file, ofstream& output_file){
    output_file << fixed << setprecision(7);
    while(true){
        for(int i=0; i<3; i++){
            Point p;
            input_file >> p.x >> p.y >> p.z;
            if(input_file.eof()){
                return;
            }
            Matrix transformed_p = transformMatrix * pointToMatrix(p);
            p.x = transformed_p.get(0,0)/transformed_p.get(3,0);
            p.y = transformed_p.get(1,0)/transformed_p.get(3,0);
            p.z = transformed_p.get(2,0)/transformed_p.get(3,0);
    
            output_file << p.x << " " << p.y << " " << p.z << endl;
        }
        output_file << endl;
    }
}

void config_setup(ifstream& config_file){
    config_file >> screenWidth >> screenHeight
                >> screenLeft
                >> screenBottom
                >> zFront
                >> zRear;

    dx = -2*screenLeft / screenWidth;
    dy = -2*screenBottom / screenHeight;
    left_x = screenLeft + (dx/2);
    right_x = -left_x;
    top_y = -(screenBottom) - (dy/2);
    bottom_y = -top_y;

    cout<< "dx = " << dx << endl;
    cout<< "dy = " << dy << endl;
    cout<< "left_x = " << left_x << endl;
    cout<< "right_x = " << right_x << endl;
    cout<< "top_y = " << top_y << endl;
    cout<< "bottom_y = " << bottom_y << endl;

    zBuffer = new double*[screenHeight];
    for(int i=0; i<screenHeight; i++){
        zBuffer[i] = new double[screenWidth];
        for(int j=0; j<screenWidth; j++){
            zBuffer[i][j] = zRear;  
        }
    }

    frameBuffer = new Color*[screenHeight];
    for(int i=0; i<screenHeight; i++){
        frameBuffer[i] = new Color[screenWidth];
        for(int j=0; j<screenWidth; j++){
            frameBuffer[i][j] = Color(); // Initialize with black color(default)
        }
    }
}


void readTriangles(ifstream& stage3_file){
    triangles.clear();
    while(true){
        double x1, y1, z1, x2, y2, z2, x3, y3, z3;
        stage3_file >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
        if(stage3_file.eof()){
            break;
        }
        Point p1(x1,y1,z1);
        Point p2(x2,y2,z2);
        Point p3(x3,y3,z3);
        triangles.push_back(Triangle(p1,p2,p3));
    }
    stage3_file.close();
}

double getArea(double x1, double y1, double x2, double y2, double x3, double y3){
    return fabs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)) / 2.0);
}

bool isPointInsideTriangle(double p_x, double p_y, const Triangle& t){
    double area_ABC = getArea(t.p1.x, t.p1.y, t.p2.x, t.p2.y, t.p3.x, t.p3.y);
    double area_PAB = getArea(p_x, p_y, t.p1.x, t.p1.y, t.p2.x, t.p2.y);
    double area_PBC = getArea(p_x, p_y, t.p2.x, t.p2.y, t.p3.x, t.p3.y);
    double area_PCA = getArea(p_x, p_y, t.p3.x, t.p3.y, t.p1.x, t.p1.y);

    return fabs(area_ABC - (area_PAB + area_PBC + area_PCA)) < epsilon;
}

void constructBuffers(){
    for(const Triangle& t : triangles){
        double min_x = min({t.p1.x, t.p2.x, t.p3.x});
        double max_x = max({t.p1.x, t.p2.x, t.p3.x});
        double min_y = min({t.p1.y, t.p2.y, t.p3.y});
        double max_y = max({t.p1.y, t.p2.y, t.p3.y});

        int start_x = max(0, (int)((min_x - left_x) / dx));
        int end_x = min((screenWidth - 1), (int)round((max_x - left_x) / dx));

        int start_y = max(0, (int)((top_y - max_y) / dy));
        int end_y = min((screenHeight - 1), (int)round((top_y - min_y) / dy));

        for(int row = start_y; row <= end_y; row++){
            for(int col = start_x; col <= end_x; col++){
                double p_x = left_x + col * dx;
                double p_y = top_y - row * dy;
                if(isPointInsideTriangle(p_x, p_y, t)){
                    double z = (t.p1.z + t.p2.z + t.p3.z) / 3.0; // Average Z value for simplicity
                    if(z < zBuffer[row][col] && z >= zFront){
                        zBuffer[row][col] = z;
                        frameBuffer[row][col] = t.color; // Set the color of the triangle
                    }
                }
            }
        }
    }
}

void drawImage(const string& filepath){
    bitmap_image image(screenWidth,screenHeight);
    for(int row = 0; row < screenHeight; row++){
        for(int col = 0; col < screenWidth; col++){
            // col -> x coord; row -> y coord
            Color c = frameBuffer[row][col];
            image.set_pixel(col, row, c.r, c.g, c.b);
        }
    }
    image.save_image(filepath);
    cout << "Image saved to " << filepath << endl;

    // Free memory
    for(int i=0; i< screenHeight; i++){
        delete[] frameBuffer[i];
    }
    delete[] frameBuffer;
}

void writeZBufferToFile(const string& filepath){
    ofstream zBuffer_file(filepath);
    if(!zBuffer_file.is_open()){
        cerr << "Unable to open zBuffer file for writing" << endl;
        return;
    }
    zBuffer_file << fixed << setprecision(6);
    for(int row = 0; row < screenHeight; row++){
        for(int col = 0; col < screenWidth; col++){
            if(zBuffer[row][col] != zRear){
                zBuffer_file << zBuffer[row][col] << "\t";
            }
        }
        zBuffer_file << endl;
    }
    zBuffer_file.close();
    cout << "Z-buffer written to " << filepath << endl;

    for(int i=0; i<screenHeight; i++){
        delete[] zBuffer[i];
    }
    delete[] zBuffer;
}

int main(int argc, char** argv){
    int choice;
    cout << "Enter test case number(1-4): ";
    cin >> choice;
    
    while(choice < 1 || choice >4){
        cout << "Invalid choice. Please enter a number between 1 and 4: ";
        cin >> choice;
    }

    string scene_file_path = "./test_cases/" + to_string(choice) + "/scene.txt";
    string config_file_path = "./test_cases/" + to_string(choice) + "/config.txt";
    
    filesystem::create_directories("./output/"+to_string(choice));
    
    string stage1_file_path = "./output/" + to_string(choice) + "/stage1.txt";
    string stage2_file_path = "./output/" + to_string(choice) + "/stage2.txt";
    string stage3_file_path = "./output/" + to_string(choice) + "/stage3.txt";

    string zBuffer_file_path = "./output/" + to_string(choice) + "/zBuffer.txt";
    string image_file_path = "./output/" + to_string(choice) + "/out.bmp";

    ifstream scene_file(scene_file_path);
    ofstream stage1_file_output(stage1_file_path);
    
    if(!scene_file.is_open() || !stage1_file_output.is_open()){   
        cerr << "Unable to open required files for stage 1" << endl;
        exit(1);
    }
    
    // Stage 1

    basic_setup(scene_file);
    perform_transformations(scene_file, stage1_file_output);
    
    stage1_file_output.close();
    scene_file.close();
    
    cout << "Stage 1 completed. Output written to "<< stage1_file_path << endl;

    // Stage 2

    ifstream stage1_file_input(stage1_file_path);
    ofstream stage2_file_output(stage2_file_path);

    if(!stage1_file_input.is_open() || !stage2_file_output.is_open()){   
        cerr << "Unable to open required files for stage 2" << endl;
        exit(1);
    }

    Matrix V = getViewTransformationMatrix();

    transformPoints(V, stage1_file_input, stage2_file_output);
    
    stage1_file_input.close();
    stage2_file_output.close();
    
    cout << "Stage 2 completed. Output written to "<< stage2_file_path << endl;

    // Stage 3

    ifstream stage2_file_input(stage2_file_path);
    ofstream stage3_file_output(stage3_file_path);

    if(!stage2_file_input.is_open() || !stage3_file_output.is_open()){   
        cerr << "Unable to open required files for stage 3" << endl;
        exit(1);
    }

    Matrix P = getProjectionTransformationMatrix();

    transformPoints(P, stage2_file_input, stage3_file_output);
    stage2_file_input.close();
    stage3_file_output.close();
    cout << "Stage 3 completed. Output written to "<< stage3_file_path << endl;


    // Stage 4

    ifstream config_file(config_file_path);
    if(!config_file.is_open()){
        cerr << "Unable to open config file" << endl;
        exit(1);
    }
    config_setup(config_file);
    config_file.close();

    ifstream stage3_file_input(stage3_file_path);
    if(!stage3_file_input.is_open()){
        cerr << "Unable to open stage 3 file" << endl;
        exit(1);
    }
    readTriangles(stage3_file_input);

    constructBuffers();
    drawImage(image_file_path);
    writeZBufferToFile(zBuffer_file_path);
    cout << "Stage 4 completed." << endl;

    return 0;
}