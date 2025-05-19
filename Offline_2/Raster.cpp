#include<bits/stdc++.h>

#define pi 3.1416
#define deg2rad pi/180.0

#define epsilon 1e-4f

using namespace std;

float properValue(float value){
    float abs_value = fabs(value);
    if(abs_value < epsilon){
        return 0;
    }
    float ceil_value = ceil(abs_value);
    float delta = ceil_value - abs_value;
    if(delta < epsilon){
        return value >=0 ? ceil_value : -ceil_value;
    }
    return value;
}

class Point{
    public:
    float x,y,z;
    Point(){
        x = 0;
        y = 0;
        z = 0;
    }
    Point(float x, float y, float z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
    // this . other
    float operator*(const Point& other) const {
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
    Point operator*(float scalar) const {
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
        float length = sqrt((x * x) + (y * y) + (z * z));
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

class Triangle{
    public:
    Point p1, p2, p3;
    Triangle(Point p1, Point p2, Point p3){
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
    }
};

class Matrix{
    public:
    int rows;
    int cols;
    bool push_matrix = false;
    vector<vector<float>> m;
    Matrix(int rows, int cols){
        this->rows = rows;
        this->cols = cols;
        m.resize(rows, vector<float>(cols, 0.0f));
    }
    Matrix(vector<vector<float>> m){
        this->rows = m.size();
        this->cols = m[0].size();
        m.resize(rows, vector<float>(cols, 0.0f));
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                this->m[i][j] = m[i][j];
            }
        }
    }
    float get(int i, int j){
        return m[i][j];
    }
    void set(int i, int j, float value){
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

Point eye, look;
Vector up;
float fovY, aspectRatio, zNear, zFar;


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
                >> fovY >> aspectRatio >> zNear >> zFar;

}

Matrix getTranslationMatrix(float tx, float ty, float tz){
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

Matrix getScalingMatrix(float sx, float sy, float sz){
    Matrix scaling(4,4);
    scaling.set(0,0,sx);
    scaling.set(1,1,sy);
    scaling.set(2,2,sz);
    scaling.set(3,3,1);
    return scaling;
}

Vector rotateUsingRodrigues(Vector p, float angle, Vector& axis){
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

Matrix getRotationMatrix(float angle, Vector& axis){
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
                float x,y,z;
                scene_file >> x >> y >> z;
                Point p = Point(x,y,z);
                Matrix transformed_p = matrix_stack.top() * pointToMatrix(p);
                p.x = properValue(transformed_p.get(0,0));
                p.y = properValue(transformed_p.get(1,0));
                p.z = properValue(transformed_p.get(2,0));
                stage1_file << p.x << " " << p.y << " " << p.z << endl;
            }
            stage1_file << endl;
        }
        else if (command == "translate"){
            float tx, ty, tz;
            scene_file >> tx >> ty >> tz;

            Matrix translation_matrix = getTranslationMatrix(tx, ty, tz);
            matrix_stack.push(matrix_stack.top() * translation_matrix);
        }
        else if (command == "scale"){
            float sx, sy, sz;
            scene_file >> sx >> sy >> sz;

            Matrix scaling_matrix = getScalingMatrix(sx, sy, sz);
            matrix_stack.push(matrix_stack.top() * scaling_matrix);
        }
        else if (command == "rotate"){
            float angle, ax, ay, az;
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

Matrix getViewTransformationMatrix(Vector& l, Vector& r, Vector& u){
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
            float value = j == 0 ? v.x : (j == 1 ? v.y : v.z);
            R.set(i,j,value);
        }
    }
    R.set(3,3,1);

    Matrix V = R * T;
    return V;
}

void transformPoints(Matrix& V, ifstream& stage1_file, ofstream& stage2_file){
    stage2_file << fixed << setprecision(7);
    while(true){
        for(int i=0; i<3; i++){
            Point p;
            stage1_file >> p.x >> p.y >> p.z;
            if(stage1_file.eof()){
                return;
            }
            Matrix transformed_p = V * pointToMatrix(p);
            p.x = properValue(transformed_p.get(0,0));
            p.y = properValue(transformed_p.get(1,0));
            p.z = properValue(transformed_p.get(2,0));
    
            stage2_file << p.x << " " << p.y << " " << p.z << endl;
        }
        stage2_file << endl;
    }
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
    
    Vector l = look - eye;
    l.normalize();
    Vector r = l ^ up;
    r.normalize();
    Vector u = r ^ l;
    u.normalize();

    ifstream stage1_file_input(stage1_file_path);
    ofstream stage2_file_output(stage2_file_path);

    if(!stage1_file_input.is_open() || !stage2_file_output.is_open()){   
        cerr << "Unable to open required files for stage 2" << endl;
        exit(1);
    }

    Matrix V = getViewTransformationMatrix(l, r, u);

    transformPoints(V, stage1_file_input, stage2_file_output);
    stage1_file_input.close();
    stage2_file_output.close();
    cout << "Stage 2 completed. Output written to "<< stage2_file_path << endl;

    return 0;
}