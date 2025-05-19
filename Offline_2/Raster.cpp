#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<GL/glut.h>
#include <vector>
#include<stack>
#include<filesystem>

using namespace std;

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
};

class Perspective{
    public:
    float fovY, aspectRatio, zNear, zFar;
    Perspective(){
        fovY = 0;
        aspectRatio = 0;
        zNear = 0;
        zFar = 0;
    }
    Perspective(float fovY, float aspectRatio, float zNear, float zFar){
        this->fovY = fovY;
        this->aspectRatio = aspectRatio;
        this->zNear = zNear;
        this->zFar = zFar;
    }
};


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
                for (int k = 0; k < 4; ++k) {
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

Point eye, look, up;
Perspective perspective;
vector<Triangle> triangles;


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
                >> perspective.fovY >> perspective.aspectRatio >> perspective.zNear >> perspective.zFar;

}

void perform_transformations(ifstream& scene_file, ofstream& stage1_file){
    stack<Matrix> matrix_stack;
    matrix_stack.push(getIdentityMatrix());
    while (true){
        string command;
        getline(scene_file, command);
        stringstream ss(command);
        ss >> command;
        if(command == "triangle"){
            cout<< "Reading triangle" << endl;
            for(int i=0;i<3;i++){
                float x,y,z;
                scene_file >> x >> y >> z;
                Point p = Point(x,y,z);
                Matrix transformed_p = matrix_stack.top() * pointToMatrix(p);
                p.x = transformed_p.get(0,0);
                p.y = transformed_p.get(1,0);
                p.z = transformed_p.get(2,0);
                stage1_file << p.x << " " << p.y << " " << p.z << endl;
            }
        }
        else if(command == "end"){
            break;
        }

    }
    stage1_file.close();
    cout << "Stage 1 completed. Output written to stage1.txt" << endl;
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
    filesystem::create_directories("./output");
    ifstream scene_file(scene_file_path);
    ofstream stage1_file("./output/stage1.txt");
    if(scene_file.is_open()){
        if(!stage1_file.is_open()){
            cerr << "Unable to open stage1 file for writing." << endl;
            return 1;
        }
        basic_setup(scene_file);
        perform_transformations(scene_file, stage1_file);
        
        scene_file.close();
    }
    else{
        cerr << "Unable to open scene file: " << scene_file_path << endl;
    }

    return 0;
}