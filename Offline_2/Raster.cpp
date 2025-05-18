#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<GL/glut.h>
#include <vector>

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


Point eye, look, up;
Perspective perspective;
vector<Triangle> triangles;



void basic_setup(ifstream& scene_file){
    string line;
    
    getline(scene_file, line);
    stringstream ss(line);
    ss >> eye.x >> eye.y >> eye.z;
    
    getline(scene_file, line);
    ss.clear();
    ss.str(line);
    ss >> look.x >> look.y >> look.z;
    
    getline(scene_file, line);
    ss.clear();
    ss.str(line);
    ss >> up.x >> up.y >> up.z;
    
    getline(scene_file, line);
    ss.clear();
    ss.str(line);
    ss >> perspective.fovY >> perspective.aspectRatio >> perspective.zNear >> perspective.zFar;

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

    ifstream scene_file(scene_file_path);

    if(scene_file.is_open()){
        basic_setup(scene_file);
        


        scene_file.close();

    }
    else{
        cerr << "Unable to open scene file: " << scene_file_path << endl;
    }


    return 0;
}