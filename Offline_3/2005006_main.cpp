#include "2005006_classes.hpp"
#include <bits/stdc++.h>
#include<GL/glut.h>

using namespace std;

vector<Object> objects;
vector<Light> lights;

string scene_file = "scene.txt";
int recursion_levels;
int image_width;
int image_height;

void Sphere :: draw(){
    // lalala
}
void Triangle :: draw(){
    // lalala
}
void Floor :: draw(){
    // lalala
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

int main(){
    load_data();
    return 0;
}