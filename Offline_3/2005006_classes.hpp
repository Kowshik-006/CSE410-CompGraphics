#include<bits/stdc++.h>
#include<GL/glut.h>
#define deg2rad M_PI/180.0

enum rotation_type {yaw,pitch,roll};
class Vector{
    public:
    double x,y,z;
    Vector() : x(0), y(0), z(0) {}
    Vector(double x, double y, double z) : x(x), y(y), z(z) {}
    // this . other
    double operator*(const Vector& other) const {
        return (x * other.x) + (y * other.y) + (z * other.z);
    }
    // this x other
    Vector operator^(const Vector& other) const {
        return Vector(
            (y * other.z) - (other.y * z),
            (other.x * z) - (x * other.z),
            (x * other.y) - (other.x * y)
        );
    }
    // this * scalar
    Vector operator*(double scalar) const {
        return Vector(x * scalar, y * scalar, z * scalar);
    }
    // this + other
    Vector operator+(const Vector& other) const {
        return Vector(x + other.x, y + other.y, z + other.z);
    }

    Vector operator-(const Vector& other) const {
        return Vector(x - other.x, y - other.y, z - other.z);
    }
    Vector normalize(){
        double length = sqrt((x * x) + (y * y) + (z * z));
        if(length == 0){
            return *this;
        }
        return Vector(x/length, y/length, z/length);
    }
};
class Color{
    public:
    double r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(double r, double g, double b) : r(r), g(g), b(b) {}
};
class CoEfficients{
    public: 
    double ambient, diffuse, specular, reflection;
    CoEfficients() : ambient(0), diffuse(0), specular(0), reflection(0) {}
    CoEfficients(double ambient, double diffuse, double specular, double reflection) 
        : ambient(ambient), diffuse(diffuse), specular(specular), reflection(reflection) {}
};
class Object{
    public:
    Vector reference_point;
    double height, width, length;
    Color color;
    CoEfficients coefficients;
    int shine;
    Object(){}
    Object(Vector ref_point, double height, double width, double length) : 
        reference_point(ref_point), height(height), width(width), length(length){}
    virtual void draw() {}   
    void setColor(Color color) {
        this->color = color;
    }
    void setShine(int shine) {
        this->shine = shine;
    }
    void setCoefficients(CoEfficients coefficients) {
        this->coefficients = coefficients;
    }
};


class Sphere : public Object {
    public: 
    Sphere(Vector center, double radius){
        reference_point = center;
        length = radius;
    }
    void draw() override;
};

class Triangle : public Object {
    public:
    Vector a, b, c;
    Triangle(Vector a, Vector b, Vector c) : a(a), b(b), c(c) {}
    void draw() override;
};

class Light{
    public:
    Vector light_position;
    Color color;
    Light(Vector position, Color color) : light_position(position), color(color){}
};

class PointLight : public Light{
    public:
    PointLight(Vector position, Color color) : Light(position, color){}
};

class SpotLight : public Light {
    public:
    Vector light_direction;
    double cutoff_angle;
    SpotLight(Vector position, Color color, Vector direction, double cutoff_angle) 
        : Light(position, color), light_direction(direction), cutoff_angle(cutoff_angle){}
};

class Floor : public Object {
    Floor(double floor_width, double tile_width){
        reference_point = Vector(-floor_width/2, -floor_width/2, 0);
        length = tile_width;
    }
    void draw();
};
