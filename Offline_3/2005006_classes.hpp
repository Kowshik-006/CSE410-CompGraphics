#include<bits/stdc++.h>

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
    double length() const {
        return sqrt((x * x) + (y * y) + (z * z));
    }
};

class Matrix{
    public:
    double m[3][3];
    Matrix(Vector column1, Vector column2, Vector column3){
        m[0][0] = column1.x; m[0][1] = column2.x; m[0][2] = column3.x;
        m[1][0] = column1.y; m[1][1] = column2.y; m[1][2] = column3.y;
        m[2][0] = column1.z; m[2][1] = column2.z; m[2][2] = column3.z;
    }
    double determinant() const {
        return  m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) 
               -m[0][1] * (m[1][0] * m[2][2] - m[2][0] * m[1][2]) 
               +m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
    }
};
class Color{
    public:
    double r, g, b;
    Color() : r(0), g(0), b(0) {}
    Color(double r, double g, double b) : r(r), g(g), b(b) {}
    Color operator*(double scalar) const {
        return Color(r * scalar, g * scalar, b * scalar);
    }
    Color operator+(const Color& other) const {
        return Color(r + other.r, g + other.g, b + other.b);
    }
    Color operator*(const Color& color) const {
        return Color(r * color.r, g * color.g, b * color.b);
    }
};

class CoEfficients{
    public: 
    double ambient, diffuse, specular, reflection;
    CoEfficients() : ambient(0), diffuse(0), specular(0), reflection(0) {}
    CoEfficients(double ambient, double diffuse, double specular, double reflection) 
        : ambient(ambient), diffuse(diffuse), specular(specular), reflection(reflection) {}
};

class Ray{
    public: 
    Vector start;
    Vector direction;
    Ray(Vector start, Vector direction) : start(start), direction(direction.normalize()){}
};

class Object{
    private:
    double A, B, C, D, E, F, G, H, I, J; // General object coefficients
    bool valid(const Ray& ray, double t);
    public:
    Vector reference_point;
    double height, width, length;
    Color color;
    CoEfficients coefficients;
    int shine;
    Object(){}
    Object(Vector ref_point, double length, double width, double height) : 
        reference_point(ref_point), length(length), width(width), height(height){}
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
    void setGeneralCoefficients(double A, double B, double C, double D, 
                            double E, double F, double G, double H, double I, double J){
        this->A = A; this->B = B; this->C = C; this->D = D;
        this->E = E; this->F = F; this->G = G; this->H = H;
        this->I = I; this->J = J;
    }
    virtual double intersect(const Ray& ray, Color& color, int level);
    virtual Vector getNormal(const Vector& p);
};


class Sphere : public Object {
    public: 
    Sphere(Vector center, double radius){
        reference_point = center;
        length = radius;
    }
    void draw() override;
    double intersect(const Ray& ray, Color& color, int level) override;
    Vector getNormal(const Vector& p) override;
};

class Triangle : public Object {
    public:
    Vector a, b, c;
    Triangle(Vector a, Vector b, Vector c) : a(a), b(b), c(c) {}
    void draw() override;
    double intersect(const Ray& ray, Color& color, int level) override;
    Vector getNormal(const Vector& p) override;
};

class Light{
    public:
    Vector light_position;
    Color color;
    Light(Vector position, Color color) : light_position(position), color(color){}
    virtual ~Light() {}
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
    public:
    Floor(double floor_width, double tile_width){
        reference_point = Vector(-floor_width/2, -floor_width/2, 0);
        length = tile_width;
    }
    void draw() override;
    double intersect(const Ray& ray, Color& color, int level) override;
    Color getColor(Vector p);
    Vector getNormal(const Vector& p) override;
};


