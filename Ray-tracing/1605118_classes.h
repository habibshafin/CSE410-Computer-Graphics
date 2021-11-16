#include<stdio.h>
#include<iostream>
#include<bits/stdc++.h>

#define pi (2*acos(0.0))
using namespace std;


struct point{
    double x,y,z;
};

class Vector3D{
public:
    double x,y,z;
    Vector3D() {};
    Vector3D(double a, double b, double c){
        x=a;
        y=b;
        z=c;
    }
    Vector3D Minus(Vector3D b){
        Vector3D res(x-b.x, y-b.y, z-b.z);
        return res;
    }
    double dotMultiply(Vector3D b){
        return x*b.x + y*b.y + z*b.z;
    }
    void Normalise(){
        double val = sqrt( x*x + y*y + z*z);
        x = x/val;
        y = y/val;
        z = z/val;
    }
    Vector3D add( Vector3D b){
        Vector3D res( x+b.x, y+b.y, z+b.z );
        return res;
    }
    Vector3D multiply(double val){
        Vector3D res( x*val, y*val , z*val );
        return res;
    }
    Vector3D cross(Vector3D b){
    Vector3D res( (y*b.z - z*b.y),( z*b.x - x*b.z ),( x*b.y - y*b.x ) );
    return res;
    }
};

class Ray{
public:
    Vector3D start;
    Vector3D dir;
    Ray(){};
    Ray(Vector3D s, Vector3D d){
        start.x = s.x; start.y = s.y; start.z = s.z;
        dir.x = d.x; dir.y = d.y; dir.z = d.z;
    }
    void print(){
        cout<<"start : "<<start.x<<" , "<<start.y<<" , "<<start.z<<endl;
        cout<<"dir : "<<dir.x<<" , "<<dir.y<<" , "<<dir.z<<endl;
    }
};

class Object{
public:
    Vector3D reference_point;
    double height,width,length;
    double color[3];
    double coEfficients[4];
    int shine;
    Object(){};
    virtual void draw(){};
    virtual double getMinT(Ray *r){};
    virtual double intersect(Ray *r, double *color, int level){
        return -1.0;
    }
    void setReferencePoint(double x , double y, double z){
        reference_point.x = x;
        reference_point.y = y;
        reference_point.z = z;
    }
    void setShine(int sh){
        shine = sh;
    }
    void setColor(double r ,double g, double b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }
    void setCoefficients(double a, double b, double c, double d){
        coEfficients[0] = a;
        coEfficients[1] = b;
        coEfficients[2] = c;
        coEfficients[3] = d;
    }
};

class Light{
public:
    Vector3D reference_point;
    double color[3];
    Light(){};
    void setReferencePoint(double x , double y, double z){
        reference_point.x = x;
        reference_point.y = y;
        reference_point.z = z;
    }
    void setColor(double r ,double g, double b){
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }
    void draw(){
        glPushMatrix();
            glTranslated(reference_point.x,reference_point.y, reference_point.z);
            glColor3f(color[0],color[1],color[2]);
            glutSolidSphere(2,5,5);
        glPopMatrix();
    }
};

extern point pos,u,r,l;
extern vector<Object*> objList;
extern vector<Light*> lights;



class Sphere : public Object{
public:
    Sphere(){};
    Sphere(Vector3D center, double rad){
        setReferencePoint(center.x, center.y, center.z);
        length = rad;
    }
    void draw(){
        /*Vector3D points[100][100];
        double radius = length;
        int stacks = 50;
        int slices = 50;
        int i,j;
        double h,r;
        //generate points
        for(i=0;i<=stacks;i++){
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++){
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }*/
        glPushMatrix();
        glTranslated(reference_point.x,reference_point.y, reference_point.z);
        glColor3f(color[0],color[1],color[2]);

        //draw quads using generated points
        /*for(i=0;i<stacks;i++){
            //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
            glColor3f(color[0],color[1],color[2]);
            for(j=0;j<slices;j++){
			glBegin(GL_QUADS);{
            //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
        }*/
        glutSolidSphere(length,50,50);
        glPopMatrix();
    }

    double getMinT(Ray *r){
        Vector3D Ro = r->start.Minus(reference_point);
        Vector3D Rd = r->dir;
        //Rd.Normalise();

        double a = Rd.dotMultiply(Rd);
        double b = 2.0 * Ro.dotMultiply(Rd);
        double c = Ro.dotMultiply(Ro) - length*length;

        double descriminant = (b*b) - (4.0*a*c);
        //cout<<"a:b:c:desc"<<a<<","<<b<<","<<c<<" , "<<descriminant<<endl;

        if(descriminant<0){
            return -1;
        }
        else{
            double t = min(( -b + sqrt(descriminant) )/2,( -b - sqrt(descriminant) )/2);
            t = t<=0 ? max(( -b + sqrt(descriminant) )/2,( -b - sqrt(descriminant) )/2) : t;
            return t;
        }
    }
    double intersect(Ray *r, double *color, int level){
        //return getMinT(r);
        double t = getMinT(r);
        if(t>=0){
            color[0] = this->color[0]*coEfficients[0];
            color[1] = this->color[1]*coEfficients[0];
            color[2] = this->color[2]*coEfficients[0];
        }
        if(level==0) return t;

        Vector3D intersectionPoint = r->start.add(r->dir.multiply(t));
        Vector3D normal = intersectionPoint.Minus(reference_point);
        normal.Normalise();

        for(int i =0; i<lights.size(); i++){
            Vector3D lightDir = lights.at(i)->reference_point.Minus(intersectionPoint);
            lightDir.Normalise();
            Ray *l = new Ray(lights.at(i)->reference_point, lightDir);

            double distance = sqrt(lightDir.x*lightDir.x + lightDir.y*lightDir.y + lightDir.z*lightDir.z);
            int j;
            for( j=0; j<objList.size(); j++){
                double dist = objList.at(j)->getMinT(l);
                if(dist>0 && dist<distance)
                    break;
            }
            if ( j!=objList.size())
                continue;

            double lambertValue = normal.dotMultiply(l->dir);
            Vector3D reflectedRay = normal.multiply((-2)*r->dir.dotMultiply(normal)).add(r->dir);
            double phongValue = reflectedRay.dotMultiply(r->dir);

            if(lambertValue < 0)
                lambertValue = 0;
            else if(lambertValue>1)
                lambertValue = 1;
            if(phongValue < 0)
                phongValue = 0;
            else if(phongValue>1)
                phongValue = 1;

            double mul = (lambertValue*coEfficients[1]+pow(phongValue,shine)*coEfficients[2]);
            color[0] = color[0] + lights.at(i)->color[0]*mul;
            color[1] = color[1] + lights.at(i)->color[1]*mul;
            color[2] = color[2] + lights.at(i)->color[2]*mul;

        }
    }
};

class Triangle : public Object{
public:
    Vector3D vtxA, vtxB, vtxC;
    Triangle(){};
    Triangle(Vector3D vtx_a, Vector3D vtx_b, Vector3D vtx_c ){
        vtxA.x = vtx_a.x; vtxA.y = vtx_a.y; vtxA.z = vtx_a.z;
        vtxB.x = vtx_b.x; vtxB.y = vtx_b.y; vtxB.z = vtx_b.z;
        vtxC.x = vtx_c.x; vtxC.y = vtx_c.y; vtxC.z = vtx_c.z;
    }
    void draw(){
        //cout<<"draw : triangle"<<endl;
        glBegin(GL_TRIANGLES);{
            glColor3f(color[0], color[1], color[2]);
            glVertex3f(vtxA.x, vtxA.y, vtxA.z);
            glVertex3f(vtxB.x, vtxB.y, vtxB.z);
            glVertex3f(vtxC.x, vtxC.y, vtxC.z);
        }glEnd();
    }
    double getMinT(Ray *r){
        //cout<<"get min t called "<<endl;
        Vector3D normal = vtxB.Minus(vtxA).cross(vtxC.Minus(vtxA));
        normal.Normalise();
        //cout<<"normal : "<<normal.x<<", "<<normal.y<<" , "<<normal.z<<endl;
        double den = normal.dotMultiply(r->dir);
        //cout<<"den1 : "<<den<<endl;
        if(abs(den) < 0.00001){
            return -1;
        }
        if(den<0){
            normal = normal.multiply(-1);
            den = normal.dotMultiply(r->dir);
        }
        //cout<<"den 2 : "<<den<<endl;
        //double t =  (normal.dotMultiply(r->start) + normal.dotMultiply(vtxA)) / den;
        double t =  (normal.dotMultiply(vtxA.Minus(r->start)) ) / den;
        //cout<<" t : "<<t<<endl;
        if(t<0)
            return -1;

        Vector3D n = vtxB.Minus(vtxA).cross(vtxC.Minus(vtxA));
        Vector3D edge0 = vtxB.Minus(vtxA);
        Vector3D edge1 = vtxC.Minus(vtxB);
        Vector3D edge2 = vtxA.Minus(vtxC);
        Vector3D P = r->start.add(r->dir.multiply(t));
        //cout<<" P: "<<P.x<<" , "<<P.y<<" , "<<P.z<<endl;
        Vector3D C0 =  P.Minus(vtxA);
        Vector3D C1 = P.Minus(vtxB);
        Vector3D C2 = P.Minus(vtxC);

        if(n.dotMultiply(edge0.cross(C0))<0 || n.dotMultiply(edge1.cross(C1))<0 || n.dotMultiply(edge2.cross(C2))<0.00001)
        {
            return -1;
        }
        else
            return t;
    }
    double intersect(Ray *r, double *color, int level){
        double t = getMinT(r);
        if(t>=0){
            color[0] = this->color[0]*coEfficients[0];
            color[1] = this->color[1]*coEfficients[0];
            color[2] = this->color[2]*coEfficients[0];
        }
        return t;
    }
};

class General : public Object{
public:
    double a,b,c,d,e,f,g,h,i,j;
    General(){};
    General(double k,double l,double m,double n,double o,double p,double q,double r,double s,double t){
        a=k;b=l;c=m;d=n;e=o;f=p;g=q;h=r;i=s;j=t;
    }
    void draw(){
        //cout<<"draw : general : do nothing "<<endl;
    }
    double getMinT(Ray *r){
        double A = a*(r->dir.x*r->dir.x) + b*(r->dir.y*r->dir.y) + c*(r->dir.z*r->dir.z);
        A = A + d*r->dir.x*r->dir.y + e*r->dir.x*r->dir.z + f*r->dir.y*r->dir.z;
        double B = 2*a*r->start.x*r->dir.x + 2*b*r->start.y*r->dir.y + 2*c*r->start.z*r->dir.z;
        B = B + d*(r->start.x*r->dir.y + r->start.y*r->dir.x);
        B = B + e*(r->start.x*r->dir.z + r->start.z*r->dir.x);
        B = B + f*(r->start.y*r->dir.z + r->start.z*r->dir.y);
        B = B + g*r->dir.x + h*r->dir.y + i*r->dir.z;
        double C = a*(r->start.x*r->start.x) + b*(r->start.y*r->start.y) + c*(r->start.z*r->start.z);
        C = C + d*r->start.x*r->start.y + e*r->start.x*r->start.z + f*r->start.y*r->start.z;
        C = C + g*r->start.x + h*r->start.y + i*r->start.z + j;

        if(abs(A) < 0.001){
            return -C/B;
        }
        double discriminant = B*B - 4.0*A*C;
        if(discriminant<0.00001){
            return -1;
        }
        double t1 = (-B + sqrt(discriminant))/(2*A);
        double t2 = (-B - sqrt(discriminant))/(2*A);
        double t = min(t1,t2);
        Vector3D intersectionPoint = r->start.add(r->dir.multiply(t));
        if(t>0){
            if(length!=0){
                if((intersectionPoint.x-reference_point.x)>(length/2) || (intersectionPoint.x-reference_point.x) < (length/2))
                    t = -1;
            }
            if(width!=0){
                if((intersectionPoint.y-reference_point.y)>(width/2) || (intersectionPoint.y-reference_point.y) < (width/2))
                    t = -1;
            }
            if(height!=0){
                if((intersectionPoint.z-reference_point.z)>height || (intersectionPoint.z-reference_point.z)<0)
                    t = -1;
            }
        }
        if(t<=0){
            t = max(t1,t2);
        }

        intersectionPoint = r->start.add(r->dir.multiply(t));
        if(length!=0){
            if((intersectionPoint.x-reference_point.x)>(length/2) || (intersectionPoint.x-reference_point.x) < (length/2))
                t = -1;
        }
        if(width!=0){
            if((intersectionPoint.y-reference_point.y)>(width/2) || (intersectionPoint.y-reference_point.y) < (width/2))
                t = -1;
        }
        if(height!=0){
            if((intersectionPoint.z-reference_point.z)>height || (intersectionPoint.z-reference_point.z)<0)
                t = -1;
        }

        return t;
    }
    double intersect(Ray *r, double * color, int level){
        double t = getMinT(r);
        if(t>=0){
            color[0] = this->color[0]*coEfficients[0];
            color[1] = this->color[1]*coEfficients[0];
            color[2] = this->color[2]*coEfficients[0];
        }
        return t;
    }
};

class Floor : public Object{
public:
    Floor(double floorWidth, double tilewidth){
        setReferencePoint(-floorWidth/2,-floorWidth/2, 0);
        length = tilewidth;
    }
    void draw(){
        int tilesCount = (-reference_point.x*2)/length;

        for(int i=0; i<tilesCount; i++){
            for(int j=0; j<tilesCount; j++){
                if(i%2==0 && j%2==0 || i%2==1 && j%2==1)
                    glColor3f(1,1,1);
                else
                    glColor3f(0,0,0);
                glBegin(GL_QUADS);{
                    glVertex3f(reference_point.x+i*length+length ,reference_point.y+j*length+length,0);
                    glVertex3f(reference_point.x+i*length+length ,reference_point.y+j*length,0);
                    glVertex3f(reference_point.x+i*length ,reference_point.y+j*length,0);
                    glVertex3f(reference_point.x+i*length ,reference_point.y+j*length+length,0);
                }glEnd();
            }
        }
    }
    double getMinT(Ray *r){
        Vector3D normal(0,0,1);
        double den = normal.dotMultiply(r->dir);
        if(abs(den)<0.000001){
            return -1;
        }
        double t = - (normal.dotMultiply(r->start))/den;
        return t;
    }
    double intersect(Ray *r, double *color, int level){
        double t = getMinT(r);
        if(t>0){
            Vector3D intersectionPoint = r->start.add(r->dir.multiply(t));
            int axis_x = abs(intersectionPoint.x-reference_point.x)/length;
            int axis_y = abs(intersectionPoint.y-reference_point.y)/length;
            if(abs(intersectionPoint.x) <= -reference_point.x && abs(intersectionPoint.y)<=-reference_point.y){
            if((axis_x%2==0 && axis_y%2==0) ||(axis_x%2==1 && axis_y%2==1)){
                color[0] = 1; color[1] = 1; color[2] = 1;
            }else{
                color[0] = 0; color[1] = 0; color[2] = 0;
            }
            }
        }
        return t;
    }
};


