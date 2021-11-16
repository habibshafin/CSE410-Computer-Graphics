#include<iostream>
#include<fstream>
#include<vector>
#include<stack>
#include<cmath>
#include <iomanip>
#include <bits/stdc++.h>
#include "bitmap_image.hpp"

using namespace std;

struct Point{
    double x;
    double y;
    double z;
};
double pointDistance(Point a , Point b){
    return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y)+(a.z-b.z)*(a.z-b.z));
}
struct Point eye,look,up;
double fovY, aspectRatio, near, far;
typedef vector<vector<double> > matrixArray;
int pushN = 0;
vector<int> pushCount;

class Matrix{
public:
    matrixArray arr;
    Matrix(int r, int c){
        arr = vector<vector<double> >(r);
        for(int i=0; i<r; i++)
            arr[i] = vector<double>(c);
    }
    void print(){
        for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                cout<<fixed<<setprecision(7)<<arr[i][j]<<"   ";
            }
            cout<<endl;
        }
    }
    Matrix MatrixMultiplication(Matrix a){
        Matrix res(4,4);
        double val;
        for(int i=0; i<4; i++){
            for(int c=0; c<4; c++){
                val = 0;
                for(int j=0; j<4; j++){
                    val = val + arr[i][j] * a.arr[j][c];
                }
                res.arr[i][c] = val;
            }
        }
        /*for(int i=0; i<4; i++){
            for(int j=0; j<4; j++){
                cout<<res.arr[i][j]<<"  ";
            }cout<<endl;
        }*/
        return res;
    }
};

Matrix getIdentityMatrix(){
    Matrix m(4,4);
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            if(i==j){
                m.arr[i][j] = 1;
            }

            else
                m.arr[i][j] = 0;
        }
    }

    return m;
}
stack<Matrix> transMatrix;

class Triangle{
public:
    Point points[3];
    int color[3];
    void print(ofstream &out){
        for(int i=0; i<3; i++){
            out<<fixed<<setprecision(7)<<points[i].x<<"  "<<points[i].y<<"  "<<points[i].z<<endl;
        }
        out<<endl;
    }
};
vector<Triangle> TriangleVector;

Point PointTransMultiply(Matrix trans, Point p){
    Matrix pMatrix(4,1);
    pMatrix.arr[0][0] = p.x;
    pMatrix.arr[1][0] = p.y;
    pMatrix.arr[2][0] = p.z;
    pMatrix.arr[3][0] = 1;
    vector<double> tArr(4);
    for(int i=0;i<4; i++){
        double val = 0 ;
        for(int j=0;j<4; j++){
            val = val + trans.arr[i][j]*pMatrix.arr[j][0];
        }
        tArr[i] = val;
    }
    for(int i=0;i<3;i++){
        tArr[i] = tArr[i]/tArr[3];
    }
    Point ret = {tArr[0],tArr[1],tArr[2]};
    return ret;
}

Matrix TranslateFunction(Matrix trans, double x, double y, double z){
    Matrix translate = getIdentityMatrix();
    translate.arr[0][3] = x;
    translate.arr[1][3] = y;
    translate.arr[2][3] = z;
    return trans.MatrixMultiplication(translate);
}

Matrix ScaleFunctions(Matrix trans, double x, double y, double z){
    Matrix scale = getIdentityMatrix();
    scale.arr[0][0] = x;
    scale.arr[1][1] = y;
    scale.arr[2][2] = z;
    return trans.MatrixMultiplication(scale);
}
Point multiply(Point p, double val){
    Point res = { p.x*val, p.y*val, p.z*val};
    return res;
}
double dotMultiply(Point a, Point b ){
    return a.x*b.x + a.y*b.y + a.z*b.z ;
}
Point crossMultiply(Point a, Point b){
    Point res = { a.y*b.z-a.z*b.y , -a.x*b.z+a.z*b.x, a.x*b.y-a.y*b.x };
    return res;
}
Point normaliseVector(Point p){
    double val = sqrt(p.x*p.x + p.y*p.y + p.z*p.z);
    Point res = { p.x/val,p.y/val, p.z/val };
    return res;
}
Point Rodrigues(Point axis, Point a, double angle ){
    //cout<<angle<<endl;
    Point p1 = multiply(axis,cos(angle*M_PI/180));
    //cout<<"p1 : "<<p1.x<<"  "<<p1.y<<"  "<<p1.z<<endl;
    double dotP = dotMultiply(a,axis);
    Point p2 = multiply(a,dotP*(1-cos(angle*M_PI/180)));
    //cout<<"p2 : "<<p2.x<<"  "<<p2.y<<"  "<<p2.z<<endl;
    Point p3 = multiply(crossMultiply(a,axis),sin(angle*M_PI/180));
    //cout<<"p3 : "<<p3.x<<"  "<<p3.y<<"  "<<p3.z<<endl;
    Point res = { p1.x+p2.x+p3.x, p1.y+p2.y+p3.y, p1.z+p2.z+p3.z};
    return res;
}
Matrix rotateFunction(Matrix trans, Point a, double angle){
    Matrix rotMat = getIdentityMatrix();
    Point normalisedA = normaliseVector(a);
    //cout<<"normalised A : "<< normalisedA.x<<"  "<<normalisedA.y<<" "<<normalisedA.z <<endl;
    Point c[3];
    c[0] = Rodrigues({1,0,0},normalisedA, angle);
    c[1] = Rodrigues({0,1,0},normalisedA, angle);
    c[2] = Rodrigues({0,0,1},normalisedA, angle);
    /*cout<<"c1 : "<<c[0].x<<"   "<<c[0].y<<"  "<<c[0].z<<endl;
    cout<<"c2 : "<<c[1].x<<"   "<<c[1].y<<"  "<<c[1].z<<endl;
    cout<<"c3 : "<<c[2].x<<"   "<<c[2].y<<"  "<<c[2].z<<endl;*/
    for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            double val;
            if(i==0) val = c[j].x;
            else if(i==1) val = c[j].y;
            else val = c[j].z;
            rotMat.arr[i][j] = val;
        }
    }
    /*cout<<endl<<"In rotation : "<<endl;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            cout<<rotMat.arr[i][j]<<"  ";
        }
        cout<<endl;
    }*/

    return trans.MatrixMultiplication(rotMat);
}

bool comparePointY(Point p1, Point p2){
    return (p1.y < p2.y);
}

void assignColors(){
    for(int i=0; i<TriangleVector.size(); i++){
        for(int j=0;j<3; j++){
            TriangleVector[i].color[j] = rand()%256;
        }
    }
}

int main(){
    //ifstream ifile("config.txt", ios::in);
    Matrix identity = getIdentityMatrix();
    transMatrix.push(identity);
    pushCount.push_back(0);

    ifstream sc;
    sc.open("scene4.txt");
    ofstream out("stage1.txt");
    double x,y,z;
    sc>>x>>y>>z;
    eye = {x, y, z};
    sc>>x>>y>>z;
    look = {x, y, z};
    sc>>x>>y>>z;
    up = {x, y, z};
    sc>>fovY>>aspectRatio>>near>> far;

    string command;
    int counter = 0;

    ///section for testing purpose
    /*Point u = {2,1,-3};
    Point v = {0,4,5};
    Point r = crossMultiply(u,v);
    cout<<r.x<<" " <<r.y<<" "<<r.z<<endl;
    */

    while(true){
        /*cout<<"translation matrix :"<<endl;
        transMatrix.top().print();
        cout<<endl<<endl;*/

        sc>>command;
        if (command=="triangle"){
            Triangle t;
            //cout<<endl;transMatrix.top().print();
                //cout<<endl;
            for(int i=0; i<3; i++){
                sc>>x>>y>>z;
                t.points[i] = PointTransMultiply(transMatrix.top(),{x,y,z});
                out<<fixed<<setprecision(7)<<t.points[i].x<<"  "<<t.points[i].y<<"  "<<t.points[i].z<<endl;
            }
            TriangleVector.push_back(t);
            out<<endl;
        }else if (command=="translate"){
            sc>>x>>y>>z;
            transMatrix.push(TranslateFunction(transMatrix.top(),x,y,z));
            pushCount[pushN]++;
            //cout<<endl<<"After translation"<<endl;
        }
        else if (command=="scale"){
            sc>>x>>y>>z;
            transMatrix.push(ScaleFunctions(transMatrix.top(),x,y,z));
            pushCount[pushN]++;
            //cout<<endl<<"After scaling"<<endl;
        }else if (command=="rotate"){
            double rotAngle;
            sc>>rotAngle;
            sc>>x>>y>>z;
            transMatrix.push(rotateFunction(transMatrix.top(), {x,y,z}, rotAngle));
            pushCount[pushN]++;
        }
        else if (command=="push"){
            pushN++;
            pushCount.push_back(0);
        }else if (command=="pop"){
            for(int i=0; i<pushCount[pushN]; i++){
                transMatrix.pop();
            }
            pushCount.pop_back();
            pushN--;
        }
        else if (command=="end"){
            break;
        }

    }
    sc.close();
    out.close();
    ///VIEW TRANSFORMATION
    Point l = normaliseVector({look.x-eye.x,look.y-eye.y,look.z-eye.z});
    Point r = normaliseVector(crossMultiply(l,up));
    Point u = crossMultiply(r,l);
    //cout<<l.x<<" "<<l.y<<" "<<l.z<<endl;

    Matrix T = getIdentityMatrix();
    T.arr[0][3] = -eye.x;
    T.arr[1][3] = -eye.y;
    T.arr[2][3] = -eye.z;
    Matrix R = getIdentityMatrix();
    R.arr[0][0] = r.x;
    R.arr[0][1] = r.y;
    R.arr[0][2] = r.z;
    R.arr[1][0] = u.x;
    R.arr[1][1] = u.y;
    R.arr[1][2] = u.z;
    R.arr[2][0] = -l.x;
    R.arr[2][1] = -l.y;
    R.arr[2][2] = -l.z;

    Matrix V = R.MatrixMultiplication(T);

    ofstream out2("stage2.txt");
    for(int i=0; i<TriangleVector.size(); i++){
        for(int j = 0; j<3; j++){
            TriangleVector[i].points[j] = PointTransMultiply(V,TriangleVector[i].points[j]);
        }
        TriangleVector[i].print(out2);
    }
    out2.close();

    ///PROJECTION TRANSFORMATION
    double fovX = fovY * aspectRatio;
    double t = near * tan((fovY/2)*M_PI/180);
    double r2 = near * tan((fovX/2)*M_PI/180);
    Matrix P = getIdentityMatrix();
    P.arr[0][0] = near/r2;
    P.arr[1][1] = near/t;
    P.arr[2][2] = -(far+near)/(far-near);
    P.arr[2][3] = -(2*far*near)/(far-near);
    P.arr[3][2] = -1;
    P.arr[3][3] = 0;
    ofstream out3("stage3.txt");
    for(int i=0; i<TriangleVector.size(); i++){
        for(int j = 0; j<3; j++){
            TriangleVector[i].points[j] = PointTransMultiply(P,TriangleVector[i].points[j]);
        }
        TriangleVector[i].print(out3);
    }
    out3.close();
    ifstream conf;
    ofstream out4("z_buffer.txt");

    conf.open("config.txt");
    int screenWidth, screenHeight;
    double leftlimitX,bottomlimitY,frontZ,rearZ;
    conf>>screenWidth>>screenHeight>>leftlimitX>>bottomlimitY>>frontZ>>rearZ;
    double rightlimitX = -leftlimitX;
    double toplimitY = -bottomlimitY;
    double dx = (rightlimitX - leftlimitX) / screenWidth;
    double dy = (toplimitY - bottomlimitY ) / screenHeight;
    double topY = toplimitY - dy;
    double leftX = leftlimitX + dx;

    //double zbuffer[screenHeight][screenWidth];
    double** zbuffer = new double*[screenWidth];
    Point ** pixelinfo = new Point*[screenWidth];
    for(int c = 0; c < screenWidth; c++){
        zbuffer[c] = new double[screenHeight];
        pixelinfo[c] = new Point[screenHeight];
    }
    assignColors();



    for(int i=0; i<screenHeight; i++){
        for(int j =0; j<screenWidth; j++){
            zbuffer[i][j] = rearZ;
            pixelinfo[i][j].x = 0;
            pixelinfo[i][j].y = 0;
            pixelinfo[i][j].z = 0;
        }
    }

    ///Z-buffer algorithm
    for(int i=0; i<TriangleVector.size(); i++){
        double topScanline = min(toplimitY,max(TriangleVector[i].points[0].y,max(TriangleVector[i].points[1].y,TriangleVector[i].points[2].y)));
        double bottomScanline = max(bottomlimitY,min(TriangleVector[i].points[0].y,min(TriangleVector[i].points[1].y,TriangleVector[i].points[2].y)));
        //cout<<"top and bottom"<<topScanline<<"   "<<bottomScanline<<"   ";

        int rowstart = ceil((toplimitY-topScanline)/dy);
        int rowfinish = ceil((toplimitY-bottomScanline)/dy);
        //cout<<rowstart<<"  "<<rowfinish<<endl;
        for(int j=rowstart; j<=rowfinish; j++){
            Point D,E,F;
            D.y = toplimitY- j*dy;
            E.y = D.y;
            F.y = D.y;
            Triangle t = TriangleVector[i];
            D.x = ( (D.y-t.points[0].y)/(t.points[1].y-t.points[0].y) ) * (t.points[1].x-t.points[0].x) + t.points[0].x;
            E.x = ( (E.y-t.points[1].y)/(t.points[2].y-t.points[1].y) ) * (t.points[2].x-t.points[1].x) + t.points[1].x;
            F.x = ( (F.y-t.points[2].y)/(t.points[0].y-t.points[2].y) ) * (t.points[0].x-t.points[2].x) + t.points[2].x;

            D.z = ( (D.y-t.points[0].y)/(t.points[1].y-t.points[0].y) ) * (t.points[1].z-t.points[0].z) + t.points[0].z;
            E.z = ( (E.y-t.points[1].y)/(t.points[2].y-t.points[1].y) ) * (t.points[2].z-t.points[1].z) + t.points[1].z;
            F.z = ( (F.y-t.points[2].y)/(t.points[0].y-t.points[2].y) ) * (t.points[0].z-t.points[2].z) + t.points[2].z;

            vector<Point> correctPs;
            /*cout<<(pointDistance(t.points[0],D)+pointDistance(t.points[1],D))<<"   "<< pointDistance(t.points[0],t.points[1])<<endl;
            cout<<(pointDistance(t.points[1],E)+pointDistance(t.points[2],E))<<"   "<< pointDistance(t.points[1],t.points[2])<<endl;
            cout<<(pointDistance(t.points[2],F)+pointDistance(t.points[0],F))<<"   "<< pointDistance(t.points[2],t.points[0])<<endl;
            cout<<endl;*/
            if(abs((pointDistance(t.points[0],D)+pointDistance(t.points[1],D)) - pointDistance(t.points[0],t.points[1]))<0.001)
                correctPs.push_back(D);
            if(abs((pointDistance(t.points[1],E)+pointDistance(t.points[2],E)) - pointDistance(t.points[1],t.points[2]))<0.001)
                correctPs.push_back(E);
            if(abs((pointDistance(t.points[2],F)+pointDistance(t.points[0],F)) - pointDistance(t.points[2],t.points[0]))<0.001)
                correctPs.push_back(F);

            /*if(correctPs.size()==1)
                cout<<"at vertex point"<<endl;
            else if(correctPs.size()==2)
                cout<<"two intersecting points"<<endl;
            else cout<<"jhamela"<<endl;
            /*cout<<D.x<<"  "<<D.y<<"  "<<D.z<<endl;
            cout<<E.x<<"  "<<E.y<<"  "<<E.z<<endl;
            cout<<F.x<<"  "<<F.y<<"  "<<F.z<<endl;*/
            /*cout<<endl;
            for(int l=0; l<correctPs.size(); l++)
                cout<<correctPs[l].x<<"  "<<correctPs[l].y<<endl;
            cout<<endl;*/
            int colstart,colfinish;
            double colLeft,colRight,startZ,finishZ,dz;
            if(correctPs.size()==2){
                //colLeft = max(min(correctPs[0].x,correctPs[1].x), leftlimitX);
                //colRight = min( max(correctPs[0].x,correctPs[1].x),rightlimitX);
                if(correctPs[0].x<correctPs[1].x){
                    colLeft = max(leftlimitX,correctPs[0].x);
                    colRight = min(rightlimitX,correctPs[1].x);
                    dz = (correctPs[1].z-correctPs[0].z)/(correctPs[1].x-correctPs[0].x);
                    startZ = correctPs[0].z + ((colLeft - correctPs[0].x)/dx) * dz;
                    finishZ = correctPs[0].z + ((colRight - correctPs[0].x)/dx) * dz;
                }else{
                    colLeft = max(leftlimitX,correctPs[1].x);
                    colRight = min(rightlimitX,correctPs[0].x);
                    dz = (correctPs[0].z-correctPs[1].z)/(correctPs[0].x-correctPs[1].x);
                    startZ = correctPs[1].z + ((colLeft - correctPs[1].x)/dx) * dz;
                    finishZ = correctPs[1].z + ((colRight - correctPs[1].x)/dx) * dz;
                }
                colstart = round((colLeft-leftlimitX)/dx);
                colfinish = round((colRight-leftlimitX)/dx);
                if(colstart<0)
                    colstart = 0;
                if(colfinish>=screenWidth)
                    colfinish = screenWidth-1;
            }
            //cout<<"dz : "<<dz<<endl;
            //cout<<"colstart :"<<colstart<<"  colfinish : "<<colfinish<<endl;
            //cout<<"startZ :"<<startZ<<"  finishZ : "<<finishZ<<endl;*/
            double dz1 = (finishZ-startZ)/(colfinish-colstart);
            double cnt = 0;
            for(int c=colstart; c<=colfinish; c++){
                double zval = startZ + dz1*cnt;
                //cout<<zval<<endl;
                if(zval < zbuffer[j][c]){
                    zbuffer[j][c] = zval;
                    pixelinfo[j][c].x = TriangleVector[i].color[0];
                    pixelinfo[j][c].y = TriangleVector[i].color[1];
                    pixelinfo[j][c].z = TriangleVector[i].color[2];
                }
                cnt = cnt + 1;
            }
        }

    }

    for(int i=0;i<500; i++){
        for(int j=0; j<500; j++){
            if(zbuffer[i][j]<2)
                out4<<zbuffer[i][j]<<"\t";
        }
        out4<<endl;
    }

    bitmap_image image(screenWidth,screenHeight);

    for(int h=0; h<screenHeight; h++){
        for(int w=0; w<screenWidth; w++){
            image.set_pixel(h,w,pixelinfo[w][h].x,pixelinfo[w][h].y,pixelinfo[w][h].z);
        }
    }

    image.save_image("output.bmp");;
    free(zbuffer);
    free(pixelinfo);

    return 0;
}

