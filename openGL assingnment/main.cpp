#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <glut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double cylinderAngle;
double cylinderAngleL;
double cylinderAngleZ;
double cylinderRotAngle;

struct point
{
	double x,y,z;
};

struct point pos = {100,100, 0};
struct point u = {0,0,1};
struct point r = {-1/sqrt(2),1/sqrt(2) ,0};
struct point l = {-1/sqrt(2),-1/sqrt(2) ,0};

struct point crossMultiply(struct point p1, struct point p2){
    struct point res;
    res.x = p1.y * p2.z - p1.z * p2.y;
    res.y = p1.z * p2.x - p1.x * p2.z;
    res.z = p1.x * p2.y - p1.y * p2.x;
    return res;
};

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);
		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		//glColor3f(1,1,1);
		for(j=0;j<slices;j++)
		{
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
	}
}

void drawHorFSphere(double radius,int slices,int stacks, int choice){
    struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].y=h;
			points[i][j].z=r*cos(((double)j/(double)slices)*2*pi);
		}
	}
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		glColor3f(1,0,0);
		for(j=0;j<slices;j++)
		{
		    if(j%2)
                glColor3f(1,1,1);
            else
                glColor3f(0,0,0);
		    glBegin(GL_QUADS);{
                if(choice==1 || choice==3){
                //right hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

                }
                if (choice==2 || choice==3){
                    //left hemisphere
                    glVertex3f(points[i][j].x,-points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,-points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,-points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,-points[i+1][j].y,points[i+1][j].z);
                }
			}glEnd();
		}
	}
}

void drawSS()
{
    glColor3f(1,0,0);
    drawSquare(20);

    glRotatef(angle,0,0,1);
    glTranslatef(110,0,0);
    glRotatef(2*angle,0,0,1);
    glColor3f(0,1,0);
    drawSquare(15);

    glPushMatrix();
    {
        glRotatef(angle,0,0,1);
        glTranslatef(60,0,0);
        glRotatef(2*angle,0,0,1);
        glColor3f(0,0,1);
        drawSquare(10);
    }
    glPopMatrix();

    glRotatef(3*angle,0,0,1);
    glTranslatef(40,0,0);
    glRotatef(4*angle,0,0,1);
    glColor3f(1,1,0);
    drawSquare(5);
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
		    struct point uxr;
		    uxr = crossMultiply(u,r);
		    r.x = r.x*cos((5.0*pi/180.0)) + uxr.x * sin((5.0*pi/180.0));
		    r.y = r.y*cos((5.0*pi/180.0)) + uxr.y * sin((5.0*pi/180.0));
		    r.z = r.z*cos((5.0*pi/180.0)) + uxr.z * sin((5.0*pi/180.0));
            struct point uxl;
            uxl = crossMultiply(u,l);
			l.x = l.x*cos((5.0*pi/180.0)) + uxl.x * sin((5.0*pi/180.0));
		    l.y = l.y*cos((5.0*pi/180.0)) + uxl.y * sin((5.0*pi/180.0));
		    l.z = l.z*cos((5.0*pi/180.0)) + uxl.z * sin((5.0*pi/180.0));
            break;
        case '2':
		    //struct point uxr;
		    uxr = crossMultiply(u,r);
		    r.x = r.x*cos((-5.0*pi/180.0)) + uxr.x * sin((-5.0*pi/180.0));
		    r.y = r.y*cos((-5.0*pi/180.0)) + uxr.y * sin((-5.0*pi/180.0));
		    r.z = r.z*cos((-5.0*pi/180.0)) + uxr.z * sin((-5.0*pi/180.0));
            //struct point uxl;
            uxl = crossMultiply(u,l);
			l.x = l.x*cos((-5.0*pi/180.0)) + uxl.x * sin((-5.0*pi/180.0));
		    l.y = l.y*cos((-5.0*pi/180.0)) + uxl.y * sin((-5.0*pi/180.0));
		    l.z = l.z*cos((-5.0*pi/180.0)) + uxl.z * sin((-5.0*pi/180.0));
            break;
        case '3':
		    struct point rxu;
		    rxu = crossMultiply(r,u);
		    u.x = u.x*cos((5.0*pi/180.0)) + rxu.x * sin((5.0*pi/180.0));
		    u.y = u.y*cos((5.0*pi/180.0)) + rxu.y * sin((5.0*pi/180.0));
		    u.z = u.z*cos((5.0*pi/180.0)) + rxu.z * sin((5.0*pi/180.0));
            struct point rxl;
            rxl = crossMultiply(r,l);
			l.x = l.x*cos((5.0*pi/180.0)) + rxl.x * sin((5.0*pi/180.0));
		    l.y = l.y*cos((5.0*pi/180.0)) + rxl.y * sin((5.0*pi/180.0));
		    l.z = l.z*cos((5.0*pi/180.0)) + rxl.z * sin((5.0*pi/180.0));
            break;
        case '4':
		    //struct point rxu;
		    rxu = crossMultiply(r,u);
		    u.x = u.x*cos((5.0*pi/180.0)) + rxu.x * sin((-5.0*pi/180.0));
		    u.y = u.y*cos((5.0*pi/180.0)) + rxu.y * sin((-5.0*pi/180.0));
		    u.z = u.z*cos((5.0*pi/180.0)) + rxu.z * sin((-5.0*pi/180.0));
            //struct point rxl;
            rxl = crossMultiply(r,l);
			l.x = l.x*cos((5.0*pi/180.0)) + rxl.x * sin((-5.0*pi/180.0));
		    l.y = l.y*cos((5.0*pi/180.0)) + rxl.y * sin((-5.0*pi/180.0));
		    l.z = l.z*cos((5.0*pi/180.0)) + rxl.z * sin((-5.0*pi/180.0));
            break;
        case '5':
		    struct point lxu;
		    lxu = crossMultiply(l,u);
		    u.x = u.x*cos((5.0*pi/180.0)) + lxu.x * sin((5.0*pi/180.0));
		    u.y = u.y*cos((5.0*pi/180.0)) + lxu.y * sin((5.0*pi/180.0));
		    u.z = u.z*cos((5.0*pi/180.0)) + lxu.z * sin((5.0*pi/180.0));
            struct point lxr;
		    lxr = crossMultiply(l,r);
		    r.x = r.x*cos((5.0*pi/180.0)) + lxr.x * sin((5.0*pi/180.0));
		    r.y = r.y*cos((5.0*pi/180.0)) + lxr.y * sin((5.0*pi/180.0));
		    r.z = r.z*cos((5.0*pi/180.0)) + lxr.z * sin((5.0*pi/180.0));
            break;
        case '6':
		    //struct point lxu;
		    lxu = crossMultiply(l,u);
		    u.x = u.x*cos((5.0*pi/180.0)) + lxu.x * sin((-5.0*pi/180.0));
		    u.y = u.y*cos((5.0*pi/180.0)) + lxu.y * sin((-5.0*pi/180.0));
		    u.z = u.z*cos((5.0*pi/180.0)) + lxu.z * sin((-5.0*pi/180.0));
            //struct point lxr;
		    lxr = crossMultiply(l,r);
		    r.x = r.x*cos((5.0*pi/180.0)) + lxr.x * sin((-5.0*pi/180.0));
		    r.y = r.y*cos((5.0*pi/180.0)) + lxr.y * sin((-5.0*pi/180.0));
		    r.z = r.z*cos((5.0*pi/180.0)) + lxr.z * sin((-5.0*pi/180.0));
            break;

        case 'q':
            if(cylinderAngle<45.0)
                cylinderAngle = cylinderAngle + 5.0;
            break;
        case 'w':
            if(cylinderAngle>-45)
                cylinderAngle = cylinderAngle - 5.0;
            break;
        case 'e':
            if(cylinderAngleL<45.0)
                cylinderAngleL = cylinderAngleL + 5.0;
            break;
        case 'r':
            if(cylinderAngleL>-45)
                cylinderAngleL = cylinderAngleL - 5.0;
            break;
        case 'a':
            if(cylinderAngleZ<45.0)
                cylinderAngleZ = cylinderAngleZ + 5.0;
            break;
        case 's':
            if(cylinderAngleZ>-45)
                cylinderAngleZ = cylinderAngleZ - 5.0;
            break;
        case 'd':
            if(cylinderRotAngle<45.0)
                cylinderRotAngle = cylinderRotAngle + 5.0;
            break;
        case 'f':
            if(cylinderRotAngle>-45)
                cylinderRotAngle = cylinderRotAngle - 5.0;
            break;
		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_UP:		//down arrow key
			//cameraHeight -= 3.0;
			pos.x = pos.x + 3 * l.x;
			pos.y = pos.y + 3 * l.y;
			pos.z = pos.z + 3 * l.z;
			break;
		case GLUT_KEY_DOWN:		// up arrow key
			//cameraHeight += 3.0;
			pos.x = pos.x - 3 * l.x;
			pos.y = pos.y - 3 * l.y;
			pos.z = pos.z - 3 * l.z;
			break;

		case GLUT_KEY_RIGHT:
			//cameraAngle += 0.03;
			pos.x = pos.x + 3 * r.x;
			pos.y = pos.y + 3 * r.y;
			pos.z = pos.z + 3 * r.z;
			break;
		case GLUT_KEY_LEFT:
			//cameraAngle -= 0.03;
			pos.x = pos.x - 3 * r.x;
			pos.y = pos.y - 3 * r.y;
			pos.z = pos.z - 3 * r.z;
			break;

		case GLUT_KEY_PAGE_UP:
			pos.x = pos.x + 3 * u.x;
			pos.y = pos.y + 3 * u.y;
			pos.z = pos.z + 3 * u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
			pos.x = pos.x - 3 * u.x;
			pos.y = pos.y - 3 * u.y;
			pos.z = pos.z - 3 * u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void printMinisquare(){
    glColor3f(1,0,0);
    drawSquare(200);
    drawAxes();
    glutPostRedisplay();
}
void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
		    if (state==GLUT_DOWN){
                printf("clicked");
		        glutIdleFunc(printMinisquare);

		    }
            /*
		    glTranslatef(-10,0,0);
            glColor3f(1,0,0);
            drawSquare(100);*/
            //drawAxes();
        break;

		case GLUT_RIGHT_BUTTON:
			glTranslatef(0,290,0);
            glRotatef(90,1,0,0);
		    drawSquare(100);
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}
void drawCylinder(double radius,int slices,int stacks, int length)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<stacks;i++)
	{
		//h=radius*sin(((double)i/(double)stacks)*(pi/2));
		h = i*(length/stacks);
		//r=radius*cos(((double)i/(double)stacks)*(pi/2));
		r = radius;
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].y=h;
			points[i][j].z=r*cos(((double)j/(double)slices)*2*pi);
		}
	}
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		glColor3f(1,0,0);
		for(j=0;j<slices;j++)
		{
		    if(j%2)
                glColor3f(1,1,1);
            else
                glColor3f(0,0,0);
		    glBegin(GL_QUADS);{
                //right hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

                    //left hemisphere
                    glVertex3f(points[i][j].x,-points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,-points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,-points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,-points[i+1][j].y,points[i+1][j].z);
                }
			}glEnd();
		}
	}

void drawOuterCylinder(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<stacks/3;i++)
	{
		h=radius*sin(3*((double)i/(double)stacks)*(pi/2));
		r=radius + ((double)i/4.0)* radius*sin((3*(double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=h;
			points[i][j].z=r*sin(((double)j/(double)slices)*2*pi);
		}
	}
	for(i=0;i<stacks/4;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		//glColor3f(1,0,0);
		for(j=0;j<slices;j++)
		{
		    if(j%2==0)
                glColor3f(1,1,1);
            else
                glColor3f(0,0,0);
		    glBegin(GL_QUADS);{
                //right hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}


void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

	//now give three info
	//1. where is the camera (viewer)?
	//2. where is the camera looking?
	//3. Which direction is the camera's UP direction?

	//gluLookAt(100,100,100,	0,0,0,	0,0,1);
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);;


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	//drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);
    //glColor3f(1,0,0);
	//drawSphere(30,20,10);
	glPushMatrix();{
        glRotatef(cylinderAngle,0,0,1);
	drawHorFSphere(30,60,20,2);//drawing the left-half sphere

    glRotatef(cylinderAngleL,1,0,0);
	drawHorFSphere(30,60,20,1);//drawing the right-half sphere

	//cylinder start
	glTranslatef(0,40,0);
	glRotatef(cylinderAngleZ,1,0,0);
	glRotatef(cylinderRotAngle,0,1,0);
	drawHorFSphere(10,60,20,2);//cylinder half sphere
    glTranslatef(0,45,0);
    drawCylinder(10,60,25, 50);//cylinder cylinder
    glTranslatef(0,48,0);
    drawOuterCylinder(10,60,20);
	}
	glPopMatrix();

    glPushMatrix();{
        glTranslatef(0,300,0);
        glRotatef(90,1,0,0);
        glColor3f(.7,0.75,0.7);
        drawSquare(100);
    }
    glPopMatrix();
    //glPushMatrix();
    //glTranslatef(0,-2,0);

    //drawAxes();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
	cylinderAngle = 0;
	cylinderAngleL = 0;
	cylinderAngleZ = 0;
	cylinderRotAngle = 0;

	//clear the screen
	glClearColor(0,0,0,0);

	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);

	//initialize the matrix
	glLoadIdentity();

	//give PERSPECTIVE parameters
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
