#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
#include "viewer/Arcball.h"                           /*  Arc Ball  Interface         */
#include "bmp/RgbImage.h"
#include "MyMesh.h"
#include <fstream>



using namespace MeshLib;

/* window width and height */
int win_width, win_height;
int gButton;
int startx, starty;
int shadeFlag = 1;
bool showMesh = true;
bool showUV = false;
bool showAxis = false;
int colorMode = 0;

/* rotation quaternion and translation vector for the object */
CQrot       ObjRot(0, 0, 1, 0);
CPoint      ObjTrans(0, 0, 0);

/* global mesh */
CMyMesh mesh;
/* arcball object */
CArcball arcball;

int textureFlag = 2;
/* texture id and image */
GLuint texName;
RgbImage image;
bool hasTexture = false;


void _disk_harmonic_map( const std::string& output) {

   CHarmonicMapper<CMyMesh> mapper( &mesh );

   mapper._iterative_map(0.0000001);

   mesh.write_m(output.c_str());
}

//copy frame buffer to an image
/*! save frame buffer to an image "snap_k.bmp"
*/
void readFrameBuffer()
{
    static int id = 0;

    GLfloat* buffer = new GLfloat[win_width * win_height * 3];
    assert(buffer);
    glReadBuffer(GL_FRONT_LEFT);
    glReadPixels(0, 0, win_width, win_height, GL_RGB, GL_FLOAT, buffer);

    RgbImage image(win_height, win_width);

    for (int i = 0; i < win_height; i++)
        for (int j = 0; j < win_width; j++)
        {
            float r = buffer[(i * win_width + j) * 3 + 0];
            float g = buffer[(i * win_width + j) * 3 + 1];
            float b = buffer[(i * win_width + j) * 3 + 2];

            image.SetRgbPixelf(i, j, r, g, b);
        }
    delete[]buffer;

    char name[256];
    std::ostringstream os(name);
    os << "snape_" << id++ << ".bmp";
    image.WriteBmpFile(os.str().c_str());

}

/*! initialize bitmap image texture */
void initializeBmpTexture()
{
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glGenTextures(1, &texName);
    glBindTexture(GL_TEXTURE_2D, texName);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,   GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    //	int ImageWidth  = image.GetNumRows();
    //	int ImageHeight = image.GetNumCols();
    int ImageWidth = image.GetNumCols();
    int ImageHeight = image.GetNumRows();
    GLubyte* ptr = (GLubyte*)image.ImageData();

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB,
        ImageWidth,
        ImageHeight,
        0,
        GL_RGB,
        GL_UNSIGNED_BYTE,
        ptr);

    if (textureFlag == 1)
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
    else if (textureFlag == 2)
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glEnable(GL_TEXTURE_2D);
}

/*! setup the object, transform from the world to the object coordinate system */
void setupObject(void)
{
    double rot[16];

    glTranslated(ObjTrans[0], ObjTrans[1], ObjTrans[2]);
    ObjRot.convert(rot);
    glMultMatrixd((GLdouble*)rot);
}

/*! the eye is always fixed at world z = +5 */
void setupEye(void) {
    glLoadIdentity();
    gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);
}

/*! setup light */
void setupLight()
{
    GLfloat lightOnePosition[4] = { 0, 0, 1, 0 };
    GLfloat lightTwoPosition[4] = { 0, 0, -1, 0 };
    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
    glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);
}

void renderStrokeFontString(
    float x,
    float y,
    float z,
    void* font,
    char* string) {

    char* c;
    glPushMatrix();
    glTranslatef(x, y, z);
    glScalef(1 / 1800.0, 1 / 1800.0, 1 / 1800.0);

    for (c = string; *c != '\0'; c++) {
        glutStrokeCharacter(font, *c);
    }

    glPopMatrix();
}

/*! draw axis */
void drawAxis()
{
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);

    glLineWidth(1.5);
    //x axis
    glColor3f(1.0, 0.0, 0.0);	//red
    glBegin(GL_LINES);
    glVertex3d(0, 0, 0);
    glVertex3d(1, 0, 0);
    glEnd();
    renderStrokeFontString(1, 0, 0, GLUT_STROKE_MONO_ROMAN, "X");

    //y axis
    glColor3f(0.0, 1.0, 0);		//green
    glBegin(GL_LINES);
    glVertex3d(0, 0, 0);
    glVertex3d(0, 1, 0);
    glEnd();
    renderStrokeFontString(0, 1, 0, GLUT_STROKE_MONO_ROMAN, "Y");

    //z axis
    glColor3f(0.0, 0.0, 1.0);	//blue
    glBegin(GL_LINES);
    glVertex3d(0, 0, 0);
    glVertex3d(0, 0, 1);
    glEnd();
    renderStrokeFontString(0, 0, 1, GLUT_STROKE_MONO_ROMAN, "Z");

    glLineWidth(1.0);
}

/*! draw mesh */
void drawMesh()
{
    glEnable(GL_LIGHTING);
    if (hasTexture)
        glBindTexture(GL_TEXTURE_2D, texName);

    for (CMyMesh::MeshFaceIterator fiter(&mesh); !fiter.end(); ++fiter)
    {
        glBegin(GL_POLYGON);
        CMyFace* pf = *fiter;
        for (CMyMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
        {
            CMyVertex* v = *fviter;
            CPoint& p = v->point();
            CPoint2& uv = v->uv();
            CPoint& rgb = v->rgb();
            CPoint n;
            switch (shadeFlag)
            {
            case 0:
                n = pf->normal();
                break;
            case 1:
                n = v->normal();
                break;
            }
            glNormal3d(n[0], n[1], n[2]);
            glTexCoord2d(uv[0], uv[1]);
            switch (colorMode) {

            case 0:
                glColor3f(rgb[0], rgb[1], rgb[2]);
                break;
            case 1:
                glColor3f(p[0], p[1], p[2]);
                break;
            case 2:
                glColor3f(n[0], n[1], n[2]);
                break;
            default:
                break;
            }
            glVertex3d(p[0], p[1], p[2]);
        }
        glEnd();
    }
}

/*! draw uv */
void drawUv()
{
    //glEnable(GL_LIGHTING);
    glDisable(GL_LIGHTING);
    if (hasTexture)
        glBindTexture(GL_TEXTURE_2D, texName);

    for (CMyMesh::MeshFaceIterator fiter(&mesh); !fiter.end(); ++fiter)
    {
        glBegin(GL_POLYGON);
        CMyFace* pf = *fiter;
        for (CMyMesh::FaceVertexIterator fviter(pf); !fviter.end(); ++fviter)
        {
            CMyVertex* v = *fviter;
            CPoint& p = v->point();
            CPoint2& uv = v->uv();
            CPoint& rgb = v->rgb();
            CPoint n;
            switch (shadeFlag)
            {
            case 0:
                n = pf->normal();
                break;
            case 1:
                n = v->normal();
                break;
            }
            glNormal3d(n[0], n[1], n[2]);
            glTexCoord2d(uv[0], uv[1]);
            switch (colorMode) {
                
                case 0:
                    glColor3f(rgb[0], rgb[1], rgb[2]);
                    break;
                case 1:
                    glColor3f(p[0], p[1], p[2]);
                    break;
                case 2:
                    glColor3f(n[0], n[1], n[2]);
                    break;
                default:
                    break;
            }
            glVertex3d(uv[0], uv[1], 0);
            
        }
        glEnd();
    }
}

void drawSharpEdges()
{
    glDisable(GL_LIGHTING);

    glLineWidth(1.5);
    glColor3f(1, 0, 0);
    glBegin(GL_LINES);
    for (CMyMesh::MeshEdgeIterator eiter(&mesh); !eiter.end(); ++eiter)
    {
        CMyEdge* pE = *eiter;
        if (pE->sharp() == true)
        {
            CMyVertex* p0 = mesh.edgeVertex1(pE);
            CMyVertex* p1 = mesh.edgeVertex2(pE);
            glColor3f(1.0f, 0.0f, 0.0f);
            glVertex3f(p0->point()[0], p0->point()[1], p0->point()[2]);
            glVertex3f(p1->point()[0], p1->point()[1], p1->point()[2]);
        }
    }
    glEnd();
    glLineWidth(1.0);
}

/*! display call back function
*/
void display()
{
    /* clear frame buffer */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    setupLight();
    /* transform from the eye coordinate system to the world system */
    setupEye();
    glPushMatrix();
    /* transform from the world to the ojbect coordinate system */
    setupObject();

    /* draw the axis */
    if (showAxis)
        drawAxis();

    /* draw sharp edges */
    drawSharpEdges();

    /* draw the mesh */
    switch (textureFlag)
    {
    case 0:
        glDisable(GL_TEXTURE_2D);
        break;
    case 1:
        glEnable(GL_TEXTURE_2D);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        break;
    case 2:
        glEnable(GL_TEXTURE_2D);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        break;
    }
    if (showMesh)
        drawMesh();
    if (showUV)
        drawUv();

    glPopMatrix();
    glutSwapBuffers();
}

/*! Called when a "resize" event is received by the window. */
void reshape(int w, int h)
{
    float ar;
    //std::cout << "w:" << w << "\th:" << h << std::endl;
    win_width = w;
    win_height = h;

    ar = (float)(w) / h;
    glViewport(0, 0, w, h);               /* Set Viewport */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // magic imageing commands
    gluPerspective(40.0, /* field of view in degrees */
        ar, /* aspect ratio */
        0.1, /* Z near */
        100.0 /* Z far */);

    glMatrixMode(GL_MODELVIEW);

    glutPostRedisplay();
}

/*! helper function to remind the user about commands, hot keys */
void help()
{
    printf("0  -  Show the coordinate axis\n");
    printf("1  -  Show the original mesh\n");
    printf("2  -  Show the uv parametrization\n");
    printf("w  -  Wireframe Display\n");
    printf("f  -  Flat Shading \n");
    printf("s  -  Smooth Shading\n");
    printf("?  -  Help Information\n");
    printf("esc - quit\n");
}

void specialKey(GLint key, GLint x, GLint y)
{
    glutPostRedisplay();
}

/*! Keyboard call back function */
void keyBoard(unsigned char key, int x, int y)
{
    switch (key)
    {
    case '0':
        showAxis = !showAxis;
        break;
    case '1':
        showMesh = !showMesh;
        break;
    case '2':
        showUV = !showUV;
        break;
    case 'c':
        colorMode = (colorMode + 1) % 3;
        break;
    case 'f':
        //Flat Shading
        glPolygonMode(GL_FRONT, GL_FILL);
        shadeFlag = 0;
        break;
    case 's':
        //Smooth Shading
        glPolygonMode(GL_FRONT, GL_FILL);
        shadeFlag = 1;
        break;
    case 'w':
        //Wireframe mode
        glPolygonMode(GL_FRONT, GL_LINE);
        break;
    case 't':
        textureFlag = (textureFlag + 1) % 3;
        break;
    case 'o':
        readFrameBuffer();
        break;
    case '?':
        help();
        break;
    case 27:
        exit(0);
        break;
    }
    glutPostRedisplay();
}

/*! setup GL states */
void setupGLstate() {
    GLfloat lightOneColor[] = { 1, 1, 1, 1.0 };
    GLfloat globalAmb[] = { .1, .1, .1, 1 };
    GLfloat lightOnePosition[] = { .0, 0.0, 1.0, 1.0 };
    GLfloat lightTwoPosition[] = { .0, 0.0, -1.0, 1.0 };

    glEnable(GL_CULL_FACE);
    glFrontFace(GL_CCW);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.35, 0.53, 0.70, 0);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_LIGHT1);
    glEnable(GL_LIGHT2);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);

    glLightfv(GL_LIGHT1, GL_DIFFUSE, lightOneColor);
    glLightfv(GL_LIGHT2, GL_DIFFUSE, lightOneColor);
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, globalAmb);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    glLightfv(GL_LIGHT1, GL_POSITION, lightOnePosition);
    glLightfv(GL_LIGHT2, GL_POSITION, lightTwoPosition);

    const GLfloat specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 64.0f);

    GLfloat mat_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f };
    GLfloat mat_diffuse[] = { 0.01f, 0.01f, 0.01f, 1.0f };
    GLfloat mat_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat mat_shininess[] = { 32 };

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
}

/*! mouse click call back function */
void  mouseClick(int button, int state, int x, int y) {
    /* set up an arcball around the Eye's center
    switch y coordinates to right handed system  */

    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        gButton = GLUT_LEFT_BUTTON;
        arcball = CArcball(win_width, win_height, x - win_width / 2, win_height - y - win_height / 2);
    }

    if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
        startx = x;
        starty = y;
        gButton = GLUT_MIDDLE_BUTTON;
    }

    if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) {
        startx = x;
        starty = y;
        gButton = GLUT_RIGHT_BUTTON;
    }
    return;
}

/*! mouse motion call back function */
void mouseMove(int x, int y)
{
    CPoint trans;
    CQrot  rot;

    /* rotation, call arcball */
    if (gButton == GLUT_LEFT_BUTTON)
    {
        rot = arcball.update(x - win_width / 2, win_height - y - win_height / 2);
        ObjRot = rot * ObjRot;
        glutPostRedisplay();
    }

    /*xy translation */
    if (gButton == GLUT_MIDDLE_BUTTON)
    {
        double scale = 10. / win_height;
        trans = CPoint(scale * (x - startx), scale * (starty - y), 0);
        startx = x;
        starty = y;
        ObjTrans = ObjTrans + trans;
        glutPostRedisplay();
    }

    /* zoom in and out */
    if (gButton == GLUT_RIGHT_BUTTON) {
        double scale = 10. / win_height;
        trans = CPoint(0, 0, scale * (starty - y));
        startx = x;
        starty = y;
        ObjTrans = ObjTrans + trans;
        glutPostRedisplay();
    }

}


/*! Normalize mesh
* \param pMesh the input mesh
*/
void normalizeMesh(CMyMesh* pMesh)
{
    CPoint s(0, 0, 0);
    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        s = s + v->point();
    }
    s = s / pMesh->numVertices();

    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        CPoint p = v->point();
        p = p - s;
        v->point() = p;
    }

    double d = 0;
    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        CPoint p = v->point();
        for (int k = 0; k < 3; k++)
        {
            d = (d > fabs(p[k])) ? d : fabs(p[k]);
        }
    }

    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        CPoint p = v->point();
        p = p / d;
        v->point() = p;
    }
};

/*! Compute the face normal and vertex normal
* \param pMesh the input mesh
*/
void computeNormal(CMyMesh* pMesh)
{
    for (CMyMesh::MeshVertexIterator viter(pMesh); !viter.end(); ++viter)
    {
        CMyVertex* v = *viter;
        CPoint n(0, 0, 0);
        for (CMyMesh::VertexFaceIterator vfiter(v); !vfiter.end(); ++vfiter)
        {
            CMyFace* pF = *vfiter;

            CPoint p[3];
            CHalfEdge* he = pF->halfedge();
            for (int k = 0; k < 3; k++)
            {
                p[k] = he->target()->point();
                he = he->he_next();
            }

            CPoint fn = (p[1] - p[0]) ^ (p[2] - p[0]);
            pF->normal() = fn / fn.norm();
            n += fn;
        }

        n = n / n.norm();
        v->normal() = n;
    }
};

void initOpenGL(int argc, char* argv[])
{
    if (hasTexture)
        image.LoadBmpFile(argv[2]);

    /* glut stuff */
    glutInit(&argc, argv);                /* Initialize GLUT */
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(600, 600);
    glutCreateWindow("Mesh Viewer");	  /* Create window with given title */
    glViewport(0, 0, 600, 600);

    glutDisplayFunc(display);             /* Set-up callback functions */
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseMove);
    glutKeyboardFunc(keyBoard);
    glutSpecialFunc(&specialKey);
    setupGLstate();

    if (hasTexture)
        initializeBmpTexture();

    glutMainLoop();                       /* Start GLUT event-processing loop */
}


/*! main function for viewer
*/
int main(int argc, char* argv[])
{
    std::string buf;
    std::ifstream ifs;
    ifs.open(argv[2], std::ios::in);
    /*
    for (int i = 0; getline(ifs, buf); i++) {
        cornerID[i] = std::stoi(buf);
    }
    if (argc < 2)
    {
        std::cout << "Usage: input.m [texture.bmp]" << std::endl;
        return -1;
    }*/

    if (argc > 2)
        hasTexture = true;


    std::string mesh_name(argv[1]);
    if (strutil::endsWith(mesh_name, ".obj"))
    {
        mesh.read_obj(mesh_name.c_str());
    }
    if (strutil::endsWith(mesh_name, ".m"))
    {
        mesh.read_m(mesh_name.c_str());
    }
    if (strutil::endsWith(mesh_name, ".off"))
    {
        mesh.read_off(mesh_name.c_str());
    }

    normalizeMesh(&mesh);
    computeNormal(&mesh);
    _disk_harmonic_map(argv[2]);
    mesh.outputMeshInfo();
    mesh.testIterator();

    mesh.geometricImage();
    initOpenGL(argc, argv);
    return 0;
}

