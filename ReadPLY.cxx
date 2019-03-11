#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkPLYReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkImageData.h"
#include "vtkObjectFactory.h"
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataReader.h>

#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
//#include <vtkNormals.h>

#include <vector>
#include <iostream>
#include <time.h>

using namespace std;

class Triangle 
{
    public:
        double X[3];
        double Y[3];
        double Z[3];
        double fieldValue[3];

        double normals[3][3];
};

std::vector<Triangle>
GetTriangles(std::string inputFilename) {

    vtkPLYReader *reader = vtkPLYReader::New();
    reader->SetFileName ( inputFilename.c_str() );
    cerr << "Reading" << endl;
    reader->Update();
    cerr << "Done reading" << endl;
    if (reader->GetOutput()->GetNumberOfCells() == 0) {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }

    vtkPolyData *pd = reader->GetOutput();
        cerr << "Hit 53" << endl;

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
        cerr << "Hit 57" << endl;
    vtkCellArray *cells = pd->GetPolys();
    
        cerr << "Hit 61" << endl;
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    #if VTK_MAJOR_VERSION <= 5
    normalGenerator->SetInput(pd);
    #else
    normalGenerator->SetInputData(pd);
    #endif
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOff();
    normalGenerator->Update();
    /*
    // Optional settings
    normalGenerator->SetFeatureAngle(0.1);
    normalGenerator->SetSplitting(1);
    normalGenerator->SetConsistency(0);
    normalGenerator->SetAutoOrientNormals(0);
    normalGenerator->SetComputePointNormals(1);
    normalGenerator->SetComputeCellNormals(0);
    normalGenerator->SetFlipNormals(0);
    normalGenerator->SetNonManifoldTraversal(1);
    */

    pd = normalGenerator->GetOutput();
        cerr << "Hit 63" << endl;

    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);

    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal(); cells->GetNextCell(npts, ptIds); idx++) {
        if (npts != 3) {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        
        tris[idx].fieldValue[0] = 1.0;  //TODO FIX THIS;
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        tris[idx].fieldValue[1] = 1.0;  //TODO FIX THIS;
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
        tris[idx].fieldValue[2] = 1.0;  //TODO FIX THIS;
    }

    cerr << "Finished getting" << endl;

    return tris;

}

class vtkBunnyMapper : public vtkOpenGLPolyDataMapper
{
    public:
        GLuint displayList;
        bool initialized;
        std::vector<Triangle> tris;
    public:
        static vtkBunnyMapper *New();

        vtkBunnyMapper()
        {
            initialized = false;
            tris = GetTriangles("bunny/reconstruction/bun_zipper.ply");
        }
        
        void RemoveVTKOpenGLStateSideEffects()
        {
     float Info[4] = { 0, 0, 0, 1 };
     glLightModelfv(GL_LIGHT_MODEL_AMBIENT, Info);
     float ambient[4] = { 1,1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
     float diffuse[4] = { 1, 1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse);
     float specular[4] = { 1, 1, 1, 1.0 };
     glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
   }


   void SetupLight(void)
   {
       glEnable(GL_LIGHTING);
       glEnable(GL_LIGHT0);
       GLfloat diffuse0[4] = { 0.8, 0.8, 0.8, 1 };
       GLfloat ambient0[4] = { 0.2, 0.2, 0.2, 1 };
       GLfloat specular0[4] = { 0.0, 0.0, 0.0, 1 };
       GLfloat pos0[4] = { 1, 2, 3, 0 };
       glLightfv(GL_LIGHT0, GL_POSITION, pos0);
       glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse0);
       glLightfv(GL_LIGHT0, GL_AMBIENT, ambient0);
       glLightfv(GL_LIGHT0, GL_SPECULAR, specular0);
       glDisable(GL_LIGHT1);
       glDisable(GL_LIGHT2);
       glDisable(GL_LIGHT3);
       glDisable(GL_LIGHT5);
       glDisable(GL_LIGHT6);
       glDisable(GL_LIGHT7);
   }

    virtual void RenderPiece(vtkRenderer *ren, vtkActor *act) {
        RemoveVTKOpenGLStateSideEffects();
        SetupLight();

        glEnable(GL_COLOR_MATERIAL);
        glDisable(GL_TEXTURE_1D);
        
        // Write all the triangles
        glBegin(GL_TRIANGLES);
        {   
            unsigned char *buff;
            time_t timer;
            time(&timer);
            
            for (Triangle t: tris) {
                for (int i = 0; i < 3; i++) {
                    glColor3ub(timer%255, timer%255, timer%255); //TODO Revise this
                    glNormal3f(t.normals[i][0], t.normals[i][1], t.normals[i][2]);
                    glVertex3f(t.X[i], t.Y[i], t.Z[i]);
                }       
            }       
        }     
        glEnd();

    }
};

vtkStandardNewMacro(vtkBunnyMapper);

int main ( int argc, char *argv[] )
{
  if(argc != 2)
    {
    std::cout << "Usage: " << argv[0] << "  Filename(.ply)" << std::endl;
    return EXIT_FAILURE;
    }

  std::string inputFilename = argv[1];

  vtkSmartPointer<vtkSphereSource> sphere =
      vtkSmartPointer<vtkSphereSource>::New();
  sphere->SetThetaResolution(100);
  sphere->SetPhiResolution(50);

  // Visualize
  vtkSmartPointer<vtkBunnyMapper> mapper =
    vtkSmartPointer<vtkBunnyMapper>::New();
  mapper->SetInputConnection(sphere->GetOutputPort());

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();

  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);

  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  renderer->SetBackground(0.1804,0.5451,0.3412); // Sea green

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
