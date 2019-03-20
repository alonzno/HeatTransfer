//#if !(__has_feature(cxx_variadic_templates))
//#define _LIBCPP_HAS_NO_VARIADICS
//#endif

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
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkObjectBase.h>

#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkCommand.h>

#include <vector>
#include <iostream>
#include <stdlib.h>
//#include <unistd.h>
#include <time.h>

#include "sim.hxx"

#define HT_SIZE 256

using namespace std;

void movePt(double *pt) {
    pt[0] = (pt[0] - (-0.09468989819288253784)) / (0.06100910156965255737 - (-0.09468989819288253784)) * (HT_SIZE-1);
    pt[1] = (pt[1] - (0.03298740088939666748)) / (0.187321007251739502 - (0.03298740088939666748)) * (HT_SIZE-1);
    pt[2] = (pt[2] - (-0.06187359988689422607)) / (0.05879969894886016846 - (-0.06187359988689422607)) * (HT_SIZE-1);
}

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
    
//    double minX = INFINITY, maxX = -INFINITY, minY = INFINITY, maxY = -INFINITY, minZ = INFINITY, maxZ = -INFINITY;
    int idx;
    for (idx = 0, cells->InitTraversal(); cells->GetNextCell(npts, ptIds); idx++) {
        if (npts != 3) {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        
//        movePt(pt);
        
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        
//        if (pt[0] < minX) minX = pt[0];
//        if (pt[1] < minY) minY = pt[1];
//        if (pt[2] < minZ) minZ = pt[2];
//        if (pt[0] > maxX) maxX = pt[0];
//        if (pt[1] > maxY) maxY = pt[1];
//        if (pt[2] > maxZ) maxZ = pt[2];
        
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        
        tris[idx].fieldValue[0] = 1.0;  //TODO FIX THIS;
        pt = pts->GetPoint(ptIds[1]);
        
//        movePt(pt);
        
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        
//        if (pt[0] < minX) minX = pt[0];
//        if (pt[1] < minY) minY = pt[1];
//        if (pt[2] < minZ) minZ = pt[2];
//        if (pt[0] > maxX) maxX = pt[0];
//        if (pt[1] > maxY) maxY = pt[1];
//        if (pt[2] > maxZ) maxZ = pt[2];
        
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        tris[idx].fieldValue[1] = 1.0;  //TODO FIX THIS;
        pt = pts->GetPoint(ptIds[2]);
        
        //movePt(pt);
        
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        
//        if (pt[0] < minX) minX = pt[0];
//        if (pt[1] < minY) minY = pt[1];
//        if (pt[2] < minZ) minZ = pt[2];
//        if (pt[0] > maxX) maxX = pt[0];
//        if (pt[1] > maxY) maxY = pt[1];
//        if (pt[2] > maxZ) maxZ = pt[2];
        
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
        tris[idx].fieldValue[2] = 1.0;  //TODO FIX THIS;
    }
    
//    cerr << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << minX << " " << maxX << endl;
//    cerr << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << minY << " " << maxY << endl;
//    cerr << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << minZ << " " << maxZ << endl;


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
     float INFINITYo[4] = { 0, 0, 0, 1 };
     glLightModelfv(GL_LIGHT_MODEL_AMBIENT, INFINITYo);
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

        
        cerr << "Rendering" << endl;

        glEnable(GL_COLOR_MATERIAL);
        glDisable(GL_TEXTURE_1D);
        
        // Write all the triangles
        glBegin(GL_TRIANGLES);
        {   
            unsigned char *buff;
//            time_t timer;
//            time(&timer);
//            //int rv = std::rand();
//
//            cerr << "Time: " << timer << endl;
//            cerr << "RV: " << rv << endl;
            
            for (Triangle t: tris) {
                for (int i = 0; i < 3; i++) {
                    //cerr << t.fieldValue[i] << endl;
                    glColor3d(t.fieldValue[i]/100, t.fieldValue[i]/100, t.fieldValue[i]/100); //TODO Revise thi
                    glNormal3f(t.normals[i][0], t.normals[i][1], t.normals[i][2]);
                    glVertex3f(t.X[i], t.Y[i], t.Z[i]);
                }       
            }       
        }     
        glEnd();

    }
};

vtkStandardNewMacro(vtkBunnyMapper);

HeatTransfer ht(256, 0.0);

class vtkTimerCallback1 : public vtkCommand 
{
    public:
        static vtkTimerCallback1 *New() {
            vtkTimerCallback1 *cb = new vtkTimerCallback1;
            cb->TimerCount = 0;
            return cb;
        }
    
        void UpdateMapper(vtkBunnyMapper *mapper) {
            float a;
            double *d = (double *) malloc(sizeof(double) * 3);
//            cerr << "Crash 2" << endl;
            for (int j = 0; j < mapper -> tris.size(); j++) {
                //cerr << "Crash 3" << endl;

                for (int i = 0; i < 3; i++) {
//                    cerr << "Crash 4" << endl;
                    d[0] = mapper -> tris[j].X[i] + 0;
                    d[1] = mapper -> tris[j].Y[i] + 0;
                    d[2] = mapper -> tris[j].Z[i] + 0;
                    
//                    cerr << d[0] << " "
//                    << d[1] << " " << d[2] << endl;

                    movePt(d);
//                    cerr << d[0] << " "
//                    << d[1] << " " << d[2] << endl;
                    
                    a = ht.getHeat(d[0], d[1], d[2]);
                    if (isnan(a)) {
                        exit(1);
                    }
//                    cerr << a << endl;
//                    cerr << "Crash 5" << endl;

                    mapper -> tris[j].fieldValue[i] = a;
                }
            }
            free(d);
        }

        virtual void Execute(vtkObject *caller, unsigned long eventId, void *vtkNotUsed(callData)) {
            if (vtkCommand::TimerEvent == eventId) {
                ++(this->TimerCount);
            }
            std::cout << this->TimerCount << std::endl;
            ht.Advance();
            
            

            vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::SafeDownCast(caller);
            
            // Advance the simulation
            
            // Replace the tempeeratures at all points
            
            // Update sim
            //cerr << "Crash 0" << endl;
            //sleep(100000);
            UpdateMapper( (vtkBunnyMapper *) actor -> GetMapper());
            //cerr << "Crash 1" << endl;

            actor -> GetMapper() -> Modified();
            actor -> GetMapper() -> Update();
            iren -> GetRenderWindow() -> Render();

        }
  
    private:
        int TimerCount;
    public:
        vtkActor *actor;

};

//vtkStandardNewMacro(vtkBunnyMapper);

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

  //sphere -> SetCenter(0.0, 0.0, 0.0);
  sphere->Update();


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

  //renderWindow->SetFullScreen(true);
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); 
  renderWindowInteractor->SetInteractorStyle(style);

  renderWindowInteractor->Initialize();

  vtkSmartPointer<vtkTimerCallback1> cb =
      vtkSmartPointer<vtkTimerCallback1>::New();
  cb->actor = actor;
    
    vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
    renderer->SetActiveCamera(camera);
    
  renderWindowInteractor -> AddObserver(vtkCommand::TimerEvent, cb);
    
  int timerId = renderWindowInteractor -> CreateRepeatingTimer(500);
  std::cout << "timerId: " << timerId << std::endl;

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
