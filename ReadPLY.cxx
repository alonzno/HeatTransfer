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

#include <vtkFloatArray.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkCommand.h>

#include <vector>
#include <iostream>
#include <stdlib.h>
#include <time.h>

#include "sim.hxx"

#define HT_SIZE 128

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

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    #if VTK_MAJOR_VERSION <= 5
    normalGenerator->SetInput(pd);
    #else
    normalGenerator->SetInputData(pd);
    #endif
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOff();
    normalGenerator->Update();

    pd = normalGenerator->GetOutput();

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
        
        tris[idx].fieldValue[0] = 100.0;  //TODO FIX THIS;
        pt = pts->GetPoint(ptIds[1]);
        
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        tris[idx].fieldValue[1] = 100.0;  //TODO FIX THIS;
        pt = pts->GetPoint(ptIds[2]);
        
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
        tris[idx].fieldValue[2] = 100.0;  //TODO FIX THIS;
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

    void SetUpTexture(){
	GLubyte Texture3[15] = {
	130, 78, 160,
	55, 181, 74,
	66, 190, 216,
	245, 237, 46,
	238, 56, 35,
	};
        glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB, 5, 0, GL_RGB, GL_UNSIGNED_BYTE, Texture3);
        glEnable(GL_COLOR_MATERIAL);
	glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        initialized = true;

}
    virtual void RenderPiece(vtkRenderer *ren, vtkActor *act) {
        RemoveVTKOpenGLStateSideEffects();
        SetupLight();

        
        cerr << "Rendering" << endl;

	if (!initialized){
		SetUpTexture();
	}
	glEnable(GL_TEXTURE_1D);
	float ambient[3] = {1, 1, 1};
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambient);
        glBegin(GL_TRIANGLES);
        {   
            
            for (Triangle t: tris) {
                for (int i = 0; i < 3; i++) {                    
                    glNormal3f(t.normals[i][0], t.normals[i][1], t.normals[i][2]);
		    glTexCoord1d((t.fieldValue[i]+100.0)/200.0); 
                    glVertex3f(t.X[i], t.Y[i], t.Z[i]);
                }       
            }       
        }     
        glEnd();
	glDisable(GL_TEXTURE_1D);

    }
};

vtkStandardNewMacro(vtkBunnyMapper);

HeatTransfer ht(HT_SIZE, 0.0);

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
            for (int j = 0; j < mapper -> tris.size(); j++) {

                for (int i = 0; i < 3; i++) {
                    d[0] = mapper -> tris[j].X[i] + 0;
                    d[1] = mapper -> tris[j].Y[i] + 0;
                    d[2] = mapper -> tris[j].Z[i] + 0;
                    
                    movePt(d);
                    
                    a = ht.getHeat(d[0], d[1], d[2]);
                    if (isnan(a)) {
                        exit(1);
                    }

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
            
            UpdateMapper( (vtkBunnyMapper *) actor -> GetMapper());

            actor -> GetMapper() -> Modified();
            actor -> GetMapper() -> Update();
            iren -> GetRenderWindow() -> Render();

        }
  
    private:
        int TimerCount;
    public:
        vtkActor *actor;

};

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

  sphere->Update();

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

  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New(); 
  renderWindowInteractor->SetInteractorStyle(style);

  renderWindowInteractor->Initialize();

  vtkSmartPointer<vtkTimerCallback1> cb =
      vtkSmartPointer<vtkTimerCallback1>::New();
  cb->actor = actor;
    
    vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
    renderer->SetActiveCamera(camera);
    
  renderWindowInteractor -> AddObserver(vtkCommand::TimerEvent, cb);
    
  int timerId = renderWindowInteractor -> CreateRepeatingTimer(100);
  std::cout << "timerId: " << timerId << std::endl;

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
