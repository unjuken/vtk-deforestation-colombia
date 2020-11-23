#include <vtkSmartPointer.h>

#include <vtkActor.h>
#include <vtkDelaunay2D.h>
#include <vtkLookupTable.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPNGReader.h>
#include "vtkExtractSurface.h"
#include <vtkImageDataGeometryFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkColorTransferFunction.h>
#include <vtkProbeFilter.h>
#include "vtkWarpScalar.h"

// For compatibility with new VTK generic data arrays
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif


int main(int argc, char* argv[])
{
    std::string path = "/Users/aramirez/Documents/SciViz/Proyecto/ColoredElevationMap/COL Height Map (ASTER 30m).png";
    // Read the image
    vtkSmartPointer<vtkPNGReader> reader = vtkSmartPointer<vtkPNGReader>::New();
    reader->SetFileName(path.c_str());
    reader->Update();


    //geometry filter
    vtkSmartPointer<vtkImageDataGeometryFilter> geometry = vtkSmartPointer<vtkImageDataGeometryFilter>::New();
    geometry->SetInputConnection(reader->GetOutputPort());
    geometry->Update();

    vtkSmartPointer<vtkWarpScalar> warp = vtkSmartPointer<vtkWarpScalar>::New();
    warp->SetInputConnection(geometry->GetOutputPort());
    warp->SetScaleFactor(0.0009743095565169255);
    warp->Update();


  // Triangulate the grid points
  /*vtkSmartPointer<vtkDelaunay2D> delaunay =
    vtkSmartPointer<vtkDelaunay2D>::New();
  delaunay->SetInputData(warp->GetOutput());
  delaunay->Update();*/

    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter =vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputConnection(warp->GetOutputPort());
    smoothFilter->SetNumberOfIterations(15);
    smoothFilter->SetRelaxationFactor(0.001);
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();

    // Update normals on newly smoothed polydata
    vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalGenerator->SetInputConnection(smoothFilter->GetOutputPort());
    normalGenerator->ComputePointNormalsOn();
    normalGenerator->ComputeCellNormalsOn();
    normalGenerator->Update();

  vtkPolyData* outputPolyData = normalGenerator->GetOutput();

  double bounds[6];
  outputPolyData->GetBounds(bounds);

  // Find min and max z
  double minz = bounds[4];
  double maxz = bounds[5];

  std::cout << "minz: " << minz << std::endl;
  std::cout << "maxz: " << maxz << std::endl;

  // Create the color map
  vtkSmartPointer<vtkLookupTable> colorLookupTable =
    vtkSmartPointer<vtkLookupTable>::New();
  //colorLookupTable->SetTableRange(minz, maxz);
  //colorLookupTable->Build();
    auto tableSize = (maxz-minz)*10;
  //MakeLUTFromCTF(, colorLookupTable);

    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->SetColorSpaceToDiverging();
    ctf->AddRGBPoint(0.0, 0.865003, 0.865003, 0.865003);
    ctf->AddRGBPoint(tableSize/50, 0.333333, 0.478431, 0.164706);
    ctf->AddRGBPoint(tableSize, 0.705882, 0.0156863, 0.14902);

  // Generate the colors for each point based on the color map
  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  std::cout << "There are " << outputPolyData->GetNumberOfPoints()
            << " points." << std::endl;

  for(int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
  {
    double p[3];
    outputPolyData->GetPoint(i,p);

    double dcolor[3];

    ctf->GetColor(p[2]*10, dcolor);
    /*std::cout << "dcolor: "
              << dcolor[0] << " "
              << dcolor[1] << " "
              << dcolor[2] << std::endl;*/
    unsigned char color[3];
    for(unsigned int j = 0; j < 3; j++)
    {
      color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
    }

    if(p[2] < 0.1)
    {
        color[0] = 9;
        color[1] = 46;
        color[2] = 93;
    }

    colors->InsertNextTupleValue(color);
  }

  outputPolyData->GetPointData()->SetScalars(colors);


  // Create a mapper and actor
  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData(outputPolyData);

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // Create a renderer, render window, and interactor
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  // Add the actor to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(.1, .2, .3);

  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
