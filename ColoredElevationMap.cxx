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
#include <vtkXMLPolyDataWriter.h>
#include <vtkPNGReader.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkColorTransferFunction.h>
#include <vtkFieldDataToAttributeDataFilter.h>
#include <vtkMaskPoints.h>
#include <vtkGlyph3D.h>
#include <vtkImageActor.h>
#include <vtkImageMathematics.h>
#include <vtkImageShiftScale.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageGradient.h>
#include <vtkImageRGBToHSV.h>
#include <vtkImageData.h>
#include <vtkImageReader2Factory.h>
#include <vtkNamedColors.h>
#include <vtkArrowSource.h>
#include <vtkSliderRepresentation3D.h>
#include <vtkSliderWidget.h>
#include <vtkMarchingCubes.h>
#include "vtkWarpScalar.h"
#include "vtkDoubleArray.h"
#include "vtkImageMapper3D.h"

// For compatibility with new VTK generic data arrays
#ifdef vtkGenericDataArray_h
#define InsertNextTupleValue InsertNextTypedTuple
#endif

// The callback does the work.
// The callback keeps a pointer to the sphere whose resolution is
// controlled. After constructing the callback, the program sets the
// SphereSource of the callback to
// the object to be controlled.
class vtkSliderCallback : public vtkCommand
{
public:
    static vtkSliderCallback *New()
    {
        return new vtkSliderCallback;
    }
    virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
        vtkSliderWidget *sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
        for (int i = 0; i < 19; ++i) {
            renderer->RemoveActor(deforestationActors[i]);
        }
        int year = static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue();
        renderer->AddActor(deforestationActors[year]);
    }
    vtkSliderCallback():renderer(0) {}
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkActor> deforestationActors [19];
};

vtkSmartPointer<vtkActor> GetElevationActor(){
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

    return actor;
}

vtkSmartPointer<vtkActor> GetGradientActor(int year){
    auto namedColors = vtkSmartPointer<vtkNamedColors>::New();
    namedColors->SetColor("Bkg", 0.2, 0.3, 0.6);

    vtkSmartPointer<vtkImageData> originalImage;
    vtkSmartPointer<vtkImageData> image;

    int onRatio = 1;
    double scaleFactor = 2.0;
    std::string fileName = "/Users/aramirez/Documents/SciViz/Proyecto/COL Forest Change/b_" + std::to_string(year) + ".1080.png";
    // Read an image
    auto readerFactory = vtkSmartPointer<vtkImageReader2Factory>::New();
    vtkSmartPointer<vtkImageReader2> reader;
    reader.TakeReference(readerFactory->CreateImageReader2(fileName.c_str()));
    reader->SetFileName(fileName.c_str());

    // Convert to HSV and extract the Value
    auto hsvFilter = vtkSmartPointer<vtkImageRGBToHSV>::New();
    hsvFilter->SetInputConnection(reader->GetOutputPort());

    auto extractValue = vtkSmartPointer<vtkImageExtractComponents>::New();
    extractValue->SetInputConnection(hsvFilter->GetOutputPort());
    extractValue->SetComponents(2);
    extractValue->Update();

    image = reader->GetOutput();
    originalImage = reader->GetOutput();

    // Use 1% of the points
    onRatio = image->GetPointData()->GetScalars()->GetNumberOfTuples() /
              (image->GetPointData()->GetScalars()->GetNumberOfTuples() * 0.003);
    scaleFactor = 1.0;

    // Compute the gradient of the Value
    auto gradientFilter = vtkSmartPointer<vtkImageGradient>::New();
    gradientFilter->SetInputData(image);
    gradientFilter->SetDimensionality(2);
    gradientFilter->Update();

    // Extract the x component of the gradient
    auto extractXFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
    extractXFilter->SetComponents(0);
    extractXFilter->SetInputConnection(gradientFilter->GetOutputPort());

    double xRange[2];

    extractXFilter->Update();
    extractXFilter->GetOutput()->GetPointData()->GetScalars()->GetRange(xRange);

    // Gradient could be negative, so take the absolute value
    auto imageAbsX = vtkSmartPointer<vtkImageMathematics>::New();
    imageAbsX->SetOperationToAbsoluteValue();
    imageAbsX->SetInputConnection(extractXFilter->GetOutputPort());

    // Scale the output (0,255)
    auto shiftScaleX = vtkSmartPointer<vtkImageShiftScale>::New();
    shiftScaleX->SetOutputScalarTypeToUnsignedChar();
    shiftScaleX->SetScale(255 / xRange[1]);
    shiftScaleX->SetInputConnection(imageAbsX->GetOutputPort());

    // Extract the y component of the gradient
    auto extractYFilter = vtkSmartPointer<vtkImageExtractComponents>::New();
    extractYFilter->SetComponents(1);
    extractYFilter->SetInputConnection(gradientFilter->GetOutputPort());

    double yRange[2];
    extractYFilter->Update();
    extractYFilter->GetOutput()->GetPointData()->GetScalars()->GetRange(yRange);

    // Gradient could be negative, so take the absolute value
    auto imageAbsY = vtkSmartPointer<vtkImageMathematics>::New();
    imageAbsY->SetOperationToAbsoluteValue();
    imageAbsY->SetInputConnection(extractYFilter->GetOutputPort());

    // Scale the output (0,255)
    auto shiftScaleY = vtkSmartPointer<vtkImageShiftScale>::New();
    shiftScaleY->SetOutputScalarTypeToUnsignedChar();
    shiftScaleY->SetScale(255 / yRange[1]);
    shiftScaleY->SetInputConnection(imageAbsY->GetOutputPort());

    // Create the Glyphs for the gradient
    auto arrowSource = vtkSmartPointer<vtkArrowSource>::New();

    // The gradient is 2D but Glyph3D needs 3D vectors. Add a 0 z-component
    // Also, ImageGradient generates a 2-component scalar for the
    // gradient, but Glyph3D needs normals or vectors
    auto zeroes = vtkSmartPointer<vtkDoubleArray>::New();
    zeroes->SetNumberOfComponents(1);
    zeroes->SetName("Zero");
    zeroes->SetNumberOfTuples(gradientFilter->GetOutput()
                                      ->GetPointData()
                                      ->GetScalars()
                                      ->GetNumberOfTuples());
    zeroes->FillComponent(0, 0);
    gradientFilter->GetOutput()->GetPointData()->AddArray(zeroes);

    std::string scalarName =
            gradientFilter->GetOutput()->GetPointData()->GetScalars()->GetName();


    auto scalarsToVectors =
            vtkSmartPointer<vtkFieldDataToAttributeDataFilter>::New();
    scalarsToVectors->SetInputConnection(gradientFilter->GetOutputPort());
    scalarsToVectors->SetInputFieldToPointDataField();
    scalarsToVectors->SetOutputAttributeDataToPointData();
    scalarsToVectors->SetVectorComponent(0, scalarName.c_str(), 0);
    scalarsToVectors->SetVectorComponent(1, scalarName.c_str(), 1);
    scalarsToVectors->SetVectorComponent(2, "Zero", 0);

    // Select a small percentage of the gradients
    auto maskPoints = vtkSmartPointer<vtkMaskPoints>::New();
    maskPoints->SetInputConnection(scalarsToVectors->GetOutputPort());
    maskPoints->RandomModeOn();
    maskPoints->SetOnRatio(onRatio);

    auto vectorGradientGlyph = vtkSmartPointer<vtkGlyph3D>::New();
    vectorGradientGlyph->SetSourceConnection(arrowSource->GetOutputPort());
    vectorGradientGlyph->SetInputConnection(maskPoints->GetOutputPort());
    vectorGradientGlyph->SetScaleModeToScaleByVector();
    vectorGradientGlyph->SetVectorModeToUseVector();
    vectorGradientGlyph->SetScaleFactor(scaleFactor);

    // Visualize

    // (xmin, ymin, xmax, ymax)

    auto vectorGradientMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vectorGradientMapper->SetInputConnection(
            vectorGradientGlyph->GetOutputPort());
    vectorGradientMapper->ScalarVisibilityOff();

    auto vectorGradientActor = vtkSmartPointer<vtkActor>::New();
    vectorGradientActor->SetMapper(vectorGradientMapper);
    vectorGradientActor->GetProperty()->SetColor(
            namedColors->GetColor3d("red").GetData());
    vectorGradientActor->SetPosition(0,0,50);

    return vectorGradientActor;
}

void EnableSlider(vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor){


    //vtkSmartPointer<vtkSliderCallback> callback = vtkSmartPointer<vtkSliderCallback>::New();
    //callback->ModelSource = flesh;

    //sliderWidget->AddObserver(vtkCommand::InteractionEvent,callback);
}

int main(int argc, char* argv[])
{

    // Create a renderer, render window, and interactor
    vtkSmartPointer<vtkRenderer> renderer =
            vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renderWindow =
            vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(1500, 1500);
    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    renderWindowInteractor->SetRenderWindow(renderWindow);


    vtkSmartPointer<vtkActor> elevationActor = GetElevationActor();
    vtkSmartPointer<vtkActor> deforestationActors [19] = {};

    for (int i = 1; i < 20; ++i) {
        deforestationActors[i-1] = GetGradientActor(2000 + i);

    }

  // Add the actor to the scene
  renderer->AddActor(deforestationActors[0]);
  renderer->AddActor(elevationActor);
  renderer->SetBackground(.1, .2, .3);

    vtkSmartPointer<vtkSliderRepresentation3D> sliderRep = vtkSmartPointer<vtkSliderRepresentation3D>::New();
    sliderRep->SetMinimumValue(1);
    sliderRep->SetMaximumValue(19);
    sliderRep->SetValue(1);
    sliderRep->SetTitleText("Year");
    sliderRep->SetTitleHeight(0.03);
    sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToWorld();
    sliderRep->GetPoint1Coordinate()->SetValue(0,-300,0);
    sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToWorld();
    sliderRep->GetPoint2Coordinate()->SetValue(1080,-300,0);
    sliderRep->SetSliderLength(0.05);
    sliderRep->SetSliderWidth(0.05);
    sliderRep->SetEndCapLength(0.05);

    vtkSmartPointer<vtkSliderWidget> sliderWidget = vtkSmartPointer<vtkSliderWidget>::New();
    sliderWidget->SetInteractor(renderWindowInteractor);
    sliderWidget->SetRepresentation(sliderRep);
    sliderWidget->SetAnimationModeToAnimate();
    sliderWidget->On();

  // Render and interact
  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}




