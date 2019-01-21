/*=========================================================================

			Project 3 for CIS 410 (W18)
			ApplyBlueHotColorMap, ApplyDifferenceColorMap,
			ApplyHSVColorMap, EvaluateFieldAtLocation, and 
			Interpolation functions implemented by 
			Jacob Brown 1/24/2018

===========================================================================*/
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSetReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkFloatArray.h>

// ****************************************************************************
//  Function: GetNumberOfPoints
//
//  Arguments:
//     dims: an array of size 3 with the number of points in X, Y, and Z.
//           2D data sets would have Z=1
//
//  Returns:  the number of points in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfPoints(const int *dims)
{
    // 3D
    //return dims[0]*dims[1]*dims[2];
    // 2D
    return dims[0]*dims[1];
}

// ****************************************************************************
//  Function: GetNumberOfCells
//
//  Arguments:
//
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the number of cells in a rectilinear mesh
//
// ****************************************************************************

int GetNumberOfCells(const int *dims)
{
    // 3D
    //return (dims[0]-1)*(dims[1]-1)*(dims[2]-1);
    // 2D
    return (dims[0]-1)*(dims[1]-1);
}


// ****************************************************************************
//  Function: GetPointIndex
//
//  Arguments:
//      idx:  the logical index of a point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1]
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the point index
//
// ****************************************************************************

int GetPointIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*dims[0]*dims[1]+idx[1]*dims[0]+idx[0];
    // 2D
    return idx[1]*dims[0]+idx[0];
}


// ****************************************************************************
//  Function: GetCellIndex
//
//  Arguments:
//      idx:  the logical index of a cell.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  the cell index
//
// ****************************************************************************

int GetCellIndex(const int *idx, const int *dims)
{
    // 3D
    //return idx[2]*(dims[0]-1)*(dims[1]-1)+idx[1]*(dims[0]-1)+idx[0];
    // 2D
    return idx[1]*(dims[0]-1)+idx[0];
}

// ****************************************************************************
//  Function: GetLogicalPointIndex
//
//  Arguments:
//      idx (output):  the logical index of the point.
//              0 <= idx[0] < dims[0]
//              1 <= idx[1] < dims[1] 
//              2 <= idx[2] < dims[2] (or always 0 if 2D)
//      pointId:  a number between 0 and (GetNumberOfPoints(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalPointIndex(int *idx, int pointId, const int *dims)
{
    // 3D
    // idx[0] = pointId%dim[0];
    // idx[1] = (pointId/dims[0])%dims[1];
    // idx[2] = pointId/(dims[0]*dims[1]);

    // 2D
    idx[0] = pointId%dims[0];
    idx[1] = pointId/dims[0];
}


// ****************************************************************************
//  Function: GetLogicalCellIndex
//
//  Arguments:
//      idx (output):  the logical index of the cell index.
//              0 <= idx[0] < dims[0]-1
//              1 <= idx[1] < dims[1]-1 
//              2 <= idx[2] < dims[2]-1 (or always 0 if 2D)
//      cellId:  a number between 0 and (GetNumberOfCells(dims)-1).
//      dims: an array of size 3 with the number of points in X, Y, and Z.
//            2D data sets would have Z=1
//
//  Returns:  None (argument idx is output)
//
// ****************************************************************************

void GetLogicalCellIndex(int *idx, int cellId, const int *dims)
{
    // 3D
    // idx[0] = cellId%(dims[0]-1);
    // idx[1] = (cellId/(dims[0]-1))%(dims[1]-1);
    // idx[2] = cellId/((dims[0]-1)*(dims[1]-1));

    // 2D
    idx[0] = cellId%(dims[0]-1);
    idx[1] = cellId/(dims[0]-1);
}

// Interpolation Function
float Interpolation(float value, float a, float b, float fa, float fb) {
	return (b - a) != 0 ? fa + ((value - a) / (b - a)) * (fb - fa) : 0;
}

// ****************************************************************************
//  Function: EvaluateFieldAtLocation
//
//  Arguments:
//     pt: a two-dimensional location
//     dims: an array of size two.  
//              The first number is the size of the array in argument X, 
//              the second the size of Y.
//     X: an array (size is specified by dims).  
//              This contains the X locations of a rectilinear mesh.
//     Y: an array (size is specified by dims).  
//              This contains the Y locations of a rectilinear mesh.
//     F: a scalar field defined on the mesh.  Its size is dims[0]*dims[1].
//
//   Returns: the interpolated field value. 0 if the location is out of bounds.
//
// ****************************************************************************

float EvaluateFieldAtLocation(const float *pt, const int *dims, const float *X, 
                              const float *Y, const float *F){
	int bottomLeftPointGet[2];
	int bottomRightPointGet[2];
	int topLeftPointGet[2];
	int topRightPointGet[2];

	// Checks to see if values are outside arrays before searching
	if (pt[0] < X[0] || pt[0] >= X[dims[0] - 1] || pt[1] < Y[0] || pt[1] >= Y[dims[1] - 1])
		return 0;

	// Run through X array, seeing if its between two values
	for (int i = 0; i < dims[0] - 1; i++) {
		if (X[i] <= pt[0] && X[i + 1] > pt[0]) {
			bottomLeftPointGet[0] = i;
			break;
		}
	}

	// Run through Y-array, seeing if its between two values
	for (int i = 0; i < dims[1] - 1; i++) {
		if (Y[i] <= pt[1] && Y[i + 1] > pt[1]) {
			bottomLeftPointGet[1] = i;
			break;
		}
	}

	// Uses logical point index to find the other indices of the cell
	bottomRightPointGet[0] = bottomLeftPointGet[0] + 1;
	bottomRightPointGet[1] = bottomLeftPointGet[1];

	topLeftPointGet[0] = bottomLeftPointGet[0];
	topLeftPointGet[1] = bottomLeftPointGet[1] + 1;

	topRightPointGet[0] = bottomLeftPointGet[0] + 1;
	topRightPointGet[1] = bottomLeftPointGet[1] + 1;

	// Uses indices to find all four F-values of cell
	float bottom_left_f = F[GetPointIndex(bottomLeftPointGet, dims)];
	float bottom_right_f = F[GetPointIndex(bottomRightPointGet, dims)];
	float top_left_f = F[GetPointIndex(topLeftPointGet, dims)];
	float top_right_f = F[GetPointIndex(topRightPointGet, dims)];

	// Interpolates between bottom left & right and top left & right
	float bottom_interp = Interpolation(pt[0], X[bottomLeftPointGet[0]], X[bottomRightPointGet[0]], bottom_left_f, bottom_right_f);
	float top_interp = Interpolation(pt[0], X[topLeftPointGet[0]], X[topRightPointGet[0]], top_left_f, top_right_f);

	// Returns final interpolation
	return Interpolation(pt[1], Y[bottomLeftPointGet[1]], Y[topLeftPointGet[1]], bottom_interp, top_interp);
}


void
WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
    full_filename += ".png";
    vtkPNGWriter *writer = vtkPNGWriter::New();
    writer->SetInputData(img);
    writer->SetFileName(full_filename.c_str());
    writer->Write();
    writer->Delete();
}

vtkImageData *
NewImage(int width, int height)
{
    vtkImageData *image = vtkImageData::New();
    image->SetDimensions(width, height, 1);
    //image->SetWholeExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetUpdateExtent(0, width-1, 0, height-1, 0, 0);
    //image->SetNumberOfScalarComponents(3);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
    //image->AllocateScalars();

    return image;
}

// ****************************************************************************
//  Function: ApplyBlueHotColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using the blue 
//     hot color map.
//
//     The blue hot color map has:
//        F=0: (0,0,128) 
//        F=1: (255,255,255) 
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************

void
ApplyBlueHotColorMap(float F, unsigned char *RGB){

	// Interpolates scalar value on blue hot color map
	RGB[0] = (int)Interpolation(F, 0, 1, 0, 255);
	RGB[1] = (int)Interpolation(F, 0, 1, 0, 255);
	RGB[2] = (int)Interpolation(F, 0, 1, 128, 255);
}

// ****************************************************************************
//  Function: ApplyDifferenceColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using a divergent colormap
//
//     The divergent color map has:
//        F=0: (0,0,128) 
//        F=0.5: (255,255,255) 
//        F=1: (128, 0, 0)
//       and smooth interpolation in between
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyDifferenceColorMap(float F, unsigned char *RGB){

	// Interpolates with different values depending
	// scalar (F == 0.5 included just for completeness)
	if (F < 0.5) {
		RGB[0] = (int)Interpolation(F, 0.5, 0, 255, 0);
		RGB[1] = (int)Interpolation(F, 0.5, 0, 255, 0);
		RGB[2] = (int)Interpolation(F, 0.5, 0, 255, 128);
	} else if (F > 0.5) {
		RGB[0] = (int)Interpolation(F, 1, 0.5, 128, 255);
		RGB[1] = (int)Interpolation(F, 1, 0.5, 0, 255);
		RGB[2] = (int)Interpolation(F, 1, 0.5, 0, 255);
	} else {
		RGB[0] = 255;
		RGB[1] = 255;
		RGB[2] = 255;
	}
}

// ****************************************************************************
//  Function: ApplyBHSVColorMap
//
//  Purpose: 
//     Maps a normalized scalar value F (0<=F<=1) to a color using an HSV rainbow colormap
//
//     The rainbow colormap uses a saturation =1.0, value = 1.0, 
//     and interpolates hue from 0 to 360 degrees 
//
//  Arguments:
//       F (input):     a scalar value between 0 and 1
//       RGB (output):  the location to store the color
//      
// ****************************************************************************
void
ApplyHSVColorMap(float F, unsigned char *RGB){

	// Interpolates hue based on scalar input
	float hue = Interpolation(F, 0, 1, 0, 360) / 60.0f;

	int i = floor(hue);
	float f = hue - i;

	// p is excluded because it'll never be anything
	// other than zero
	float q = 1.0f - f;
	float t = 1.0f - (1.0f - f);

	// Switches based on hue location in colormap
	switch (i) {
		case (0): 
			RGB[0] = 255;
			RGB[1] = (int)(t * 255);
			RGB[2] = 0;
			break;

		case(1):
			RGB[0] = (int)(q * 255);
			RGB[1] = 255;
			RGB[2] = 0;
			break;

		case(2): 
			RGB[0] = 0;
			RGB[1] = 255;
			RGB[2] = (int)(t * 255);
			break;

		case(3):
			RGB[0] = 0;
			RGB[1] = (int)(q * 255);
			RGB[2] = 255;
			break;

		case(4):
			RGB[0] = (int)(t * 255);
			RGB[1] = 0;
			RGB[2] = 255;
			break;

		case(5):
			RGB[0] = 255;
			RGB[1] = 0;
			RGB[2] = (int)(q * 255);
			break;
	}
}

int main()
{
    int  i, j;

    vtkDataSetReader *rdr = vtkDataSetReader::New();
    rdr->SetFileName("proj3_data.vtk");
    rdr->Update();

    int dims[3];
    vtkRectilinearGrid *rgrid = (vtkRectilinearGrid *) rdr->GetOutput();
    rgrid->GetDimensions(dims);

    float *X = (float *) rgrid->GetXCoordinates()->GetVoidPointer(0);
    float *Y = (float *) rgrid->GetYCoordinates()->GetVoidPointer(0);
    float *F = (float *) rgrid->GetPointData()->GetScalars()->GetVoidPointer(0);
    
    int nx = 500;
    int ny = 500;

    vtkImageData *images[3];
    unsigned char *buffer[3];
    for (i = 0 ; i < 3 ; i++)
    {
        images[i] = NewImage(nx, ny);
        buffer[i] = (unsigned char *) images[i]->GetScalarPointer(0,0,0);
    }

    for (i = 0 ; i < 3*nx*ny ; i++)
        for (j = 0 ; j < 3 ; j++)
            buffer[j][i] = 0;

	// Change back to nx, ny
    for (i = 0 ; i < nx; i++)
        for (j = 0 ; j < ny; j++)
        {
            // ITERATE OVER PIXELS
			float pt[2];
			pt[0] = Interpolation(i, 0, nx - 1, -9, 9);
			pt[1] = Interpolation(j, 0, ny - 1, -9, 9);

			float f = EvaluateFieldAtLocation(pt, dims, X, Y, F);
			float normalizedF = (f - 1.2) / 3.82;

            // I TAKE OVER HERE
            int offset = 3*(j*nx+i);
            ApplyBlueHotColorMap(normalizedF, buffer[0]+offset);
			ApplyDifferenceColorMap(normalizedF, buffer[1] + offset);
			ApplyHSVColorMap(normalizedF, buffer[2] + offset);
        }

    WriteImage(images[0], "bluehot");
    WriteImage(images[1], "difference");
    WriteImage(images[2], "hsv");
}
