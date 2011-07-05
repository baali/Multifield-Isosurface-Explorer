/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: Cone.cxx,v $

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vtkFloatArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStructuredPoints.h"
#include "vtkPointData.h"
#include "vtkStructuredPointsReader.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkHexahedron.h"
#include "vtkDataSetAttributes.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkByteSwap.h"
#include "vtkCellArray.h"
#include "vtkCellType.h"

#include "Util/vector.h"
#include "Util/matrix.h"
#include "Util/projective.h"

#include <fstream>
#include <sys/time.h>
#include <sstream>
#include <iomanip>

//Our OpenCL Particle Systemclass
#include "cll.h"

#define BINS 100
CL* example;

using namespace Eigen;
using namespace std;

float *minmax;
int sargc = 0;
int gotpoint = 0;
char *argsv[] = { "qhull" };

int mark = 0;

float bins[BINS + 10];
float binst[BINS + 10];
float increment;

ofstream op;

void
WriteKappa (char *filename)
{
  FILE *fp = fopen (filename, "w");
  double fCurrent = minmax[0];
  double df = (minmax[1] - minmax[0]) / (double) BINS;
  if (df > 0)
    {
      for (int i = 0; i < BINS; ++i)
	{
	  fprintf (fp, "%f, %f\n", fCurrent, bins[i]);
	  fCurrent += df;
	}
    }
  fclose (fp);

  fp = fopen ("statistics-kappa.csv", "w");
  fCurrent = minmax[0];
  df = (minmax[1] - minmax[0]) / (double) BINS;
  if (df > 0)
    {
      for (int i = 0; i < BINS; ++i)
	{
	  fprintf (fp, "%f, %f\n", fCurrent, binst[i]);
	  fCurrent += df;
	}
    }
  fclose (fp);

}

#define MODESINGLE "0"
#define MODETWO "1"
#define MODETHREE "2"
#define MODEDERIVED "3"
#define MAX_SOURCE_SIZE (0x100000)

int
main (int argc, char *argv[])
{
  //QCoreApplication app(argc, argv);

  char buf[32];
  if (argc < 2)
    {
      printf ("usage : kappa3d mode\n");
      exit (1);
    }
  if (!strcmp (MODETWO, argv[1]))
    {
      if (argc < 6)
	{
	  printf ("Usage : kappa3d %s fname gname input output\n", MODETWO);
	  exit (1);
	}
      vtkUnstructuredGridReader *reader = vtkUnstructuredGridReader::New ();
      vtkDataArray *fArray = NULL;
      vtkDataArray *gArray = NULL;
      vtkDataArray *hArray = NULL;
      reader->SetFileName (argv[4]);
      reader->ReadAllScalarsOn ();
      // printf("Reading VTK %s source: Done\n", argv[4]);
      vtkUnstructuredGrid *uGrid;
      uGrid = reader->GetOutput ();
      uGrid->Update ();
      int numFunction = reader->GetNumberOfScalarsInFile ();
      vtkPoints *points = uGrid->GetPoints ();

      fArray = uGrid->GetPointData ()->GetScalars (argv[2]);
      gArray = uGrid->GetPointData ()->GetScalars (argv[3]);
      // printf("Getting scalar values of %s and %s: Done\n", argv[2], argv[3]);
      int numCells = uGrid->GetNumberOfCells ();
      int numPoints = uGrid->GetNumberOfPoints();
      printf("Numcells %d\n", numCells);
      printf("NumPoints %d\n", numPoints);

      // vtkStructuredPointsReader *reader = vtkStructuredPointsReader::New ();
      // vtkDataArray *fArray = NULL;
      // vtkDataArray *gArray = NULL;
      // vtkDataArray *hArray = NULL;
      // reader->SetFileName (argv[4]);
      // reader->ReadAllScalarsOn ();
      // // vtkUnstructuredGrid *uGrid;
      // vtkStructuredPoints *uGrid;
      // uGrid = reader->GetOutput ();
      // uGrid->Update ();
      // int numFunction = reader->GetNumberOfScalarsInFile ();
      // // vtkStructuredUGrid *uGrid = uGrid->GetUGrid ();

      // fArray = uGrid->GetPointData ()->GetScalars (argv[2]);
      // gArray = uGrid->GetPointData ()->GetScalars (argv[3]);

      vtkDataArray *xArray = uGrid->GetPointData ()->GetScalars ("x");
      vtkDataArray *yArray = uGrid->GetPointData ()->GetScalars ("y");
      vtkDataArray *zArray = uGrid->GetPointData ()->GetScalars ("z");
      double *xlim = xArray->GetRange();
      double *ylim = yArray->GetRange();
      double *zlim = zArray->GetRange();
      int XDIM, YDIM, ZDIM;
      // XDIM = atoi(argv[6]);
      // YDIM = atoi(argv[7]);
      // ZDIM = atoi(argv[8]);
      XDIM = xlim[1] + 1;
      YDIM = ylim[1] + 1;
      ZDIM = zlim[1] + 1;
      printf("Coords %d, %d, %d\n", XDIM, YDIM, ZDIM );
      printf("cells from points %d\n", (XDIM - 1)*(YDIM - 1)*(ZDIM - 1));
      // return 0;

      // int ids[8];
      // int idtet[4];
 
      int kappaFlag = 0;

      if (gArray != NULL)
	  kappaFlag = 1;

      for (int i = 0; i < BINS + 10; ++i)
	{
	  bins[i] = 0;
	  binst[i] = 0;
	}
      // printf("Initialization of bins: Done\n");
      minmax = new float[2];
      double *dminmax;
      dminmax = fArray->GetRange ();
      minmax[0] = (float)dminmax[0];
      minmax[1] = (float)dminmax[1];
      // printf("Mimmax values %f, %f \n", minmax[0], minmax[1]);
      increment = (minmax[1] - minmax[0]) / (float) BINS;
      // printf("increment is %f\n", increment);
      double *pt;
      timeval tim;
      double t1, t2, t3;

      // OpenCL Stuff
      // printf("Hello, OpenCL\n");

      gettimeofday(&tim, NULL);
      t1=tim.tv_sec+(tim.tv_usec/1000000.0);

      CL example;

      // load and build our CL program from the file
      std::ifstream file("../pointData.cl");
      std::string source_str(
			     std::istreambuf_iterator<char>(file),
			     (std::istreambuf_iterator<char>()));
      example.loadProgram(source_str);

      gettimeofday(&tim, NULL);
      t2 = tim.tv_sec+(tim.tv_usec/1000000.0);
      printf("%.6lf seconds for building kernel\n", t2-t1);

      double totalTime = 0;
      double tInitial, tFinal;
      int p = 0;
      int step = 12;
      // numPoints = XDIM * YDIM * step;
      // numPoints = 8;
      // numCells = 1;
      numCells = ((XDIM - 1)*(YDIM - 1)*(step - 1));
      Point *point = new Point[numPoints];
      float (*binS)[110] = new float[numCells][110];
      float (*binsT)[110] = new float[numCells][110];

      // int ids[8];
      // ids[0] = 171;
      // ids[1] = ids[0]+1;
      // ids[2] = ids[1] + (XDIM) * (YDIM);
      // ids[3] = ids[0] + (XDIM) * (YDIM);
      // ids[4] = ids[0] + (XDIM);
      // ids[5] = ids[1] + (XDIM);
      // ids[6] = ids[2] + (XDIM);
      // ids[7] = ids[3] + (XDIM);

      if (increment > 0)
	{
	  // Reading all the points data.
	  for (int ptIndex = 0; ptIndex < numPoints; ptIndex++)
	    {
	      pt = uGrid->GetPoint (ptIndex);
	      // pt = uGrid->GetPoint (ids[ptIndex]);
	      // printf("%d\n", ids[ptIndex]);
	      // printf("[%f, %f, %f]\n", pt[0], pt[1], pt[2]);
	      point[ptIndex].v = Vec4(pt[0], pt[1], pt[2], 1);
	      if (fArray->GetDataType () == VTK_FLOAT)
		{
		  // point[ptIndex].f = point[ptIndex].h = fArray->GetTuple1 (ids[ptIndex]);
		  point[ptIndex].f = fArray->GetTuple1 (ptIndex);
                  point[ptIndex].h = fArray->GetTuple1 (ptIndex);
		  if (gArray != NULL)
		    // point[ptIndex].g = gArray->GetTuple1 (ids[ptIndex]);
		    point[ptIndex].g = gArray->GetTuple1 (ptIndex);
		  // printf("%f, %f\n", point[ptIndex].f, point[ptIndex].g);
		}	      	      
	    }
	  // return 0;
	  gettimeofday(&tim, NULL);
	  t3 = tim.tv_sec+(tim.tv_usec/1000000.0);
	  printf("Time spent to read all points %f\n", (t3-t2));
	  int batchCount = 0;
	  // Loop for covering batches
	  while ( batchCount <= ZDIM/step)
	  // while ( batchCount < 5)
	    {
	      numCells = ((XDIM - 1)*(YDIM - 1)*(step - 1));
	      // printf("Number of cells in this case %d\n", numCells);
	      example.loadPoints(point, 
				 kappaFlag, minmax, increment, 
				 batchCount, numCells, 
				 XDIM, YDIM, ZDIM);
	      example.loadArguments();
	      example.pointKernel(binS, binsT, bins, binst, 
				  numCells);
	      batchCount += 1;
	    }
	  // Handling rest of cells
	  numCells = ((XDIM - 1)*(YDIM - 1)*(ZDIM - 1)) - 
	    (batchCount * ((XDIM - 1)*(YDIM - 1)*(step - 1)));
	  if (numCells != 0)
	    {
	      printf("Number of cells in this case %d\n", numCells);
	      example.loadPoints(point, 
	  			 kappaFlag, minmax, increment, 
	  			 batchCount, numCells, 
	  			 XDIM, YDIM, ZDIM);
	      example.loadArguments();
	      example.pointKernel(binS, binsT, bins, binst, 
	  			  numCells);
	    }
	  // double t4 = tim.tv_sec+(tim.tv_usec/1000000.0);
	  // printf("Time on GPU for points: %f\n", (t4 - t3));

	  WriteKappa (argv[5]);
	  example.timing();
	}
    }
}
