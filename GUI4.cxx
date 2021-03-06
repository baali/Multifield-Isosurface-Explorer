/*=========================================================================

  Program:   Visualization Toolkit
  Module:    GUI4.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/*=========================================================================

  Copyright 2004 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/

/*========================================================================
 For general information about using VTK and Qt, see:
 http://www.trolltech.com/products/3rdparty/vtksupport.html
=========================================================================*/

/*
Feature List
* Set Axis names.
* Set title of Chart.
* Get MouseClick working for Charts.
* If possible try threading for handing visualization and calculation part
  separately.
 */

#include "GUI4.h"

#include <QMenu>

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkCommand.h"

#include "vtkContourFilter.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkInteractorStyle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkHexahedron.h"
#include "vtkDataSetAttributes.h"
#include "vtkUnstructuredGridReader.h"
#include "vtkFloatArray.h"
#include "vtkByteSwap.h"
#include "vtkCellArray.h"
#include "vtkCellType.h"
#include "vtkEventQtSlotConnect.h"

#include "Util/vector.h"
#include "Util/matrix.h"
#include "Util/projective.h"

#include "vtkTDxInteractorStyleCamera.h"
#include "vtkTDxInteractorStyleSettings.h"
#include "QVTKInteractor.h"
#include "vtkChartXY.h"

#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkDelimitedTextReader.h>
#include <vtkDoubleArray.h>
#include <vtkTable.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkPlot.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include "vtkObjectFactory.h"
#include <vtkAxis.h>
#include <vtkVariant.h>
#include <vtkPlotPoints.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <stdio.h>
#include <sys/time.h>
#include <fstream>
#include "cll.h"

vtkStandardNewMacro(vtkNewChart);
GUI4::GUI4()
{
  this->setupUi(this);
  ureader = vtkUnstructuredGridReader::New();
  uGrid = vtkUnstructuredGrid::New();
  chart = vtkNewChart::New();

  //Number of bins
  BINS = 100;
  bins = new float[BINS + 10];
  binst = new float[BINS + 10];

  //clears combobox
  connect(actionOpen, SIGNAL(triggered()), comboBox, SLOT(clear()));
  connect(actionOpen, SIGNAL(triggered()), comboBox_2, SLOT(clear()));

  connect(actionOpen, SIGNAL(triggered()), this, SLOT(OpenFile()));
  //Connecting slider with Slot to update IsoValue
  connect(horizontalSlider,SIGNAL(sliderReleased()),this,SLOT(SetIsoValue()));
  //Updating lineEdit
  connect(horizontalSlider,SIGNAL(valueChanged(int)), this, SLOT(SetLineEdit(int)));

  //Connecting Combobox selection with Combobox_2
  connect(comboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(DisableButton(int)));
  connect(comboBox_2,SIGNAL(currentIndexChanged(int)), this, SLOT(DisableButton(int)));

  //Updating horizontalSlider based on updating comboBox
  connect(comboBox,SIGNAL(currentIndexChanged(int)), this, SLOT(UpdateSlider(int)));

  //Connecting calculate button with OpenCL code
  connect(pushButton, SIGNAL(clicked()), this, SLOT(CalculateKappa()));

  //Connecting lineEdit text change with slider and isosurface
  connect(lineEdit, SIGNAL(returnPressed()), this, SLOT(UpdateLineEdit()));

  // create a window to make it stereo capable and give it to QVTKWidget
  vtkRenderWindow* renwin = vtkRenderWindow::New();
  // vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  // vtkActor* actor = vtkActor::New();
  
  // Set up the view
  view = vtkContextView::New();

  // Ren1 = view->GetRenderer();
  // create a window to make it stereo capable and give it to QVTKWidget
  renwin = vtkRenderWindow::New();
  renwin->StereoCapableWindowOn();

  qVTK2->SetRenderWindow(renwin);
  renwin->Delete();

  // QVTKInteractor *iren2=qVTK2->GetInteractor();

  // add a renderer
  Ren2 = vtkRenderer::New();

  // Setting up Contour filter
  contours = vtkContourFilter::New();
  // contours->UseScalarTreeOn();

  textActor = vtkTextActor::New();
  textActor->GetTextProperty()->SetFontSize ( 36 );
  textActor->SetPosition2 ( 10, 40 );
  textActor->GetTextProperty()->SetColor ( 1.0,0.0,0.0 );

  Connections = vtkEventQtSlotConnect::New();
  Connections->Connect(qVTK1->GetRenderWindow()->GetInteractor(),
		       vtkCommand::LeftButtonPressEvent,
		       this,
		       SLOT(UpdateCoords()));

  // Building the kernel
  timeval tim;
  double t1, t2;
  gettimeofday(&tim, NULL);
  t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  
  std::ifstream file("../pointData.cl");
  std::string source_str(
			 std::istreambuf_iterator<char>(file),
			 (std::istreambuf_iterator<char>()));
  example.loadProgram(source_str);
  gettimeofday(&tim, NULL);
  t2 = tim.tv_sec+(tim.tv_usec/1000000.0);
  printf("%.6lf seconds for building kernel\n", t2-t1);
}

GUI4::~GUI4()
{
  // Ren1->Delete();
  Ren2->Delete();
}

void GUI4::UpdateCoords()
{
  vtkNewChart* newchart = (vtkNewChart*)chart;
  // std::cout<<newchart->chartPos.X()<<" "<<newchart->chartPos.Y()<<endl;
  contours->SetValue(0, newchart->chartPos.X());
  horizontalSlider->setValue(newchart->chartPos.X());
  QString str;
  str.sprintf("%f", newchart->chartPos.X());
  lineEdit->setText(str);
  
  if( chart->GetNumberOfPlots() == 5)
    chart->RemovePlot(4);

  vtkSmartPointer<vtkTable> stable =
    vtkSmartPointer<vtkTable>::New();
  vtkSmartPointer<vtkFloatArray> arrX =
    vtkSmartPointer<vtkFloatArray>::New();
  arrX->SetName("X Axis");
  stable->AddColumn(arrX);
  vtkSmartPointer<vtkFloatArray> arrC =
    vtkSmartPointer<vtkFloatArray>::New();
  arrC->SetName("Cosine");
  stable->AddColumn(arrC);
  int numPoints = 1;
  stable->SetNumberOfRows(numPoints);
  stable->SetValue(0, 0, newchart->chartPos.X());
  stable->SetValue(0, 1, newchart->chartPos.Y());
  vtkPlot *Spoints = chart->AddPlot(vtkChart::POINTS);
  Spoints->SetLabel("");
  Spoints->SetInput(stable, 0, 1);
  Spoints->SetColor(0, 0, 0, 255);
  vtkPlotPoints::SafeDownCast(Spoints)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
  vtkPlotPoints::SafeDownCast(Spoints)->SetMarkerSize(12);

  qVTK2->update();
}

void GUI4::OpenFile()
{
  fileName = QFileDialog::getOpenFileName(this, tr("Open Dataset"), tr("VTK Files (*.vtk)")).toStdString();
  if(fileName == "")
    {      
      printf("Cant open the file.\n");
      return;
    }
  ureader->SetFileName(fileName.c_str());
  ureader->ReadAllScalarsOn ();
  uGrid = ureader->GetOutput();
  uGrid->Update ();
  // Setting up parameters for combobox
  for(int i = 0; i < ureader->GetNumberOfScalarsInFile(); i++)
    {
      comboBox->addItem(ureader->GetScalarsNameInFile(i));
      comboBox_2->addItem(ureader->GetScalarsNameInFile(i));
    }
  comboBox_2->setCurrentIndex(1);

  pushButton->setEnabled(TRUE);

  vtkDataArray *fArray = NULL;
  double *dminmax;
  // This all sucks big time :(
  std::cout<<"Setting filename done"<<endl;
  // Setting up vtk reader and Contour filters
  std::string scalarName = comboBox->itemText(0).toStdString();
  fArray = uGrid->GetPointData ()->GetScalars (scalarName.c_str());
  dminmax = fArray->GetRange ();

  // Setting up slider parameters
  horizontalSlider->setRange(dminmax[0], dminmax[1]);
  // horizontalSlider->setTickInterval((dminmax[1] - dminmax[0])/100);
  // horizontalSlider->setSingleStep((dminmax[1] - dminmax[0])/100);
  horizontalSlider->setValue((dminmax[1] + dminmax[0])/2);
  // std::cout<<dminmax[1]<<" "<<dminmax[0]<<" "<<(dminmax[1] + dminmax[0])/2<<endl;
  //setting initial value of TextLabel
  QString str;
  str.sprintf("%d", horizontalSlider->value());
  lineEdit->setText(str);

  contours->SetInput(ureader->GetOutput());
  contours->SetValue(0, (dminmax[1] + dminmax[0])/2);
  vtkPolyDataMapper *contMapper = vtkPolyDataMapper::New();
  contMapper->SetInput(contours->GetOutput());
  contMapper->SetScalarRange(0.0, 1.2);

  // create an actor for the contours
  contActor = vtkActor::New();
  contActor->SetMapper(contMapper);

  textActor->SetInput ( scalarName.c_str());

  Ren2->AddViewProp(contActor);
  Ren2->AddActor2D ( textActor );
  Ren2->SetBackground(1,1,1);
  contActor->Delete();
  contMapper->Delete();
  qVTK2->GetRenderWindow()->AddRenderer(Ren2);
  qVTK2->update();

  //Updating the Plot with blank.
  if(view->GetScene()->GetNumberOfItems() != 0)
    {
      view->GetScene()->ClearItems();
      // qVTK1->SetRenderWindow(view->GetRenderWindow());
      qVTK1->update();
    }
}

void GUI4::DisableButton(int index)
{
  if(comboBox->currentIndex() == comboBox_2->currentIndex())
    pushButton->setEnabled(FALSE);
  else
    pushButton->setEnabled(TRUE);
}

void GUI4::UpdateSlider(int index)
{
  vtkDataArray *fArray = NULL;
  double *dminmax;
  // This all sucks big time :(
  std::string scalarName = comboBox->currentText().toStdString();
  fArray = uGrid->GetPointData ()->GetScalars (scalarName.c_str());
  dminmax = fArray->GetRange ();

  // Setting up slider parameters
  horizontalSlider->setRange(dminmax[0], dminmax[1]);
  // horizontalSlider->setTickInterval((dminmax[1] - dminmax[0])/100);
  // horizontalSlider->setSingleStep((dminmax[1] - dminmax[0])/100);
  horizontalSlider->setValue((dminmax[1] + dminmax[0])/2);
  // std::cout<<dminmax[1]<<" "<<dminmax[0]<<" "<<(dminmax[1] + dminmax[0])/2<<endl;
  //setting initial value of TextLabel
  QString str;
  str.sprintf("%d", horizontalSlider->value());
  lineEdit->setText(str);
  //Updating the contour filter based on comboBox

  ureader->SetScalarsName (scalarName.c_str());
  contours->SetInput(ureader->GetOutput());
  contours->SetValue(0, (dminmax[1] + dminmax[0])/2);
  vtkPolyDataMapper *contMapper = vtkPolyDataMapper::New();
  contMapper->SetInput(contours->GetOutput());
  contMapper->SetScalarRange(0.0, 1.2);

  // // create an actor for the contours
  contActor = vtkActor::New();
  contActor->SetMapper(contMapper);
  textActor->SetInput ( scalarName.c_str());
  Ren2->Clear();
  Ren2->AddViewProp(contActor);
  Ren2->AddActor2D ( textActor );
  Ren2->SetBackground(1,1,1);
  // contActor->Delete();
  // contMapper->Delete();
  qVTK2->GetRenderWindow()->AddRenderer(Ren2);
  qVTK2->update();
  if(chart->GetNumberOfPlots() != 0)
    {
      chart->ClearPlots();       
      qVTK1->update();
    }
}

void GUI4::UpdateLineEdit()
{
  QString str;
  str = lineEdit->text();
  contours->SetValue(0, str.toFloat());
  horizontalSlider->setValue(str.toFloat());
  // this updating of slider changes value in line edit too :P
  lineEdit->setText(str);
  qVTK2->update();
}

void GUI4::SetLineEdit(int value)
{
  QString str;
  str.sprintf("%d", value);
  lineEdit->setText(str);
}

void GUI4::SetIsoValue()
{
  // Have to make it work in a way that Renderer can update iso-surface automatically.
  int isoValue = horizontalSlider->value();
  contours->SetValue(0, isoValue);
  qVTK2->update();
}

void GUI4::WriteKappa (char *filename)
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

  view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

  std::string inputFilename = "kappa.csv";
  vtkSmartPointer<vtkDelimitedTextReader> reader =
    vtkSmartPointer<vtkDelimitedTextReader>::New();
  reader->SetFileName(filename);
  reader->DetectNumericColumnsOn();
  reader->SetFieldDelimiterCharacters(",");
  reader->Update();
  std::cout << "done Reading file"<< endl;
  vtkTable* table = reader->GetOutput();
  if(chart->GetNumberOfPlots() != 0)
    {
      chart->ClearPlots();       
      qVTK1->update();
    }
  int numRows = table->GetNumberOfRows();
  vtkSmartPointer<vtkFloatArray> varDensity =
    vtkSmartPointer<vtkFloatArray>::New();
  varDensity->SetNumberOfValues(100);
  varDensity->SetName("Variation Density");
  table->AddColumn(varDensity);
  float max = 0;
  for(int row = 0; row < numRows; row++)
    {
      float kappa = table->GetValue(row, 1).ToFloat();
      if ( kappa > max)
  	max = kappa;
    }
  for(int row = 0; row < numRows; row++)
    {
      float kappa = table->GetValue(row, 1).ToFloat();
      table->SetValue(row, 2, (kappa/max));
    }
  table->Update();

  view->GetScene()->AddItem(chart);
  line = chart->AddPlot(vtkChart::LINE);
  line->SetInput(table, 0, 2);
  line->SetColor(255, 0, 0, 255);
  line->SetWidth(2.0);

  vtkPlot *points = chart->AddPlot(vtkChart::POINTS);
  points->SetInput(table, 0, 2);
  points->SetColor(255, 0, 0, 255);
  points->SetWidth(1.0);
  points->SetLabel("");
  vtkPlotPoints::SafeDownCast(points)->SetMarkerStyle(vtkPlotPoints::CIRCLE);

  //Reading iso-statistics file and adding it 
  // Plots
  vtkSmartPointer<vtkDelimitedTextReader> statsReader =
    vtkSmartPointer<vtkDelimitedTextReader>::New();

  inputFilename = "statistics-kappa.csv";
  statsReader->SetFileName(inputFilename.c_str());
  statsReader->DetectNumericColumnsOn();
  statsReader->SetFieldDelimiterCharacters(",");
  statsReader->Update();
  // std::cout << "done Second file"<< endl;
  vtkTable* statTable = statsReader->GetOutput();
  vtkSmartPointer<vtkFloatArray> isoStats =
    vtkSmartPointer<vtkFloatArray>::New();
  isoStats->SetNumberOfValues(100);
  isoStats->SetName("Isosurface statistics");
  statTable->AddColumn(isoStats);
  max = 0;
  for(int row = 0; row < numRows; row++)
    {
      float kappa = statTable->GetValue(row, 1).ToFloat();
      if ( kappa > max)
  	max = kappa;
    }
  for(int row = 0; row < numRows; row++)
    {
      float kappa = statTable->GetValue(row, 1).ToFloat();
      statTable->SetValue(row, 2, (kappa/max));
    }
  statTable->Update();

  view->GetScene()->AddItem(chart);
  line = chart->AddPlot(vtkChart::LINE);
  line->SetInput(statTable, 0, 2);
  line->SetColor(0, 255, 0, 255);
  line->SetWidth(2.0);
  chart->SetShowLegend(1);

  vtkPlot *statPoints = chart->AddPlot(vtkChart::POINTS);
  statPoints->SetInput(statTable, 0, 2);
  statPoints->SetColor(0, 255, 0, 255);
  statPoints->SetWidth(1.0);
  statPoints->SetLabel("");
  vtkPlotPoints::SafeDownCast(statPoints)->SetMarkerStyle(vtkPlotPoints::CIRCLE);

  // Setting Axis labels(Figure out how to get greek symbols)
  // X-Axis
  vtkAxis* axis = chart->GetAxis(1);
  axis->SetTitle("isovalues ("+
		 comboBox->currentText().toStdString()+")");
  // axis->SetNumberOfTicks(10);
  axis->SetNotation(0);
  axis->Delete();    
  // Y-Axis
  axis = chart->GetAxis(0);
  axis->SetNotation(0);
  axis->SetMaximum(1.1);
  axis->SetBehavior(1);
  axis->SetTitle("psi("+comboBox->currentText().toStdString()+
		 ",{"+comboBox->currentText().toStdString()+
		 ","+comboBox_2->currentText().toStdString()+
		 "},r)");
  // axis->SetNotation(0);

  // vtkSmartPointer<vtkLegendBoxActor> legend = 
  //   vtkSmartPointer<vtkLegendBoxActor>::New();
  // legend->SetNumberOfEntries(2);
  // double color[3] = {1,0,0};
  // legend->SetEntry(0, static_cast<vtkPolyData *> (NULL), "variation density", color);
  // vtkSmartPointer<vtkLineSource> legendSource = 
  //   vtkSmartPointer<vtkLineSource>::New();
  // color = {0, 1, 0};
  // vtkSmartPointer<vtkPolyData> legendLine = legendSource->GetOutput();
  // legend->SetEntry(1, static_cast<vtkPolyData *> (NULL), "Isovalues Statistics", color);
  // view->GetRenderer()->
  view->SetInteractor(qVTK1->GetInteractor());
  qVTK1->SetRenderWindow(view->GetRenderWindow());
  qVTK1->update();
}

void GUI4::CalculateKappa()
{
  vtkDataArray *fArray = NULL;
  vtkDataArray *gArray = NULL;
  vtkDataArray *hArray = NULL;
  // vtkPoints *points = uGrid->GetPoints ();
  std::string scalarF = comboBox->currentText().toStdString();
  std::string scalarG = comboBox_2->currentText().toStdString();;
  fArray = uGrid->GetPointData ()->GetScalars (scalarF.c_str());
  gArray = uGrid->GetPointData ()->GetScalars (scalarG.c_str());
  
  int numCells = uGrid->GetNumberOfCells ();
  printf("Number of cells %d\n", numCells);
  int numPoints = uGrid->GetNumberOfPoints();

  // int ids[8];
  // int idtet[4];
  // Vec4 (*v)[8] = new Vec4[100000][8];
  // float (*fs)[8] = new float[100000][8];
  // float (*gs)[8] = new float[100000][8];
  // float (*hs)[8] = new float[100000][8];
  int kappaFlag;
  // float (*binS)[110] = new float[100000][110];
  // float (*binsT)[110] = new float[100000][110];

  // Getting X, Y, Z dimensions of dataset
  // There should be a way by which we can specify this
  // apart from embedding x, y, z in dataset.
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

  timeval tim;
  double t1, t2, t3;

  if (gArray != NULL)
    kappaFlag = 1;
  
  for (int i = 0; i < BINS + 10; ++i)
    {
      bins[i] = 0;
      binst[i] = 0;
    }

  minmax = new float[2];
  double *dminmax;
  dminmax = fArray->GetRange ();
  minmax[0] = (float)dminmax[0];
  minmax[1] = (float)dminmax[1];
  // printf("Mimmax values %f, %f \n", minmax[0], minmax[1]);
  increment = (minmax[1] - minmax[0]) / (float) BINS;
  // printf("increment is %f\n", increment);
  double *pt;

  // double totalTime = 0;
  // double tInitial, tFinal;

  // int p = 0;
  // int step = 12;
  Point *point = new Point[numPoints];  
  // numCells = ((XDIM - 1)*(YDIM - 1)*(step - 1));
  // Reduce step size in case of bigger dimension
  // while (numCells > 300000)
  //   {
  //     step /= 2;
  //     numCells = ((XDIM - 1)*(YDIM - 1)*(step - 1));
  //     if ( step <= 1)
  // 	{
  // 	  printf("Too small step size.\n");
  //         step = 4;
  // 	  // return;
  // 	}
  //   }
  // float (*binS)[110] = new float[numCells][110];
  // float (*binsT)[110] = new float[numCells][110];

  gettimeofday(&tim, NULL);      
  t2 = tim.tv_sec+(tim.tv_usec/1000000.0);
  if (increment > 0)
    {
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
      numCells = 100000;
      // while ( batchCount <= ZDIM/step)
      while (batchCount < ((XDIM - 1)*(YDIM - 1)*(ZDIM - 1)/numCells))
	// while ( batchCount < 5)
	{
	  // numCells = ((XDIM - 1)*(YDIM - 1)*(step - 1));
	  // printf("Number of cells in this case %d\n", numCells);
	  example.loadPoints(point, 
			     kappaFlag, minmax, increment, 
			     batchCount, numCells, 
			     XDIM, YDIM, ZDIM);
	  example.loadArguments();
	  // example.pointKernel(binS, binsT, bins, binst, 
	  // 		      numCells);
	  example.pointKernel(bins, binst, numCells);
	  batchCount += 1;
	}
      // Handling rest of cells
      // numCells = ((XDIM - 1)*(YDIM - 1)*(ZDIM - 1)) - 
      //   (batchCount * ((XDIM - 1)*(YDIM - 1)*(step - 1)));
      numCells  = ((XDIM - 1)*(YDIM - 1)*(ZDIM - 1)) - 
        (batchCount * 100000);
      if (numCells != 0)
	{
	  // printf("Number of cells in this case %d\n", numCells);
	  example.loadPoints(point, 
			     kappaFlag, minmax, increment, 
			     batchCount, numCells, 
			     XDIM, YDIM, ZDIM);
	  example.loadArguments();
	  // example.pointKernel(binS, binsT, bins, binst, 
	  // 		      numCells);
	  example.pointKernel(bins, binst, numCells);
	}
      WriteKappa ("kappa.csv");
      example.timing();
    }

  // Isosurface of first scalar field.
  // contours->SetInput(ureader->GetOutput());
  // contours->SetValue(0, (dminmax[1] + dminmax[0])/2);
  // vtkPolyDataMapper *contMapper = vtkPolyDataMapper::New();
  // contMapper->SetInput(contours->GetOutput());
  // contMapper->SetScalarRange(0.0, 1.2);

  // // create an actor for the contours
  // contActor = vtkActor::New();
  // contActor->SetMapper(contMapper);
  // Ren2->AddViewProp(contActor);
  // Ren2->SetBackground(1,1,1);
  // contActor->Delete();
  // contMapper->Delete();
  // qVTK2->GetRenderWindow()->AddRenderer(Ren2);
  // qVTK2->update();
}
