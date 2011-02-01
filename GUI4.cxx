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

#include "GUI4.h"

#include <QMenu>

#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkCommand.h"
#include "vtkEventQtSlotConnect.h"
#include "vtkSphereSource.h"

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

#include "vtkTDxInteractorStyleCamera.h"
#include "vtkTDxInteractorStyleSettings.h"
#include "QVTKInteractor.h"

#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkDelimitedTextReader.h>
#include <vtkDoubleArray.h>
#include <vtkTable.h>
#include <vtkPointData.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>

GUI4::GUI4()
{
  this->setupUi(this);
  
  // Setting up vtk reader and Contour filters
  vtkUnstructuredGridReader *ureader = vtkUnstructuredGridReader::New ();
  ureader->SetFileName ("/media/sda6/isabel/01-250x250.vtk");
  ureader->ReadAllScalarsOn ();
  vtkUnstructuredGrid *uGrid;
  uGrid = ureader->GetOutput ();
  uGrid->Update ();
  vtkDataArray *fArray = NULL;
  double *dminmax;
  fArray = uGrid->GetPointData ()->GetScalars ("Pf");
  dminmax = fArray->GetRange ();

  // Setting up slider parameters
  horizontalSlider->setRange(dminmax[0], dminmax[1]);
  horizontalSlider->setTickInterval((dminmax[1] - dminmax[0])/100);
  horizontalSlider->setSingleStep((dminmax[1] - dminmax[0])/100);

  //Connecting slider with Slot to update IsoValue
  connect(horizontalSlider,SIGNAL(sliderReleased()),this,SLOT(SetIsoValue()));

  // create a window to make it stereo capable and give it to QVTKWidget
  vtkRenderWindow* renwin = vtkRenderWindow::New();
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  vtkActor* actor = vtkActor::New();
  
  // Set up the view
  vtkSmartPointer<vtkContextView> view = 
    vtkSmartPointer<vtkContextView>::New();
  view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);
  view->SetInteractor(qVTK1->GetInteractor());
  // Ren1 = view->GetRenderer();

  std::string inputFilename = "He+100-op_csv";
  vtkSmartPointer<vtkDelimitedTextReader> reader =
    vtkSmartPointer<vtkDelimitedTextReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->DetectNumericColumnsOn();
  reader->SetFieldDelimiterCharacters(",");
  reader->Update();
  std::cout << "done Reading file"<< endl;
  vtkTable* table = reader->GetOutput();

  // Add multiple line plots, setting the colors etc
  vtkSmartPointer<vtkChartXY> chart = 
    vtkSmartPointer<vtkChartXY>::New();
  view->GetScene()->AddItem(chart);
  vtkPlot *line = chart->AddPlot(vtkChart::LINE);
  line->SetInput(table, 0, 1);
  line->SetColor(255, 0, 0, 255);
  line->SetWidth(1.0);

  qVTK1->SetRenderWindow(view->GetRenderWindow());

  // create a window to make it stereo capable and give it to QVTKWidget
  renwin = vtkRenderWindow::New();
  renwin->StereoCapableWindowOn();

  qVTK2->SetUseTDx(true);
  qVTK2->SetRenderWindow(renwin);
  renwin->Delete();

  QVTKInteractor *iren2=qVTK2->GetInteractor();
  
  // add a renderer
  Ren2 = vtkRenderer::New();
  qVTK2->GetRenderWindow()->AddRenderer(Ren2);
  
  // Setting up Contour filter
  contours = vtkContourFilter::New();
  contours->SetInput(ureader->GetOutput());
  contours->SetValue(1, (dminmax[1] - dminmax[0])/2);
  vtkPolyDataMapper *contMapper = vtkPolyDataMapper::New();
  contMapper->SetInput(contours->GetOutput());
  contMapper->SetScalarRange(0.0, 1.2);

  // create an actor for the contours
  contActor = vtkActor::New();
  contActor->SetMapper(contMapper);
  Ren2->AddViewProp(contActor);
  Ren2->SetBackground(1,1,1);
  contActor->Delete();
  contMapper->Delete();

}

GUI4::~GUI4()
{
  // Ren1->Delete();
  Ren2->Delete();
}

void GUI4::SetIsoValue()
{
  // Have to make it work in a way that Renderer can update iso-surface automatically.
  int isoValue = horizontalSlider->value();
  contours->SetValue(0, isoValue);
  contActor->GetMapper()->GetInput()->Update();
  Ren2->AddViewProp(contActor);
  std::cout<<contours->GetValue(0)<<endl;
}
