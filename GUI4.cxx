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
  
  // create a window to make it stereo capable and give it to QVTKWidget
  vtkRenderWindow* renwin = vtkRenderWindow::New();
  vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
  vtkActor* actor = vtkActor::New();
  
  // renwin->StereoCapableWindowOn();

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
  line->SetColor(0, 255, 0, 255);
  line->SetWidth(1.0);

  // QVTKInteractor *iren1=qVTK1->GetInteractor();
  // iren1->SetRenderWindow(view->GetRenderWindow());
  // std::cout<<"Set Render window"<< endl;
  // iren1->Initialize();
  // iren1->Start();
  qVTK1->SetRenderWindow(view->GetRenderWindow());

  // create a window to make it stereo capable and give it to QVTKWidget
  renwin = vtkRenderWindow::New();
  renwin->StereoCapableWindowOn();

  qVTK2->SetUseTDx(true);
  qVTK2->SetRenderWindow(renwin);
  renwin->Delete();

  QVTKInteractor *iren2=qVTK2->GetInteractor();
  vtkInteractorStyle *s2=
    static_cast<vtkInteractorStyle *>(iren2->GetInteractorStyle());
  
  const double angleSensitivity=0.02;
  const double translationSensitivity=0.001;

  vtkTDxInteractorStyle *t2=s2->GetTDxStyle();
  t2->GetSettings()->SetAngleSensitivity(angleSensitivity);
  t2->GetSettings()->SetTranslationXSensitivity(translationSensitivity);
  t2->GetSettings()->SetTranslationYSensitivity(translationSensitivity);
  t2->GetSettings()->SetTranslationZSensitivity(translationSensitivity);

  // add a renderer
  Ren2 = vtkRenderer::New();
  qVTK2->GetRenderWindow()->AddRenderer(Ren2);
  
  // put sphere in other window
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
  contours = vtkContourFilter::New();
  contours->SetInput(ureader->GetOutput());
  contours->SetValue(1, (dminmax[1] - dminmax[0])/2);
  vtkPolyDataMapper *contMapper = vtkPolyDataMapper::New();
  contMapper->SetInput(contours->GetOutput());
  contMapper->SetScalarRange(0.0, 1.2);
  // create an actor for the contours
  vtkActor *contActor = vtkActor::New();
  contActor->SetMapper(contMapper);
  Ren2->AddViewProp(contActor);
  Ren2->SetBackground(1,1,1);
  contActor->Delete();
  contMapper->Delete();

  Connections = vtkEventQtSlotConnect::New();

  // Connections->Connect(qVTK1->GetRenderWindow()->GetInteractor(), 
  //                      vtkCommand::EnterEvent,
  //                      radio1, 
  //                      SLOT(animateClick()));
  
  // // connect window enter event to radio button slot
  // Connections->Connect(qVTK2->GetRenderWindow()->GetInteractor(), 
  //                      vtkCommand::EnterEvent,
  //                      radio2, 
  //                      SLOT(animateClick()));
  
  // update coords as we move through the window
  Connections->Connect(qVTK1->GetRenderWindow()->GetInteractor(),
                       vtkCommand::MouseMoveEvent,
                       this,
                       SLOT(updateCoords(vtkObject*)));
  
  // update coords as we move through the window
  Connections->Connect(qVTK2->GetRenderWindow()->GetInteractor(),
                       vtkCommand::MouseMoveEvent,
                       this,
                       SLOT(updateCoords(vtkObject*)));
  // Connection of slider widget with ISO surface value.
  
  // connect(this->horizontalSlider, SIGNAL(valueChanged(int)), this, SLOT(SetIsoValue(int)));

  Connections->PrintSelf(cout, vtkIndent());
}

GUI4::~GUI4()
{
  // Ren1->Delete();
  Ren2->Delete();

  Connections->Delete();
}

void GUI4::SetIsoValue(int isoValue)
{
    contours->SetValue(1, isoValue);
}

void GUI4::updateCoords(vtkObject* obj)
{
  // get interactor
  vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(obj);
  // get event position
  int event_pos[2];
  iren->GetEventPosition(event_pos);
  // update label
  QString str;
  str.sprintf("x=%d : y=%d", event_pos[0], event_pos[1]);
  coord->setText(str);
}

void GUI4::popup(vtkObject * obj, unsigned long, 
           void * client_data, void *,
           vtkCommand * command)
{
  // A note about context menus in Qt and the QVTKWidget
  // You may find it easy to just do context menus on right button up,
  // due to the event proxy mechanism in place.
  
  // That usually works, except in some cases.
  // One case is where you capture context menu events that 
  // child windows don't process.  You could end up with a second 
  // context menu after the first one.

  // See QVTKWidget::ContextMenuEvent enum which was added after the 
  // writing of this example.

  // get interactor
  vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(obj);
  // consume event so the interactor style doesn't get it
  command->AbortFlagOn();
  // get popup menu
  QMenu* popupMenu = static_cast<QMenu*>(client_data);
  // get event location
  int* sz = iren->GetSize();
  int* position = iren->GetEventPosition();
  // remember to flip y
  QPoint pt = QPoint(position[0], sz[1]-position[1]);
  // map to global
  QPoint global_pt = popupMenu->parentWidget()->mapToGlobal(pt);
  // show popup menu at global point
  popupMenu->popup(global_pt);
}

void GUI4::color1(QAction* color)
{
  // if(color->text() == "Background White")
  //   Ren1->SetBackground(1,1,1);
  // else if(color->text() == "Background Black")
  //   Ren1->SetBackground(0,0,0);
  // else if(color->text() == "Stereo Rendering")
  // {
  //   Ren1->GetRenderWindow()->SetStereoRender(!Ren1->GetRenderWindow()->GetStereoRender());
  // }
  qVTK1->update();
}

void GUI4::color2(QAction* color)
{
  if(color->text() == "Background White")
    this->Ren2->SetBackground(1,1,1);
  else if(color->text() == "Background Black")
    this->Ren2->SetBackground(0,0,0);
  else if(color->text() == "Stereo Rendering")
  {
    this->Ren2->GetRenderWindow()->SetStereoRender(!this->Ren2->GetRenderWindow()->GetStereoRender());
  }
  qVTK2->update();
}

