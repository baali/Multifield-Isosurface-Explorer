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

/*========================================================================
 !!! WARNING for those who want to contribute code to this file.
 !!! If you use a commercial edition of Qt, you can modify this code.
 !!! If you use an open source version of Qt, you are free to modify
 !!! and use this code within the guidelines of the GPL license.
 !!! Unfortunately, you cannot contribute the changes back into this
 !!! file.  Doing so creates a conflict between the GPL and BSD-like VTK
 !!! license.
=========================================================================*/

#ifndef _GUI_h
#define _GUI_h

#include <QMainWindow>
#include "ui_GUI4.h"
#include <qslider.h>
#include <qlabel.h>
#include <qcombobox.h>
#include <qpushbutton.h>
#include <qfiledialog.h>
#include "vtkChartXY.h"
#include "vtkChartLegend.h"
#include "vtkTooltipItem.h"
#include "cll.h"

class vtkRenderer;
class vtkEventQtSlotConnect;
class vtkObject;
class vtkCommand;
class vtkContourFilter;
class vtkActor;
class vtkUnstructuredGridReader;
class vtkUnstructuredGrid;
class vtkContextView;
class vtkPlot;
class vtkTextActor;
//class vtkChartXY;

class VTK_CHARTS_EXPORT vtkNewChart: public vtkChartXY
{
 public:
  vtkVector2f chartPos;
  vtkTypeMacro(vtkNewChart, vtkChartXY);
  static vtkNewChart* New();
  void SetTooltipInfo(const vtkContextMouseEvent& mouse,
		      const vtkVector2f &plotPos,
		      vtkIdType seriesIndex, vtkPlot* plot, vtkIdType segmentIndex)
  {
    /* chartPos.X = plotPos.X(); */
    /* chartPos.Y = plotPos.Y(); */
    chartPos = vtkVector2f(plotPos.X(), plotPos.Y());
    /* std::cout<<plotPos.X()<<" "<<plotPos.Y()<<std::endl; */
    vtkChartXY::SetTooltipInfo(mouse, plotPos, seriesIndex, plot, segmentIndex);
  }
 protected:
  vtkNewChart(){
  }
};

class GUI4 : public QMainWindow, public Ui::GUI
{
  Q_OBJECT
public:
  GUI4();
  ~GUI4();
  int BINS;
  float *bins;
  float *binst;
  float *minmax;
  float increment;
  //CL example;

  void WriteKappa(char*);
      
public slots:
  void OpenFile();
  void SetIsoValue();
  void SetLineEdit(int);
  void DisableButton(int);
  void UpdateSlider(int);
  void CalculateKappa();
  void UpdateCoords();
  void UpdateLineEdit();

protected:
  vtkRenderer* Ren1;
  vtkRenderer* Ren2;
  vtkContourFilter* contours;
  vtkTextActor* textActor;
  vtkActor *contActor;
  std::string fileName;

  //load and build our CL program from the file
  CL example;

  vtkUnstructuredGridReader *ureader;
  vtkUnstructuredGrid *uGrid;
  vtkContextView *view;
  vtkChartXY *chart;
  // NewChart *chartNew;
  vtkPlot *line;
  vtkEventQtSlotConnect* Connections;
  //name more sensible
  /* CL* example; */
};

#endif // _GUI_h

