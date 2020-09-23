#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QSlider>
#include <QSpinBox>
#include <QTableView>
#include <QScatterSeries>
#include <QTableWidget>

#include <QtCharts>

/*
#include <QVTKOpenGLNativeWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkDataSetMapper.h>
#include <vtkCamera.h>
#include <vtkProp3DCollection.h>
#include <vtkNamedColors.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkProperty.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkLabeledDataMapper.h>
#include <vtkIdFilter.h>
#include <vtkCellCenters.h>
#include <vtkActor2D.h>
#include <vtkProperty2D.h>
#include <vtkPointData.h>
*/

#include <vector>

#include "model.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void showEvent( QShowEvent* event ) override;
    void UpdateGUI();

    std::vector<Model*> models;
    Model *selectedModel;


    QGridLayout *gridLayout;

    // toolbar
    QSlider *slider;                // in toolbar
    QSpinBox *spin;                 // in toolbar


    // plots
    QChartView *chartViewGraphs, *chartViewCircles;

    QLineSeries *series_tangential, *series_normal, *series_mohr;
    QScatterSeries *series_selected, *series_mohr_selected;

    QChart *chart_line, *chart_line_mohr;
/*
    // vtk
    QVTKOpenGLNativeWidget *qt_vtk_widget;
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkNamedColors> colors;

    vtkNew<vtkPoints> points;

    vtkNew<vtkUnstructuredGrid> ugrid;
    vtkNew<vtkDataSetMapper> dataSetMapper;
    vtkNew<vtkActor> actor_mesh;

    // arrows
    vtkNew<vtkPoints> points_origin;
    vtkNew<vtkUnstructuredGrid> ugrid_arrow_origin;
    vtkNew<vtkDoubleArray> arrowCoords;
    vtkNew<vtkArrowSource> arrowSource;
    vtkNew<vtkGlyph3D> glyph3D;
    vtkNew<vtkPolyDataMapper> mapper_arrows;
    vtkNew<vtkActor> actor_arrows;

    // sector numbers
    vtkNew<vtkUnstructuredGridGeometryFilter> geometryFilter;
    vtkNew<vtkLabeledDataMapper> labledDataMapper;
    vtkNew<vtkActor2D> actor_labels;
    vtkNew<vtkIdFilter> idfilter;
    vtkNew<vtkCellCenters> cellCenters;
*/
    // info
    QTableWidget *table;
    QTableWidgetItem *twi_selectedIdx, *twi_angle, *twi_normal_traction, *twi_tangential_traction;

private slots:
    void sliderValueChanged(int val);
    void spinValueChanged(int val);


private:
    Ui::MainWindow *ui;


};
#endif // MAINWINDOW_H
