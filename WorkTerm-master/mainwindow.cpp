#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <cfloat>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    Model *model = new Model;
    model->isBoundary = false;
    model->AddSector(-0.0804524,-0.126599,7.1294e-06,0.0694121,-0.132974,6.35204e-06,4.62408,1.44387,-0.0311261);
    model->AddSector(0.0694121,-0.132974,6.35204e-06,0.149865,-0.0063741,1.42622e-05,3.5138,1.1411,-0.21346);
    model->AddSector(0.149865,-0.0063741,1.42622e-05,0.0804524,0.126599,7.43639e-06,3.55001,1.10216,0.218794);
    model->AddSector(0.0804524,0.126599,7.43639e-06,-0.0694121,0.132974,6.11511e-06,4.62,1.46422,-0.0410333);
    model->AddSector(-0.0694121,0.132974,6.11511e-06,-0.149865,0.00637403,1.39815e-05,3.58378,1.18685,-0.207684);
    model->AddSector(-0.149865,0.00637403,1.39815e-05,-0.0804524,-0.126599,7.1294e-06,3.60604,1.10554,0.212817);
    model->getMax(0,2*M_PI);
    models.push_back(model);

    model = new Model;
    model->isBoundary = true;
    model->AddSector(-0.0867755,0.05456,0.000800575,-0.085338,-0.0808328,-0.00291115,18.8618,5.24858,0.604117);
    model->AddSector(-0.085338,-0.0808328,-0.00291115,0.0112856,-0.0864482,-0.00227359,10.1852,-0.517818,2.41207);
    model->AddSector(0.0112856,-0.0864482,-0.00227359,0.0867755,-0.05456,-0.000736006,1.68372,-1.84484,3.50632);
    model->getMax(0,2*M_PI);
    models.push_back(model);

    model = new Model;
    model->isBoundary = true;
    model->AddSector(0.0871197,-0.0953093,-0.00147892,0.104189,0.00711211,-0.000556284,9.46015,1.12215,-2.48601);
    model->AddSector(0.104189,0.00711211,-0.000556284,0.0421382,0.0464292,0.000203371,3.49938,-0.69808,0.171079);
    model->getMax(0,2*M_PI);
    models.push_back(model);

    model = new Model;
    model->isBoundary = true;
    model->AddSector(0.140231,0.0389385,2.62128e-06,0.0433272,0.11605,-4.62367e-06,0.55391,0.647449,-0.0864676);
    model->AddSector(0.0433272,0.11605,-4.62367e-06,-0.0928865,0.0981676,-6.09364e-06,1.26464,0.957015,-0.260301);
    model->AddSector(-0.0928865,0.0981676,-6.09364e-06,-0.140188,-0.0389793,2.08664e-06,1.50313,1.34882,-0.21918);
    model->getMax(0,2*M_PI);
    models.push_back(model);

    model = new Model;
    model->isBoundary = false;
    model->AddSector(-0.121666,-0.0413276,7.31515e-06,-0.0404978,-0.137792,8.16867e-06,1.32592,1.34946,-0.305451);
    model->AddSector(-0.0404978,-0.137792,8.16867e-06,0.0832074,-0.100172,6.08423e-06,1.50142,1.34715,-0.218812);
    model->AddSector(0.0832074,-0.100172,6.08423e-06,0.114863,0.0179242,1.46605e-06,1.26295,0.95555,-0.260039);
    model->AddSector(0.114863,0.0179242,1.46605e-06,0.0806698,0.136044,-9.70057e-08,1.66774,1.73367,-0.514309);
    model->AddSector(0.0806698,0.136044,-9.70057e-08,-0.025141,0.132821,2.93461e-06,1.2932,1.63841,-0.433953);
    model->AddSector(-0.025141,0.132821,2.93461e-06,-0.103108,0.0621824,5.95843e-06,1.21137,1.41362,-0.326353);
    model->AddSector(-0.103108,0.0621824,5.95843e-06,-0.121666,-0.0413276,7.31515e-06,1.15086,1.24152,-0.317631);
    model->getMax(0,2*M_PI);
    models.push_back(model);


    selectedModel = models[0];

    gridLayout = new QGridLayout;
    //vbox1 = new QVBoxLayout;

    // charts
    chartViewGraphs = new QChartView;
    chartViewCircles = new QChartView;

    series_tangential = new QLineSeries;
    series_tangential->setName("tangential");

    series_normal = new QLineSeries;
    series_normal->setName("normal");

    series_selected = new QScatterSeries;

    series_mohr = new QLineSeries;
    series_mohr->setName("normal vs tangential traction");

    series_mohr_selected = new QScatterSeries;

    chart_line = new QChart;
    chart_line->addSeries(series_tangential);
    chart_line->addSeries(series_normal);
    chart_line->addSeries(series_selected);

    chart_line_mohr = new QChart;
    chart_line_mohr->addSeries(series_mohr);
    chart_line_mohr->addSeries(series_mohr_selected);

    chartViewGraphs->setChart(chart_line);
    chartViewCircles->setChart(chart_line_mohr);

    chart_line_mohr->createDefaultAxes();
    chart_line->createDefaultAxes();





    // summary
    table = new QTableWidget;
    table->setRowCount(4);
    table->setColumnCount(2);

    QTableWidgetItem *row0 = new QTableWidgetItem("index");
    QTableWidgetItem *row1 = new QTableWidgetItem("angle");
    QTableWidgetItem *row2 = new QTableWidgetItem("normal traction");
    QTableWidgetItem *row3 = new QTableWidgetItem("tangential traction");

    twi_selectedIdx = new QTableWidgetItem;

    twi_angle = new QTableWidgetItem;
    twi_normal_traction = new QTableWidgetItem;
    twi_tangential_traction = new QTableWidgetItem;

    table->setItem(0,1, twi_selectedIdx);
    table->setItem(1,1, twi_angle);
    table->setItem(2,1, twi_normal_traction);
    table->setItem(3,1, twi_tangential_traction);
    table->setItem(0,0, row0);
    table->setItem(1,0, row1);
    table->setItem(2,0, row2);
    table->setItem(3,0, row3);


    // toopbar
    slider = new QSlider(Qt::Horizontal);
    slider->setRange(0,Model::number_of_ponts-1);
    ui->toolBar->addWidget(slider);

    spin = new QSpinBox();
    spin->setRange(0,models.size()-1);
    ui->toolBar->addWidget(spin);

    QWidget *w = new QWidget();
    w->setLayout(gridLayout);
    setCentralWidget(w);
/*
    // vtk
    points_origin->InsertPoint(0, 0,0,0);
    qt_vtk_widget = new QVTKOpenGLNativeWidget();
    qt_vtk_widget->setRenderWindow(renderWindow);
    renderer->SetBackground(colors->GetColor3d("White").GetData());
    renderer->AddActor(actor_mesh);
    renderer->AddActor(actor_arrows);
    renderer->AddActor(actor_labels);
    renderWindow->AddRenderer(renderer);
    gridLayout->addWidget(qt_vtk_widget, 1, 0);
*/
    gridLayout->addWidget(chartViewGraphs, 0, 0);
    gridLayout->addWidget(chartViewCircles, 0, 1);
    gridLayout->addWidget(table, 1, 1);

    connect(slider, SIGNAL(valueChanged(int)), this, SLOT(sliderValueChanged(int)));
    connect(spin, SIGNAL(valueChanged(int)), this, SLOT(spinValueChanged(int)));

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::showEvent( QShowEvent*)
{
    // for testing
    UpdateGUI();
}

void MainWindow::sliderValueChanged(int)
{
    UpdateGUI();
}



void MainWindow::UpdateGUI()
{
    // TABLE
    int selected_idx = slider->value();
    Model::Result &res = selectedModel->results[selected_idx];

    QString strIdx = QString::number(selected_idx);
    twi_selectedIdx->setText(strIdx);

    twi_angle->setText(QString::number(res.angle_fwd*180/M_PI));
    twi_normal_traction->setText(QString::number(res.trac_normal));
    twi_tangential_traction->setText(QString::number(res.trac_tangential));

    // PLOTS
    double ymin = DBL_MAX, ymax = -DBL_MAX;
    double xmin = DBL_MAX, xmax = -DBL_MAX;
    series_mohr->clear();
    series_normal->clear();
    series_tangential->clear();


    for(unsigned i=0;i<Model::number_of_ponts;i++)
    {
        Model::Result &res = selectedModel->results[i];
        double val_x = res.trac_normal;
        double val_y = res.trac_tangential;

        xmin = std::min(xmin, val_x);
        xmax = std::max(xmax, val_x);
        ymin = std::min(ymin, val_y);
        ymax = std::max(ymax, val_y);
        series_mohr->append(val_x,val_y);
    }
    double span_y = ymax-ymin;
    double span_x = xmax-xmin;
    double span = std::max(span_y, span_x)*0.65;
    double avgx = (xmax+xmin)/2;

    QList<QAbstractAxis*> axes = chart_line_mohr->axes();
    axes[0]->setRange(avgx-span, avgx+span);
    axes[1]->setRange(-span, span);


    xmin = ymin = DBL_MAX;
    xmax = ymax = -DBL_MAX;
    for(unsigned i=0;i<Model::number_of_ponts;i++)
    {
        Model::Result &res = selectedModel->results[i];
        double val_x = res.angle_fwd;
        xmin = std::min(xmin, val_x);
        xmax = std::max(xmax, val_x);
        ymin = std::min(ymin, res.trac_tangential);
        ymax = std::max(ymax, res.trac_tangential);
        ymin = std::min(ymin, res.trac_normal);
        ymax = std::max(ymax, res.trac_normal);
        series_normal->append(val_x,res.trac_normal);
        series_tangential->append(val_x,res.trac_tangential);
    }
    span_y = ymax-ymin;
    span_x = xmax-xmin;
    span = std::max(span_y, span_x)*0.65;
    avgx = (xmax+xmin)/2;

    axes = chart_line->axes();
    axes[0]->setRange(0, selectedModel->max_angle);
    axes[1]->setRange(ymin-span_y*0.1, ymax+span_y*0.1);


    // plots - selections
    series_selected->clear();
    series_mohr_selected->clear();
    series_selected->append(res.angle_fwd, res.trac_normal);
    series_mohr_selected->append(res.trac_normal, res.trac_tangential);

    // VTK
    int nSectors = selectedModel->fan.size();
/*
    points->SetNumberOfPoints(nSectors*3);
    vtkIdType pts2[3];
    ugrid->Reset();
    for(int i=0;i<nSectors;i++)
    {
        Model::Sector &s = selectedModel->fan[i];
        points->SetPoint(i*3+0, 0, 0, 0);
        points->SetPoint(i*3+1, s.u.x(), s.u.y(), 0);
        points->SetPoint(i*3+2, s.v.x(), s.v.y(), 0);

        pts2[0] = i*3+0;
        pts2[1] = i*3+1;
        pts2[2] = i*3+2;
        ugrid->InsertNextCell(VTK_TRIANGLE, 3, pts2);
    }

    points->Modified();
    ugrid->SetPoints(points);

    dataSetMapper->SetInputData(ugrid);
    actor_mesh->SetMapper(dataSetMapper);
    actor_mesh->GetProperty()->SetColor(colors->GetColor3d("Seashell").GetData());
    actor_mesh->GetProperty()->SetLineWidth(0.5);
    actor_mesh->GetProperty()->EdgeVisibilityOn();

    // sector numbers
    geometryFilter->SetInputData(ugrid);
    geometryFilter->Update();

    idfilter->SetInputConnection(geometryFilter->GetOutputPort());
    idfilter->PointIdsOff();
    idfilter->CellIdsOn();
    idfilter->FieldDataOff();

    cellCenters->SetInputConnection(idfilter->GetOutputPort());

    labledDataMapper->SetInputConnection(cellCenters->GetOutputPort());
    labledDataMapper->SetLabelModeToLabelScalars();

    actor_labels->SetMapper(labledDataMapper);
    actor_labels->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());

    // arrow
    ugrid_arrow_origin->Reset();
    pts2[0] = 0;
    ugrid_arrow_origin->InsertNextCell(VTK_VERTEX, 1, pts2);
    ugrid_arrow_origin->SetPoints(points_origin);

    arrowCoords->Reset();
    arrowCoords->SetNumberOfComponents(3);
    arrowCoords->SetName("arrowCoords");
    arrowCoords->SetNumberOfTuples(1);

    double d[3]={};
    d[0] = res.tn.x();
    d[1] = res.tn.y();
    arrowCoords->SetTuple(0, d);
    arrowCoords->Modified();

    ugrid_arrow_origin->GetPointData()->AddArray(arrowCoords);
    ugrid_arrow_origin->GetPointData()->SetActiveVectors("arrowCoords");

    arrowSource->SetShaftRadius(0.02);
    arrowSource->SetTipRadius(0.03);
    arrowSource->SetTipLength(0.06);
    glyph3D->SetSourceConnection(arrowSource->GetOutputPort());
    glyph3D->SetVectorModeToUseVector();
    glyph3D->SetInputData(ugrid_arrow_origin);
    glyph3D->OrientOn();
    glyph3D->ScalingOn();
//    glyph3D->ScalingOff();
//    glyph3D->SetColorModeToColorByVector();
//    glyph3D->SetScaleModeToScaleByVector();
    glyph3D->SetScaleFactor(0.1);
    mapper_arrows->SetInputConnection(glyph3D->GetOutputPort());
    actor_arrows->SetMapper(mapper_arrows);

    renderWindow->Render();
*/


}
void MainWindow::spinValueChanged(int val)
{
    selectedModel = models[val];
    UpdateGUI();
}

