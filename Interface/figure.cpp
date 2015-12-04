#include "figure.h"
#include "types.h"

Figure::Figure(QCustomPlot *figure, figure_type fig_type)

{
    // To be sure, that the Program only try to access the momory if there are data
    computed = false;

    // First set global pointer to Figure Handle
    fig = figure;
    type = fig_type;


    clear_figure();
    if (type==CONSTELLATION_RX1||type==CONSTELLATION_RX2||type==CONSTELLATION_TX1||type==CONSTELLATION_TX2)
        basic_type = CONSTELLATION;

    else if (type==TIME_SIGNAL_RX1||type==TIME_SIGNAL_RX2||type==TIME_SIGNAL_TX1||type==TIME_SIGNAL_TX2)
        basic_type = TIME_SIGNAL;

    else
        basic_type = NONE;

    replot();


}

Figure::~Figure()
{
    clear_figure();
}


void Figure::replot()
{
    clear_figure();
    if (computed == false)
        return;

    switch(basic_type){
        case CONSTELLATION: set_figure_constellation(); plot_constellation(); break;
        case TIME_SIGNAL: set_figure_time_signal(); plot_time_signal(); break;
    }

    fig->replot();
}

// Clear Plots
void Figure::clear_figure()
{
    fig->clearGraphs();
}

// Set Figure to configured Parameter
void Figure::set_figure_constellation()
{
    char buffer[20];
    // give the axes some labels:
    fig->xAxis->setLabel("real");
    fig->yAxis->setLabel("imag");
    fig->legend->setVisible(true);

    // Set graph 0
    fig->addGraph();
    switch(type){
        case CONSTELLATION_RX1: sprintf(buffer, "RX1"); break;
        case CONSTELLATION_RX2: sprintf(buffer, "RX2"); break;
        case CONSTELLATION_TX1: sprintf(buffer, "TX1"); break;
        case CONSTELLATION_TX2: sprintf(buffer, "TX2"); break;
        default: sprintf(buffer, "No Signal");
    }
    fig->graph(0)->setName(buffer);
    fig->graph(0)->setPen(QPen(Qt::blue));

    // Line Settings:
    fig->graph(0)->setLineStyle(QCPGraph::lsNone);
    fig->graph(0)->setScatterStyle(QCPScatterStyle::ssCircle);

}

void Figure::set_figure_time_signal()
{
    char buffer[20];
    // give the axes some labels:
    fig->xAxis->setLabel("Time/sek");
    fig->yAxis->setLabel("Magnitude");
    fig->legend->setVisible(true);

    // Set graph 0
    fig->addGraph();
    switch(type){
        case TIME_SIGNAL_RX1: sprintf(buffer, "RX1 Real"); break;
        case TIME_SIGNAL_RX2: sprintf(buffer, "RX2 Real"); break;
        case TIME_SIGNAL_TX1: sprintf(buffer, "TX1 Real"); break;
        case TIME_SIGNAL_TX2: sprintf(buffer, "TX2 Real"); break;
        default: sprintf(buffer, "No Signal");
    }
    fig->graph(0)->setName(buffer);
    fig->graph(0)->setPen(QPen(Qt::blue));


    // Set graph 1
    fig->addGraph();
    switch(type){
        case TIME_SIGNAL_RX1: sprintf(buffer, "RX1 Imag"); break;
        case TIME_SIGNAL_RX2: sprintf(buffer, "RX2 Imag"); break;
        case TIME_SIGNAL_TX1: sprintf(buffer, "TX1 Imag"); break;
        case TIME_SIGNAL_TX2: sprintf(buffer, "TX2 Imag"); break;
        default: sprintf(buffer, "No Signal");
    }
    fig->graph(1)->setName(buffer);
    fig->graph(1)->setPen(QPen(Qt::red));
}

// Plot Data to a figure on UI
void Figure::plot_constellation()
{
    // assign data to graph:
    fig->graph(0)->setData(data->real, data->imag);

    // set axes ranges, so we see all data:
    float min, max;
    min = *std::min_element(data->real.begin(), data->real.end());
    max = *std::max_element(data->real.begin(), data->real.end());
    fig->xAxis->setRange(-ceil(-min), ceil(max));
    min = *std::min_element(data->imag.begin(), data->imag.end());
    max = *std::max_element(data->imag.begin(), data->imag.end());
    fig->yAxis->setRange(-ceil(-min), ceil(max));
}

void Figure::plot_time_signal()
{
    // Compute Time:
    unsigned int sample_nbr = data->real.size();

    QVector<double> time;

    for (unsigned int i=0; i<sample_nbr; i++)
        time.push_back(((float)i)/parameter->fs);

    // assign data to graph:
    fig->graph(0)->setData(time, data->real);
    fig->graph(1)->setData(time, data->imag);

    // set axes ranges, so we see all data:
    float min, max;
    min = *std::min_element(time.begin(), time.end());
    max = *std::max_element(time.begin(), time.end());
    fig->xAxis->setRange(-ceil(-min), max);
    min = *std::min_element(data->imag.begin(), data->imag.end());
    max = *std::max_element(data->imag.begin(), data->imag.end());
    fig->yAxis->setRange(-ceil(-min), ceil(max));
}

