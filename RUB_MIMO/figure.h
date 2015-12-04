#ifndef FIGURES_H
#define FIGURES_H

#include "mainwindow.h"

class Figure : public MainWindow
{
public:
    explicit Figure(QCustomPlot *figure, figure_type fig_type = NONE);
    ~Figure();

    void replot();

    bool computed;
    QCorrelationData *data;

private:

    // Clear Plots
    void clear_figure();

    // Set Figure
    void set_figure_constellation();
    void set_figure_time_signal();

    // Plot a Signal to a figure(customPlot)
    void plot_constellation();
    void plot_time_signal();

    QCustomPlot *fig;
    figure_type type, basic_type;

};



#endif // FIGURES_H

