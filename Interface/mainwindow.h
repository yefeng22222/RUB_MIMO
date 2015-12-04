#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qcustomplot.h>
#include "uhd/device.hpp"
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include "fstream"
#include "types.h"
#include "usrp_device.h"

#define CONFIG_FILE "dev_config.json"

class usrp_device;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    /// Functions:
    MainWindow(QCustomPlot *figure, figure_type fig_type);
    MainWindow(QWidget *parent = 0, bool createplots = false);
    ~MainWindow();

    // TODO: What does the function
    Parameter * get_parameter();

    // Get the devices that are connected (LAN, USB)
    void * update_attached_device_list();

    // Update the Vector, with the Devices from File
    void * update_saved_device_list();

    // Saves all Devices to Config File
    void * save_config();


    /// Variables
    Parameter * parameter;

private:

/// Variables:
    Ui::MainWindow *ui;

    // List with Devices
    std::vector<usrp_device> devices;
    // The selected Devs, valid range: 0 - nbr_of_devs-1
    unsigned int current_dev_rx, current_dev_tx;



/// Functions:

    // Execute a Program and get the response
    std::string execute(const char* cmd);

    // Print some Content to log window
    void print_to_log(std::string content);

    // Read a Binary File and copy it to a Vector
    // ATTENTION: The vector contains real and imag floats: <Real> <Imag> <Real> .....
    QCorrelationData *read_binary_file(char * filename);

    // Read all Binary Files and save to ram
    int get_data();

    // Plot all Figures
    void replot_figures();

    // Set Standard Variables to Formularfields
    void set_standard_variables();

    // Check wether File exists
    bool is_file_exist(const char *fileName);

    // Reconfigure Devs for executation
    void reconfig_devs();

private slots:

    // Defines what happens when Buttons are clicked
    void on_run_button_clicked();

    // Object in TX Dev Identification in GUI changed
    void on_tx_dev_id_currentIndexChanged(int index);

    // Object in RX Dev Identification in GUI changed
    void on_rx_dev_id_currentIndexChanged(int index);
};




#endif // MAINWINDOW_H
