#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "figure.h"

Figure *const_tx_1, *const_rx_1, *const_tx_2, *const_rx_2, *time_tx_1, *time_tx_2, *time_rx_1, *time_rx_2;

using namespace std;
MainWindow::MainWindow(QCustomPlot *figure, figure_type fig_type)
{
// empty constructor, preserves some issues
}


MainWindow::MainWindow(QWidget *parent, bool createplots) :
      QMainWindow(parent),
      ui(new Ui::MainWindow)
{
    // Chose "true" only once and only in the main.cpp,
    // else the program crashes, because of some cild classes who try
    // to create new plots
    if (createplots)
    {
        ui->setupUi(this);
        const_tx_1 = new Figure(ui->const_tx_1, CONSTELLATION_TX1);
        const_tx_2 = new Figure(ui->const_tx_2, CONSTELLATION_TX2);
        const_rx_1 = new Figure(ui->const_rx_1, CONSTELLATION_RX1);
        const_rx_2 = new Figure(ui->const_rx_2, CONSTELLATION_RX2);
        time_tx_1 = new Figure(ui->time_tx_1, TIME_SIGNAL_TX1);
        time_tx_2 = new Figure(ui->time_tx_2, TIME_SIGNAL_TX2);
        time_rx_1 = new Figure(ui->time_rx_1, TIME_SIGNAL_RX1);
        time_rx_2 = new Figure(ui->time_rx_2, TIME_SIGNAL_RX2);
        parameter = new Parameter;


        update_saved_device_list();
        update_attached_device_list();
        set_standard_variables();
        save_config();

    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

Parameter * MainWindow::get_parameter()
{
    return parameter;
}


// Update the Vector with all connected URSP devices
void * MainWindow::update_attached_device_list()
{
    uhd::device_addr_t hint;        // an empty hint to discovers all devices
    uhd::device_addrs_t dev_addrs = uhd::device::find(hint);

    int avaliable;
    std::string temp;

    /// Only for devleopment,to put the string from uhd::find in a file
    /// --------------------------------------------------------------------------------------
    //fstream dev_file;
    //dev_file.open("dev_output.txt", ios::out);
    //std::string dev_string;
    //dev_string.assign("Only for Devlopment");
    /// --------------------------------------------------------------------------------------


    for(unsigned i=0; i<dev_addrs.size(); i++)
    {
        // Check wether the Dev is already available
        avaliable = -1;
        temp.assign(dev_addrs.at(i).to_string());

        /// Only for devleopment
        /// --------------------------------------------------------------------------------------
        //dev_string.append("\n");
        //dev_string.append(temp);
        /// --------------------------------------------------------------------------------------


        for (unsigned p=0; p<devices.size(); p++)
        {
            // Compare the serial of each device in Vector with the serial in the UHD::String
            // If found, update Adress and Name
            avaliable = devices.at(p).serial.compare(temp.substr(temp.find("serial=")+7));
            if (avaliable == 0)
            {devices.at(p).construct_from_uhd(dev_addrs.at(i)); break;}
        }

        // Else add the Device to the end of the vector
        if (avaliable!= 0)
            devices.push_back(dev_addrs.at(i));
    }
    /// Only for devleopment
    /// --------------------------------------------------------------------------------------
    //dev_file << dev_string;
    //dev_file.close();
    /// --------------------------------------------------------------------------------------
}

// Update the Vector, with the Devices from File
void * MainWindow::update_saved_device_list()
{
    if (is_file_exist(CONFIG_FILE) == true)
    {
        ifstream file;
        file.open(CONFIG_FILE);

        char buffer[20];
        Json::Value root;
        Json::Reader reader;

        reader.parse(file, root);

        int nbr = root.get("nbr_of_devs", 0).asInt();

        for (int i=0; i<nbr; i++)
        {
            sprintf(buffer, "Device_%d", i);
            devices.push_back(root.get(buffer, "NULL"));
        }
        file.close();
    }
}

// Saves all Devices to Config File
void * MainWindow::save_config()
{
    fstream file;
    file.open(CONFIG_FILE, ios::out);
    char buffer[20];

    Json::Value root;
    int nbr = devices.size();

    for (int i=0; i<nbr; i++)
    {
        sprintf(buffer, "Device_%d", i);
        root[buffer] = devices.at(i).to_json_string();
    }
    root["nbr_of_devs"] = nbr;

    file << root;
    file.close();
}


/// Private:

// Execute a Program and get the response
std::string MainWindow::execute(const char* cmd)
{
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
    if(fgets(buffer, 128, pipe) != NULL)
        result += buffer;
    }
    pclose(pipe);
    return result;
}

// Show Strings in the log Field on UI
void MainWindow::print_to_log(std::string content)
{
    ui->log->append(QString::fromStdString(content));
}

// Read a Binary File and copy it to a Vector
// ATTENTION: The vector contains real and imag floats: <Real> <Imag> <Real> .....
QCorrelationData * MainWindow::read_binary_file(char *filename)
{
    FILE *f = fopen(filename, "rb");
    if (f == NULL)
    {
        std::string buffer("\n----ERROR: File '");
        buffer.append(filename);
        buffer.append("' can't be open----\n");
        print_to_log(buffer);
        QCorrelationData *null_vector; return null_vector;
    }

    // Get File Size
    fseek(f, 0, SEEK_END);
    unsigned long lSize = ftell(f);
    rewind(f);

    // allocate Memory (vecor is dynamical)
    std::vector<float> buffer;
    float v;
    for (unsigned long i=0; i<lSize/4; i++)
    {
        fread((void*)&v, sizeof(float), 1, f);
        buffer.push_back(v);
    }

    int size = buffer.size()/2;
    QCorrelationData *data = new QCorrelationData;

    // Convert Signal Vector to two QVectors (Real and Imaginary)

    for (int i=0; i<size; i++)
    {
       data->real.push_back(buffer[2*i]);
       data->imag.push_back(buffer[2*i+1]);
    }
    return data;
}

// Get all Data from Files
int MainWindow::get_data()
{
    ///TODO:
    /// Use the correct Filenames an Paths
    const_tx_1->data = read_binary_file("/tmp/tx_data1.dat");
    const_tx_2->data = read_binary_file("/tmp/tx_data2.dat");
    const_rx_1->data = read_binary_file("/tmp/rx_data1.dat");
    const_rx_2->data = read_binary_file("/tmp/rx_data2.dat");
    time_tx_1->data = read_binary_file("/tmp/tx_sig1.dat");
    time_tx_2->data = read_binary_file("/tmp/tx_sig2.dat");
    time_rx_1->data = read_binary_file("/tmp/rx_sig1.dat");
    time_rx_2->data = read_binary_file("/tmp/rx_sig2.dat");
    return 0;
}

// Plot all the figures
void MainWindow::replot_figures()
{
    const_rx_1->replot();
    const_rx_2->replot();
    const_tx_1->replot();
    const_tx_2->replot();
    time_rx_1->replot();
    time_rx_2->replot();
    time_tx_1->replot();
    time_tx_2->replot();
}

// Set standard variables in UI
void MainWindow::set_standard_variables()
{
    // No Device Avaliable
    if (devices.size()<1)
    {
        /// TODO: Fehler Fenster "No Dev avaliable" erstellen
        return;
    }

    // List with powers of 2
    QStringList list;
    list.append(QString::number(0));
    list.append(QString::number(1));
    for (int i=0; i<=10; i++)
        list.append(QString::number((2 << i)));

    // Create a List with devices
    QStringList dev_list;
    for (unsigned i=0; i<devices.size(); i++)
        dev_list.append(QString(devices.at(i).id.c_str()));


    ui->tx_dev_id->addItems(dev_list);
    ui->rx_dev_id->addItems(dev_list);

    ui->fc->setValidator(new QDoubleValidator(0, 1000000, 2, this) );
    ui->fs->setValidator(new QDoubleValidator(0, 1000000, 2, this) );
    ui->cyclic_prefix_len->setValidator(new QIntValidator(0, 200, this) );
    ui->tx_gain->setValidator(new QDoubleValidator(0, 200, 3, this) );
    ui->rx_gain->setValidator(new QDoubleValidator(0, 200, 3, this) );
    ui->training_sequences->setValidator(new QDoubleValidator(0, 20, 3, this) );

    ui->no_subcarriers->addItems(list);
    ui->no_subcarriers->setCurrentText("2048");
    ui->no_nullcarriers->addItems(list);
    ui->no_nullcarriers->setCurrentText("0");

    ui->cyclic_prefix_len->setText("5");
    ui->tx_gain->setText("5");
    ui->rx_gain->setText("5");
    ui->fc->setText("1000");
    ui->fs->setText("10000");
    ui->training_sequences->setText("5");
}

// Check wether file exist
bool MainWindow::is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

// Reconfig the Devs for execution
void MainWindow::reconfig_devs()
{
    // Reconfig TX Device
    devices.at(current_dev_tx).center_freq = ui->fc->text().toDouble();
    devices.at(current_dev_tx).sampling_rate = ui->fs->text().toDouble();
    devices.at(current_dev_tx).modulation_type = ui->modulation_type->currentIndex();
    devices.at(current_dev_tx).nbr_of_nullcarriers = ui->no_nullcarriers->currentText().toInt();
    devices.at(current_dev_tx).nbr_of_subcarriers = ui->no_subcarriers->currentText().toInt();
    devices.at(current_dev_tx).tx_gain = ui->tx_gain->text().toDouble();
    devices.at(current_dev_tx).rx_gain = ui->rx_gain->text().toDouble();
    devices.at(current_dev_tx).prefix_length = ui->cyclic_prefix_len->text().toDouble();
    devices.at(current_dev_tx).training_sequences = ui->training_sequences->text().toInt();

    // Reconfig RX Device
    devices.at(current_dev_rx).center_freq = ui->fc->text().toDouble();
    devices.at(current_dev_rx).sampling_rate = ui->fs->text().toDouble();
    devices.at(current_dev_rx).modulation_type = ui->modulation_type->currentIndex();
    devices.at(current_dev_rx).nbr_of_nullcarriers = ui->no_nullcarriers->currentText().toInt();
    devices.at(current_dev_rx).nbr_of_subcarriers = ui->no_subcarriers->currentText().toInt();
    devices.at(current_dev_rx).tx_gain = ui->tx_gain->text().toDouble();
    devices.at(current_dev_rx).rx_gain = ui->rx_gain->text().toDouble();
    devices.at(current_dev_rx).prefix_length = ui->cyclic_prefix_len->text().toDouble();
    devices.at(current_dev_rx).training_sequences = ui->training_sequences->text().toInt();
}

/// Private Slots:

// Defines what Happens when the Run Button in the UI is clicked
void MainWindow::on_run_button_clicked()
{
    print_to_log("Processing...\n");
    reconfig_devs();

    ///TODO:
    /// Execute real Communitcation
    /// actually only a Dummy
    //std::string result = execute("cp -r /mnt/Sascha/*.dat /tmp/");
    //print_to_log(result);
    //get_data();

    const_rx_1->computed = true;
    const_rx_2->computed = true;
    const_tx_1->computed = true;
    const_tx_2->computed = true;
    time_rx_1->computed = true; time_rx_1 ->parameter = parameter;
    time_rx_2->computed = true; time_rx_2 ->parameter = parameter;
    time_tx_1->computed = true; time_tx_1 ->parameter = parameter;
    time_tx_2->computed = true; time_tx_2 ->parameter = parameter;

    replot_figures();

}


// Object in TX Dev Identification in GUI changed
void MainWindow::on_tx_dev_id_currentIndexChanged(int index)
{

}

// Object in RX Dev Identification in GUI changed
void MainWindow::on_rx_dev_id_currentIndexChanged(int index)
{

}
