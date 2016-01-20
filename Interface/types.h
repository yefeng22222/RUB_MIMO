#include "mainwindow.h"
#include <sstream>

#ifndef TYPES_H
#define TYPES_H


// Figure Types
#define CONSTELLATION 9
#define TIME_SIGNAL 10
#define CONSTELLATION_TX1 11
#define CONSTELLATION_TX2 12
#define CONSTELLATION_RX1 13
#define CONSTELLATION_RX2 14
#define TIME_SIGNAL_TX1 15
#define TIME_SIGNAL_TX2 16
#define TIME_SIGNAL_RX1 17
#define TIME_SIGNAL_RX2 18
#define NONE 19

// Communication modes
#define MODE_SISO 0
#define MODE_RX_DIVERSITY 1
#define MODE_RX_ZF 2
#define MODE_RX_BEAMFORMING 3
#define MODE_TX_BEAMFORMING 4

// Device Type
#define TX_DEV 0
#define RX_DEV 1


// Macro for converting to String
#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()



// Type for the Complex Signals
typedef struct {
    QVector<double> real;
    QVector<double> imag;
} QCorrelationData;

typedef int figure_type;

typedef struct {
    std::string mode;
    int fc;
    int fs;
    std::string tx_dev_id;
    std::string rx_dev_id;
} Parameter;


#endif // TYPES_H
