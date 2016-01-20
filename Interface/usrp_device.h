#ifndef USRP_DEVICE_H
#define USRP_DEVICE_H

#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <json/json.h>
#include "types.h"


#define MOD_QUAM4 1
#define MOD_QUAM16 2
#define MOD_QUAM32 3
#define MOD_QUAM64 4


// This class gives a object, that represent a usrp_device
class usrp_device
{
public:

    // Constructor for attached Device
    usrp_device(uhd::device_addr_t dev_adress);

    // Constructor for saved Devices
    usrp_device(Json::Value dev);

    // Devices Parameter
    std::string addr;
    std::string serial;
    std::string id;
    std::string dev_type;
    std::string product;
    double tx_gain;
    double rx_gain;
    double center_freq;
    double sampling_rate;
    unsigned int nbr_of_subcarriers;
    unsigned int nbr_of_nullcarriers;
    double prefix_length;
    double training_sequences;
    std::string subdev_spec;

    unsigned int modulation_type;
    // TX OR RX
    unsigned int mode;
    bool attached;


    /// Functions:
    ///

    // Construct a Device from uhd address
    void construct_from_uhd (uhd::device_addr_t dev_address);

    // Returns Parameter in a JSON String
    Json::Value to_json_string();

    // Configure and start UHD Dev for Communication


private:

    // UHD Handler to work with the device
    uhd::device_addr_t dev_addr;

    /// Functions

    // Get the Standart Subdevice Specifications
    char * get_standard_subdev();

};

#endif // USRP_DEVICE_H
