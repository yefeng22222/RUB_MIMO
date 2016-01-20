#include "usrp_device.h"
#include "json/json.h"

// Constructor for attached Devices
usrp_device::usrp_device(uhd::device_addr_t dev_adress)
{
    usrp_device::construct_from_uhd(dev_adress);
}

// Constructor for saved Devices in Json Config
usrp_device::usrp_device(Json::Value dev)
{
    serial = dev.get("Serial", "NULL").asString();
    id = dev.get("ID", "NULL").asString();
    addr = dev.get("Adress", "NULL").asString();
    dev_type = dev.get("Type", "NULL").asString();
    product = dev.get("Product", "NULL").asString();
    tx_gain = dev.get("TX Gain", 0).asDouble();
    rx_gain = dev.get("RX Gain", 0).asDouble();
    center_freq = dev.get("Center Freq.", 0).asDouble();
    sampling_rate = dev.get("Samp. Rate", 0).asDouble();
    nbr_of_subcarriers = dev.get("Number of Subcarriers", 0).asUInt();
    nbr_of_nullcarriers = dev.get("Number of Nullcarriers", 0).asUInt();
    prefix_length = dev.get("Prefix Length", 0).asDouble();
    training_sequences = dev.get("Training Sequences", 0).asDouble();
    subdev_spec.assign(dev.get("Subdevice Specifications", get_standard_subdev()).asString());
}

// Returns a formattet String in JSON Format, for the Device
Json::Value usrp_device::to_json_string()
{
    Json::Value subroot;
    subroot["ID"] = id;
    subroot["Serial"] = serial;
    subroot["Adress"] = addr;
    subroot["Type"] = dev_type;
    subroot["Product"] = product;
    subroot["TX Gain"] = tx_gain;
    subroot["RX Gain"] = rx_gain;
    subroot["Center Freq."] = center_freq;
    subroot["Samp. Rate"] = sampling_rate;
    subroot["Number of Subcarriers"] = nbr_of_subcarriers;
    subroot["Number of Nullcarriers"] = nbr_of_nullcarriers;
    subroot["Prefix Length"] = prefix_length;
    subroot["Training Sequences"] = training_sequences;
    subroot["Subdevice Specifications"] = subdev_spec;


    return subroot;
}

// Construct a Device from uhd address
void usrp_device::construct_from_uhd (uhd::device_addr_t dev_address)
{
    attached = true;
    dev_addr = dev_address;
    std::string buffer;

    unsigned int pos1 = 0, pos2 = 0;
    std::string temp = dev_addr.to_string();

    // Write the Parameter from uhd object to local variables for easier accessability
    pos1 = temp.find("type=")+5;
    pos2 = temp.find(",", pos1);
    dev_type = temp.substr(pos1, pos2-pos1);


    pos1 = temp.find("name=")+5;
    pos2 = temp.find(",", pos1);
    id = temp.substr(pos1, pos2-pos1);

    pos1 = temp.find("serial=")+7;
    pos2 = temp.find(",", pos1);
    serial = temp.substr(pos1, pos2-pos1);

    pos1 = temp.find("addr=")+5;
    if (pos1 != std::string::npos+5)
    {
        pos2 = temp.find(",", pos1);
        addr = temp.substr(pos1, pos2-pos1);
    }
    else  // Seem to be an USB Device
    {
        addr.assign("serial=");
        addr.append(serial);
    }

    pos1 = temp.find("product=")+8;
    if (pos1 != std::string::npos+8)
    {
        pos2 = temp.find(",", pos1);
        product = temp.substr(pos1, pos2-pos1);
    }
    else
        product.assign("No Info");
}

/// --------------------------------------------------------------------------------------------------
/// Private

// Get the Standart Subdevice Specifications
char * usrp_device::get_standard_subdev()
{
    if (mode == RX_DEV && (dev_type.compare("b210") == 0 || dev_type.compare("B210") == 0))
        return "A:A A:B";

    else if (mode == TX_DEV && (dev_type.compare("b210") == 0 || dev_type.compare("B210") == 0))
        return "A:B A:A";

    else if (dev_type.compare("x300") == 0 || dev_type.compare("X300") == 0)
        return "A:0 B:0";

    else if (dev_type.compare("n200") == 0 || dev_type.compare("N200") == 0)
        return "A:B A:A";

    else
        return "A:0 B:0";
}
