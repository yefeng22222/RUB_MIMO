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
}

// Reconstructor, for a device, that is attached and saved in config
void usrp_device::reconstruct(uhd::device_addr_t dev_adress)
{
    usrp_device::construct_from_uhd(dev_adress);
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

    return subroot;
}


/// --------------------------------------------------------------------------------------------------
/// Private

// Construct a Device from uhd address
void usrp_device::construct_from_uhd (uhd::device_addr_t dev_address)
{
    attached = true;
    dev_addr = dev_address;

    unsigned int pos1 = 0, pos2 = 0;
    std::string temp = dev_addr.to_string();

    // Write the Parameter from uhd object to local variables for easier accessability
    pos1 = temp.find("type=")+5;
    pos2 = temp.find(",", pos1);
    dev_type = temp.substr(pos1, pos2-pos1);

    pos1 = temp.find("addr=")+5;
    if (pos1 != std::string::npos+5)
    {
        pos2 = temp.find(",", pos1);
        addr = temp.substr(pos1, pos2-pos1);
    }
    else  // Seem to be an USB Device
    {
        addr.assign("USB");
    }

    pos1 = temp.find("name=")+5;
    pos2 = temp.find(",", pos1);
    id = temp.substr(pos1, pos2-pos1);

    pos1 = temp.find("serial=")+7;
    pos2 = temp.find(",", pos1);
    serial = temp.substr(pos1, pos2-pos1);

    pos1 = temp.find("product=")+8;
    if (pos1 != std::string::npos+8)
    {
        pos2 = temp.find(",", pos1);
        product = temp.substr(pos1, pos2-pos1);
    }
    else
        product.assign("No Info");
}
