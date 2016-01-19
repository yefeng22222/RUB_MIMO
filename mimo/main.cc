/*
 * Copyright (c) 2014, 2015 Manu T S
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <uhd/usrp/multi_usrp.hpp>
#include <iostream>
#include <uhd/utils/thread_priority.hpp>
#include <uhd/utils/safe_main.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <ctime>
#include <liquid/liquid.h>
#include <volk/volk.h>
#include <volk/volk_malloc.h>
#include <unistd.h>
#include <gnuradio/gr_complex.h>
#include "framing.h"
#include "config.h"

namespace po = boost::program_options;

/// CHANGED: Generally Changed:
/// replaced NUM_SUBCARRIERS with d_options.M
/// added d_options as Global Variable




// global variablds
unsigned int pid;
unsigned int num_frames_detected;
unsigned int num_valid_headers_received;
unsigned int num_valid_bytes_received;
std::vector<unsigned int> num_valid_bytes;
unsigned int num_valid_packets_received;
time_t tx_begin, tx_end, rx_begin, rx_end;
unsigned long int tx_sig_len;
std::vector<gr_complex *> tx_sig;
std::vector<gr_complex *> rx_sig;
std::vector<unsigned int *> tx_data;
std::vector<unsigned int *> rx_data;
unsigned int received_sample_counter;
unsigned int num_data_carriers,
             num_pilot_carriers,
             num_occupied_carriers,
             num_null_carriers;

int _callback(unsigned char *  _header,
              int              _header_valid,
              unsigned char *  _payload,
              unsigned int     _payload_len,
              int              _payload_valid,
              framesyncstats_s _stats,
              void *           _userdata)
{
  // update global counters
  num_frames_detected++;

  if (_header_valid)
    num_valid_headers_received++;

  if (_payload_valid) {
    num_valid_packets_received++;
    num_valid_bytes_received += _payload_len;
  }
  return 0;
}
// rxmd
uhd::rx_metadata_t rxmd;

// thread synchronization
bool start_ac = false;
bool stop_rx  = false;
bool start_mimo = false;
pthread_mutex_t mutex_txrx;
pthread_cond_t condition_ac;
pthread_cond_t condition_mimo;

// CSI matrix
gr_complex *** G;
// precoding matrix
gr_complex *** W;

// Compute precoding vector given CSI
void * design_mimo_precoder ()
{
  printf("******* design mimo precoder *********\n");
  return NULL;
}

void * callback (std::vector<gr_complex *> x, unsigned int occupied_carriers) {
  num_frames_detected++;
  if(num_valid_packets_received == PID_MAX)
    return NULL;
  num_valid_packets_received++;
  for(unsigned int stream = 0; stream < NUM_STREAMS; stream++) {
    memmove(rx_sig[stream] + received_sample_counter,
	    x[stream],
	    sizeof(gr_complex)*occupied_carriers);
  }
  
#if DEBUG_PRINT_VERBOSE
  printf("******* callback %4u ********\n",
         num_valid_packets_received);
  for(unsigned int i = 0; i < occupied_carriers; i++) {
    printf("%3u,\t", i);
    for(unsigned stream = 0; stream < NUM_STREAMS; stream++) {
      printf("%3.6f + %3.6fi,\t\t",
             //      real(rx_sig[stream][occupied_carriers + i]),
             //      imag(rx_sig[stream][occupied_carriers + i]));
             real(x[stream][i]),
             imag(x[stream][i]));
    }
    printf("\n");
  }
#endif

  received_sample_counter += occupied_carriers;
  return NULL;
}

// structure to hold command line inputs
typedef struct
{
  double cent_freq;         // center frequency of transmission
  double samp_rate;         // usrp samping rate
  float dsp_gain;           // dsp gain
  double txgain;            // tx frontend gain
  double rxgain;            // rx frontend gain
  unsigned int M;           // number of subcarriers
  unsigned int cp_len;      // length of cyclic prefix
  bool verbose;             // verbosity of the app


  // ADDED: Dev Adress Specifications
  std::string rx_adress;
  std::string tx_adress;

  std::string rx_subdev_spec;
  std::string tx_subdev_spec;

} options;

// structure to hold usrp configurations
typedef struct
{
  double cent_freq;               // center frequency
  double samp_rate;               // sampling ratea
  float rf_gain;                  // sampling ratea
  clock_source_type clock_source; // clock source TODO List the sources
  time_source_type time_source;   // time source TODO List the sources 

} usrp_config;

options d_options;

// function to read commandline options
// argc: argument count
// argv: argument list
// parameter options: pointer to a struct_options object
void read_options(int       argc,
                  char  **  argv,
                  options * d_options)
{
  //set the operating parameters
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",
     "help message")
    ("freq,f",
     po::value<double>(&(d_options->cent_freq))->default_value(
                                            CENTER_FREQUENCY),
     "RF center frequency in Hz")
    ("rate,r",
     po::value<double>(&(d_options->samp_rate))->default_value(
                                            SAMPLING_RATE),
     "USRP Sampling rate")
    ("dsp_gain",
     po::value<float>(&(d_options->dsp_gain))->default_value(
                                          BASEBAND_GAIN),
     "TX DSP gain")
    ("tx_gain",
     po::value<double>(&(d_options->txgain))->default_value(
                                         TX_FRONTEND_GAIN),
     "TX Front end gain")
    ("rx_gain",
     po::value<double>(&(d_options->rxgain))->default_value(
                                         RX_FRONTEND_GAIN),
     "RX Front end gain")
    ("num_subcarriers",
     po::value<unsigned int>(&(d_options->M))->default_value(
                                          NUM_SUBCARRIERS),
     "Number of OFDM subcarriers")
    ("cp_len",
     po::value<unsigned int>(&(d_options->cp_len))->default_value(
                                               CP_LENGTH),
     "Cyclic Prefix Length")

     // ADDED: DEV Adress Specifications
    ("rx_addr", po::value<std::string>(&(d_options->rx_adress))->default_value(B210_RX), "RX Device Adress")
    ("tx_addr", po::value<std::string>(&(d_options->tx_adress))->default_value(B210_TX), "TX Device Adress")

    // ADDED: Subdev_specs  --- WARNING: Changed code at Line
    ("tx_subdev", po::value<std::string>(&(d_options->tx_subdev_spec))->default_value(B210_SUBDEV_SPEC_TX), "TX Subdevice Specs")
    ("rx_subdev", po::value<std::string>(&(d_options->rx_subdev_spec))->default_value(B210_SUBDEV_SPEC_RX), "RX Subdevice Specs")

    ("verbose,v",
     "Verbose")
    ("quite,q",
     "Quite")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  //print the help message
  if (vm.count("help")) {
    std::cout << boost::format("ofdmtxrx %s") % desc << std::endl;
    std::cout
      << std::endl
      << "Basic OFDM P2P Link\n"
      << std::endl;
    exit(0);
  }

  //check sanity of options
  if (vm.count("verbose") && vm.count("quite")) {
    std::cout << "Conflicting Options Verbose and Quite."
      << " Please use only one of those."
      << std::endl;
    exit(0);
  }

  if (vm.count("verbose"))
    d_options->verbose = true;
}

// initialize USRP
// parameter u: sptr to USRP
// parameter configs: values for freq, gain etc
// parameter is_tx: set to true if USRP is transmitting
// TODO set device time to 0 before tuning
void init_usrp(uhd::usrp::multi_usrp::sptr u,
               usrp_config * configs,
               bool is_tx)
{
  size_t num_chans;     // number of channels in USRP
  // set clock source
  switch(configs->clock_source) {
    case(CLOCK_SOURCE_NONE):
      break;
    case(CLOCK_SOURCE_EXTERNAL):
      u->set_clock_source("external");
      break;
    default:
      std::cout << "Clock source not known. Exiting\n";
      exit(1);
  }
  // set time source
  switch(configs->time_source) {
    case(TIME_SOURCE_NONE):
      break;
    case(TIME_SOURCE_EXTERNAL):
      u->set_time_source("external");
      break;
    default:
      std::cout << "Time source not known. Exiting\n";
      exit(1);
  }

  if(is_tx) {
    // set subdev specs
    u->set_tx_subdev_spec(uhd::usrp::subdev_spec_t(d_options.tx_subdev_spec.c_str()), uhd::usrp::multi_usrp::ALL_MBOARDS);
///---------------------  OLD CODE  --------------------------
//    if(u->get_mboard_name() == "X300") {
//      u->set_tx_subdev_spec(
//          uhd::usrp::subdev_spec_t(X300_SUBDEV_SPEC),
//          uhd::usrp::multi_usrp::ALL_MBOARDS);
//    }
//    else if((u->get_mboard_name() == "N200") ||
//            (u->get_mboard_name() == "N200r4")){
//      u->set_tx_subdev_spec(
//          uhd::usrp::subdev_spec_t(N200_SUBDEV_SPEC),
//          uhd::usrp::multi_usrp::ALL_MBOARDS);
//    }
//    else if(u->get_mboard_name() == "B210") {
//      u->set_tx_subdev_spec(
//          uhd::usrp::subdev_spec_t(B210_SUBDEV_SPEC_TX),
//          uhd::usrp::multi_usrp::ALL_MBOARDS);
//    }
//    else {
//      std::cout << "TX Motherboard not compatible\n"
//		<< "Subdevice specification for "
//                << u->get_mboard_name()
//                << " not known. Exiting\n";
//      exit(1);
//    }
///---------------------------------------------------------------
    num_chans = u->get_tx_num_channels();
    // set freq, gain and antenna
    // TODO pass antenna as a parameter
    for (size_t chan = 0; chan < num_chans; chan++) {
      u->set_tx_rate(configs->samp_rate, chan);
      u->set_tx_gain(configs->rf_gain, chan);
      u->set_tx_antenna("TX/RX", chan);
    }
    // tune the channels simultaneously to obtain
    // same phase on all the channels
    uhd::time_spec_t cmd_time = u->get_time_now() + 
                                uhd::time_spec_t(0.1);
    u->set_command_time(cmd_time);
    uhd::tune_request_t tx_tune_request(configs->cent_freq);
    for (size_t chan = 0; chan < num_chans; chan++) 
      u->set_tx_freq(tx_tune_request, chan);
    u->clear_command_time();
  }

  // IS RX
  else {
    // set subdev specs
    u->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(d_options.tx_subdev_spec.c_str()), uhd::usrp::multi_usrp::ALL_MBOARDS);
///---------------------  OLD CODE  --------------------------
//    if(u->get_mboard_name() == "X300") {
//      u->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(
//            X300_SUBDEV_SPEC),
//            uhd::usrp::multi_usrp::ALL_MBOARDS);
//    }
//    else if((u->get_mboard_name() == "N200") ||
//            (u->get_mboard_name() == "N200r4")) {
//      u->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(
//            N200_SUBDEV_SPEC),
//            uhd::usrp::multi_usrp::ALL_MBOARDS);
//    }
//    else if(u->get_mboard_name() == "B210") {
//      u->set_rx_subdev_spec(uhd::usrp::subdev_spec_t(
//            B210_SUBDEV_SPEC_RX),
//            uhd::usrp::multi_usrp::ALL_MBOARDS);
//    }
//    else {
//      std::cout << "RX Motherboard not compatible\n"
//                << "Subdevice specification for "
//                << u->get_mboard_name()
//                << " not known. Exiting\n";
//      exit(1);
//    }
///-----------------------------------------------------------
    num_chans = u->get_rx_num_channels();
    // set freq, gain and antenna
    // TODO pass antenna as a parameter
    for (size_t chan = 0; chan < num_chans; chan++) {
      u->set_rx_rate(configs->samp_rate, chan);
      u->set_rx_gain(configs->rf_gain, chan);
      u->set_rx_antenna("TX/RX", chan);
    }
    // tune the channels simultaneously to obtain
    // same phase on all the channels
    uhd::time_spec_t cmd_time = u->get_time_now() + 
                                uhd::time_spec_t(0.1);
    u->set_command_time(cmd_time);
    uhd::tune_request_t rx_tune_request(configs->cent_freq);
    for (size_t chan = 0; chan < num_chans; chan++) 
      u->set_rx_freq(rx_tune_request, chan);
    u->clear_command_time();
  }
}

// tx beamforming setup. CSI feedback involved.
#if TX_BEAMFORMING
// structure to load data to tx_thread
typedef struct
{
  uhd::usrp::multi_usrp::sptr * tx;
  framegen * fg;
} tx_thread_data;

// structure to load data to rx_thread
typedef struct
{
  uhd::usrp::multi_usrp::sptr * rx;
  framesync * fs;
} rx_thread_data;

// rx_worker thread
void * rx_worker (void * _data)
{
  rx_thread_data * data = (rx_thread_data *)_data;
  // rx buffer length
  unsigned int rx_buffer_len;
  // number of channels available on USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples read in each call of recv
  unsigned int num_samples_read;
  // usrp buffer
  std::vector<std::complex<float> *> rx_buffer;
  // log files
  std::vector<FILE *> fp;
  // state returned from framesync
  framesync_states_t sync_state = STATE_SEEK_PLATEAU;

  num_channels = (*(data->rx))->get_rx_num_channels();
  assert(num_channels == NUM_STREAMS);
  for(size_t chan = 0; chan < num_channels; chan++) {
    channels.push_back(chan);
    if(LOG) {
      fp.push_back(fopen(boost::str(
        boost::format("%srx%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
    }
  }

  // initialize rx streamer
  uhd::stream_args_t rx_stream_args(CPU, WIRE);
  rx_stream_args.channels = channels;
  uhd::rx_streamer::sptr rx_stream = 
    (*(data->rx))->get_rx_stream(rx_stream_args);

  // initialilze rx_buffer_len, allocate memory
  rx_buffer_len = rx_stream->get_max_num_samps();
  for(size_t chan = 0; chan < num_channels; chan++) {
    // FIXME is volk malloc required here??
    rx_buffer.push_back((std::complex<float> *) malloc 
      (sizeof(std::complex<float>)*rx_buffer_len));
  }

  // start receiving in 0.1 seconds
  uhd::stream_cmd_t stream_cmd(
      uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
  stream_cmd.stream_now = false;
  stream_cmd.time_spec = (*(data->rx))->get_time_now() 
                       + uhd::time_spec_t(0.1);
  rx_stream->issue_stream_cmd(stream_cmd);
  // lock mutex and signal tx to start
  std::cout << "Rx locking the mutex and signalling tx\n";
  pthread_mutex_lock(&mutex_txrx);
  start_ac = true;
  pthread_cond_signal(&condition_ac);
  pthread_mutex_unlock(&mutex_txrx);

  rx_begin = time(NULL);
  // NOTE: timout should be larger if rx_buffer_length
  // is larger.
  float timeout = 0.2;
  // FIXME on what condition does rx stop streaming?
  // TODO receive sigal from tx to stop
  bool break_loop = false;
  while(1)
  {
    num_samples_read = rx_stream->recv(rx_buffer,
                                       rx_buffer_len,
                                       rxmd,
                                       timeout);
    if(rxmd.error_code)
      break; // TODO error code handling to be done in main

    // logging
    if(LOG) {
      for(size_t chan = 0; chan < num_channels; chan++) {
        assert(fwrite(rx_buffer[chan],
                      sizeof(std::complex<float>),
                      num_samples_read,
                      fp[chan]) == num_samples_read);
      }
    }

    // process the received samples through framesync
    sync_state = (data->fs)->execute(rx_buffer, num_samples_read);
    if(sync_state == STATE_WAIT) {
      // stop rx streaming
      stream_cmd.stream_mode = 
        uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
      rx_stream->issue_stream_cmd(stream_cmd);
      // compute channel
      (data->fs)->estimate_channel();
      (data->fs)->get_G(G);
      design_mimo_precoder();
      // start streaming
      stream_cmd.stream_mode = 
        uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS;
      stream_cmd.stream_now = false;
      stream_cmd.time_spec = (*(data->rx))->get_time_now() 
                           + uhd::time_spec_t(0.01);
      rx_stream->issue_stream_cmd(stream_cmd);
      // lock mutex and signal tx to start zf
      std::cout << "Rx signalling tx to start zf\n";
      pthread_mutex_lock(&mutex_txrx);
      start_mimo = true;
      pthread_cond_signal(&condition_mimo);
      pthread_mutex_unlock(&mutex_txrx);
    }

    // check of stop_rx is set
    pthread_mutex_lock(&mutex_txrx);
    if(stop_rx) {
      std::cout << "stop_rx is set\n";
      break_loop = true;
    }
    pthread_mutex_unlock(&mutex_txrx);
    if(break_loop)
      break;
  }

  rx_end = time(NULL);
  // stop rx streaming
  stream_cmd.stream_mode = 
    uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
  rx_stream->issue_stream_cmd(stream_cmd);

  // free memory
  for(size_t chan = 0; chan < num_channels; chan++) {
    free(rx_buffer[chan]);
  }
  // close files
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++)
      fclose(fp[chan]);
  }
  std::cout << "Exiting rx thread\n";
  pthread_exit(NULL);
}

// tx_worker thread
void * tx_worker (void * _data)
{
  tx_thread_data * data = (tx_thread_data *)_data;
  // tx_buffer_len = length of sync word
  unsigned int tx_buffer_len = 
    (CP_LENGTH + d_options.M)*
    (NUM_STREAMS*NUM_ACCESS_CODES + 1);
  // number of channels available with USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples sent in each call of send
  unsigned int num_samples_sent;
  // number of samples read from framegen
  unsigned int num_samples_read;
  // usrp buffer
  std::vector<std::complex<float> *> tx_buffer;
  // framegen buffer
  std::vector<std::complex<float> *> fg_buffer;
  // tx multiplier gain
  std::vector<std::complex<float> > bb_gain;
  // log files
  std::vector<FILE *> fp;

  num_channels = (*(data->tx))->get_tx_num_channels();
  assert(num_channels == NUM_STREAMS);
  size_t volk_alignment = volk_get_alignment();
  std::complex<float> Z(0.0, 0.0);
  for(size_t chan = 0; chan < num_channels; chan++) {
    channels.push_back(chan);
    tx_buffer.push_back((std::complex<float> *) volk_malloc 
      (sizeof(std::complex<float>)*tx_buffer_len,
       volk_alignment));
    fg_buffer.push_back((std::complex<float> *) volk_malloc 
      (sizeof(std::complex<float>)*tx_buffer_len,
       volk_alignment));
    std::fill(tx_buffer[chan], tx_buffer[chan] + tx_buffer_len,
              Z);
    if(LOG) {
      fp.push_back(fopen(boost::str(
        boost::format("%stx%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
    }
    bb_gain.push_back(Z);
  }

  // reset baseband gains
  bb_gain[0] = BASEBAND_GAIN;
  bb_gain[1] = BASEBAND_GAIN;

  // initialize tx streamer
  uhd::stream_args_t tx_stream_args(CPU, WIRE);
  tx_stream_args.channels = channels;
  uhd::tx_streamer::sptr tx_stream = 
    (*(data->tx))->get_tx_stream(tx_stream_args);

  /*
   * Lock mutex and wait for signal.  
   * Note that the pthread_cond_wait routine will 
   * automatically and atomically unlock mutex while it waits.
   * Also, note that if start_tx is true before this 
   * routine is run by the tx thread, the loop will be 
   * skipped to prevent pthread_cond_wait from never returning. 
   */
  pthread_mutex_lock(&mutex_txrx);
  while (!(start_ac)) {
    // blocks tx thread
    std::cout << "Waiting for signal from rx\n";
    pthread_cond_wait(&condition_ac, &mutex_txrx);
    std::cout << "Sigal received from rx thread\n";
  }
  pthread_mutex_unlock(&mutex_txrx);

  // begin tx
  uhd::tx_metadata_t txmd;
  txmd.start_of_burst = true;
  txmd.end_of_burst = false;
  txmd.has_time_spec = true;
  txmd.time_spec = (*(data->tx))->get_time_now()
                 + uhd::time_spec_t(0.1);

  tx_begin = time(NULL);

  // transmit zeros to flush out the buffers
  num_samples_sent = tx_stream->send(tx_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(tx_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }

  num_samples_read = (data->fg)->write_sync_words(fg_buffer);
  assert(num_samples_read <= tx_buffer_len);
  // apply baseband gain
  // FIXME see if this multiply const can be done in place.
  // see if it is benificial to do this in separate
  // thread for each channel.
  for(size_t chan = 0; chan < num_channels; chan++) {
    volk_32fc_s32fc_multiply_32fc(tx_buffer[chan],
                                  fg_buffer[chan],
                                  bb_gain[chan],
                                  num_samples_read);
    // fill fg_buffer with zeros
    std::fill(fg_buffer[chan], fg_buffer[chan] + tx_buffer_len,
              Z);
  }
  // set tx metadata
  txmd.start_of_burst = false;
  txmd.has_time_spec = false;

  // send sync words
  num_samples_sent = tx_stream->send(tx_buffer,
                                     num_samples_read,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == num_samples_read);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(tx_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }
  // stop streaming
  txmd.end_of_burst = true;
  num_samples_sent = tx_stream->send(fg_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(fg_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }
  
  pthread_mutex_lock(&mutex_txrx);
  while (!(start_mimo)) {
    // blocks tx thread
    std::cout << "Waiting for signal from rx to begin zf\n";
    pthread_cond_wait(&condition_mimo, &mutex_txrx);
    std::cout << "Begin zf sigal received from rx thread\n";
  }
  pthread_mutex_unlock(&mutex_txrx);

  // update the mimo beamformer
  (data->fg)->set_W(W);
  //************* MIMO ******************//
  std::cout << "***** begin mimo transmission ******\n";
  txmd.start_of_burst = true;
  txmd.end_of_burst = false;
  txmd.has_time_spec = true;
  txmd.time_spec = (*(data->tx))->get_time_now()
                 + uhd::time_spec_t(0.1);
  num_samples_sent = tx_stream->send(fg_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(fg_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }
  txmd.start_of_burst = false;
  txmd.has_time_spec = false;

  for(pid = 0; pid < PID_MAX; pid++) {
    num_samples_read = (data->fg)->write_mimo_packet(fg_buffer);
    assert(num_samples_read <= tx_buffer_len);
    // apply baseband gain
    // FIXME see if this multiply const can be done in place.
    // see if it is benificial to do this in separate
    // thread for each channel.
    for(size_t chan = 0; chan < num_channels; chan++) {
      volk_32fc_s32fc_multiply_32fc(tx_buffer[chan],
                                    fg_buffer[chan],
                                    bb_gain[chan],
                                    num_samples_read);
    }
    num_samples_sent = tx_stream->send(tx_buffer,
                                       num_samples_read,
                                       txmd,
                                       1.0);
    assert(num_samples_sent == num_samples_read);
    // logging
    if(LOG) {
      for(size_t chan = 0; chan < num_channels; chan++) {
        assert(fwrite(tx_buffer[chan],
                      sizeof(std::complex<float>),
                      num_samples_sent,
                      fp[chan]) == num_samples_sent);
      }
    }
  }
  // stop streaming
  for(size_t chan = 0; chan < num_channels; chan++)
    std::fill(tx_buffer[chan], tx_buffer[chan] + tx_buffer_len, Z);
  txmd.end_of_burst = true;
  num_samples_sent = tx_stream->send(tx_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(tx_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }
  
  // sleep for some time and signal rx to stop
  tx_end = time(NULL);
  usleep(1000000);
  pthread_mutex_lock(&mutex_txrx);
  stop_rx = true;
  pthread_mutex_unlock(&mutex_txrx);

  // free memory
  for(size_t chan = 0; chan < num_channels; chan++) {
    volk_free(tx_buffer[chan]);
    volk_free(fg_buffer[chan]);
  }
  // close files
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++)
      fclose(fp[chan]);
  }
  // return
  std::cout << "Exiting tx thread\n";
  pthread_exit(NULL);
}

// rx beamforming setup. no CSI feedback involved.
#else
// structure to load data to tx_thread
typedef struct
{
  uhd::usrp::multi_usrp::sptr * tx;
  rx_beamforming::framegen * fg;
} tx_thread_data;

// structure to load data to rx_thread
typedef struct
{
  uhd::usrp::multi_usrp::sptr * rx;
  rx_beamforming::framesync * fs;
} rx_thread_data;

// rx_worker thread
void * rx_worker (void * _data)
{
  rx_thread_data * data = (rx_thread_data *)_data;
  // rx buffer length
  unsigned int rx_buffer_len;
  // number of channels available on USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples read in each call of recv
  unsigned int num_samples_read;
  // total number of samples read
  unsigned long int num_accumulated_samples;
  // usrp buffer
  std::vector<std::complex<float> *> rx_buffer;
  // log files
  std::vector<FILE *> fp;

  num_channels = (*(data->rx))->get_rx_num_channels();
  assert(num_channels == NUM_STREAMS);
  for(size_t chan = 0; chan < num_channels; chan++) {
    channels.push_back(chan);
    fp.push_back(fopen(boost::str(
      boost::format("%srx%d.dat") % LOG_DIR % (chan + 1)).c_str(),
      "wb"));
  }

  // initialize rx streamer
  uhd::stream_args_t rx_stream_args(CPU, WIRE);
  rx_stream_args.channels = channels;
  uhd::rx_streamer::sptr rx_stream = 
    (*(data->rx))->get_rx_stream(rx_stream_args);

  // initialilze rx_buffer_len, allocate memory
  rx_buffer_len = rx_stream->get_max_num_samps();
  for(size_t chan = 0; chan < num_channels; chan++) {
    // FIXME is volk malloc required here??
    rx_buffer.push_back((std::complex<float> *) malloc 
      (sizeof(std::complex<float>)*rx_buffer_len));
  }

  // start receiving in 0.1 seconds
  uhd::stream_cmd_t stream_cmd(
      uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
  stream_cmd.stream_now = false;
  stream_cmd.time_spec = (*(data->rx))->get_time_now() 
                       + uhd::time_spec_t(0.1);
  rx_stream->issue_stream_cmd(stream_cmd);
  // lock mutex and signal tx to start
  std::cout << "Rx locking the mutex and signalling tx\n";
  pthread_mutex_lock(&mutex_txrx);
  start_ac = true;
  pthread_cond_signal(&condition_ac);
  pthread_mutex_unlock(&mutex_txrx);

  rx_begin = time(NULL);
  // NOTE: timout should be larger if rx_buffer_length
  // is larger.
  float timeout = 0.2;
  // FIXME on what condition does rx stop streaming?
  // TODO receive sigal from tx to stop
  bool continue_loop = true;
  num_accumulated_samples = 0;
  while(continue_loop)
  {
    num_samples_read = rx_stream->recv(rx_buffer,
                                       rx_buffer_len,
                                       rxmd,
                                       timeout);
    if(rxmd.error_code)
      break; // TODO error code handling to be done in main

    // logging
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(rx_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_read,
                    fp[chan]) == num_samples_read);
    }

    num_accumulated_samples += num_samples_read;

    // check of stop_rx is set
    pthread_mutex_lock(&mutex_txrx);
    if(stop_rx) {
      std::cout << "stop_rx is set\n";
      continue_loop = false;
    }
    pthread_mutex_unlock(&mutex_txrx);
  }

  rx_end = time(NULL);
  // stop rx streaming
  stream_cmd.stream_mode = 
    uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
  rx_stream->issue_stream_cmd(stream_cmd);

  // reallocate memory and reopen files
  for(size_t chan = 0; chan < num_channels; chan++) {
    free(rx_buffer[chan]);
    rx_buffer[chan] = (gr_complex *) malloc 
      (sizeof(gr_complex)*num_accumulated_samples);
    assert(freopen(NULL, "rb", fp[chan]));
    assert(fseek(fp[chan], 0, SEEK_SET) == 0);
    assert(num_accumulated_samples == 
           fread(rx_buffer[chan],
                 sizeof(gr_complex),
                 num_accumulated_samples,
                 fp[chan]));
  }
  // process the received samples through framesync.
  // use computed CSI to invert the channel at the RX.
  // no CSI feedback involved. offline decoding
  (data->fs)->execute(rx_buffer, num_accumulated_samples);

  // close files
  for(size_t chan = 0; chan < num_channels; chan++)
    fclose(fp[chan]);
  // free memory
  for(size_t chan = 0; chan < num_channels; chan++) {
    free(rx_buffer[chan]);
  }

  std::cout << "Exiting rx thread\n";
  pthread_exit(NULL);
}

// tx_worker thread
void * tx_worker (void * _data)
{
  tx_thread_data * data = (tx_thread_data *)_data;
  // tx_buffer_len = length of sync word
  unsigned int tx_buffer_len = 
    (CP_LENGTH + d_options.M)*
    (NUM_STREAMS*NUM_ACCESS_CODES + 1);
  // number of channels available with USRP
  unsigned int num_channels;
  // channel vector
  std::vector<size_t> channels;
  // number of samples sent in each call of send
  unsigned int num_samples_sent;
  // number of samples read from framegen
  unsigned int num_samples_read;
  // usrp buffer
  std::vector<std::complex<float> *> tx_buffer;
  // framegen buffer
  std::vector<std::complex<float> *> fg_buffer;
  // tx multiplier gain
  std::vector<std::complex<float> > bb_gain;
  // log files
  std::vector<FILE *> fp;
  // framegen input buffer
  std::vector<gr_complex *> fg_in_buffer;
  // num of frequency symbols transmitted
  unsigned int num_symbols_transmitted = 0;

  num_channels = (*(data->tx))->get_tx_num_channels();
  assert(num_channels == NUM_STREAMS);
  size_t volk_alignment = volk_get_alignment();
  std::complex<float> Z(0.0, 0.0);
  for(size_t chan = 0; chan < num_channels; chan++) {
    channels.push_back(chan);
    tx_buffer.push_back((std::complex<float> *) volk_malloc 
      (sizeof(std::complex<float>)*tx_buffer_len,
       volk_alignment));
    fg_buffer.push_back((std::complex<float> *) volk_malloc 
      (sizeof(std::complex<float>)*tx_buffer_len,
       volk_alignment));
    fg_in_buffer.push_back((std::complex<float> *) malloc 
      (sizeof(std::complex<float>)*num_occupied_carriers));
    std::fill(tx_buffer[chan], tx_buffer[chan] + tx_buffer_len,
              Z);
    if(LOG) {
      fp.push_back(fopen(boost::str(
        boost::format("%stx%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
    }
    bb_gain.push_back(Z);
  }

  // reset baseband gains
  bb_gain[0] = BASEBAND_GAIN;
  bb_gain[1] = BASEBAND_GAIN;

  // initialize tx streamer
  uhd::stream_args_t tx_stream_args(CPU, WIRE);
  tx_stream_args.channels = channels;
  uhd::tx_streamer::sptr tx_stream = 
    (*(data->tx))->get_tx_stream(tx_stream_args);

  /*
   * Lock mutex and wait for signal.  
   * Note that the pthread_cond_wait routine will 
   * automatically and atomically unlock mutex while it waits.
   * Also, note that if start_tx is true before this 
   * routine is run by the tx thread, the loop will be 
   * skipped to prevent pthread_cond_wait from never returning. 
   */
  pthread_mutex_lock(&mutex_txrx);
  while (!(start_ac)) {
    // blocks tx thread
    std::cout << "Waiting for signal from rx\n";
    pthread_cond_wait(&condition_ac, &mutex_txrx);
    std::cout << "Sigal received from rx thread\n";
  }
  pthread_mutex_unlock(&mutex_txrx);

  // begin tx
  uhd::tx_metadata_t txmd;
  txmd.start_of_burst = true;
  txmd.end_of_burst = false;
  txmd.has_time_spec = true;
  txmd.time_spec = (*(data->tx))->get_time_now()
                 + uhd::time_spec_t(0.1);

  tx_begin = time(NULL);

  // transmit zeros to flush out the buffers
  num_samples_sent = tx_stream->send(tx_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(tx_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }

  num_samples_read = (data->fg)->write_sync_words(fg_buffer);
  assert(num_samples_read <= tx_buffer_len);
  // apply baseband gain
  // FIXME see if this multiply const can be done in place.
  // see if it is benificial to do this in separate
  // thread for each channel.
  for(size_t chan = 0; chan < num_channels; chan++) {
    volk_32fc_s32fc_multiply_32fc(tx_buffer[chan],
                                  fg_buffer[chan],
                                  bb_gain[chan],
                                  num_samples_read);
  }
  // set tx metadata
  txmd.start_of_burst = false;
  txmd.has_time_spec = false;

  // send sync words
  num_samples_sent = tx_stream->send(tx_buffer,
                                     num_samples_read,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == num_samples_read);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(tx_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }

  for(pid = 0;
      ((pid < PID_MAX) && 
       (num_symbols_transmitted + num_occupied_carriers <= tx_sig_len));
      pid++) {
    for(unsigned int chan = 0; chan < num_channels; chan++) {
      memmove(fg_in_buffer[chan],
              tx_sig[chan] + num_symbols_transmitted,
              sizeof(gr_complex)*num_occupied_carriers);
    }
    num_symbols_transmitted += num_occupied_carriers;
    num_samples_read = 
      (data->fg)->assemble_mimo_packet(fg_buffer, fg_in_buffer);
    assert(num_samples_read <= tx_buffer_len);
    // apply baseband gain
    // FIXME see if this multiply const can be done in place.
    // see if it is benificial to do this in separate

    // thread for each channel.
    for(size_t chan = 0; chan < num_channels; chan++) {
      volk_32fc_s32fc_multiply_32fc(tx_buffer[chan],
                                    fg_buffer[chan],
                                    bb_gain[chan],
                                    num_samples_read);
    }
    num_samples_sent = tx_stream->send(tx_buffer,
                                       num_samples_read,
                                       txmd,
                                       1.0);
    assert(num_samples_sent == num_samples_read);
    // logging
    if(LOG) {
      for(size_t chan = 0; chan < num_channels; chan++) {
        assert(fwrite(tx_buffer[chan],
                      sizeof(std::complex<float>),
                      num_samples_sent,
                      fp[chan]) == num_samples_sent);
      }
    }
  }
  // stop streaming
  for(size_t chan = 0; chan < num_channels; chan++)
    std::fill(tx_buffer[chan], tx_buffer[chan] + tx_buffer_len, Z);
  txmd.end_of_burst = true;
  num_samples_sent = tx_stream->send(tx_buffer,
                                     tx_buffer_len,
                                     txmd,
                                     1.0);
  assert(num_samples_sent == tx_buffer_len);
  // logging
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++) {
      assert(fwrite(tx_buffer[chan],
                    sizeof(std::complex<float>),
                    num_samples_sent,
                    fp[chan]) == num_samples_sent);
    }
  }
  
  // sleep for some time and signal rx to stop
  tx_end = time(NULL);
  usleep(1000000);
  pthread_mutex_lock(&mutex_txrx);
  stop_rx = true;
  pthread_mutex_unlock(&mutex_txrx);

  // free memory
  for(size_t chan = 0; chan < num_channels; chan++) {
    volk_free(tx_buffer[chan]);
    volk_free(fg_buffer[chan]);
    free(fg_in_buffer[chan]);
  }
  // close files
  if(LOG) {
    for(size_t chan = 0; chan < num_channels; chan++)
      fclose(fp[chan]);
  }
  // return
  std::cout << "Exiting tx thread\n";
  pthread_exit(NULL);
}
#endif

// main
int UHD_SAFE_MAIN(int argc, char **argv)
{
  uhd::set_thread_priority_safe();

  // read command line options

  read_options(argc, argv, &d_options);

  // CHANGED: Use Options instead of Hardcoded Definings
  unsigned int M = d_options.M;
  unsigned int cp_len = d_options.cp_len;


  unsigned int num_streams = NUM_STREAMS;


  // allocate space for CSI and W
  gr_complex I(1.0, 0.0);

  gr_complex _G_[d_options.M][NUM_STREAMS][NUM_STREAMS];
  gr_complex _W_[d_options.M][NUM_STREAMS][NUM_STREAMS];
  G = (gr_complex ***)_G_;
  W = (gr_complex ***)_W_;

  unsigned char * p;
  msequence ms_S0;
  std::vector<msequence> ms_S1;
  std::vector<modem> mod;
  std::vector<modem> dem;
  std::vector<FILE *> tx_data_fp;
  std::vector<FILE *> tx_sig_fp;
  std::vector<FILE *> rx_data_fp;
  std::vector<FILE *> rx_sig_fp;

  srand(time(NULL));

  ms_S1.resize(num_streams);
  p = (unsigned char *) malloc (sizeof(unsigned char)*M);
  ofdmframe_init_default_sctype(p, M);
  ofdmframe_validate_sctype(p, M, &num_null_carriers,
                            &num_pilot_carriers,
                            &num_data_carriers);
  num_occupied_carriers = num_pilot_carriers + num_data_carriers;
  tx_sig_len = num_occupied_carriers*PID_MAX;
  for(unsigned int chan = 0; chan < num_streams; chan++) {
    num_valid_bytes.push_back(0);
    mod.push_back(modem_create(MODEM_SCHEME));
    dem.push_back(modem_create(MODEM_SCHEME));
    tx_sig.push_back((gr_complex *) malloc 
                     (sizeof(gr_complex)*tx_sig_len));
    rx_sig.push_back((gr_complex *) malloc 
                     (sizeof(gr_complex)*PID_MAX*(M + cp_len)));
    tx_data.push_back((unsigned int *) malloc 
                      (sizeof(unsigned int)*tx_sig_len));
    rx_data.push_back((unsigned int *) malloc 
                      (sizeof(unsigned int)*tx_sig_len));
    std::fill(tx_data[chan], tx_data[chan] + tx_sig_len, 0);
    std::fill(tx_sig[chan], tx_sig[chan] + tx_sig_len, gr_complex(0.0, 0.0));
#if SISO
    if(chan == SISO_TX) {
      for(unsigned int sample = 0; sample < tx_sig_len; sample++) {
        tx_data[chan][sample] = rand() % ARITY;
        modem_modulate(mod[chan], tx_data[chan][sample], tx_sig[chan] + sample);
      }
    }
#else
#if SAME_SIGNAL_ON_ALL_TX
    if(chan == 0) {
      for(unsigned int sample = 0; sample < tx_sig_len; sample++) {
        tx_data[chan][sample] = rand() % ARITY;
        modem_modulate(mod[chan], tx_data[chan][sample], tx_sig[chan] + sample);
      }
    }
    else {
      memmove(tx_data[chan], tx_data[0], sizeof(unsigned int)*tx_sig_len);
      memmove(tx_sig[chan], tx_sig[0], sizeof(gr_complex)*tx_sig_len);
    }
#else
    for(unsigned int sample = 0; sample < tx_sig_len; sample++) {
      tx_data[chan][sample] = rand() % ARITY;
      modem_modulate(mod[chan], tx_data[chan][sample], tx_sig[chan] + sample);
    }
#endif
#endif
  }
#if LOG
  for(unsigned int chan = 0; chan < num_streams; chan++) {
      tx_data_fp.push_back(fopen(boost::str(
        boost::format("%stx_data%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
      tx_sig_fp.push_back(fopen(boost::str(
        boost::format("%stx_sig%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
      rx_data_fp.push_back(fopen(boost::str(
        boost::format("%srx_data%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
      rx_sig_fp.push_back(fopen(boost::str(
        boost::format("%srx_sig%d.dat") % LOG_DIR % (chan + 1)).c_str(),
        "wb"));
      fwrite(tx_data[chan], sizeof(unsigned int), tx_sig_len,
             tx_data_fp[chan]);
      fwrite(tx_sig[chan], sizeof(gr_complex), tx_sig_len,
             tx_sig_fp[chan]);
      fclose(tx_data_fp[chan]);
      fclose(tx_sig_fp[chan]);
  }
#endif

  // NOTE: initialize the sequence generator for other streams as well
  // in case the number of streams is more that 2. Here we assume that
  // we only have 2 streams.
  ms_S0    = msequence_create(LFSR_SMALL_LENGTH, LFSR_SMALL_0_GEN_POLY, 1);
  ms_S1[0] = msequence_create(LFSR_LARGE_LENGTH, LFSR_LARGE_0_GEN_POLY, 1);
  ms_S1[1] = msequence_create(LFSR_LARGE_LENGTH, LFSR_LARGE_1_GEN_POLY, 1);
#if TX_BEAMFORMING
  framegen fg(M,
              cp_len,
              num_streams,
              NUM_ACCESS_CODES,
              p,
              ms_S0,
              ms_S1);
  fg.set_W(W);
  msequence_reset(ms_S0);
  msequence_reset(ms_S1[0]);
  msequence_reset(ms_S1[1]);
  framesync fs(M,
               cp_len,
               num_streams,
               NUM_ACCESS_CODES,
               p,
               ms_S0,
               ms_S1,
               callback);
#else
  rx_beamforming::
  framegen fg(M,
              cp_len,
              num_streams,
              NUM_ACCESS_CODES,
              p,
              ms_S0,
              ms_S1);
  msequence_reset(ms_S0);
  msequence_reset(ms_S1[0]);
  msequence_reset(ms_S1[1]);
  rx_beamforming::
  framesync fs(M,
               cp_len,
               num_streams,
               NUM_ACCESS_CODES,
               p,
               ms_S0,
               ms_S1,
               callback);
#if SISO
  fs.set_siso_tx(SISO_TX);
  fs.set_siso_rx(SISO_RX);
#endif
#endif
  // destroy msequence ojbects.
  msequence_destroy(ms_S0);
  msequence_destroy(ms_S1[0]);
  msequence_destroy(ms_S1[1]);

  usrp_config * tx_config;
  usrp_config * rx_config;
  tx_thread_data * tx_thr_data;
  rx_thread_data * rx_thr_data;
  uhd::usrp::multi_usrp::sptr tx;
  uhd::usrp::multi_usrp::sptr rx;

  num_frames_detected         = 0;
  num_valid_headers_received  = 0;
  num_valid_bytes_received    = 0;
  num_valid_packets_received  = 0;
  received_sample_counter     = 0;

  pthread_t tx_thread, rx_thread;

  /// CHANGED:
  /// removed Hardcoded Adress and added options
  // initialize USRP
  rx = uhd::usrp::multi_usrp::make(
       uhd::device_addr_t(d_options.tx_adress.c_str()));
  tx = uhd::usrp::multi_usrp::make(
       uhd::device_addr_t(d_options.rx_adress.c_str()));
  tx_config = (usrp_config *) malloc 
              (sizeof(usrp_config));
  rx_config = (usrp_config *) malloc 
              (sizeof(usrp_config));
  // assign tx and rx sampling rate freq etc
  tx_config->cent_freq = d_options.cent_freq;
  rx_config->cent_freq = d_options.cent_freq;
  tx_config->samp_rate = d_options.samp_rate;
  rx_config->samp_rate = d_options.samp_rate;
  tx_config->rf_gain   = d_options.txgain;
  rx_config->rf_gain   = d_options.rxgain;
  // TODO Add command line options 
  // for clock and time source
  tx_config->clock_source = CLOCK_SOURCE;
  rx_config->clock_source = CLOCK_SOURCE;
  tx_config->time_source = TIME_SOURCE;
  rx_config->time_source = TIME_SOURCE;

  init_usrp(tx, tx_config, true);
  init_usrp(rx, rx_config, false);

  // initialize worker threads
  tx_thr_data = (tx_thread_data *) malloc (sizeof(tx_thread_data));
  rx_thr_data = (rx_thread_data *) malloc (sizeof(rx_thread_data));
  tx_thr_data->tx = &tx;
  rx_thr_data->rx = &rx;
  tx_thr_data->fg = &fg;
  rx_thr_data->fs = &fs;

  if(pthread_create(&tx_thread, NULL, tx_worker, (void *)tx_thr_data)){
    printf("Error invoking tx thread\n");
    return 1;
  }
  else
    printf("****** tx thread invoked ********\n");
  if(pthread_create(&rx_thread, NULL, rx_worker, (void *)rx_thr_data)){
    printf("Error invoking rx thread\n");
    return 1;
  }
  else
    printf("****** rx thread invoked ********\n");
  pthread_join(tx_thread, NULL);
  pthread_join(rx_thread, NULL);

  // free malloc
  free(tx_thr_data);
  free(rx_thr_data);
  free(tx_config);
  free(rx_config);

  // demodulation
#if SISO
  for(unsigned int sample = 0; sample < tx_sig_len; sample++) {
    modem_demodulate(dem[SISO_RX], rx_sig[SISO_RX][sample], rx_data[SISO_RX] + sample);
    if(rx_data[SISO_RX][sample] == tx_data[SISO_TX][sample]) {
      num_valid_bytes_received++;
    }
  }
#else
  for(unsigned int chan = 0; chan < num_streams; chan++) {
    for(unsigned int sample = 0; sample < tx_sig_len; sample++) {
      modem_demodulate(dem[chan], rx_sig[chan][sample], rx_data[chan] + sample);
      if(rx_data[chan][sample] == tx_data[chan][sample]) {
        num_valid_bytes[chan]++;
      }
    }
  }
#endif
  for(unsigned int chan = 0; chan < num_streams; chan++) {
#if LOG
    fwrite(rx_data[chan], sizeof(unsigned int), tx_sig_len,
             rx_data_fp[chan]);
    fwrite(rx_sig[chan], sizeof(gr_complex), tx_sig_len,
             rx_sig_fp[chan]);
    fclose(rx_data_fp[chan]);
    fclose(rx_sig_fp[chan]);
#endif
    modem_destroy(mod[chan]);
    modem_destroy(dem[chan]);
    free(tx_data[chan]);
    free(rx_data[chan]);
    free(tx_sig[chan]);
    free(rx_sig[chan]);
  }

  // print experiment output
  printf("    plateau width 1         : %6lu\n", fs.get_plateau_end(0) - 
                                                 fs.get_plateau_start(0) + 1);
  printf("    plateau start 1         : %6lu\n", fs.get_plateau_start(0));
  printf("    plateau end   1         : %6lu\n", fs.get_plateau_end(0));
  printf("    plateau width 2         : %6lu\n", fs.get_plateau_end(1) - 
                                                 fs.get_plateau_start(1) + 1);
  printf("    plateau start 2         : %6lu\n", fs.get_plateau_start(1));
  printf("    plateau end   2         : %6lu\n", fs.get_plateau_end(1));
  printf("    *************************\n");
  printf("    frames sync index 2     : %6lu\n", fs.get_sync_index());
  printf("    num samples processed   : %6llu\n", fs.get_num_samples_processed());
  printf("    frames transmitted      : %6u\n", pid);
  printf("    num_occupied_carriers   : %6u\n", num_occupied_carriers);
#if SISO
  printf("    symbols transmitted     : %6lu\n", tx_sig_len);
#else
  printf("    symbols transmitted     : %6lu\n", num_streams*tx_sig_len);
#endif
  printf("    frames detected         : %6u\n", num_frames_detected);
  printf("    valid headers           : %6u\n", num_valid_headers_received);
  printf("    valid packets           : %6u\n", num_valid_packets_received);
#if SISO
  printf("    valid symbols received  : %6u\n", num_valid_bytes_received);
  printf("    symbol error rate       : %1.6f%%\n",
         (float(tx_sig_len - num_valid_bytes_received)/float(tx_sig_len))*100);
#else
  for(unsigned int chan = 0; chan < num_streams; chan++) {
    printf("    valid symbols received %u: %6u\n", chan, num_valid_bytes[chan]);
    printf("    symbol error rate      %u: %1.6f%%\n", chan,
           (float(tx_sig_len - num_valid_bytes[chan])/float(tx_sig_len))*100);
  }
#endif
  printf("    tx run time             : %6lu\n", tx_end - tx_begin);
  printf("    rx run time             : %6lu\n", rx_end - rx_begin);
  printf("    bit rate                : %e bps\n",
         float(num_valid_bytes_received*8)/(rx_end - rx_begin));
  if(rxmd.error_code)
  {
    printf("    rx thread exit with     : %s\n", rxmd.strerror().c_str());
  }
  return 0;
}
