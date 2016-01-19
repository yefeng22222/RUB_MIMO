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

#ifndef FRAMING_H
#define FRAMING_H

#include <vector>
#include <fftw3.h>
#include <gnuradio/gr_complex.h>
#include <complex>
#include <liquid/liquid.h>
#include <boost/format.hpp>
#include "config.h"

// callback
typedef void * (* mimo_callback) (std::vector<gr_complex *>,
                                  unsigned int occupied_carriers);

  // receiver state
typedef enum {
  STATE_SEEK_PLATEAU = 0,
  STATE_SAVE_ACCESS_CODES,
  STATE_WAIT,
  STATE_MIMO
} framesync_states_t;

namespace rx_beamforming {
  class framegen {
   private:
    // number of subcarriers
    unsigned int M;
    // cyclic prefix length
    unsigned int cp_len;
    // symbol length
    unsigned int symbol_len;
    // number of data streams
    unsigned int num_streams;
    // number of access codes
    unsigned int num_access_codes;
    // subcarrier allocation
    unsigned char * p;
    // number of null subcarriers
    unsigned int M_null;
    // number of pilot subcarriers
    unsigned int M_pilot;
    // number of data subcarriers
    unsigned int M_data;
    // frequency domain buffer
    std::vector<gr_complex *> X;
    // time domain buffer
    std::vector<gr_complex *> x;
    // scaling factors
    gr_complex g_data;
    // transform object
    std::vector<fftwf_plan> ifft;
    gr_complex dft_normalizer;
    // PLCP short
    gr_complex * S0;
    gr_complex * s0;
    // PLCP long
    std::vector<gr_complex *> S1;
    std::vector<gr_complex *> s1;
  
    // volk buffer for volk operations
    std::complex<float> * volk_buff_fc1;
    std::complex<float> * volk_buff_fc2;
    std::complex<float> * volk_buff_fc3;
    float * volk_buff_f1;
  
   public:
    // constructor
    framegen(unsigned int _M,
             unsigned int _cp_len,
             unsigned int _num_streams,
             unsigned int _num_access_codes,
             unsigned char * const &_p,
             msequence const &_ms_S0,
             std::vector<msequence> const &_ms_S1);
    // destructor
    ~framegen();
    void print();
  
    unsigned int
    write_sync_words(std::vector<std::complex<float> *> tx_buff);
    unsigned int
      assemble_mimo_packet(std::vector<gr_complex *> tx_buff,
                           std::vector<gr_complex *> in_buff);
    unsigned int get_num_streams();
  };
  
  class framesync {
   private:
    // number of subcarriers
    unsigned int M;
    // M/2
    unsigned int M2;
    // cyclic prefix length
    unsigned int cp_len;
    // symbol length
    unsigned int symbol_len;
    // number of data streams
    unsigned int num_streams;
    // number of access codes
    unsigned int num_access_codes;
    // subcarrier allocation
    unsigned char * p;
    // number of null subcarriers
    unsigned int M_null;
    // number of pilot subcarriers
    unsigned int M_pilot;
    // number of data subcarriers
    unsigned int M_data;
    unsigned int M_occupied;
    // counter
    unsigned int sample_counter;
    // frequency domain buffer
    std::vector<gr_complex *> X;
    // time domain buffer
    std::vector<gr_complex *> x;
    // ofdm symbol
    std::vector<gr_complex *> symbol;
    // CSI
    std::vector<std::vector<std::vector<gr_complex> > > G;
    // Receive beamformer
    std::vector<std::vector<std::vector<gr_complex> > > W;
    // transform object
    std::vector<fftwf_plan> fft;
    gr_complex dft_normalizer;
    // PLCP short
    gr_complex * S0;
    gr_complex * s0;
    // PLCP long
    std::vector<gr_complex *> S1;
    std::vector<gr_complex *> s1;
    // normalizer
    float * normalize_gain;

    // callback
    mimo_callback callback;
  
    // volk buffer for volk operations
    std::complex<float> * volk_buff_fc1;
    std::complex<float> * volk_buff_fc2;
    std::complex<float> * volk_buff_fc3;
    float * volk_buff_f1;

    unsigned int siso_tx;
    unsigned int siso_rx;
  
    void increment_state();
  
    // sync 
    unsigned long int sync_index;
    unsigned long long int num_samples_processed;
    void execute_sc_sync(gr_complex _x[]);
    void execute_save_access_codes(gr_complex _x[]);
    framesync_states_t state;
    void execute_mimo_decode(gr_complex _x[]);
    void execute_siso_decode(gr_complex _x[]);
    float execute_sc_sync(gr_complex _x, unsigned int stream);

    // Schmidl & Cox
    std::vector<wdelaycf> delay_M2;
    std::vector<windowcf> access_code_buffer;
    unsigned int access_code_buffer_len;
    unsigned int tx_sig_len;
    std::vector<firfilt_crcf> crosscorrelator;
    std::vector<firfilt_rrrf> normalizer;
    std::vector<unsigned long int> plateau_start;
    std::vector<unsigned long int> plateau_end;
    std::vector<bool> in_plateau;
    std::vector<FILE *> f_sc_debug_out;
   public:
    // constructor
    framesync(unsigned int _M,
              unsigned int _cp_len,
              unsigned int _num_streams,
              unsigned int _num_access_codes,
              unsigned char * const &_p,
              msequence const &_ms_S0,
              std::vector<msequence> const &_ms_S1,
              mimo_callback _callback);
    // destructor
    ~framesync();
    void print();
    unsigned long int get_sync_index();
    std::vector<std::vector<std::vector<gr_complex> > > get_G();
    unsigned long long int get_num_samples_processed();
    framesync_states_t execute(std::vector<gr_complex *> const &in_buff,
                             unsigned int num_samples);
    unsigned long int get_plateau_start(unsigned int stream);
    unsigned long int get_plateau_end(unsigned int stream);
    void reset();
    void estimate_channel();
    void compute_receive_beamformer();

    void set_siso_tx(unsigned int _tx);
    void set_siso_rx(unsigned int _rx);
  };
} // namespace rx_beamforming

namespace tx_beamforming {
}

// initialize default subcarrier allocation
//  _M      :   number of subcarriers
//  _p      :   output subcarrier allocation vector, [size: _M x 1]
//
// key: '.' (null), 'P' (pilot), '+' (data)
// .+++P+++++++P.........P+++++++P+++
//
void ofdmframe_init_default_sctype(unsigned char *_p, unsigned int _M);

// validate subcarrier type (count number of null, pilot, and data
// subcarriers in the allocation)
//  _p          :   subcarrier allocation array, [size: _M x 1]
//  _M_null     :   output number of null subcarriers
//  _M_pilot    :   output number of pilot subcarriers
//  _M_data     :   output number of data subcarriers
void ofdmframe_validate_sctype(const unsigned char * _p,
                               unsigned int _M,
                               unsigned int * _M_null,
                               unsigned int * _M_pilot,
                               unsigned int * _M_data);

// print subcarrier allocation to screen
//
// key: '.' (null), 'P' (pilot), '+' (data)
// .+++P+++++++P.........P+++++++P+++
//
void ofdmframe_print_sctype(const unsigned char * _p, unsigned int _M);

// generate short sequence symbols
//  _p      :   subcarrier allocation array
//  _S0     :   output symbol (freq)
//  _s0     :   output symbol (time)
//  _M_S0   :   total number of enabled subcarriers in S0
//  ms      :   pn_sequence generator to generate S0
void ofdmframe_init_S0(const unsigned char * _p,
                       unsigned int _M,
                       std::complex<float> * _S0,
                       std::complex<float> * _s0,
                       msequence ms);

// generate long sequence symbols
//  _p      :   subcarrier allocation array
//  _S1     :   output symbol (freq)
//  _s1     :   output symbol (time)
//  _M_S1   :   total number of enabled subcarriers in S1
//  ms      :   pn_sequence generator to generate S1
void ofdmframe_init_S1(const unsigned char * _p,
                       unsigned int _M,
                       unsigned int _num_access_codes,
                       std::complex<float> * _S1,
                       std::complex<float> * _s1,
                       msequence ms);

gr_complex liquid_cexpjf(float theta);
float cabsf(gr_complex z);
float cargf(gr_complex z);
float fabsf(float x);
gr_complex conjf(gr_complex z);

// currently only for 2 X 2 matrix
float invert(std::vector<std::vector<gr_complex> > &W,
             std::vector<std::vector<gr_complex> > const &G);

#endif // FRAMING_H
