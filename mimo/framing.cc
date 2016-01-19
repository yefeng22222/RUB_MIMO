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

#include "framing.h"
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <volk/volk.h>
#include <volk/volk_malloc.h>
#include <string.h>

float sqrt2 = sqrtf(2.0);

gr_complex CHANNEL_I[2][2] =
{
  {gr_complex(1.0, 0.0),
   gr_complex(0.0, 0.0)},
  {gr_complex(0.0, 0.0),
   gr_complex(1.0, 0.0)}
};
gr_complex BPSK_CONSTELLATION[] = 
{
  gr_complex(-1.0, 0.0),
  gr_complex(1.0,  0.0),
};
gr_complex QPSK_CONSTELLATION[] = 
{
  gr_complex(sqrt2, sqrt2),
  gr_complex(-sqrt2, sqrt2),
  gr_complex(-sqrt2, -sqrt2),
  gr_complex(sqrt2, -sqrt2)
};

inline gr_complex
liquid_cexpjf(float theta)
{
  return std::polar(1.0f, theta);
}

inline float
cabsf(gr_complex z)
{
  return std::abs(z);
}

inline float
cargf(gr_complex z)
{
  return std::arg(z);
}

inline float
fabsf(float x)
{
  return std::abs(x);
}

inline gr_complex
conjf(gr_complex z)
{
  return std::conj(z);
}

namespace rx_beamforming {
  framegen::framegen(unsigned int _M,
                     unsigned int _cp_len,
                     unsigned int _num_streams,
                     unsigned int _num_access_codes,
                     unsigned char * const &_p,
                     msequence const &_ms_S0,
                     std::vector<msequence> const &_ms_S1) 
  {
    // assign values for variables
    M = _M;
    cp_len = _cp_len;
    symbol_len = M + cp_len;
    num_streams = _num_streams;
    num_access_codes = _num_access_codes;

    ifft.resize(num_streams);
    X.resize(num_streams);
    x.resize(num_streams);
    S1.resize(num_streams);
    s1.resize(num_streams);
    p  = (unsigned char *) malloc (sizeof(unsigned char)*M);
    memmove(p, _p, sizeof(unsigned char)*M);
    S0 = (gr_complex *) malloc (sizeof(gr_complex)*M);
    s0 = (gr_complex *) malloc (sizeof(gr_complex)*M);
    // sanity check for subcarrier allocation
    ofdmframe_validate_sctype(p,
                              M,
                              &M_null,
                              &M_pilot,
                              &M_data);
    // initialize S0 sequence
    ofdmframe_init_S0(p,
                      M,
                      S0,
                      s0,
                      _ms_S0);
    dft_normalizer = 1.0f / sqrtf(M_pilot + M_data);
  
    // get volk alignment for volk_malloc
    size_t volk_alignment = volk_get_alignment();
    for(unsigned int i = 0; i < num_streams; i++)
    {
      // memory allocation
      // NOTE: all to be freed in destructor
      S1[i] = (gr_complex *) malloc 
              (sizeof(gr_complex)*
               M*num_access_codes);
      s1[i] = (gr_complex *) malloc 
              (sizeof(gr_complex)*
               M*num_access_codes);
      X[i] = (std::complex<float> *) fftwf_malloc
             (sizeof(std::complex<float>)*M);
      x[i] = (std::complex<float> *) fftwf_malloc
             (sizeof(std::complex<float>)*M);
      assert(volk_is_aligned(X[i]));
      assert(volk_is_aligned(x[i]));
      ifft[i] = fftwf_plan_dft_1d(M,
                                  reinterpret_cast<fftwf_complex *>(X[i]),
                                  reinterpret_cast<fftwf_complex *>(x[i]),
                                  FFTW_BACKWARD,
                                  FFTW_ESTIMATE);
      // initialize S1 sequence
      ofdmframe_init_S1(p,
                        M,
                        num_access_codes,
                        S1[i],
                        s1[i],
                        _ms_S1[i]);
    }
    volk_buff_fc1 = (std::complex<float> *) volk_malloc 
                    (sizeof(std::complex<float>)*M,
                     volk_alignment);
    volk_buff_fc2 = (std::complex<float> *) volk_malloc 
                    (sizeof(std::complex<float>)*M,
                     volk_alignment);
    volk_buff_fc3 = (std::complex<float> *) volk_malloc 
                    (sizeof(std::complex<float>)*M,
                     volk_alignment);
    volk_buff_f1 = (float *) volk_malloc 
                   (sizeof(float)*
                    (M + cp_len)*3*num_streams,
                    volk_alignment);
  }
  
  unsigned int
  framegen::get_num_streams()
  {
    return num_streams;
  }
  
  unsigned int
  framegen::write_sync_words(std::vector<std::complex<float> *> tx_buff)
  {
    assert(num_streams == tx_buff.size());
    unsigned int sample_index = 0;
    unsigned int total_count = (num_access_codes*num_streams + 1)*
                               (M + cp_len);
    // set tx_buff to 0;
    for(unsigned int stream = 0; stream < num_streams; stream++) {
      std::fill(tx_buff[stream],
                tx_buff[stream] + total_count,
                gr_complex(0.0, 0.0));
    }
    // write S0 onto ch0 and then S1 in TDMA
    memmove(tx_buff[0] + sample_index,
            s0 + M - cp_len,
            sizeof(std::complex<float>)*cp_len);
    sample_index += cp_len;
    memmove(tx_buff[0] + sample_index,
            s0,
            sizeof(std::complex<float>)*M);
    sample_index += M;
    for(unsigned int ac_id = 0; ac_id < num_access_codes; ac_id++) {
      for(unsigned int stream = 0; stream < num_streams; stream++) {
        // add cyclic prefix
        memmove(tx_buff[stream] + sample_index,
                s1[stream] + M*ac_id + M - cp_len,
                sizeof(std::complex<float>)*cp_len);
        sample_index += cp_len;
        // write ofdm symbol
        memmove(tx_buff[stream] + sample_index,
                s1[stream] + M*ac_id,
                sizeof(std::complex<float>)*M);
        sample_index += M;
      }
    }
    // total count of sync samples
    assert(sample_index == total_count);
    return sample_index;
  }
  
  unsigned int
  framegen::assemble_mimo_packet(std::vector<gr_complex *> tx_buff,
                                 std::vector<gr_complex *> in_buff)
  {
    assert(tx_buff.size() == num_streams);
    assert(in_buff.size() == num_streams);
    for(unsigned int stream = 0; stream < num_streams; stream++) {
      for(unsigned int i = 0, j = 0; i < M; i++) {
        if(p[i] == OFDMFRAME_SCTYPE_NULL)
          X[stream][i] = 0.0f;
        else
          X[stream][i] = in_buff[stream][j++];
      }
      fftwf_execute(ifft[stream]);
      volk_32fc_s32fc_multiply_32fc(volk_buff_fc1, x[stream], dft_normalizer, M);
      // cyclic prefix
      memmove(tx_buff[stream],
              volk_buff_fc1 + M - cp_len,
              sizeof(gr_complex)*cp_len);
      // ofdm symbol
      memmove(tx_buff[stream] + cp_len,
              volk_buff_fc1,
              sizeof(gr_complex)*M);
    }
    return symbol_len;
  }
  
  framegen::~framegen()
  {
    for(unsigned int i = 0; i < num_streams; i++)
    {
      free(S1[i]);
      free(s1[i]);
      fftwf_free(X[i]);
      fftwf_free(x[i]);
      fftwf_destroy_plan(ifft[i]);
    }
    free(p);
    free(S0);
    free(s0);
    volk_free(volk_buff_fc1);
    volk_free(volk_buff_fc2);
    volk_free(volk_buff_fc3);
    volk_free(volk_buff_f1);
  }
  
  void framegen::print()
  {
      printf("ofdmframegen:\n");
      printf("    num subcarriers     :   %-u\n", M);
      printf("      - NULL            :   %-u\n", M_null);
      printf("      - pilot           :   %-u\n", M_pilot);
      printf("      - data            :   %-u\n", M_data);
      printf("    cyclic prefix len   :   %-u\n", cp_len);
      printf("    ");
      ofdmframe_print_sctype(p, M);
  }
  
  framesync::framesync(unsigned int _M,
                       unsigned int _cp_len,
                       unsigned int _num_streams,
                       unsigned int _num_access_codes,
                       unsigned char * const &_p,
                       msequence const &_ms_S0,
                       std::vector<msequence> const &_ms_S1,
                       mimo_callback _callback) 
  {
    // assign values for variables
    M = _M;
    cp_len = _cp_len;
    symbol_len = M + cp_len;
    num_streams = _num_streams;
    num_access_codes = _num_access_codes;
    M2 = M / 2;
    access_code_buffer_len = symbol_len*(num_access_codes*num_streams + 4);
    tx_sig_len = PID_MAX*symbol_len;
    callback = _callback;
    dft_normalizer = sqrt(1.0/float(M));
  
    // resize vectors to num stream
    G.resize(M);
    W.resize(M);
  
    fft.resize(num_streams);
    X.resize(num_streams);
    x.resize(num_streams);
    symbol.resize(num_streams);
    S1.resize(num_streams);
    s1.resize(num_streams);
    p  = (unsigned char *) malloc (sizeof(unsigned char)*M);
    memmove(p, _p, sizeof(unsigned char)*M);

    for(unsigned int sc = 0; sc < M; sc++) {
      G[sc].resize(num_streams);
      W[sc].resize(num_streams);
      for(unsigned int rx_stream = 0; rx_stream < num_streams; rx_stream++) {
        G[sc][rx_stream].resize(num_streams);
        W[sc][rx_stream].resize(num_streams);
        for(unsigned int tx_stream = 0; tx_stream < num_streams; tx_stream++) {
          if((p[sc] != OFDMFRAME_SCTYPE_NULL) && (tx_stream == rx_stream)) {
            G[sc][rx_stream][tx_stream] = 1.0;
            W[sc][rx_stream][tx_stream] = 1.0;
          }
          else {
            G[sc][rx_stream][tx_stream] = 0.0;
            W[sc][rx_stream][tx_stream] = 0.0;
          }
        }
      }
    }

    S0 = (gr_complex *) malloc (sizeof(gr_complex)*M);
    s0 = (gr_complex *) malloc (sizeof(gr_complex)*M);
    // sanity check for subcarrier allocation
    ofdmframe_validate_sctype(p,
                              M,
                              &M_null,
                              &M_pilot,
                              &M_data);
    M_occupied = M_data + M_pilot;
    dft_normalizer = 1.0f / sqrtf(M_occupied);


    // Schmidl & Cox
    delay_M2.resize(num_streams);
    crosscorrelator.resize(num_streams);
    normalizer.resize(num_streams);
    plateau_start.resize(num_streams);
    plateau_end.resize(num_streams);
    float taps_crosscorrelator[M2];
    float taps_normalizer[M];
    for(unsigned int i = 0; i < M2; i++)
      taps_crosscorrelator[i] = -1.0f; // FIXME why not 1.0f
    for(unsigned int i = 0; i < M; i++)
      taps_normalizer[i] = 0.5f;
    in_plateau.resize(num_streams);
    access_code_buffer.resize(num_streams);

    // get volk alignment for volk_malloc
    size_t volk_alignment = volk_get_alignment();
    for(unsigned int i = 0; i < num_streams; i++)
    {
      // memory allocation
      // NOTE: all to be freed in destructor
      S1[i] = (gr_complex *) malloc 
              (sizeof(gr_complex)*
               M*num_access_codes);
      s1[i] = (gr_complex *) malloc 
              (sizeof(gr_complex)*
               M*num_access_codes);
      symbol[i] = (gr_complex *) malloc 
                  (sizeof(gr_complex)*symbol_len);
      X[i] = (std::complex<float> *) fftwf_malloc
             (sizeof(std::complex<float>)*M);
      x[i] = (std::complex<float> *) fftwf_malloc
             (sizeof(std::complex<float>)*M);
      assert(volk_is_aligned(X[i]));
      assert(volk_is_aligned(x[i]));
      fft[i] = fftwf_plan_dft_1d(M,
                                 reinterpret_cast<fftwf_complex *>(x[i]),
                                 reinterpret_cast<fftwf_complex *>(X[i]),
                                 FFTW_FORWARD,
                                 FFTW_ESTIMATE);
      // initialize S1 sequence
      ofdmframe_init_S1(p,
                        M,
                        num_access_codes,
                        S1[i],
                        s1[i],
                        _ms_S1[i]);
      // Schmidl & Cox
      delay_M2[i] = wdelaycf_create(M2);
      crosscorrelator[i] = firfilt_crcf_create(taps_crosscorrelator, M2);
      normalizer[i] = firfilt_rrrf_create(taps_normalizer, M);
      in_plateau[i] = false;
      plateau_start[i] = 0;
      plateau_end[i] = 0;
      access_code_buffer[i] = windowcf_create(access_code_buffer_len +
                                              tx_sig_len);

#if DEBUG_LOG
      f_sc_debug_out.push_back(fopen(boost::str(boost::format("%s%d.dat")
                                                % F_SC_DEBUG_OUT_PREFIX
                                                % (i + 1)).c_str(),
                                "wb"));
      if(!(f_sc_debug_out[i])) {
        printf("***** error opening file %s\n",
               boost::str(boost::format("%s%d.dat") 
                          % F_SC_DEBUG_OUT_PREFIX
                          % (i + 1)).c_str());
        exit(1);
      }
#endif
    }
    normalize_gain = (float *) volk_malloc 
                     (sizeof(float)*M_occupied,
                      volk_alignment);
    std::fill(normalize_gain,
              normalize_gain + M_occupied,
              1.0f);
    // initialize S0 sequence
    ofdmframe_init_S0(p,
                      M,
                      S0,
                      s0,
                      _ms_S0);
    volk_buff_fc1 = (std::complex<float> *) volk_malloc 
                    (sizeof(std::complex<float>)*M,
                     volk_alignment);
    volk_buff_fc2 = (std::complex<float> *) volk_malloc 
                    (sizeof(std::complex<float>)*M,
                     volk_alignment);
    volk_buff_fc3 = (std::complex<float> *) volk_malloc 
                    (sizeof(std::complex<float>)*M,
                     volk_alignment);
    volk_buff_f1 = (float *) volk_malloc 
                   (sizeof(float)*
                    (M + cp_len)*3*num_streams,
                    volk_alignment);
  
    // setting the counters
    sync_index = 0;
    num_samples_processed = 0;
    state = STATE_SEEK_PLATEAU;
    siso_tx = 0;
    siso_rx = 0;
  }
  
  void framesync::print()
  {
      printf("ofdmframegen:\n");
      printf("    num subcarriers     :   %-u\n", M);
      printf("      - NULL            :   %-u\n", M_null);
      printf("      - pilot           :   %-u\n", M_pilot);
      printf("      - data            :   %-u\n", M_data);
      printf("    cyclic prefix len   :   %-u\n", cp_len);
      printf("    ");
      ofdmframe_print_sctype(p, M);
  }
  
  unsigned long int framesync::get_sync_index()
  {
    return sync_index;
  }
  
  std::vector<std::vector<std::vector<gr_complex> > >
  framesync::get_G()
  {
    return G;
  }
  
  void framesync::reset()
  {
    state = STATE_SEEK_PLATEAU;
  }
  
  unsigned long long int framesync::get_num_samples_processed()
  {
    return num_samples_processed;
  }
  
  framesync_states_t
  framesync::execute(std::vector<gr_complex *>
                     const &in_buff,
                     unsigned int num_samples)
  {
    gr_complex _x[num_streams];
    bool break_loop = false;
  
    // TODO process a vector of inputs, rather than
    // one by one in a for loop.
    for(unsigned int i = 0; i < num_samples; i++){
      for(unsigned int stream = 0; stream < num_streams; stream++) {
        _x[stream] = in_buff[stream][i];
      }
  
      // FIXME correct frequency offset
      switch(state) {
        case STATE_SEEK_PLATEAU:
          execute_sc_sync(_x);
          break;
        case STATE_SAVE_ACCESS_CODES:
          execute_save_access_codes(_x);
          break;
        case STATE_MIMO:
          break_loop = true;
          break;
        default:
          printf("framesync: state %d not handled, exiting\n", state);
          exit(1);
      }
      num_samples_processed++;
      if(break_loop)
        break;
    }
    return state;
  }
    
  void framesync::execute_siso_decode(gr_complex _x[])
  {
    symbol[siso_rx][sample_counter] = _x[siso_rx];
    sample_counter++;
    if(sample_counter < symbol_len)
      return;
    memmove(x[siso_rx], symbol[siso_rx] + cp_len, sizeof(gr_complex)*M);
    fftwf_execute(fft[siso_rx]);
    volk_32fc_s32fc_multiply_32fc(volk_buff_fc1,
                                  X[siso_rx],
                                  dft_normalizer,
                                  M);
    memmove(X[siso_rx], volk_buff_fc1, sizeof(gr_complex)*M);
    for(unsigned int stream = 0; stream < num_streams; stream++) {
      std::fill(x[stream], x[stream] + M, gr_complex(0.0, 0.0));
    }
    unsigned int j = 0;
    for(unsigned int sc = 0; sc < M; sc++) {
      if(p[sc] != OFDMFRAME_SCTYPE_NULL) {
        x[siso_rx][j] = X[siso_rx][sc]/(G[sc][siso_rx][siso_tx]);
        j++;
      }
    }
    callback(x, M_occupied);
    sample_counter = 0;
  }

  void framesync::execute_mimo_decode(gr_complex _x[])
  {
    for(unsigned int rx_stream = 0; rx_stream < num_streams; rx_stream++) {
      symbol[rx_stream][sample_counter] = _x[rx_stream];
    }
    sample_counter++;
    if(sample_counter < symbol_len)
      return;

#if DEBUG_PRINT_VERBOSE
    printf("**** execute_mimo_decode: processing OFDM frame\n");
    for(unsigned int i = 0; i < symbol_len; i++) {
      printf("%3u,\t", i);
      for(unsigned int stream = 0; stream < num_streams; stream++) {
        printf("%3.6f + %3.6fi,\t\t",
               real(symbol[stream][i]),
               imag(symbol[stream][i]));
      }
      printf("\n");
    }
#endif

    for(unsigned int rx_stream = 0; rx_stream < num_streams; rx_stream++) {
      memmove(x[rx_stream], symbol[rx_stream] + cp_len,
              sizeof(gr_complex)*M);
      fftwf_execute(fft[rx_stream]);
      volk_32fc_s32fc_multiply_32fc(volk_buff_fc1,
                                    X[rx_stream],
                                    dft_normalizer,
                                    M);
      memmove(X[rx_stream], volk_buff_fc1, sizeof(gr_complex)*M);                                    
    }
    // FIXME works only for 2 x 2.
    // Extend for higher configurations
    unsigned int j = 0;
    for(unsigned int sc = 0; sc < M; sc++) {
      if(p[sc] != OFDMFRAME_SCTYPE_NULL) {
        x[0][j] = W[sc][0][0]*X[0][sc] +
                  W[sc][0][1]*X[1][sc];
        x[1][j] = W[sc][1][0]*X[0][sc] +
                  W[sc][1][1]*X[1][sc];
        j++;     
      }
    }
    volk_32fc_32f_multiply_32fc(X[0],
                                x[0],
                                normalize_gain,
                                M_occupied);
    volk_32fc_32f_multiply_32fc(X[1],
                                x[1],
                                normalize_gain,
                                M_occupied);
    callback(X, M_occupied);
    sample_counter = 0;
  }
  
  void framesync::execute_sc_sync(gr_complex _x[])
  {
    float y[num_streams];
    bool proceed = true;
    for(unsigned int stream = 0; stream < num_streams; stream++) {
      windowcf_push(access_code_buffer[stream], _x[stream]);
      y[stream] = execute_sc_sync(_x[stream], stream);
#if DEBUG_LOG
      fwrite(&y[stream], sizeof(float), 1, f_sc_debug_out[stream]);
#endif
      if(y[stream] > PLATEAU_THREASHOLD)
      {
        if(in_plateau[stream])
          plateau_end[stream] = num_samples_processed;
        else {
          in_plateau[stream] = true;
          plateau_start[stream] = num_samples_processed;
          plateau_end[stream] = num_samples_processed;
        }
      }
      else
        in_plateau[stream] = false;
      proceed = proceed &&
                (plateau_end[stream] - plateau_start[stream] > cp_len) &&
                (in_plateau[stream]);
    }
    if(proceed) {
      for(unsigned int stream = 0; stream < num_streams; stream++)
        sync_index += plateau_start[stream];
      sync_index /= num_streams;
      printf("***** proceed to save access codes ***** \n");
      state = STATE_SAVE_ACCESS_CODES;
    }
  }

  float framesync::execute_sc_sync(gr_complex _x, unsigned int stream)
  {
    gr_complex x;
    wdelaycf_read(delay_M2[stream], &x);
    wdelaycf_push(delay_M2[stream], _x);
    firfilt_crcf_push(crosscorrelator[stream], conj(x)*_x);
    firfilt_crcf_execute(crosscorrelator[stream], &x);
    float z = real(_x)*real(_x) + imag(_x)*imag(_x);
    firfilt_rrrf_push(normalizer[stream], z);
    firfilt_rrrf_execute(normalizer[stream], &z);
    return ((real(x)*real(x) + imag(x)*imag(x))/(z*z));
  }
  
  void framesync::execute_save_access_codes(gr_complex _x[])
  {
    if(num_samples_processed - sync_index < 
       tx_sig_len + access_code_buffer_len - symbol_len) {
      for(unsigned int chan = 0; chan < num_streams; chan++) {
        windowcf_push(access_code_buffer[chan], _x[chan]);
      }
      return;
    }
    printf("**** access codes saved *****\n");
    estimate_channel();
    state = STATE_MIMO;
  }
  
  void framesync::estimate_channel() {
    unsigned int _ac_id;
    unsigned int max_ac_id = num_streams*num_access_codes;
    gr_complex xyz;
    gr_complex * buffer_ptr[num_streams];
    float max_corr [num_streams][max_ac_id];
    unsigned int corr_indices [num_streams][max_ac_id];
    float * corr_values [num_streams][max_ac_id];
    float * s0_corr[num_streams];
    float max_s0_corr[num_streams];
    unsigned int s0_corr_index[num_streams];
    FILE * corr_files [num_streams][max_ac_id];
    FILE * corr_file_s0[num_streams];
    for(unsigned int chan = 0; chan < num_streams; chan++) {
      for(unsigned int ac_id = 0; ac_id < max_ac_id; ac_id++) {
        max_corr[chan][ac_id] = 0.0f;
        corr_indices[chan][ac_id] = 0;
        corr_values[chan][ac_id] = (float *) malloc 
                                   (sizeof(float)*(access_code_buffer_len - M));
	std::fill(corr_values[chan][ac_id],
		  corr_values[chan][ac_id] + access_code_buffer_len - M,
		  0.0);
#if DEBUG_LOG
        corr_files[chan][ac_id] = fopen(boost::str(boost::format("%s%d_%d.dat")
                                                   % CORR_FILE_PREFIX
                                                   % (chan + 1)
                                                   % (ac_id + 1)).c_str(),
                                        "wb");
        assert(corr_files[chan][ac_id]);
#endif
      }
      s0_corr[chan] = (float *) malloc (sizeof(float)*(access_code_buffer_len - M));
      max_s0_corr[chan] = 0.0f;
      s0_corr_index[chan] = 0;
      std::fill(s0_corr[chan],
		s0_corr[chan] + access_code_buffer_len - M,
		0.0);
#if DEBUG_LOG
      corr_file_s0[chan] = fopen(boost::str(boost::format("%s%d_0.dat")
                                            % CORR_FILE_PREFIX
                                            % (chan + 1)).c_str(),
                                 "wb");
      assert(corr_file_s0[chan]);
#endif
    }

    for(unsigned int chan = 0; chan < num_streams; chan++)
      windowcf_read(access_code_buffer[chan], &buffer_ptr[chan]);

#if USE_NEW_CHANNEL_EST
    unsigned int sample;
    for(unsigned int i = 0; i < symbol_len; i++) {
      for(unsigned int rx_stream = 0;
	  rx_stream < num_streams; rx_stream++) {
	sample = i;
	memmove(x[rx_stream], buffer_ptr[rx_stream] + sample,
		sizeof(gr_complex)*M);
	fftwf_execute(fft[rx_stream]);
        // correlate with S0
        volk_32fc_x2_conjugate_dot_prod_32fc(&xyz,
                                             X[rx_stream],
                                             S0,
                                             M);
        s0_corr[rx_stream][sample] = (real(xyz)*real(xyz) + 
				      imag(xyz)*imag(xyz))/float(M*M);
        if(s0_corr[rx_stream][sample] > max_s0_corr[rx_stream]) {
          max_s0_corr[rx_stream] = s0_corr[rx_stream][sample];
          s0_corr_index[rx_stream] = sample;
        }
	for(unsigned int code = 0; code < num_access_codes; code++) {
	  for(unsigned int tx_chan = 0; tx_chan < num_streams; tx_chan++) {
	    _ac_id = code*num_streams + tx_chan;
	    sample = i + symbol_len*(_ac_id + 1);
	    memmove(x[rx_stream], buffer_ptr[rx_stream] + sample,
		    sizeof(gr_complex)*M);
	    fftwf_execute(fft[rx_stream]);
	    volk_32fc_x2_conjugate_dot_prod_32fc(&xyz,
						 X[rx_stream],
						 (S1[tx_chan] + code*M),
						 M);
	    corr_values[rx_stream][_ac_id][sample] =
	      (real(xyz)*real(xyz) + imag(xyz)*imag(xyz))/float(M*M);
	    if(corr_values[rx_stream][_ac_id][sample] > 
	       max_corr[rx_stream][_ac_id]) {
	      corr_indices[rx_stream][_ac_id] = sample;
	      max_corr[rx_stream][_ac_id] = 
		corr_values[rx_stream][_ac_id][sample];
	    }
	  }
	}
      }
    }
#else
    // compute the access code indices
    for(unsigned int sample = 0; sample < access_code_buffer_len - M; sample++) {
      for(unsigned int chan = 0; chan < num_streams; chan++) {
        memmove(x[chan], buffer_ptr[chan] + sample, sizeof(gr_complex)*M);
        fftwf_execute(fft[chan]);
        // normalize
        volk_32fc_s32fc_multiply_32fc(volk_buff_fc1,
                                      X[chan],
                                      dft_normalizer,
                                      M);
        memmove(X[chan], volk_buff_fc1, sizeof(gr_complex)*M);
        // correlate with S0
        volk_32fc_x2_conjugate_dot_prod_32fc(&xyz,
                                             X[chan],
                                             S0,
                                             M);
        s0_corr[chan][sample] = real(xyz)*real(xyz) + imag(xyz)*imag(xyz);
        if(s0_corr[chan][sample] > max_s0_corr[chan]) {
          max_s0_corr[chan] = s0_corr[chan][sample];
          s0_corr_index[chan] = sample;
        }
        // correlate with access codes
        for(unsigned int code = 0; code < num_access_codes; code++) {
          for(unsigned int tx_chan = 0; tx_chan < num_streams; tx_chan++) {
            _ac_id = code*num_streams + tx_chan;
            volk_32fc_x2_conjugate_dot_prod_32fc(&xyz,
                                                 X[chan],
                                                 (S1[tx_chan] + code*M),
                                                 M);
            corr_values[chan][_ac_id][sample] = real(xyz)*real(xyz) + imag(xyz)*imag(xyz);
            if(corr_values[chan][_ac_id][sample] > max_corr[chan][_ac_id]) {
              corr_indices[chan][_ac_id] = sample;
              max_corr[chan][_ac_id] = corr_values[chan][_ac_id][sample];
            }
          }
        }
      }
    }
#endif
#if DEBUG_PRINT
    for(unsigned int chan = 0; chan < num_streams; chan++) {
      printf("******** s0_corr_index[%d]    : %6u\n", chan, s0_corr_index[chan]);
      for(unsigned int ac_id = 0; ac_id < max_ac_id; ac_id++) {
        //      printf("******** s0_corr_index[%d]    : %6u\n", chan, s0_corr_index[chan]);
      printf("********* corr_indices[%d][%d] : %6u\n", chan,
             ac_id,
             corr_indices[chan][ac_id]);
      }
    }
#endif

    // FIXME Not sure if this the correct way of doing things.
    // The indices of maximum correlation are "symbol_len" samples
    // from each other. It should be true theoretically. Not sure
    // if this is always true in practice.
    for(unsigned int code = 0; code < num_access_codes; code++) {
      for(unsigned int rx_stream = 0; rx_stream < num_streams; rx_stream++) {
        for(unsigned int tx_stream = 0; tx_stream < num_streams; tx_stream++) {
          _ac_id = code*num_streams + tx_stream;
          memmove(x[rx_stream],
                  buffer_ptr[rx_stream] + corr_indices[rx_stream][_ac_id],
                  sizeof(gr_complex)*M);
          fftwf_execute(fft[rx_stream]);
          for(unsigned int i = 0; i < M; i++) {
            if(p[i] != OFDMFRAME_SCTYPE_NULL)
              G[i][rx_stream][tx_stream] += X[rx_stream][i]/S1[tx_stream][code*M + i];
          }
        }
      }
    }
    
    for(unsigned int rx_stream = 0; rx_stream < num_streams; rx_stream++) {
      for(unsigned int tx_stream = 0; tx_stream < num_streams; tx_stream++) {
        for(unsigned int sc = 0; sc < M; sc++) {
          if(p[sc] != OFDMFRAME_SCTYPE_NULL)
            G[sc][rx_stream][tx_stream] *= dft_normalizer/float(num_access_codes);
        }
      }
    }
    
#if INVERT_CHANNEL
    unsigned int j = 0;
    for(unsigned int i = 0; i < M; i++) {
      if(p[i] != OFDMFRAME_SCTYPE_NULL)
        normalize_gain[j++] = invert(W[i], G[i]);
    }
#endif

#if DEBUG_PRINT
    for(unsigned int rx_stream = 0; rx_stream < num_streams; rx_stream++) {
      for(unsigned int tx_stream = 0; tx_stream < num_streams; tx_stream++) {
        printf("channel states between RX-%u and TX-%u\n",
               rx_stream, tx_stream);
        for(unsigned int i = 0; i < M; i++) {
          printf("%2.4f + %2.4fi \t\t %4.6f + %4.6f\n",
                 real(G[i][rx_stream][tx_stream]),
                 imag(G[i][rx_stream][tx_stream]),
                 real(W[i][rx_stream][tx_stream]),
                 imag(W[i][rx_stream][tx_stream]));
        }
      }
    }
    printf("**** normalize_gain\n");
    for(unsigned int j = 0; j < M_occupied; j++)
      printf("%3.6f\n", normalize_gain[j]);
#endif

    // process the mimo samples in the
    // access code buffer
    gr_complex _x[num_streams];
    sample_counter = 0;
    for(unsigned int i = corr_indices[1][max_ac_id - 1] + M;
        i < access_code_buffer_len + tx_sig_len;
        i++) {
      for(unsigned int stream = 0; stream < num_streams; stream++) {
        _x[stream] = buffer_ptr[stream][i];
      }
#if SISO
      execute_siso_decode(_x);
#else
      execute_mimo_decode(_x);
#endif
    }

    // free memory
    for(unsigned int chan = 0; chan < num_streams; chan++) {
      for(unsigned int ac_id = 0; ac_id < max_ac_id; ac_id++) {
#if DEBUG_LOG
        fwrite(corr_values[chan][ac_id], sizeof(float), access_code_buffer_len - M,
               corr_files[chan][ac_id]);
        fclose(corr_files[chan][ac_id]);
#endif
        free(corr_values[chan][ac_id]);
      }
#if DEBUG_LOG
      fwrite(s0_corr[chan], sizeof(float), access_code_buffer_len - M, corr_file_s0[chan]);
      fclose(corr_file_s0[chan]);
#endif
      free(s0_corr[chan]);
    }
  }

  void
  framesync::set_siso_tx(unsigned int _tx) {
    siso_tx = _tx;
  }

  void
  framesync::set_siso_rx(unsigned int _rx) {
    siso_rx = _rx;
  }

  void
  framesync::compute_receive_beamformer() {
  }
  
  unsigned long int
  framesync::get_plateau_start(unsigned int stream)
  {
    assert(stream < num_streams);
    return plateau_start[stream];
  }
  
  unsigned long int
  framesync::get_plateau_end(unsigned int stream)
  {
    assert(stream < num_streams);
    return plateau_end[stream];
  }
  
  framesync::~framesync()
  {
   for(unsigned int i = 0; i < num_streams; i++)
    {
      free(S1[i]);
      free(s1[i]);
      fftwf_free(X[i]);
      fftwf_free(x[i]);
      fftwf_destroy_plan(fft[i]);
      free(symbol[i]);
      // Schmidl & Cox
      wdelaycf_destroy(delay_M2[i]);
      firfilt_crcf_destroy(crosscorrelator[i]);
      firfilt_rrrf_destroy(normalizer[i]);
      windowcf_destroy(access_code_buffer[i]);
#if DEBUG_LOG
      fclose(f_sc_debug_out[i]);
#endif
    }
    free(normalize_gain);
    free(p);
    free(S0);
    free(s0);
    volk_free(volk_buff_fc1);
    volk_free(volk_buff_fc2);
    volk_free(volk_buff_fc3);
    volk_free(volk_buff_f1);
  }
} // namespace rx_beamforming

namespace tx_beamforming {
}

#if USE_ALL_CARRIERS
void ofdmframe_init_default_sctype(unsigned char * _p, unsigned int _M)
{
    for (unsigned int i=0; i<_M; i++)
        _p[i] = OFDMFRAME_SCTYPE_DATA;
}
#else
void ofdmframe_init_default_sctype(unsigned char * _p, unsigned int _M)
{
    // validate input
    if (_M < 6) {
        fprintf(stderr,"warning: ofdmframe_init_default_sctype(), less than 4 subcarriers\n");
    }

    unsigned int i;
    unsigned int M2 = _M/2;

    // compute guard band
    unsigned int G = 0;
    if (ADD_NULL_CARRIERS) {
      G = _M / 10;
      if (G < 2) G = 2;
    }

    // designate pilot spacing
    unsigned int P = (_M > 34) ? 8 : 4;
    unsigned int P2 = P/2;

    // initialize as NULL
    for (i=0; i<_M; i++)
        _p[i] = OFDMFRAME_SCTYPE_NULL;

    // upper band
    for (i=1; i<M2-G; i++) {
        if ( ((i+P2)%P) == 0 )
            _p[i] = OFDMFRAME_SCTYPE_PILOT;
        else
            _p[i] = OFDMFRAME_SCTYPE_DATA;
    }

    // lower band
    for (i=1; i<M2-G; i++) {
        unsigned int k = _M - i;
        if ( ((i+P2)%P) == 0 )
            _p[k] = OFDMFRAME_SCTYPE_PILOT;
        else
            _p[k] = OFDMFRAME_SCTYPE_DATA;
    }
}
#endif

void ofdmframe_validate_sctype(const unsigned char * _p,
                               unsigned int _M,
                               unsigned int * _M_null,
                               unsigned int * _M_pilot,
                               unsigned int * _M_data)
{
    // clear counters
    unsigned int M_null  = 0;
    unsigned int M_pilot = 0;
    unsigned int M_data  = 0;

    unsigned int i;
    for (i=0; i<_M; i++) {
        // update appropriate counters
        if (_p[i] == OFDMFRAME_SCTYPE_NULL)
            M_null++;
        else if (_p[i] == OFDMFRAME_SCTYPE_PILOT)
            M_pilot++;
        else if (_p[i] == OFDMFRAME_SCTYPE_DATA)
            M_data++;
        else {
            fprintf(stderr,"error: ofdmframe_validate_sctype(), invalid subcarrier type (%u)\n", _p[i]);
            exit(1);
        }
    }

    // set outputs
    *_M_null  = M_null;
    *_M_pilot = M_pilot;
    *_M_data  = M_data;
}

void ofdmframe_print_sctype(const unsigned char * _p, unsigned int _M)
{
    unsigned int i;

    printf("[");
    for (i=0; i<_M; i++) {
        unsigned int k = (i + _M/2) % _M;

        switch (_p[k]) {
        case OFDMFRAME_SCTYPE_NULL:     printf(".");    break;
        case OFDMFRAME_SCTYPE_PILOT:    printf("|");    break;
        case OFDMFRAME_SCTYPE_DATA:     printf("+");    break;
        default:
            fprintf(stderr,"error: ofdmframe_print_default_sctype(), invalid subcarrier type\n");
            exit(1);
        }
    }

    printf("]\n");
}

#if USE_NEW_INIT_S0
void ofdmframe_init_S0(const unsigned char * _p,
                       unsigned int _M,
                       std::complex<float> * _S0,
                       std::complex<float> * _s0,
                       msequence ms)
{
    unsigned int i;
    unsigned int s;
    unsigned int M_S0 = 0;
    fftwf_plan ifft;
    gr_complex * X = (gr_complex *) fftwf_malloc (sizeof(gr_complex)*_M);
    gr_complex * x = (gr_complex *) fftwf_malloc (sizeof(gr_complex)*_M);
    ifft = fftwf_plan_dft_1d(_M,
                             reinterpret_cast<fftwf_complex *>(X),
                             reinterpret_cast<fftwf_complex *>(x),
                             FFTW_BACKWARD,
                             FFTW_ESTIMATE);

    // short sequence
    for (i=0; i<_M; i++) {
        // generate symbol
        s = msequence_generate_symbol(ms, 1) & 0x01;

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            _S0[i] = 0.0f;
        } else {
            if ( (i%2) == 0 ) {
                // even subcarrer
                _S0[i] = s ? 1.0f : -1.0f;
                M_S0++;
            } else {
                // odd subcarrer (ignore)
                _S0[i] = 0.0f;
            }
        }
    }

    // ensure at least one subcarrier was enabled
    if (M_S0 == 0) {
        fprintf(stderr,"error: ofdmframe_init_S0(), no subcarriers enabled; check allocation\n");
        exit(1);
    }

    // run inverse fft to get time-domain sequence
    // TODO make in independent of liquid
    gr_complex dft_normalizer = sqrt(1.0/float(M_S0));

    memmove(X, _S0, sizeof(gr_complex)*_M);
    fftwf_execute(ifft);
    volk_32fc_s32fc_multiply_32fc(_s0,
                                  x,
                                  dft_normalizer,
                                  _M);
    fftwf_destroy_plan(ifft);
    fftwf_free(X);
    fftwf_free(x);
}
#else
void ofdmframe_init_S0(const unsigned char * _p,
                       unsigned int _M,
                       std::complex<float> * _S0,
                       std::complex<float> * _s0,
                       msequence ms)
{
    unsigned int i;
    unsigned int s;
    unsigned int M_S0 = 0;
    // short sequence
    for (i=0; i<_M; i++) {
        // generate symbol
        s = msequence_generate_symbol(ms, 3) & 0x01;

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            _S0[i] = 0.0f;
        } else {
            if ( (i%2) == 0 ) {
                // even subcarrer
                _S0[i] = s ? 1.0f : -1.0f;
                M_S0++;
            } else {
                // odd subcarrer (ignore)
                _S0[i] = 0.0f;
            }
        }
    }

    // ensure at least one subcarrier was enabled
    if (M_S0 == 0) {
        fprintf(stderr,"error: ofdmframe_init_S0(), no subcarriers enabled; check allocation\n");
        exit(1);
    }

    // run inverse fft to get time-domain sequence
    // TODO make in independent of liquid
    fft_run(_M, _S0, _s0, LIQUID_FFT_BACKWARD, 0);

    // normalize time-domain sequence level
    // TODO Do this with volk
    float g = 1.0f / sqrtf(M_S0);
    for (i=0; i<_M; i++)
      _s0[i] *= g;
}
#endif

#if USE_NEW_INIT_S1
#if MAKE_S1_QPSK
void ofdmframe_init_S1(const unsigned char * _p,
                       unsigned int _M,
                       unsigned int _num_access_codes,
                       std::complex<float> * _S1,
                       std::complex<float> * _s1,
                       msequence ms)
{
    unsigned int i, j;

    unsigned int s;
    unsigned int M_S1 = 0;
    fftwf_plan ifft;
    gr_complex * X = (gr_complex *) fftwf_malloc (sizeof(gr_complex)*_M);
    gr_complex * x = (gr_complex *) fftwf_malloc (sizeof(gr_complex)*_M);
    gr_complex dft_normalizer = sqrt(1.0/float(_M));
    ifft = fftwf_plan_dft_1d(_M,
                             reinterpret_cast<fftwf_complex *>(X),
                             reinterpret_cast<fftwf_complex *>(x),
                             FFTW_BACKWARD,
                             FFTW_ESTIMATE);

    // long sequence
    for(j = 0; j < _num_access_codes; j++) {
      M_S1 = 0;
      for (i=0; i<_M; i++) {
        // generate symbol
        // s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,2) & 0x11;
        assert(s < 4);

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            (_S1 + _M*j)[i] = 0.0f;
        } else {
          (_S1 + _M*j)[i] = QPSK_CONSTELLATION[s];
          M_S1++;
        }
      }
      // run inverse fft to get time-domain sequence
      memmove(X, _S1 + _M*j, sizeof(gr_complex)*_M);
      fftwf_execute(ifft);
      dft_normalizer = sqrt(1.0/float(M_S1));
      volk_32fc_s32fc_multiply_32fc(_s1 + _M*j,
                                    x,
                                    dft_normalizer,
                                    _M);
    }
    fftwf_destroy_plan(ifft);
    fftwf_free(X);
    fftwf_free(x);
}
#else
void ofdmframe_init_S1(const unsigned char * _p,
                       unsigned int _M,
                       unsigned int _num_access_codes,
                       std::complex<float> * _S1,
                       std::complex<float> * _s1,
                       msequence ms)
{
    unsigned int i, j;

    unsigned int s;
    unsigned int M_S1 = 0;
    fftwf_plan ifft;
    gr_complex * X = (gr_complex *) fftwf_malloc (sizeof(gr_complex)*_M);
    gr_complex * x = (gr_complex *) fftwf_malloc (sizeof(gr_complex)*_M);
    gr_complex dft_normalizer = sqrt(1.0/float(_M));
    ifft = fftwf_plan_dft_1d(_M,
                             reinterpret_cast<fftwf_complex *>(X),
                             reinterpret_cast<fftwf_complex *>(x),
                             FFTW_BACKWARD,
                             FFTW_ESTIMATE);

    // long sequence
    for(j = 0; j < _num_access_codes; j++) {
      for (i=0; i<_M; i++) {
        // generate symbol
        // s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,1) & 0x01;
        assert(s < 2);

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            (_S1 + _M*j)[i] = 0.0f;
        } else {
          (_S1 + _M*j)[i] = BPSK_CONSTELLATION[s];
          M_S1++;
        }
      }
      // run inverse fft to get time-domain sequence
      memmove(X, _S1 + _M*j, sizeof(gr_complex)*_M);
      fftwf_execute(ifft);
      volk_32fc_s32fc_multiply_32fc(_s1 + _M*j,
                                    x,
                                    dft_normalizer,
                                    _M);
    }
    fftwf_destroy_plan(ifft);
    fftwf_free(X);
    fftwf_free(x);
}
#endif
#else
#if MAKE_S1_QPSK
void ofdmframe_init_S1(const unsigned char * _p,
                       unsigned int _M,
                       unsigned int _num_access_codes,
                       std::complex<float> * _S1,
                       std::complex<float> * _s1,
                       msequence ms)
{
    unsigned int i, j;

    unsigned int s;
    unsigned int M_S1 = 0;

    // long sequence
    for(j = 0; j < _num_access_codes; j++) {
      for (i=0; i<_M; i++) {
        // generate symbol
        // s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,2) & 0x11;
        assert(s < 4);

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            (_S1 + _M*j)[i] = 0.0f;
        } else {
          (_S1 + _M*j)[i] = QPSK_CONSTELLATION[s];
          M_S1++;
        }
      }
      // run inverse fft to get time-domain sequence
      fft_run(_M, _S1 + _M*j, _s1 + _M*j, LIQUID_FFT_BACKWARD, 0);
    }

    // normalize time-domain sequence level
    float g = 1.0f / sqrtf(M_S1/_num_access_codes);
    for (i=0; i<_M*_num_access_codes; i++)
        _s1[i] *= g;
}
#else
void ofdmframe_init_S1(const unsigned char * _p,
                       unsigned int _M,
                       unsigned int _num_access_codes,
                       std::complex<float> * _S1,
                       std::complex<float> * _s1,
                       msequence ms)
{
    unsigned int i, j;

    unsigned int s;
    unsigned int M_S1 = 0;

    // long sequence
    for(j = 0; j < _num_access_codes; j++) {
      for (i=0; i<_M; i++) {
        // generate symbol
        // s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,3) & 0x01;
        assert(s < 2);

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            (_S1 + _M*j)[i] = 0.0f;
        } else {
          (_S1 + _M*j)[i] = BPSK_CONSTELLATION[s];
          M_S1++;
        }
      }
      // run inverse fft to get time-domain sequence
      fft_run(_M, _S1 + _M*j, _s1 + _M*j, LIQUID_FFT_BACKWARD, 0);
    }

    // normalize time-domain sequence level
    float g = 1.0f / sqrtf(M_S1/_num_access_codes);
    for (i=0; i<_M*_num_access_codes; i++)
        _s1[i] *= g;
}
#endif
#endif

float invert(std::vector<std::vector<gr_complex> >  &W,
             std::vector<std::vector<gr_complex> > const &G) {
  assert(G.size() == 2);
  assert(W.size() == 2);
  for(unsigned int rx_stream = 0; rx_stream < 2; rx_stream++) {
    assert(G[rx_stream].size() == 2);
    assert(W[rx_stream].size() == 2);
  }
  gr_complex det = G[0][0]*G[1][1] - G[0][1]*G[1][0];
#if INVERT_TO_UNITY
  gr_complex det_inv = 1.0f/det;
#else
  gr_complex det_inv = conj(det);
#endif
  W[0][0] = det_inv*G[1][1];
  W[1][1] = det_inv*G[0][0];
  W[1][0] = -det_inv*G[1][0];
  W[0][1] = -det_inv*G[0][1];
#if INVERT_TO_UNITY
  return 1.0;
#else
  return  1.0f/(real(det)*real(det) + imag(det)*imag(det));
#endif
}
