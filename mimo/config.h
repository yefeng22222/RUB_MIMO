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

#ifndef CONFIG_H
#define CONFIG_H

typedef enum
{
  CLOCK_SOURCE_NONE = 0,
  CLOCK_SOURCE_EXTERNAL
} clock_source_type;
typedef enum
{
  TIME_SOURCE_NONE = 0,
  TIME_SOURCE_EXTERNAL
} time_source_type;

extern gr_complex BPSK_CONSTELLATION[];
extern gr_complex QPSK_CONSTELLATION[];

// define default configurations
// USRP identities
#define N200_12                 "addr=134.147.118.212"
#define N200_15                 "addr=134.147.118.215"
#define X300A                   "addr=134.147.118.216"
#define X300B                   "addr=134.147.118.217"
#define B210_TX                 "serial=308F955"
#define B210_RX                 "serial=308F965"

// subdevice specifications
#define X300_SUBDEV_SPEC        "A:0 B:0"
#define N200_SUBDEV_SPEC        "A:0"
#define B210_SUBDEV_SPEC_TX     "A:B A:A"
#define B210_SUBDEV_SPEC_RX     "A:A A:B"

// tx/rx streamer configurations
#define CPU                     "fc32"
#define WIRE                    "sc16"

// RF front end configurations
#define CENTER_FREQUENCY        2450e6
#define SAMPLING_RATE           1.0e6
#define TX_FRONTEND_GAIN        67.0
#define RX_FRONTEND_GAIN        45.0
#define BASEBAND_GAIN           0.25

#define CLOCK_SOURCE            CLOCK_SOURCE_EXTERNAL
#define TIME_SOURCE             TIME_SOURCE_NONE

// OFDM configurations
#define NUM_SUBCARRIERS         2048
#define CP_LENGTH               152

// generator polynomials obtained from
// primitive_polys.pdf
#define LFSR_SMALL_LENGTH       12
#define LFSR_LARGE_LENGTH       13
#define LFSR_SMALL_0_GEN_POLY   010123
#define LFSR_SMALL_1_GEN_POLY   010151
#define LFSR_LARGE_0_GEN_POLY   020033
#define LFSR_LARGE_1_GEN_POLY   020047

// misc configurations
#define VERBOSITY               true
#define LOG                     true
#define LOG_DIR                 "/tmp/"
#define DATA_DIR                "data/"
#define TX_BEAMFORMING          0

#define DEBUG_LOG               true
#define F_SC_DEBUG_OUT_PREFIX   "/tmp/f_sc_"
#define CORR_FILE_PREFIX        "/tmp/corr_"
#define PLATEAU_THREASHOLD      0.95

#define SISO                    false
#define SISO_TX                 1
#define SISO_RX                 1
#define PID_MAX                 1000
#define DEBUG_PRINT             true
#define DEBUG_PRINT_VERBOSE     false
#define USE_ALL_CARRIERS        true
#define ADD_NULL_CARRIERS       true
#define BPSK_CONSTELLATION_SIZE 2
#define QPSK_CONSTELLATION_SIZE 4
#define CONSTELLATION           BPSK_CONSTELLATION
#define CONSTELLATION_SIZE      BPSK_CONSTELLATION_SIZE
#define MAKE_S1_QPSK            false
#define INVERT_CHANNEL          true
#define INVERT_TO_UNITY         false
#define NUM_ACCESS_CODES        20
#define NUM_MIMO_ACCESS_CODES   3
#define NUM_STREAMS             2
#define MODEM_SCHEME            LIQUID_MODEM_ARB32OPT
#define ARITY                   32
#define SAME_SIGNAL_ON_ALL_TX   false
#define NORMALIZE_DFT_TXRX      true
#define NORMALIZE_DFT_RX        true
#define USE_NEW_INIT_S0         true
#define USE_NEW_INIT_S1         true
#define USE_NEW_CHANNEL_EST     true

#endif
