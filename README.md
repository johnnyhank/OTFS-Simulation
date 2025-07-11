# OTFS (Orthogonal Time Frequency Space) Modulation System

This repository contains a MATLAB implementation of an OTFS (Orthogonal Time Frequency Space) modulation system. The system includes transmitter, channel, and receiver components, along with various OTFS variants such as RZP-OTFS, RCP-OTFS, CP-OTFS, and ZP-OTFS.

## System Overview

The system is designed with the following parameters:
- Number of Doppler bins (time slots): N = 16
- Number of delay bins (subcarriers): M = 64
- Carrier frequency: 4 GHz
- Subcarrier spacing: 15 kHz
- Block duration: T = 1/Δf
- QAM modulation size: 4 (4-QAM)

## System Components

### 1. Transmitter (`otfs_transmitter.m`)
The transmitter implements three different OTFS modulation methods:
- Method 1: Direct implementation using equations (4.17) and (4.20)
- Method 2: Using permutation matrix and Kronecker product (Eq. 4.35)
- Method 3: Alternative implementation of Method 2

Key features:
- Generates random input bits
- Performs QAM modulation
- Creates OTFS delay-Doppler frame
- Implements multiple modulation methods for flexibility

### 2. Channel (`otfs_channel.m`)
The channel module supports multiple channel models and OTFS variants:

Channel Models:
- EPA (Extended Pedestrian A)
- EVA (Extended Vehicular A)
- ETU (Extended Typical Urban)

OTFS Variants:
- Standard OTFS
- RZP-OTFS (Reduced Zero Padding)
- RCP-OTFS (Reduced Cyclic Prefix)
- CP-OTFS (Cyclic Prefix)
- ZP-OTFS (Zero Padding)

Features:
- Implements Time-Delay Line (TDL) channel model
- Supports both standard and synthetic channel parameters
- Handles various response calculation methods

### 3. Receiver (`otfs_receiver.m`)
The receiver implements multiple demodulation and detection methods:

Demodulation Methods:
- Method 1: Using equations (4.24) and (4.27)
- Method 2: Using permutation matrix (Eq. 4.35)
- Method 3: Alternative implementation of Method 2

Detection Methods:
- LMMSE in delay-Doppler domain
- LMMSE in time domain

Features:
- AWGN noise addition
- Channel estimation
- Symbol detection and bit recovery
- Condition number monitoring for numerical stability

## Performance Analysis

The system includes comprehensive performance analysis capabilities:
- BER (Bit Error Rate) vs SNR analysis
- Monte Carlo simulations
- Comparison between different OTFS variants
- Performance visualization through plots

## Usage

1. Set system parameters in `otfs_simulation.m`
2. Choose channel model and OTFS variant
3. Run the simulation
4. View performance results and BER plots

## Dependencies

- MATLAB (with Signal Processing Toolbox)
- Communication Systems Toolbox