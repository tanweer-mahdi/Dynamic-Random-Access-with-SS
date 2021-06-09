# Dynamic-Random-Access-with-SS
This repository contains all the functions wrote for Multi-user Detection in Uplink Grant-Free NOMA with Sinusoidal Sequences
Description of the codes:

1. mtc_data
This function generates data for an uplink mMTC system for given set of parameters.

2. snapdim.m
This function takes a recived signal data and create snapshots for those to be used in SUBSPACE based algorithm like ESPRIT

3. esprit_aud.m
This function takes a data provided by snapdim.m and estimates frequencies as well as UEs using ESPRIT. It also returns the estimate of variance.

4. spice_aud.m
This function takes data provided by snapdim.mat and returns active users set like esprit_aud.m. It uses a fast SPICE implementation, accelerated by FFT

5. act_detect.m
This function take a least square estimate and perform Neyman-Pearson criteria based statistical hypothesis tests on all timeslots to detect user activity. It returns an activation matrix and a revised least square estimate if required

6. act_detect_npc.m
Does the same thing as 5 but doesn not carry out the revision of least square estimate 

7. channel_esimator.m
Estimates the channels of active users. Also deliver the reliable UE set. 

8. data_detection.m
Detects the data of reliable UEs provided by 7. It only considers 4-QPSK at this point. 

9. fun_error.m
This function calculates Activity Error Rate (AER), Symbol Error Rate (SER) and Net Normalized Mean Square Error (NNMSE)

All functions with the prefix "vary" is baed on above functions. They are written for carrying out Monte Carlo simulations 

This work has been selected for publication in the book "5G and Beyond: The Futuristic IoT".
