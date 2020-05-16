# ONR_RAP

### Routine: A folder contains scripts which are used to produce files from raw data.

- *read_binary_data_all_files.m*: produces .mat posmv files from POSMV binary files.

- *read_binary.m*: produce a single posmv file. 

- *find_transmission_info_all_files_general.m*: produces .mat tx files from Scarlette files containing latitude, longitude, times, heading, velocities of the vessel when transmissions were made. The tx files are created based on hours of file creation (grouping is based on file’s names). The script requires POSMV files created by read_binary_data_all_files. A transmission time provided by a Scarlette file is matched with the closet POSMVdata point in time.
 
- *find_reception_info.m*: produces .mat rx files from HEM audio containing estimated arrival times, actual arrival times, and SNRs. 

- *find_reception_info_icListen.m*: same as above but works with icListen audio.
	
### RAP_function: A folder contains functions which are called by other scripts.

- *hyd_audio_prep/icListen_audio_prep.m*: read audio files and assign timestamps to data points

- *timing_lock.m*: A function nested in Scarlette transmission codes to fix transmission times to second 0 and second 30

- *tx_rx_extraction.m*: accesses tx and rx files, packs data, and puts out associated tx/rx data based on given date and time.

**x_corr**: contains arrival detection codes

- *ACO_cross_correlation_ideal/_icListen.m*: called by find_reception_info/_icListen. Finds direct arrivals using ideal replicas to calculate cross-correlations.

- *ACO_cross_correlation_ideal_bottom_bounce/_icListen.m*: same as above but also detects the second peak 

**ray_tracing_w_curvature**: contains multiple versions of ray tracing codes

- *ray_trace_w_earth_flattening.m*: traces a direct ray path based on a given ellipsoidal length (surface distance) and transducer depth. Called by find_reception_info/_icListen.

**azimuth**: contains azimuth.m function used by many scripts.

**inversion_function** : contains functions called by SS_inversion.m SS_inversion_icListen.m

- *obs_matrix3D_v2.m*: needs to be specified a vertical mode number. returns a row of an observation matrix G corresponding to an acoustic ray. Also returns an integrand value of a ray corresponding to the given vertical mode. 

- *data_error_matrix.m*: returns a data-data covariance matrix Cd. 


### Inversion: contains script to perform inversion procedure
For inverse result data files, go to Ikepili > RAP > Data > inversion_file

- *Simulation_2D.m*: simulates case studies in of 2D ocean (depth independence).

- *Simulation_3D.m*: simulates case studies in of 3D ocean (with EOF modes). Need to uncomment the mode representation model analysis section in.

- *SS_inversion.m*: inverts for sound speed field using HEM TTP data and calculates associated metrics.

- *SS_inversion_icListen.m*: inverts for sound speed field using icListen TTP data and calculates associated metrics.

- *Analyze_inversion.m*: used with SS_inversion script to visualize reconstructed measurement and compare SS fields of different solutions.

- *update_hyd_pos.m*: calculates a new hydrophone geodesic position using the original position and position offsets given by inverse solutions.

- *corr_len_cal.mlx*: calculates correlation length scales from travel time data and plot variograms of the datasets.

### Data
**TX_RX_Output** contains transmissiona and reception data files calcualted from raw data using scripts in the Routine folder. These are the files used for subsequent analyses in this project.

- h….ctd: CTD cast data

- ioc_station…csv: tide data

**inversion_file**: contains inverse solutions of sound speed perturbation fields. Import to Matlab and plot data using *SS_inversion(_icListen).m* and *Analyze_inversion.m*.

- *EOF_SS.mat*: 4 normalized EOF modes with corresponding scaled loadings.

- *RMS_SS.mat*:RMS sound speed profile from EOF analysis


### EOF Analysis

- *create_HOTS_EOF.m*: calculates the first 4 EOF modes for the SS at ACO
- *create_HOTS_EOF_temp_sal.m*: script calculates the first two modes for temperature and salinity at ACO
- *EOF_analysis.m*: function file that computes the eigenvectors and eigenvalues for a given dataset. Uses the function "acf.m" for autocorrelation.

### Travel_Time_Perturbation_plot

- *parameters_and_TTP_plot.m*: creates various plots of travel time perturbations of either “HEM” or “icListen” hydrophones

- *parameters_and_TTP_plot_spin.m*: creates various plots of travel time perturbations of the HEM hydrophone. Display only results from spin courses.

- *HEM_icListen_TTP_plot.m*: creates travel time perturbation plots of the HEM and the icListen and compare mutual data points between the two.

### Ray_Trace_Dzieciuch
ray tracing scripts from Dzieciuch. Used to compare RAP ray tracing with long-range ray tracing method

### Test_script
Contains individual analysis scripts and test scripts

- *check_tx.m*: plot transmission information

- *check_posmv.m*:plot POSMV information

- *check_rx.m*:plot plot acoustic reception information

- *HEM_depth from_One_way_traveltime.m*: Analysis on the ACO depth issue which was raised by a conflict between a reult from hydrographic survey and the previous estimate. Computes the ACO depth using one way acoustic travel times of transmissions made from the R/V KM when it is overhead the observatory.

- *Compare_Planar_and_Raycoord*: test ray tracing in the ray coordinate system and compre its results with the palnar ray tracing method used in RAP.

- *ray_tracing_function_withETF_test.m*: Ray tracing in various coordinate systems with Earth Falttening Transformation.

- *Compare_Sep_Oct_Nov_SS.mP*: compare sound speed profiles from CTD casts in Septemebr, October, and November 2018

- *Dopper_plot.m*: plot complex envelopes produced by cross-correlating a reception with frequency-shifted replicas.

- *geoidheight_test.m*: plot geoid height map over the ACO area.

- *PSD_HEM/PSD_icListen*: plot periodograms of acoustic receptions from HEM and icListen data

## ###### Below is folders pertaining to early analysis of the project ###########

### Signal_Arrival: contains scripts to plot signal arrival patterns

- *Signal_Arrival.m*: plots acoustic signals arrival pattern and corresponding surface ranges (Old Code since July 2018).

- *TTP_SNR_Acoustic_savefile_(HEM/icListen).m*: produces .mat files containing acoustic signals and complex envelopes (ttp_snr_sig_data_2018... .mat in ttp_snr_plot folder).

- *TTP_SNR_Acoustic_plot_(HEM/icListen).m*: plots received acoustic signals and corresponding complex envelopes of signals in a specific period. Used to have close-up look of signal receptions and complex envelopes of individual transmission.

- *TTP_SNR_Acoustic_plot_w_psd_icListen.m*: same as above, but also includes power spectra of individual receptions (used for spin courses).

### Regression

 - *Regression_Analysis_ttp.m*: implements Lasso regression on travel time perturbation data to determine non-physical factor contributions to the measurements.

### IRIG_B

- *IRIGB_decoder_test.m* – decode IRIG-B time code recorded in October 2018 (saved in Scarlette files)


### Bottom_bouncing_detection: scripts to find bottom bounces

- *Bottom_bouncing_detection_HEM/ic.m*: plots acoustic signals and complex envelopes and determine the first 2 peaks of the envelopes.

- *All_arrival_detection.m*: finds all signal arrivals which exceed the set SNR threshold

### RTX
files related to RTX GPS data analysis to verify functionality of the POS-MV

### SNRvRay_plot
visualize SNR vs surface distance


## Notes

1. *distance.m* is called in many scripts which is a function in Mapping Toolbox. This function replaces *dist.m* function which was found to produce inaccurate geodesic distances.

2. Need to install *Aerospace Toolbox* and *Mapping Toolbox*


## Work flow

1. Generating POSMV matlab file: Use read_binary.m for manually choose which POSMV binary file will be converted, or use read_binary_data_all_files.m to automatically convert all binary files in a folder to matlab data files.

2. Extracting transmission information from Scarlette outputs: Use find_transmission_info_all_files_general_.m to produce hourly tx files from Scarlette files. 

3. Downloading audios: For the HEM, use convert_raw_to_mat.m to convert raw HEM data to matlab and .wav files. For the icListen, the files are already saved in .wav extension.

4. Calculating acoustic travel times and saving reception data: Use find_reception_info_.m to produce rx files. This scripts calls for audios, tx files, and ideal replicas to calculate for reception times. 

5. Visualizing measurements: Use parameters_and_TTP_plot_.m to visualize the measurements based on rx files. Make various kinds of plots.

6. Inverting for a sound speed perturbation field: Use SS_inversion.m to solve for an inverse sound speed perturbation field. This script should be executed section-wise (if spatial filtering is not necessary, skip the spatially-filtering section). Also use Analyze_inversion.m for additional visualization.


