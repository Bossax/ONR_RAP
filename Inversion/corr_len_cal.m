% Corr_len_cal
% Find correlation length scale from measurement 
% use 1-D transects
% use 2-D plane
% compare length scales of both cases
clear
close all
% Load an inversion file: October 2018 (HEM)
cd /Users/testuser/Documents/ONR_RAP/Data/inversion_file
load('HEM_inverse_solution_Oct2018_originaldepth.mat')

% Extract ttp measurements contributed from ocean features and plot tx map
d_recov_ocean = G(:,1:end-3)*m_recov(1:end-3);

f1 = figure(1);
s1 = scatter(tx_lon,tx_lat);