close all; clear all; clc
%% LOAD DATA
path = ('sc_neighb_cell_pnts_cott.mat');
load(path)
fullData = sc_neighb_cell_pnts;

%% CALL METHOD
%% INPUT: DATA (2D xy coordinate)
%         plotflag: 1 is for ploting, 0 is for non-ploting
%  OUTPUT: centre (2D xy value)

M_final = FindCenter(fullData,1);