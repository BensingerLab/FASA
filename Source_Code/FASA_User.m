%This file is a part of FASA (Fatty Acid Source Analysis). Copyright (c) 2018, created by Moses Q. Wilks, Quan D. Zhou, and Joseph P. Argus (UCLA).
%Use by others is subject to the terms of the BSD 3-Clause License (see “license.txt”).

% User Handles:
sheets_to_model=[2,4];
cutoff=-1;
q=0.0107;
e=0.99;
assume_no_label_diffusion=0;
MC_successes=10;

% Calling FASA_Step2.m
FASA_Step2(sheets_to_model,cutoff,q,e,assume_no_label_diffusion,MC_successes)