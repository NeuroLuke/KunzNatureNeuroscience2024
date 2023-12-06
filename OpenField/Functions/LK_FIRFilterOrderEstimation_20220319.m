function filterOrder = LK_FIRFilterOrderEstimation_20220319(cfg)
%
% LK_FIRFilterOrderEstimation_20220319 computes the filter order for an FIR
% filter using the Fred Harris rule, given some settings in the input
% structure:
% - Fs = sampling rate of the filter
% - deltaF = transition band, in the same units as Fs
% - attenDB = target rejection in dB
%
% Use as: filterOrder = LK_FIRFilterOrderEstimation_20220319(cfg);
%
% Reference:
% https://dsp.stackexchange.com/questions/37646/filter-order-rule-of-thumb
%
% Lukas Kunz, 2022 

% Fred Harris rule of thumb
filterOrder = (cfg.Fs / cfg.deltaF) * (cfg.attenDB / 22);