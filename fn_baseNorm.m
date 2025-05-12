function [out] = fn_baseNorm(sig)
% fn_baseNorm performs standard z-score normalization of a signal. 
% @ Aamir Abbasi. Based on the function written by Tanuj Gulati.
%
% ----------- INPUTS
% sig = A vector of length N.
%
% ----------- OUTPUTS
% out = A normalized vector of length N. 
out=(sig-mean(sig))/std(sig);
