function error = sem(data)
%-----------------------------------------------------------------------------------------------------------------------------
% Function to get standard error
% Author@ Rohit
% INPUTS:-
%   -- data
% OUTPUTS:- 
%   -- error
%-----------------------------------------------------------------------------------------------------------------------------

error = std(data)/sqrt(length(data));

