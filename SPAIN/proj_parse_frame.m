function [proj] = proj_parse_frame(c, param, data_gapped)
% PROJ_PARSE_FRAME perfoms projection of coefficients c
% onto a multidimensional interval for the use of audio inpainting.
% 
% Input parameters
%       c               vector of input coefficients
%       param           structure containing mask and frame F
%       data_gapped     original degraded signal
%
% This projection
%       proj(z) = argmin_{u} ||z - u||_2 s.t. Au \in [b_1, b_2]
% can be evaluated as 
%       proj(z) = z-A^+(Az - proj(Az)),  where A^+ denotes pseudoinverse
% 
% The projection proj(Az) is computed by function proj_time.
% 
% Please note that this particular function works only for Parseval tight frame 
% (only in this case the pseudoinverse is identical to the analysis of the signal)


% Synthesis of the signal 
syn = frsyn(param.F, c);

% Compute proj(Az)
proj_temp = proj_time(syn, param.mask, data_gapped);

% Final projection
proj = c - frana(param.F, syn-proj_temp);

end