function [channel] = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,type)
% GENERATE_RAYLEIGH_CHANNEL(nb_taps,type)
%
% Generate pure Rayleigh channel of length nb_taps.
% 
% Arguments:
%     'delta_t' : channel delay profile in ns
%     'nb_taps': number of nb_taps in the channel.
%
%     'type': determine the normalization of the generated channel:
%         "normalized": The norm of the generated channel is exactly 1.
%         "statistical": The norm of the generated channel is 1 on average.
%         "epa": Extended pedestrian channel model
% Return:
%     'channel': column vector of length nb_taps of the generated channel.
%
% By Pierre Baudoux
% Last Update: February 2th 2019

switch type
    case "normalized"
        channel = (randn(nb_taps,1) + 1i*randn(nb_taps,1));
        norm = sum(abs(channel).^2)/nb_taps;
        channel = channel ./ sqrt(norm*nb_taps);

    case "statistical"
        channel = 1/sqrt(2*nb_taps) * ( randn(nb_taps,1) + ...
                                    1i*randn(nb_taps,1) );
    case "epa"
        tap_delay = (0:delta_t:(nb_taps-1)*delta_t);
        var_tmp   = zeros(1, length(tap_delay));        
        var_tmp([1 4 8 10 12 20 42]) = 10.^( [0 -1 -2 -3 -8 -17.2 -20.8]./10 );                     % EPA parameters (found on internet)      
        var_tmp  = var_tmp./(sum(var_tmp)).*nb_taps;                                                % normalized CIR power to nb_taps (a.u.)
        channel = (sqrt(var_tmp) .* 1/sqrt(2).*(randn(1,nb_taps) + 1i*randn(1,nb_taps))).';
        
    case "decorrelated"
        channel = 1/sqrt(2)*( randn(nb_subcar,1) + ...
                                    1i*randn(nb_subcar,1) );
%           norm = sum(abs(channel).^2);
%           channel = channel ./ sqrt(norm);
    otherwise
        error('Invalid argument: type');
end

end