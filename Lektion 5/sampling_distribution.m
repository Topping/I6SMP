function Y = sampling_distribution(x, distribution)

% *************************************************************************
% Sampling From Any Distribution:
%
%   [Y] = sampling_distribution(x, distribution)
%
%   Signal Randomly Sampled From A Given Distribution
%
%   DISTRIBUTION:
%           * Rayleigh               (DISTRIBUTION='Rayleigh')
%           * Exponential            (DISTRIBUTION='exp')
%
% *************************************************************************

% Rayleigh distribution
if strcmpi(distribution,'Rayleigh')
    Y = raylinv(x);
    figure()
    plot(Y,'.')
    xlabel('i'), ylabel('x'), title('x is Rayleigh distributed')

elseif strcmpi(distribution,'exp')
    Y = expinv(x);
    figure()
    plot(Y,'.')
    xlabel('i'), ylabel('x'), title('x is exponential distributed')    
    
else 
    error(['Unrecognised ''DISTRIBUTION'' parameter ' ...
          '(see help for a list of currently implemented distributions).']);
end

