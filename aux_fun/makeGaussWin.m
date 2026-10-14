function gw = makeGaussWin(fs, width_ms)
% makeGaussWin generates a Gaussian window with a specified FWHM (in ms)
%   fs        - sampling frequency (Hz)
%   width_ms  - desired Gaussian FWHM in milliseconds
%
%   gw        - Gaussian window (column vector)

    % convert width from ms to seconds
    w = width_ms * 1e-3;

    % convert FWHM to sigma (in samples)
    sigma = (w * fs) / (2 * sqrt(2 * log(2)));

    % define sample axis over ±3 sigma
    t = -3*sigma : 1 : 3*sigma;

    % construct the Gaussian
    gw = exp(-(t.^2) / (2 * sigma^2));
    gw = gw(:);  % column vector
end
