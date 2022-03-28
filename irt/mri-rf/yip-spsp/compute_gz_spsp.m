 function [kp, gz, kz, kf, t] = compute_gz_spsp(kp)
% function [kp,gz,kz,kf,t] = compute_gz_spsp(kp)
% Function that computes z gradient waveform accompanying the SPSP pulse,
% using parameters specified in structure kp, generated by kparameterSPSP.m.
% In units of g/cm. Currently can only generate waveforms made of trapezoids.
% Outputs:
% kp: updated trajectory parameter structure
% gz: z gradient waveform in g/cm
% kz,kf: together they specify the trajectory in SPSP k-space (kf-kz space)
%
% Chun-yu Yip, 4/1/2009

% Physical parameters
gam = 26751;                            % rad/sec/g; gyromagnetic ratio
gambar = gam / 2 / pi;                      % Hz/g

% Trapezoid area calculation
if kp.dgdtmax * (kp.T / 2) / 2 <= kp.gmax              % the lobe is triangular.
    gzarea = kp.T / 2 * kp.dgdtmax * (kp.T / 2) / 2 / 2;  % g s /cm
else                                             % this lobe is a trapezoid.
    r = kp.gmax / kp.dgdtmax;                      % s; rise time of gradient
                                                 % to plateau
    gzarea = kp.gmax * r + kp.gmax * (kp.T / 2 - 2 * r); % g s /cm
end

% Trapezoid generation
gztrap_whole = dotrap(gzarea, kp.gmax, kp.dgdtmax, kp.pointtime);
gztrap_half = dotrap(gzarea / 2, kp.gmax, kp.dgdtmax, kp.pointtime);

signz = +1;
gz = [];
for ii = 1:1:kp.Ntraps

  if ii == kp.Ntraps
   factor = 1;
  else
   factor = 1;
  end

  signz = -signz;
  gz = [gz factor * signz * gztrap_whole];          % Concatenate trapezoids
                                                % with alternating signs.
end

signz = -signz;
gz = [gz signz * gztrap_half];                    % Add refocusing lobe
gz = gz.';

t = kp.pointtime * [0:1:length(gz) - 1].';          % s; time vector

% Do "backward integral" to obtain excitation k-space trajectory
intgz = cumtrapz(gz);
kz = -gambar * kp.pointtime *  (ones(length(gz), 1) * intgz(end) - intgz); % /cm
kf = t - t(end);
kp.npnts = length(gz);              % pulse length (number of samples)
kp.pw = kp.pointtime * kp.npnts;    % s; pulse length (actual duration)

if 0 % Display gradient waveform
    figure;
    plot(t, gz);
    xlabel('time (sec)');
    ylabel('g/cm');
    title('Gz waveform');
    grid;
end
