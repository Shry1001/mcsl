
Exp2: RELATION OF CLUSTER SIZE TO SI
Exp3: Free Space and Two ray
Exp4: Small Scale fading channel prediction
Exp5: Study of multiple access technique
Exp6: PN code Generation
Exp7: Generation of Walsh Codes
Exp9: Alamouti Coding for mimo systems

EXP – 2 [ RELATION OF CLUSTER SIZE TO S/I]
clc; close all; clear all;
function calculate_si_and_cluster_size()
  fprintf('Select sectoring type:\n');
  fprintf('1. 360 degree sectoring\n');
  fprintf('2. 120 degree sectoring\n')
  fprintf('3. 60 degree sectoring\n');
  sectoring_type = input('Enter choice (1/2/3): ');
  N = 0;
  io = 0;
  n = 4;
  switch sectoring_type
    case 1
      N = 12;
      io = 6;
      S_I = calculate_si(N, io, n);
      fprintf('360 degree sectoring:\n');
      display_si(S_I);
    case 2
      N = 7;
      io = 2;
      S_I = calculate_si(N, io, n);
      fprintf('120 degree sectoring:\n');
      display_si(S_I);
    case 3
      N = input('enter the number of channels: ')
      io = input('Enter io value for 60 degree sectoring: ');
      S_I = calculate_si(N, io, n);
      fprintf('60 degree sectoring:\n');
      display_si(S_I);
      cluster_size = calculate_cluster_size(S_I);
      fprintf('Recommended cluster size:\n');
      if cluster_size <= 3
        fprintf('Take cluster size 3.\n');
      elseif cluster_size <= 4
        fprintf('Take cluster size 4.\n');
      elseif cluster_size <= 7
        fprintf('Take cluster size 7.\n');
      else
        fprintf('Cluster size is out of range.\n');
      end
    otherwise
      fprintf('Invalid choice. Please select 1, 2, or 3.\n');
  end
end
function S_I = calculate_si(N, io, n)
  S_I = (((3*N)^(1/2))^n)/io;
end
function display_si(S_I)
  fprintf('S/I (Standard Scale): %.2f\n', S_I);
  fprintf('S/I (dB Scale): %.2f dB\n', 10 * log10(S_I));
end
function cluster_size = calculate_cluster_size(S_I)
  if S_I <= 3
    cluster_size = 3;
  elseif S_I <= 4
    cluster_size = 4;
  elseif S_I <= 7
    cluster_size = 7;
  else
    cluster_size = NaN;
  end
end
calculate_si_and_cluster_size();







EXP – 3 [ STUDY OF MOBILE CELLULAR PROPOGATION MODELS]
PART A] Simulate Friis Free Space Propagation Model and Two-Ray Ground Reflection Model
clc;clear all;close all;
Pt = 50;          % Transmit power in Watts
Gt = 1;           % Transmit antenna gain (unitless)
Gr = 2;           % Receive antenna gain (unitless)
f = 900e6;        % Frequency in Hz
ht = 50;          % Transmit antenna height in meters
hr = 1.5;         % Receive antenna height in meters
c = 3e8;          % Speed of light in m/s
lambda = c / f;   % Wavelength in meters
distances_km = 10:10:100;
distances_m = distances_km * 1000;  % Convert to meters
Pr_free_space_dBW = zeros(size(distances_m));   % Initialize vectors for received power
Pr_free_space_dBm = zeros(size(distances_m));
Pr_two_ray_dBW = zeros(size(distances_m));
Pr_two_ray_dBm = zeros(size(distances_m));
for i = 1:length(distances_m)
    d = distances_m(i);
    % Free Space Model
    Pr_free_space = Pt * Gt * Gr * (lambda / (4 * pi * d))^2;
    Pr_free_space_dBW(i) = 10 * log10(Pr_free_space);
    Pr_free_space_dBm(i) = 10 * log10(Pr_free_space * 1e3);
    % Two Ray Ground Reflection Model
    d_break = 4 * pi * ht * hr / lambda; % Breakpoint distance
    if d <= d_break
        % Free space equation before breakpoint distance
        Pr_two_ray = Pr_free_space;
    else
        % Two ray ground reflection model after breakpoint distance
        Pr_two_ray = Pt * Gt * Gr * (ht * hr / d^2)^2;
    end
    Pr_two_ray_dBW(i) = 10 * log10(Pr_two_ray);
    Pr_two_ray_dBm(i) = 10 * log10(Pr_two_ray * 1e3);
end
disp('Distance (km)  Free Space (dBW)  Free Space (dBm)  Two Ray (dBW)  Two Ray (dBm)');
for i = 1:length(distances_km)
    fprintf('%10.1f %16.2f %16.2f %13.2f %13.2f\n', distances_km(i), ...
        Pr_free_space_dBW(i), Pr_free_space_dBm(i), Pr_two_ray_dBW(i), Pr_two_ray_dBm(i));
end
figure;
subplot(2, 1, 1);
plot(distances_km, Pr_free_space_dBW, '-o', 'DisplayName', 'Free Space Model');
hold on;
plot(distances_km, Pr_two_ray_dBW, '-x', 'DisplayName', 'Two Ray Ground Reflection Model');
xlabel('Distance (km)');
ylabel('Received Power (dBW)');
title('Received Power (dBW) vs Distance');
legend show;
grid on;
subplot(2, 1, 2);
plot(distances_km, Pr_free_space_dBm, '-o', 'DisplayName', 'Free Space Model');
hold on;
plot(distances_km, Pr_two_ray_dBm, '-x', 'DisplayName', 'Two Ray Ground Reflection Model');
xlabel('Distance (km)');
ylabel('Received Power (dBm)');
title('Received Power (dBm) vs Distance');
legend show;
grid on;

PART B] Simulate Free Space Propagation Model by Varying Frequency
clc;clear all;close all;
Pt = 50;          % Transmit power in Watts
Gt = 1;           % Transmit antenna gain (unitless)
Gr = 2;           % Receive antenna gain (unitless)
d = 10e3;         % Distance in meters (10 km)
c = 3e8;          % Speed of light in m/s
frequencies_MHz = 900:1:1800;    % Frequency range from 900 MHz to 1800 MHz
frequencies_Hz = frequencies_MHz * 1e6;  % Convert MHz to Hz
Pr_free_space_dBW = zeros(size(frequencies_Hz));  % Initialize vector for received power
Pr_free_space_dBm = zeros(size(frequencies_Hz));
for i = 1:length(frequencies_Hz)
    f = frequencies_Hz(i);
    lambda = c / f;  % Wavelength calculation
    % Free Space Model
    Pr_free_space = Pt * Gt * Gr * (lambda / (4 * pi * d))^2;
    Pr_free_space_dBW(i) = 10 * log10(Pr_free_space);
    Pr_free_space_dBm(i) = 10 * log10(Pr_free_space * 1e3);
end
figure;
subplot(2, 1, 1);
plot(frequencies_MHz, Pr_free_space_dBW, '-o');
xlabel('Frequency (MHz)');
ylabel('Received Power (dBW)');
title('Received Power (dBW) vs Frequency');
grid on;
subplot(2, 1, 2);
plot(frequencies_MHz, Pr_free_space_dBm, '-o');
xlabel('Frequency (MHz)');
ylabel('Received Power (dBm)');
title('Received Power (dBm) vs Frequency');
grid on;
disp('Frequency (MHz)  Free Space (dBW)  Free Space (dBm)');
for i = 1:10
    fprintf('%15.1f %16.2f %16.2f\n', frequencies_MHz(i), ...
        Pr_free_space_dBW(i), Pr_free_space_dBm(i));
end

PART C] Simulate Free Space Propagation Model by Varying Frequency
clc;clear all;close all;
c = 3e8;            % Speed of light in m/s
distance = 5e3;     % Distance to mobile (5 km in meters)
gain_dB = 2.55;     % Gain in dB
f = 900e6;          % Frequency in Hz
E_field = 10^-3;    % Electric field in V/m
ht = 50;            % Height of the transmitting antenna in meters
hr = 1.5;           % Height of the receiving antenna in meters
Pt = 50;            % Transmit power in Watts
Gt = 1;             % Transmit antenna gain (unitless)
Gr = 10^(gain_dB / 10);  % Convert dB gain to linear scale for receiver
% (a) Length and Effective Aperture
lambda = c / f;     % Wavelength calculation
L = lambda / 4;     % Quarter-wavelength monopole length
Ae = (lambda^2 / (4 * pi)) * 10^(gain_dB / 10); % Effective aperture
disp(['Length of receiving antenna: ', num2str(L), ' meters']);
disp(['Effective aperture of receiving antenna: ', num2str(Ae), ' square meters']);
% (b) Received Power using Two-Ray Model
Pr_mobile = (Pt * Gt * Gr * (ht * hr)^2) / (distance^2); % Received power in watts
Pr_mobile_dBm = 10 * log10(Pr_mobile * 1e3); % Convert to dBm
disp(['Received Power at mobile (Two-Ray Model): ', num2str(Pr_mobile), ' watts']);
disp(['Received Power at mobile (Two-Ray Model): ', num2str(Pr_mobile_dBm), ' dBm']);

 


EXP – 4 [ SMALL SCALE FADING]
PART A] BASED ON COHERENCE BANDWIDTH
clc;clear all;close all;
rms_delay_spread = input('Enter the RMS delay spread (in seconds): '); % User-defined values in seconds
fc = input('Enter the frequency correlation coefficient (e.g., 0.9, 0.8, etc.): '); % User-defined values between 0 and 1
if fc >= 0.9
    Bc = 1 / (50 * rms_delay_spread); % Coherence bandwidth for fc >= 0.9
elseif fc >= 0.5 && fc < 0.9
    Bc = 1 / (5 * rms_delay_spread); % Coherence bandwidth for 0.5 <= fc < 0.9
else
    fprintf('Frequency correlation coefficient is not in a valid range (0.5 to 1).\n');
    return; % Exit if the frequency correlation is invalid
end
fprintf('Calculated Coherence Bandwidth (Bc): %.4f Hz\n', Bc);
signal_bandwidth = input('Enter the signal bandwidth (in Hz): ');
symbol_period = input('Enter the symbol period (in seconds): ');
symbol_rate = 1 / symbol_period;  % Calculate symbol rate
if signal_bandwidth < Bc
    fprintf('The channel is Flat Fading.\n');
elseif signal_bandwidth >= Bc && signal_bandwidth < symbol_rate
    fprintf('The channel is likely Frequency Selective Fading.\n');
elseif signal_bandwidth >= symbol_rate
    fprintf('The channel is likely to experience ISI (Inter-Symbol Interference).\n');
else
    fprintf('Invalid input for signal bandwidth or symbol period.\n');
end
OUTPUT VALUES PUT: 1] AMPS: [1.387*10^-6, 0.7, 30000, 1*10^-6]
                                         2]GSM: [1.387*10^-6, 0.7, 200000, 1*10^-6]
PART B] BASED ON DOPPLER SPREAD
clc;clear all;close all;
velocity = input('Enter the velocity (in m/s): ');
frequency = input('Enter the frequency (in Hz): ');
c = 3e8;  % Speed of light in m/s
max_doppler_shift = (velocity / c) * frequency;
fprintf('Maximum Doppler Shift: %.2f Hz\n', max_doppler_shift);
max_doppler_spread = 2 * max_doppler_shift;
fprintf('Maximum Doppler Spread: %.2f Hz\n', max_doppler_spread);
coherence_time = 1 / max_doppler_spread;
fprintf('Coherence Time: %.4f seconds\n', coherence_time);
signal_bandwidth = input('Enter the signal bandwidth (in Hz): ');
symbol_period = input('Enter the symbol period (in seconds): ');
if signal_bandwidth > (1 / coherence_time)
    fprintf('The channel is Fast Fading.\n');
else
    fprintf('The channel is Slow Fading.\n');
end
fprintf('Symbol Period: %.4f seconds\n', symbol_period);
OUTPUT VALUES PUT: 1] AMPS: [19.4, 900000000, 30000, 1*10^-6]
                                         2]GSM: [19.4, 900000000, 200000, 1*10^-6]


EXP – 5 [ STUDY OF MULTPLE ACCESS TECHNIQUES IN MCS]
PART A] FDMA:  CALCULATE NUMBER OF CHANNELS
clc;clear all;close all;
total_spectrum_MHz = 12.5;    % Total spectrum allocation in MHz
guard_band_kHz = 10;          % Guard band in kHz
channel_bandwidth_kHz = 30;   % Channel bandwidth in kHz (assuming for AMPS)
guard_band_MHz = guard_band_kHz / 1000;  % Convert guard band to MHz
channel_bandwidth_MHz = channel_bandwidth_kHz / 1000;  % Convert channel bandwidth to MHz
usable_spectrum_MHz = total_spectrum_MHz - guard_band_MHz;
number_of_channels = floor(usable_spectrum_MHz / channel_bandwidth_MHz);
fprintf('Number of available channels: %d\n', number_of_channels);
B_T = total_spectrum_MHz;       % Total spectrum allocation in MHz
B_G = guard_band_MHz;           % Guard band in MHz
B_C = channel_bandwidth_MHz;    % Channel bandwidth in MHz
N = (B_T - 2*B_G) / B_C;  % Calculate the number of channels using the alternative formula
fprintf('Number of available channels using alternative calculation: %d\n', floor(N));

PART B-1] TDMA: NUMBER OF SIMULTANEOUS USERS
clc;clear all;close all;
total_bandwidth_MHz = 25;   % Total bandwidth in MHz
channel_bandwidth_kHz = 200; % Bandwidth per channel in kHz
speech_channels_per_radio_channel = 8; % Number of speech channels per radio channel
total_bandwidth_kHz = total_bandwidth_MHz * 1000;
num_radio_channels = total_bandwidth_kHz / channel_bandwidth_kHz;
total_users = num_radio_channels * speech_channels_per_radio_channel;  % Calculate the total number of simultaneous users
fprintf("Total number of simultaneous users: %d\n", total_users);


PART B-2] TDMA: BIT, SLOT, FRAME NUMERICAL
clc;clear all;close all;
data_rate_kbps = 270.833;          % Data rate in kbps
data_rate_bps = data_rate_kbps * 1000;  % Convert kbps to bps
bits_per_slot = 156.25;            % Number of bits per slot
slots_per_frame = 8;               % Number of slots per frame
time_duration_bit_sec = 1 / data_rate_bps;       %  Time duration of a bit in seconds
time_duration_bit_us = time_duration_bit_sec * 1e6;  % Convert to microseconds
time_duration_slot_sec = bits_per_slot * time_duration_bit_sec;  % Time duration of a slot in seconds
time_duration_slot_ms = time_duration_slot_sec * 1e3;            % Convert to milliseconds
time_duration_frame_sec = slots_per_frame * time_duration_slot_sec; % Time duration of a frame in seconds
time_duration_frame_ms = time_duration_frame_sec * 1e3;             % Convert to milliseconds
time_between_transmissions_ms = time_duration_frame_ms;  % Time between two transmissions (same as frame duration)
fprintf("Time duration of a bit: %.2f microseconds\n", time_duration_bit_us);
fprintf("Time duration of a slot: %.2f milliseconds\n", time_duration_slot_ms);
fprintf("Time duration of a frame: %.2f milliseconds\n", time_duration_frame_ms);
fprintf("Time between two transmissions: %.2f milliseconds\n", time_between_transmissions_ms);

PART B-3] TDMA: FRAME EFFICIENCY
clc;clear all;close all;
trailing_bits = 6;
guard_bits = 8.25;
training_bits = 26;
data_bits_per_burst = 58;
number_of_bursts = 2;
total_data_bits = number_of_bursts * data_bits_per_burst;  % Calculate total data bits
total_bits_in_time_slot = trailing_bits + guard_bits + training_bits + total_data_bits;
frame_efficiency = total_data_bits / total_bits_in_time_slot;
fprintf('Frame Efficiency: %.4f (or %.2f%%)\n', frame_efficiency, frame_efficiency * 100);





EXP – 6 [ GENERATION OF WALSH CODE]       use input as [1 0 1 0]
PART A] PN code using last 2 bits
clc; clear all; close all;
user_input = input('Enter a 4-bit array as a row vector (e.g., [1 0 0 1]): ');
if length(user_input) ~= 4
    error('Array length must be 4');
end
register = user_input(:); % Initialize the shift register with user input as a column vector
pn_code = [];  % Initialize PN code storage
num_shifts = 2^length(register) - 1;  % length of the PN code sequence (2^4 - 1 for a 4-bit LFSR)
% Generate the PN code sequence
for i = 1:num_shifts
    pn_code = [pn_code; register(end)];
    new_bit = xor(register(3), register(4)); % Feedback from last 2 bits
    register = [new_bit; register(1:end-1)]; % Shift right with feedback
end
fprintf('The length of the PN code generated from a 4-bit shift register is %d\n', num_shifts);
fprintf('The generated PN code is:\n');
disp(pn_code');

PART B] PN code using middle 2 bits
clc;clear all;close all;
arr = input('Enter a 4-bit array as [a b c d]: ');
if length(arr) ~= 4
    error('Array length must be 4');
end
pn_code = [];   % Initialize PN code and iteration counter
iterations = 0;
initial_state = arr;
% Generate the PN code by shifting elements
while true
    pn_code = [pn_code arr(1)];
    feedback = xor(arr(2), arr(3));
    arr = [feedback arr(1:end-1)];
    iterations = iterations + 1;
    if isequal(arr, initial_state)
        break;
    end
end
disp('PN code:');
disp(pn_code);
disp(['Number of iterations: ', num2str(iterations)]);

PART C] PN code using last  2 bits
clc;clear all;close all;
user_input = input('Enter a 4-bit array as a row vector (e.g., [1 0 0 1]): ');
if length(user_input) ~= 4
    error('Array length must be 4');
end
register = user_input(:);  % Initialize register with user input as a column vector
taps = [1, 2];  % Feedback taps
pn_code = [];  % Initialize PN code storage
seen_states = containers.Map('KeyType', 'char', 'ValueType', 'any');
register_to_key = @(reg) num2str(reg');  % Helper function to create a unique key from the register state
seen_states(register_to_key(register)) = true;  % Add the initial state to the seen states map
% Generate PN code sequence
while true
    pn_code = [pn_code, register(1)];
    feedback_bit = 0;
    for i = 1:length(taps)
        feedback_bit = xor(feedback_bit, register(taps(i)));
    end
    register = [feedback_bit; register(1:end-1)];
    reg_key = register_to_key(register);
    if isKey(seen_states, reg_key)
        break;  % Sequence has started repeating
    else
        seen_states(reg_key) = true;
    end
end
fprintf('The length of the PN code generated from a 4-bit shift register with taps at positions %s is %d\n', mat2str(taps), length(pn_code));
fprintf('The generated PN code is:\n');
disp(pn_code);


EXP – 7 [ GENERATION OF WALSH CODE]
PART A] GENERATE WALSH MATRIX 
clc;clear all;close all;
function H = construct_walsh_matrix(N)
    if mod(log2(N), 1) != 0
        error('Number of users N must be a power of 2');
    end                       % Check if N is a power of 2
    w = 1;                    % Base value for Walsh matrix construction
    H_base = [w, w; w, -w];   % Initial 2x2 Walsh matrix
    m = ceil(log2(N));        % Number of iterations to construct the matrix
    H = H_base;               % Start with the base Walsh matrix
    for i = 2:m
        H = kron(H, H_base);
    end                       % Expand the Walsh matrix by the Kronecker product
    fprintf('Walsh matrix (size %d x %d):\n', N, N);
    disp(H);
end
N = input('enter value of N: ');
H=construct_walsh_matrix(N);
PART B] CHECK ORTHOGONALITY
O = zeros(N, N);
for i = 1:N
    for j = 1:N
        O(i, j) = dot(H(i, :), H(j, :)) / N;
    end
end
disp('Orthogonality Matrix O:');
disp(O);



EXP – 9 [ ALAMOUTI CODING FOR MIMO]
clc;clear all;close all;
Nt = 2;  % Number of transmit antennas
Nr = 2;  % Number of receive antennas
I = eye(Nt * Nr);  % Identity matrix of size 4x4 (for correlation matrix)
CH = randn(Nt * Nr, 1) / sqrt(2) + 1i * randn(Nt * Nr, 1) / sqrt(2);  % Complex channel coefficients for a Rayleigh fading channel
H = reshape(I * CH, Nt, Nr);  % Reshape to a 2x2 matrix for the channel matrix
disp('Channel Matrix H:');
disp(H);
Data = (2 * round(rand(Nt, 1)) - 1) / sqrt(Nt);  % Data vector (values are +/- 1) random data of dimension 2x1
disp('Random Data Vector:');
disp(Data);
snr = 10;  % Signal-to-noise ratio in dB
sig = sqrt(0.5 / (10^(snr / 10)));  % Noise variance based on SNR
n = sig * (randn(Nr, Nt) + 1i * randn(Nr, Nt));  % Complex noise matrix (2x2)
disp('Noise Matrix:');
disp(n);
D1 = Data(1);  % First data symbol
D2 = Data(2);  % Second data symbol
X = [D1, -conj(D2); D2, conj(D1)];  % Alamouti encoded data matrix
disp('Alamouti Encoded Data Matrix X:');
disp(X);
R = H * X + n;  % Received signal matrix
disp('Received Signal Matrix R:');
disp(R);
S1 = (conj(H(1, :)) * R(:, 1) + conj(H(2, :)) * R(:, 2)) / (norm(H(1, :))^2);  % Combine for D1
S0 = (conj(H(1, :)) * R(:, 2) - conj(H(2, :)) * R(:, 1)) / (norm(H(2, :))^2);  % Combine for D2
S = [S1; S0];  % Combiner Output vector
disp('Combiner Output S:');
disp(S);





