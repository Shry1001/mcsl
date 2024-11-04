EXP2
SECTORING
clc;
close all;
clear all;
sectoring_type = input('Enter the type of sectoring (360, 120, or 60 degrees): ', 's');
if strcmp(sectoring_type, '360')
N = 12; io = 6; n = 4;
S_I = (sqrt(3 * N))^n / io;
fprintf('For 360-degree sectoring:\nS/I in standard scale: %.4f\nS/I in dB scale: %.2f dB\n', S_I, 10 * log10(S_I));
elseif strcmp(sectoring_type, '120')
N = 7; io = 2; n = 4;
S_I = (sqrt(3 * N))^n / io;
fprintf('For 120-degree sectoring:\nS/I in dB scale: %.2f dB\n', 10 * log10(S_I));
elseif strcmp(sectoring_type, '60')
N = 4; io = 1; n = 4;
S_I = (sqrt(3 * N))^n / io;
fprintf('For 60-degree sectoring:\nS/I in dB scale: %.2f dB\n', 10 * log10(S_I));
else360
disp('Invalid sectoring type. Please enter "360", "120", or "60".');
end

EXP3
Implementation of Free space and Two Ray Model, understanding of nature of Path Loss.
If a transmitter produces 50 watts of power, express the transmit power in units of (a) dBm, and (b) dBW. If 50 watts is applied to a unity gain antenna with a 900 MHz carrier frequency, find the received power in dBm at a free space distance of 100 m from the antenna, What is P (10 km) . Assume unity gain for the receiver antenna. 
Input these values:
Pt= 50 W, Gt=1, Gr=1, d=10 Km, f=900 MHz.
1.	Let the nature of the model be user defined.
2.	Take user defined values of Pt, f, Gt, Gr.
3.	Vary distance between Tx and Rx from 10 km to 100 km in steps of 10.
4.	Find power received using Free Space model in dBW and dBm.
5.	Declare ht=50m and hr=1.5m and find power received in dBW and dBm using Two Ray Gnd Reflection model for same range of distances.
6.	Develop understanding of PL vs frequency
Vary the frequency from 900 MHz to 1800 MHz in steps of 1MHz and study the curve obtained for received power for free space. 

% Given values
Pt = 50;  % Transmit power in Watts
Gt = 1;   % Transmit antenna gain
Gr = 2;   % Receive antenna gain
f = 900e6; % Frequency in Hz
ht = 50;  % Transmit antenna height in meters
hr = 1.5; % Receive antenna height in meters

% Speed of light
c = 3e8;

% Frequency in Hz
lambda = c / f;

% Distances in km
distances_km = 10:10:100;
distances_m = distances_km * 1000;  % Convert to meters

% Initialize vectors for received power
Pr_free_space_dBW = zeros(size(distances_m));
Pr_free_space_dBm = zeros(size(distances_m));
Pr_two_ray_dBW = zeros(size(distances_m));
Pr_two_ray_dBm = zeros(size(distances_m));

% Loop over distances
for i = 1:length(distances_m)
    d = distances_m(i);

    % Free Space Model
    Pr_free_space = Pt * Gt * Gr * (lambda / (4 * pi * d))^2;
    Pr_free_space_dBW(i) = 10 * log10(Pr_free_space);
    Pr_free_space_dBm(i) = 10 * log10(Pr_free_space * 1e3);

    % Two Ray Ground Reflection Model
    % Breakpoint distance
    d_break = 4 * pi * ht * hr / lambda;
    
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

% Display results
disp('Distance (km)  Free Space (dBW)  Free Space (dBm)  Two Ray (dBW)  Two Ray (dBm)');
for i = 1:length(distances_km)
    fprintf('%10.1f %16.2f %16.2f %13.2f %13.2f\n', distances_km(i), Pr_free_space_dBW(i), Pr_free_space_dBm(i), Pr_two_ray_dBW(i), Pr_two_ray_dBm(i));
end

% Plotting
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





Part B
Plot the graph for received power at the mobile using the 2-ray ground reflection model assuming the height of the transmitting antenna is 50 m and the receiving antenna is 1.5 m above ground. Vary the frequency from 900 MHz to 1800 MHz in steps of 1MHz and study the curve obtained for received power for two ray ground reflection model.
PART B
% Given values
Pt = 50;  % Transmit power in Watts
Gt = 1;   % Transmit antenna gain
Gr = 2;   % Receive antenna gain
d = 10e3; % Distance in meters (10 km)

% Speed of light
c = 3e8;

% Frequencies in Hz
frequencies_MHz = 900:1:1800;
frequencies_Hz = frequencies_MHz * 1e6;

% Initialize vector for received power
Pr_free_space_dBW = zeros(size(frequencies_Hz));
Pr_free_space_dBm = zeros(size(frequencies_Hz));

% Loop over frequencies
for i = 1:length(frequencies_Hz)
    f = frequencies_Hz(i);
    lambda = c / f;
    
    % Free Space Model
    Pr_free_space = Pt * Gt * Gr * (lambda / (4 * pi * d))^2;
    Pr_free_space_dBW(i) = 10 * log10(Pr_free_space);
    Pr_free_space_dBm(i) = 10 * log10(Pr_free_space * 1e3);
end

% Plotting
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

% Display results for first few frequencies for verification
disp('Frequency (MHz)  Free Space (dBW)  Free Space (dBm)');
for i = 1:10
    fprintf('%15.1f %16.2f %16.2f\n', frequencies_MHz(i), Pr_free_space_dBW(i), Pr_free_space_dBm(i));
end

Part C:
A mobile is located 5 km away from a base station and uses a vertical  monopole antenna with a gain of 2.55 dB to receive cellular radio signals. The E-field at 1 km from the transmitter is measured to be 10 ^-3 V/m.The carrier  frequency used for this system is 900 MHz. 
(a) Find the length and the effective  aperture  of the receiving antenna.
 (b) Find the received power at the mobile using the 2-ray ground reflection model assuming the height of the transmitting antenna is 50 m and the receiving antenna is 1.5 m above ground.
f = 900e6; 
c = 3e8;
Gt = 1.791; 
d = 5e3; 
d_o = 1e3;
E_o = 1e-3; 
h_t = 50; 
h_r = 1.5; 
lambda = c / f;
L = lambda / 4; 
Ae = Gt * (lambda^2) / (4 * pi);
Ep_d = 2 * E_o * d_o * (2 * h_t * h_r * pi) / (d^2 * lambda);
Pr = (Ep_d^2) * Ae / 377;
disp(['Length(L): ', num2str(L), 'meters']);
disp(['Effective aperture (Ae): ', num2str(Ae), ' square meters']);
disp(['Electric field at distance d, Ep(d): ', num2str(Ep_d), ' V/m']);
disp(['Received power (Pr): ', num2str(Pr), ' Watts']);



EXP4
a. To predict the type of small scale fading channel according to coherence bandwidth. 
b. To predict the type of small scale fading channel according to Doppler spread.
clc;
clear all;
close all;
% Part A: Coherence Bandwidth and Channel Type Prediction
% Step 1: Take user-defined values of RMS delay spread (τ_rms) and frequency correlation (fc)
tau_rms = input('Enter the RMS delay spread (in seconds): ');
fc = input('Enter the frequency correlation (fc): ');
% Validate inputs
if tau_rms <= 0
fprintf('Error: RMS delay spread must be a positive number.\n');
elseif fc < 0.5 || fc > 1
fprintf('Error: Frequency correlation (fc) must be between 0.5 and 1.\n');
else
% Step 2: Find coherence bandwidth (Bc)
if fc >= 0.9
Bc = 1 / (50 * tau_rms);
elseif fc >= 0.5 && fc < 0.9
Bc = 1 / (5 * tau_rms);
end
fprintf('The coherence bandwidth (Bc) is: %.2f Hz\n', Bc);
% Step 3: Take user-defined values of signal bandwidth (Bs) and symbol period (Ts)
Bs = input('Enter the signal bandwidth (in Hz): ');
Ts = input('Enter the symbol period (in seconds): ');
% Predict the type of channel based on coherence bandwidth and signal bandwidth
if Bs <= Bc
fprintf('The channel is a Flat Fading Channel.\n');
else
fprintf('The channel is a Frequency Selective Fading Channel.\n');
end
end
% Part B: Doppler Shift, Coherence Time, and Channel Type Prediction
% Step 1: Take user-defined velocity (v) in m/s and carrier frequency (f_c)
v = input('Enter the velocity (in m/s): ');
f_c = input('Enter the carrier frequency (in Hz): ');
% Validate inputs
if v < 0
fprintf('Error: Velocity must be a non-negative number.\n');
elseif f_c <= 0
fprintf('Error: Carrier frequency must be a positive number.\n');
else
% Step 2: Calculate maximum Doppler shift (f_d)
c = 3e8; % Speed of light in m/s
f_d = (v / c) * f_c;
fprintf('The maximum Doppler shift (f_d) is: %.2f Hz\n', f_d);
% Step 3: Calculate coherence time (Tc)
Tc = (0.423 / f_d);
fprintf('The coherence time (Tc) is: %.2f seconds\n', Tc);
% Step 4: Take user-defined values of signal bandwidth (Bs) and symbol period (Ts)
Bs = input('Enter the signal bandwidth (in Hz): ');
Ts = input('Enter the symbol period (in seconds): ');
% Predict the type of channel based on coherence time and symbol period
if Ts <= Tc
fprintf('The channel is a Slow Fading Channel.\n');
else
fprintf('The channel is a Fast Fading Channel.\n');
end
end


EXP5
Part A: The total spectrum allocation is 12.5 MHZ and the guard band allocated at the edge of the spectrum is 10KHz for a AMPS system. Determine the number of channels available for a FDMA system.
Code: clc; clear all; close all; Bt = 12.5*10e6; Bg = 10*10e3; Bc = 30*10e3; N = (Bt-2*Bg)/Bc; fprintf('Number of available channels: %d\n', floor(N));

Consider Global System for Mobile, which is a TDMA system that uses 25 MHz for the forward link, which is broken into radio channels of 200 kHz. If 8 speech channels are supported on a single radio channel, and if no guard band is assumed, find the number of simultaneous users that can be accommodated in GSM. 
Code: clc; clear all; close all; Bt = 25*10e6; Bc = 200*10e3; speech_channels = 8; N_Radio_Channels = Bt/Bc; Total_Users = N_Radio_Channels*speech_channels; fprintf("Total number of simultaneous users: %d\n", Total_Users);

ii) If GSM uses a frame structure where each frame consists of 8 time slots, and each time slot contains 156.25 bits, and data is transmitted at 270.833 kbps in the channel, find (a) the time duration of a bit, (b) the time duration of a slot, (c) the time duration of a frame, (d) how long must a user occupying a single time slot must wait between two simultaneous solutions. Code:
clc; clear all; close all; data_rate_kbps = 270.833*10e3; bits_per_slot = 156.25; slots_per_frame = 8; time_duration_of_a_bit = 1 / data_rate_kbps; time_duration_of_a_slot = bits_per_slot*time_duration_of_a_bit; time_duration_of_a_frame = slots_per_frame*time_duration_of_a_slot; wait_time_between_transmissions = time_duration_of_a_frame-time_duration_of_a_slot; fprintf('Time Duration of a Bit: %.8f seconds\n', time_duration_of_a_bit); fprintf('Time Duration of a Slot: %.8f seconds\n', time_duration_of_a_slot); fprintf('Time Duration of a Frame: %.8f seconds\n', time_duration_of_a_frame); fprintf('Wait Time Between Two Simultaneous Transmissions: %.8f seconds\n', wait_time_between_transmissions);

iii) If a normal GSM time slot consists of 6 trailing bits, 8.25 guard bits, 26 training bits, and 2 traffic bursts of 58 bits of data, find the frame efficiency. Code: clc;
clear all; close all; trailing_bits = 6; guard_bits = 8.25; training_bits = 26; data_bits_per_burst = 58; number_of_bursts = 2; total_data_bits = number_of_bursts * data_bits_per_burst; total_bits_in_time_slot = trailing_bits + guard_bits + training_bits + total_data_bits; frame_efficiency = total_data_bits / total_bits_in_time_slot; fprintf('Frame Efficiency: %.4f\n', frame_efficiency * 100);

EXP6

function exp6()
% Main script to take user input and generate PN code
disp('Enter a 4-bit binary array (e.g., [1 0 1 0]):');
array = input('', 's'); % Read user input as a string
array = str2num(array); % Convert string to numeric array
% Debugging output to check the input array
disp('Input array:');
disp(array);
% Validate input
if length(array) != 4
error('Input must be a 4-bit binary array.');
endif
% Validate if all elements are binary
if any(array < 0 | array > 1)
error('Array must contain only binary digits (0 or 1).');
endif
% Initialize PN code and shift register
pn_code = [];
shift_register = array;
% Generate PN code (length = 15 for a 4-bit shift register)
for i = 1:15
% Append the leftmost bit to the PN code
pn_code = [pn_code, shift_register(1)];
% Calculate the feedback bit using XOR of specific bits
feedback_bit = xor(shift_register(1), shift_register(2)); % Simple feedback for example
% Shift the register
shift_register = [shift_register(2:4), feedback_bit];
endfor
% Display the PN code and its length
disp('Generated PN code:');
disp(pn_code);
disp(['Length of PN code: ', num2str(length(pn_code))]);
endfunction

EXP7 
Task 1.
1. Let N be the number of users in a system, who are allotted Walsh codes.
2. Initialize a variable w=1.
Construct a 2x2 Walsh matrix H. (Hint: H=[w,w,w,~w]
3. Let the number of times this construction be repeated for N users be m. Note that m=ceil(log2(N)).
4. Ok Repeat construction of H m times to accommodate N users.
5. Display H.
CODE:
N = 16;
w = 1;
m = ceil(log2(N));
H = w;
for i = 1:m
H = [H, H; H, -H];
end
disp('Walsh Matrix H:');
disp(H);

TASK 2.
1. Generate the orthogonality test matrix ‘O’ in which the result of orthogonality is stored. For example, the matrix should look as follows for 4 users:
𝑂=[1001000000001001]
Code :
%continued below previous code
O = zeros(N, N);
for i = 1:N
for j = 1:N
O(i, j) = dot(H(i, :), H(j, :)) / N;
end
end
disp('Walsh Matrix H:');
disp(H)
disp('Orthogonality Matrix O:');
disp(O);






