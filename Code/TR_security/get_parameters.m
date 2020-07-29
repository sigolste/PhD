%% 
% Fonction that returns the communication's parameters in function of the
% investigated scenario


%% 1) COMMON PARAMETERS:

global c 

% Constant allocation
c = 3e8 ;                                                   % Light velocity     

cdr = pi / 180 ;                                            % Degrees to radians conversion factor
crd = 180 / pi ;                                            % Radians to degrees conversion factor


% EM inputs
fc = 2e9 ;                                                  % Carrier frequency
lambda_c    = c / fc ;                                      % Carrier wavelength

% Scenario inputs
nb_TX       = 1;                                            % Number of TX antennas
fd = 50e6 ;                                                 % Symbol rate (maybe it is better to fix the overall bandwidth? and then calculate the symbol rate as delta_f/Ndim change fs into frate, and define fd = fs * N
delta_f = 2 * fd ;                                          % Occupied frequency bandwidth
rolloff    = 0.01 ;                                         % Rolloff factor of RCC filter
fs = ceil( 1 + rolloff ) * delta_f ;                        % Base Band sampling frequency

nb_channels = 3;
M           = 2;                                            % Constellation size
k           = log2(M);                                      % Number of bit per symbols
nb_subcar   = 2;                                           % Number of subcarriers  
nb_taps     = 10;                                           % Channel taps - corresponding to the EPA channel model
nFFT        = nb_subcar;                                    % fft size
length_CP   = nb_taps;                                      % Size of cyclic prefix
nb_block    = 200;
delta_t     = 10;                                           % in ns, for channel delay profil





%% 3) INVESTIGATED SCENARIOS

%scenario_name of type "xparam_yparam_variableparam_constantparam"

switch scenario 
    
    % BER
    %1.
    case "ber_ebno_bor_alpha"
        nb_block    = 300;          % Need lot of blocks
        BOR         = [2 4 8];      % Back Of Rate (upsampling factor)
        EbN0        = [-5:20];      % Energy per bit to noise power spectral density ratio
        alpha       = 1;           % Percentage to the useful data (--> 1-alpha corresponds to the energy used for the AN)
    
    %2.
    case "ber_ebno_alpha_bor"
        nb_block    = 400;          % Need lot of blocks
        BOR         = [4];          % Back Of Rate (upsampling factor)
        EbN0        = [-5:20];      % Energy per bit to noise power spectral density ratio
        alpha       = [1 .95 .1];     % Percentage to the useful data (--> 1-alpha corresponds to the energy used for the AN)
    
    %3.
    case "ber_alpha_bor_ebno"
        BOR         = [2 4 8];          % Back Of Rate (upsampling factor)
        EbN0        = [15];      % Energy per bit to noise power spectral density ratio
        alpha       = 1:-0.05:0;    % Percentage to the useful data (--> 1-alpha corresponds to the energy used for the AN)
        
    %4.    
    case "ber_alpha_ebno_bor" 
        nb_block    = 400;
        BOR         = [4];          % Back Of Rate (upsampling factor)
        EbN0        = [5 10 15 20];      % Energy per bit to noise power spectral density ratio
        alpha       = 1:-0.05:0;    % Percentage to the useful data (--> 1-alpha corresponds to the energy used for the AN)
    
    %5.
    case "ber_bor_ebno_alpha"
        BOR         = [2 4 8 16];
        EbN0        = [5 10 15 20];
        alpha       = .9;
        
    %6.
    case "ber_bor_alpha_ebno"
        BOR         = [2 4 8 16];
        nb_block    = 200;
        alpha       = [1 .95 .1];
        EbN0        = [15];
        
        
        
        
        
    % SECRECY RATE
    
    %1.
    case "secrecy_bor_ebno_alpha"
        nb_block    = 20;
        BOR         = [2 4 8 16];
        EbN0        = [0 5 10 15 20];
        alpha       = .5;
        
    %2.    
    case "secrecy_bor_alpha_ebno"
        nb_block    = 1;         
        BOR         = [2 4 8 16];
        EbN0        = 20;
        alpha       = [1 .95 .5 .1]; 
        
    %3.    
    case "secrecy_alpha_bor_ebno"
%         nb_block    = 20;
%         BOR         = [2 4 8 16];
%         EbN0        = 15;
%         alpha       = 1:-0.05:0.05;
        nb_channels = 100;
        nb_block    = 1;
        BOR         = [2 4 8];
        EbN0        = 5;
        alpha       = 1:-0.01:0;
%         alpha_1       = 1:-0.05:0.7;
%         alpha_2       = 0.7:-0.002:0.5;
%         alpha_3     = 0.5:-0.05:0;
%         alpha = [alpha_1 alpha_2 alpha_3];
        
    %4.     
    case "secrecy_alpha_ebno_bor"
        nb_block    = 5;
        BOR         = [2];
        EbN0        = [0:5:25];
        alpha       = 1:-0.01:0;
%         alpha_1       = 1:-0.05:0.7;
%         alpha_2       = 0.7:-0.002:0.5;
%         alpha_3     = 0.5:-0.05:0;
%         alpha = [alpha_1 alpha_2 alpha_3];
        
    %5.
    case "secrecy_ebno_alpha_bor"
        nb_block    = 2;
        EbN0        = [-5:2:45];
        alpha       = [1 .95 .5 .1];
        BOR         = [2];
    
    %6.
    case "secrecy_ebno_bor_alpha"
%         nb_block    = 40;    
%         EbN0        = [-5:2:45];
%         BOR         = [2 4 8 16];
%         nb_taps     = 20;  
%         alpha       = .9;
        nb_block    = 2;    
        EbN0        = [-5:30];
        BOR         = [2 4];
        alpha       = .9;
    %otherwize error("Not a considered scenario")
    
    %7.
    case "secrecy_alpha_bor_ebno_match_filt1"
        nb_channels = 100;
        nb_block    = 1;
        BOR         = [2];
        EbN0        = 100;
        alpha       = 1:-0.1:0;
end




%% Dependant parameters
snr         = EbN0+10*log10(k);                             % Signal to noise Ratio (dB)
snr_db_dec  = bsxfun(@power, 10, (EbN0+10*log10(k))./10);    % Signal to noise Ratio (decimal)
nb_symb_per_block = nb_subcar./BOR;                         % Number of symbol per OFDM block
nb_symb     = nb_block.*nb_symb_per_block;
nb_bit      = nb_symb .* k ;                                % Number of bits to process