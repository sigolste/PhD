clear all;
%close all; 
set(0,'defaulttextinterpreter','latex')
set(0,'defaultAxesFontSize',32)
set(0,'defaultLineLineWidth',2)
set(0,'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultLegendInterpreter','latex')
set(0,'DefaultLineMarkerSize',15);
set(0, 'defaultFigurePosition',  [-1267  44   1256    872])

% set(0,'defaultMarkerIndices',[1:2:end])
%**************************************************************************
%
%   This script simulates an implementation of a secure Time Reversal
%   scheme in the frequency domain using OFDM. 
%
%   Frequency domain secure TR in a MISO configuration    
%   Complex Envelope (CE) simulations
%
%   OFDM to mitigate MPC effects
%
%   Artificial Noise added to secure the data transmission at Bob
%
%   Rayleigh channel model implemented
%   
%   Best global alpha found in order to maximize the secrecy capacity.
%   After that, waterfilling in order to maximize the instantaneous
%   capacity at Bob. It should not impact Eve's capacity s.t. the Secrecy
%   Capacity should be increased.
%
%
%   last update: 10 august 2020
%
%
%   by Sidney Golstein
%
%**************************************************************************




h = waitbar(0,'Simulation Progression...');

%% Parameters
% Simulation parameters
nb_run = 1000;              % number of experiments
nb_model = 3;


Q = 32;
U = [2 4 8 16];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters
EbN0_b = 15; % energy per bit over noise psd @Bob - dB
EbN0_e = [10]; % energy per bit over noise psd @Eve - dB
snr_b  = 20; %EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = 20; %EbN0_e + 10*log10(k);  % SNR @Eve

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance




% Optimal alpha decod1
alpha1_model_opt = optimalAlpha(U,snr_b,snr_e,"model1");

% Optimal alpha decod2
alpha2_model_opt = optimalAlpha(U,snr_b,snr_e,"model2");

% Optimal alpha decod2
alpha5_model_opt = optimalAlpha(U,snr_b,snr_e,"model5");

alpha_global = [alpha1_model_opt' alpha2_model_opt' alpha5_model_opt'];




alpha_opt           = zeros(Q,nb_run,length(U),nb_model);
% e_an_TX             = zeros(nb_run,length(U));
e_an_opt_TX         = zeros(nb_run,length(U),nb_model);

%
e_sym_decod1_b      = zeros(nb_run,length(U),nb_model);
e_noise_decod1_b    = zeros(nb_run,length(U),nb_model);

e_sym_decod1_opt_b      = zeros(nb_run,length(U),nb_model);
e_noise_decod1_opt_b    = zeros(nb_run,length(U),nb_model);

%
e_sym_decod1_e      = zeros(nb_run,length(U));
e_noise_decod1_e    = zeros(nb_run,length(U));
e_an_decod1_e       = zeros(nb_run,length(U));
e_denom_decod1_e    = zeros(nb_run,length(U));

e_sym_decod1_opt_e      = zeros(nb_run,length(U));
e_noise_decod1_opt_e    = zeros(nb_run,length(U));
e_an_decod1_opt_e       = zeros(nb_run,length(U));
e_denom_decod1_opt_e    = zeros(nb_run,length(U));

%
e_sym_decod2_e      = zeros(nb_run,length(U));
e_noise_decod2_e    = zeros(nb_run,length(U));
e_an_decod2_e       = zeros(nb_run,length(U));
e_denom_decod2_e    = zeros(nb_run,length(U));

e_sym_decod2_opt_e      = zeros(nb_run,length(U));
e_noise_decod2_opt_e    = zeros(nb_run,length(U));
e_an_decod2_opt_e       = zeros(nb_run,length(U));
e_denom_decod2_opt_e    = zeros(nb_run,length(U));

%
e_sym_decod5_e      = zeros(nb_run,length(U));
e_noise_decod5_e    = zeros(nb_run,length(U));
e_an_decod5_e       = zeros(nb_run,length(U));
e_denom_decod5_e    = zeros(nb_run,length(U));

e_sym_decod5_opt_e      = zeros(nb_run,length(U));
e_noise_decod5_opt_e    = zeros(nb_run,length(U));
e_an_decod5_opt_e       = zeros(nb_run,length(U));
e_denom_decod5_opt_e    = zeros(nb_run,length(U));

%
sinr1_b             = zeros(nb_run,length(U),nb_model);
sinr1_e             = zeros(nb_run,length(U));
sinr2_e             = zeros(nb_run,length(U));
sinr5_e             = zeros(nb_run,length(U));

sinr1_opt_b             = zeros(nb_run,length(U),nb_model);
sinr1_opt_e             = zeros(nb_run,length(U));
sinr2_opt_e             = zeros(nb_run,length(U));
sinr5_opt_e             = zeros(nb_run,length(U));

%
% capa1_b             = zeros(nb_run,length(U),nb_model);
% capa1_e             = zeros(nb_run,length(U));
% capa2_e             = zeros(nb_run,length(U));
% capa5_e             = zeros(nb_run,length(U));
% 
% capa1_opt_b             = zeros(nb_run,length(U),nb_model);
% capa1_opt_e             = zeros(nb_run,length(U));
% capa2_opt_e             = zeros(nb_run,length(U));
% capa5_opt_e             = zeros(nb_run,length(U));
% 
% %
% sr1             = zeros(nb_run,length(U),nb_model);
% sr2             = zeros(nb_run,length(U),nb_model);
% sr5             = zeros(nb_run,length(U));
% 
% sr_opt_1             = zeros(nb_run,length(U),nb_model);
% sr_opt_2             = zeros(nb_run,length(U));
% sr_opt_5             = zeros(nb_run,length(U));
tic

%% Mainloop
for iter = 1:nb_run
for bb =1:length(U)
for mm = 1:nb_model
msg_TX = randi( [0 1] , nb_bit(bb) , 1 ) ;                                         % Random bit data stream
sym_TX = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types, Nxnb_run

% Channel generation
% Rayleigh channel - independant subcarrier fading - unit power per
% subcarrier : CN(0,I_Q) - circulary complex random normal variable.
Hb_TX = channelRayleigh(Q, mu , sigma);  
He_TX = channelRayleigh(Q, mu , sigma);
Hb_RX = ctranspose(Hb_TX.');
He_RX = ctranspose(He_TX.');



% (De)-Spreading matrix creation
[matrix_spread,matrix_despread] = spreadingMatrix(Q,N(bb),U(bb));


%% Encoder

% Spreading + TR precoding
sym_spread      = matrix_spread*sym_TX;  % Qx1
sym_precoded    = Hb_TX*sym_spread; % Qx1, not weighted 

% AN generation
an = generateAN(Hb_RX,Q,U(bb),matrix_despread,energy(sym_precoded),"svd"); % Qx1, not weighted

    
alpha_to_opt = alpha_global(bb,mm).*ones(Q,1);                      % Coef alpha to optimize


% Non optimized
sym_precoded_TX     = sqrt(alpha_to_opt).*sym_precoded;     % weighted
an_TX               = sqrt(1-alpha_to_opt).*an;                       % weighted
e_an_TX             = energy(an_TX);
sym_transmitted     = sym_precoded_TX + an_TX;
e_sym_transmitted   = energy(sym_transmitted);


%% Decoder
decod1 = matrix_despread;                                   % despreading
decod2 = matrix_despread*Hb_RX*He_TX;                       % matched filter
decod5 = matrix_despread*He_TX;                             % Only He known by Eve


%% Optimization problem

% Problems creations
prob = optimproblem('ObjectiveSense', 'max');                               % Objective is to maximize f --> minimization of -f

% Variables to optimize
x = optimvar('x', Q,'LowerBound',0,'UpperBound',1,'Type','continuous');
   
%  Objective function
f   = @(x) energy(decod1*(diag(Hb_RX).*sqrt(x).*sym_precoded));             % Maximization of Bob SINR numerator (only SINR term depending on alpha) - sum over [Nx1]
fun = fcn2optimexpr(f,x,'OutputSize',[1,1]);
prob.Objective = fun;

% Nonlinear constraints
%sum_alpha = mean(x);                                                       % Constant total radiated energy
f2 = @(x) energy(decod1*(diag(Hb_RX).*sqrt(ones(Q,1)-x).*an));              % W lies in null space of H_bob - sum over [Nx1]
f3 = @(x) energy(sqrt(x).*sym_precoded + sqrt(ones(Q,1) - x).*an);          % Total radiated energy constant - sum over [Qx1] --> 1/U
f4 = @(x) energy(sqrt(ones(Q,1) - x).*an);                                  % AN energy constant - sum over [Qx1] ---> (1-alpha)/U

cstr_ortho  = fcn2optimexpr(f2,x,'OutputSize',[1,1]);                       % Orthogonality constraints
cstr_e_tot  = fcn2optimexpr(f3,x,'OutputSize',[1,1]);                       % Total energy constraint
cstr_e_an   = fcn2optimexpr(f4,x,'OutputSize',[1,1]);                       % AN energy constraint

%prob.Constraints.cons_alpha = sum_alpha == alpha_global(nb_bor);
prob.Constraints.cons_energy_AN     = cstr_e_an     == e_an_TX;%(1-alpha_global)/U(bb);
prob.Constraints.cons_energy_total  = cstr_e_tot    == e_sym_transmitted; %1/U(bb);
prob.Constraints.cons_orthog        = cstr_ortho    <= 1e-7;

% Initial vector
x0.x = alpha_to_opt;


% Problem options
options = optimoptions('fmincon','Algorithm', 'interior-point','MaxIterations',2500,'MaxFunctionEvaluations',300000);
% options.Display = 'iter-detailed';                                         % Display the problem status at each iteration
options.StepTolerance = 1e-6;
options.FunctionTolerance = 1e-7;
options.ConstraintTolerance = 1e-7;
%options.UseParallel;

% show(prob)                                                                % Show the problem, i.e., function to optimization, constraints, initial point
% Solve the problem
[sol] = solve(prob,x0,'Options',options);
% disp(sol.x)                                                               % Display the solution
alpha_opt(:,iter,bb,mm) = sol.x;                                               % Optimal coef alpha for a given data block, channel and BOR value


% END OF OPTIMIZATION PROBLEM
%% 



% Optimized transmitted signals
sym_precoded_opt_TX     = sqrt(sol.x).*sym_precoded;     % weighted
an_opt_TX               = sqrt(1-sol.x).*an;                       % weighted
e_an_opt_TX(iter,bb,mm) = energy(an_opt_TX);
sym_opt_transmitted     = sym_precoded_opt_TX + an_opt_TX;


%% Receiver
% Whole stream - non optimized 
sym_RX_b = Hb_RX*sym_transmitted;
sym_RX_e = He_RX*sym_transmitted;

% Whole stream - optimized 
sym_opt_RX_b = Hb_RX*sym_opt_transmitted;
sym_opt_RX_e = He_RX*sym_opt_transmitted;


% Useful symbol - non optimized
sym_b = Hb_RX*sym_precoded_TX; % Qx1
sym_e = He_RX*sym_precoded_TX;

% Useful symbol - optimized
sym_opt_b = Hb_RX*sym_precoded_opt_TX; % Qx1
sym_opt_e = He_RX*sym_precoded_opt_TX;


% Noise symbol - non optimized 
[noise_b, ~ ] = addNoise(sym_b , snr_b, energy(sym_precoded_TX+an_TX));     %addNoise(sym_b , snr_b, energy(sym_RX_b)); 
[noise_e, e_noise_e ] = addNoise(sym_e , snr_e, energy(sym_precoded_TX+an_TX)); % addNoise(sym_e , snr_e, energy(sym_RX_e));   

% Noise symbol - optimized 
[noise_opt_b, ~ ] = addNoise(sym_opt_b , snr_b, energy(sym_precoded_opt_TX+an_opt_TX));     %addNoise(sym_b , snr_b, energy(sym_RX_b)); 
[noise_opt_e, e_noise_opt_e ] = addNoise(sym_opt_e , snr_e, energy(sym_precoded_opt_TX+an_opt_TX)); % addNoise(sym_e , snr_e, energy(sym_RX_e));   


% AN symbol - non optimized
an_e = He_RX*an_TX; % Only @Eve since no AN effect after decod1 @Bob

% AN symbol - optimized
an_opt_e = He_RX*an_opt_TX; % Only @Eve since no AN effect after decod1 @Bob

%% Decoding
sym_decod1_b = decod1*sym_b;
noise_decod1_b = decod1*noise_b;
sym_decod1_opt_b = decod1*sym_opt_b;
noise_decod1_opt_b = decod1*noise_opt_b;

% @ Bob - non optimized
e_sym_decod1_b(iter,bb,mm)     = energy(sym_decod1_b);
e_noise_decod1_b(iter,bb,mm)   = energy(noise_decod1_b);

% @ Bob - optimized
e_sym_decod1_opt_b(iter,bb,mm)     = energy(sym_decod1_opt_b);
e_noise_decod1_opt_b(iter,bb,mm)   = energy(noise_decod1_opt_b);

% instantaneous SINRs - non optimized 
sinr1_b(iter,bb,mm) = e_sym_decod1_b(iter,bb,mm)/e_noise_decod1_b(iter,bb,mm);
sinr1_opt_b(iter,bb,mm) = e_sym_decod1_opt_b(iter,bb,mm)/e_noise_decod1_opt_b(iter,bb,mm) ;
% instantaneous SINRs - optimized 
% If optimization doesn't succeed: bob SINR = non optimized SINR
if e_sym_decod1_opt_b(iter,bb)/e_noise_decod1_opt_b(iter,bb,mm) >= sinr1_b(iter,bb,mm)
    sinr1_opt_b(iter,bb,mm) = e_sym_decod1_opt_b(iter,bb,mm)/e_noise_decod1_opt_b(iter,bb,mm) ;
else
    sinr1_opt_b(iter,bb,mm) = sinr1_b(iter,bb,mm);
end


if mm ==1 % Eve's decoder: Same as Bob
    
% @ Eve : decod1 RX sequence - Non optimized
sym_decod1_e = decod1*sym_e;
noise_decod1_e = decod1*noise_e;
an_decod1_e = decod1*an_e;

% @ Eve : decod 1 RX sequence - Optimized
sym_decod1_opt_e = decod1*sym_opt_e;
noise_decod1_opt_e = decod1*noise_opt_e;
an_decod1_opt_e = decod1*an_opt_e;

% @ Eve : decod 1 - non optimized
e_sym_decod1_e(iter,bb)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,bb)   = energy(noise_decod1_e);
e_an_decod1_e(iter,bb)      = energy(an_decod1_e);
e_denom_decod1_e(iter,bb)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 1 - optimized
e_sym_decod1_opt_e(iter,bb)     = energy(sym_decod1_opt_e);
e_noise_decod1_opt_e(iter,bb)   = energy(noise_decod1_opt_e);
e_an_decod1_opt_e(iter,bb)      = energy(an_decod1_opt_e);
e_denom_decod1_opt_e(iter,bb)   = energy(noise_decod1_opt_e + an_decod1_opt_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve: decod 1 - sinr non optimized
sinr1_e(iter,bb) = e_sym_decod1_e(iter,bb)/e_denom_decod1_e(iter,bb);

% @ Eve: decod 1 - sinr optimized
sinr1_opt_e(iter,bb) = e_sym_decod1_opt_e(iter,bb)/e_denom_decod1_opt_e(iter,bb);


elseif mm == 2 % Eve's decoder : Matched filtering
    
% @ Eve : decod2 RX sequence - Non optimized
sym_decod2_e = decod2*sym_e;
noise_decod2_e = decod2*noise_e;
an_decod2_e = decod2*an_e;

% @ Eve : decod2 RX sequence - optimized
sym_decod2_opt_e = decod2*sym_opt_e;
noise_decod2_opt_e = decod2*noise_opt_e;
an_decod2_opt_e = decod2*an_opt_e;

% @ Eve : decod 2 - non optimized
e_sym_decod2_e(iter,bb)     = energy(sym_decod2_e);
e_noise_decod2_e(iter,bb)   = energy(noise_decod2_e);
e_an_decod2_e(iter,bb)      = energy(an_decod2_e);
e_denom_decod2_e(iter,bb)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 2 - optimized
e_sym_decod2_opt_e(iter,bb)     = energy(sym_decod2_opt_e);
e_noise_decod2_opt_e(iter,bb)   = energy(noise_decod2_opt_e);
e_an_decod2_opt_e(iter,bb)      = energy(an_decod2_opt_e);
e_denom_decod2_opt_e(iter,bb)   = energy(noise_decod2_opt_e + an_decod2_opt_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve: decod 2 - sinr non optimized
sinr2_e(iter,bb) = e_sym_decod2_e(iter,bb)/e_denom_decod2_e(iter,bb);

% @ Eve: decod 2 - sinr optimized
sinr2_opt_e(iter,bb) = e_sym_decod2_opt_e(iter,bb)/e_denom_decod2_opt_e(iter,bb);

    
else % Eve's decoder : Own channel knowledge
    
% @ Eve : decod5 RX sequence - Non optimized
sym_decod5_e = decod5*sym_e;
noise_decod5_e = decod5*noise_e;
an_decod5_e = decod5*an_e;

% @ Eve : decod5 RX sequence - Optimized
sym_decod5_opt_e = decod5*sym_opt_e;
noise_decod5_opt_e = decod5*noise_opt_e;
an_decod5_opt_e = decod5*an_opt_e;

% @ Eve : decod 5 - non optimized
e_sym_decod5_e(iter,bb)     = energy(sym_decod5_e);
e_noise_decod5_e(iter,bb)   = energy(noise_decod5_e);
e_an_decod5_e(iter,bb)      = energy(an_decod5_e);
e_denom_decod5_e(iter,bb)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 5 - optimized
e_sym_decod5_opt_e(iter,bb)     = energy(sym_decod5_opt_e);
e_noise_decod5_opt_e(iter,bb)   = energy(noise_decod5_opt_e);
e_an_decod5_opt_e(iter,bb)      = energy(an_decod5_opt_e);
e_denom_decod5_opt_e(iter,bb)   = energy(noise_decod5_opt_e + an_decod5_opt_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve: decod 5 - sinr non optimized
sinr5_e(iter,bb) = e_sym_decod5_e(iter,bb)/e_denom_decod5_e(iter,bb);

% @ Eve: decod 5 - sinr optimized
sinr5_opt_e(iter,bb) = e_sym_decod5_opt_e(iter,bb)/e_denom_decod5_opt_e(iter,bb);

end


end
end
waitbar(iter / nb_run)
end
toc
%% POST PROCESSING
capa1_b = capacity(sinr1_b);
capa1_opt_b = capacity(sinr1_opt_b);

capa1_e = capacity(sinr1_e);
capa1_opt_e = capacity(sinr1_opt_e);

capa2_e = capacity(sinr2_e);
capa2_opt_e = capacity(sinr2_opt_e);

capa5_e = capacity(sinr5_e);
capa5_opt_e = capacity(sinr5_opt_e);


% instantaneous Secrecy capacity - non optimized
sr1 = secrecyCapacity(sinr1_b(:,:,1),sinr1_e);
sr2 = secrecyCapacity(sinr1_b(:,:,2),sinr2_e);
sr5 = secrecyCapacity(sinr1_b(:,:,3),sinr5_e);

% instantaneous Secrecy capacity - optimized
sr_opt_1 = secrecyCapacity(sinr1_opt_b(:,:,1),sinr1_opt_e);
sr_opt_2 = secrecyCapacity(sinr1_opt_b(:,:,2),sinr2_opt_e);
sr_opt_5 = secrecyCapacity(sinr1_opt_b(:,:,3),sinr5_opt_e);


% SR gain
sr_gain_1 = sr_opt_1 - sr1;
sr_gain_2 = sr_opt_2 - sr2;
sr_gain_5 = sr_opt_5 - sr5;

sr_gain_avg_1 = squeeze(mean(sr_gain_1,1));
sr_gain_avg_2 = squeeze(mean(sr_gain_2,1));
sr_gain_avg_5 = squeeze(mean(sr_gain_5,1));

% ergodic Secrecy capacity - non optimzed
sr1_avg = squeeze(mean(sr1,1));
sr2_avg = squeeze(mean(sr2,1));
sr5_avg = squeeze(mean(sr5,1));

% ergodic Secrecy capacity - optimzed
sr1_opt_avg = squeeze(mean(sr_opt_1,1));
sr2_opt_avg = squeeze(mean(sr_opt_2,1));
sr5_opt_avg = squeeze(mean(sr_opt_5,1));

figure;
plot(U,sr1_avg,'Marker','o','color','b'); hold on;
plot(U,sr1_opt_avg,'Marker','o','color','r'); hold on;
plot(U,sr2_avg,'Marker','square','color','m'); hold on;
plot(U,sr2_opt_avg,'Marker','square','color','g'); hold on;
plot(U,sr5_avg,'Marker','diamond','color','c'); hold on;
plot(U,sr5_opt_avg,'Marker','diamond','color','y'); hold on;
box on; grid on;
xlabel('BOR')
ylabel('Secrecy rate (bit/channel use)')
legend('Same decoder - non optimized' , 'Same decoder - optimized' , ...
    'Matched filtering - non optimized' , 'Matched filtering - optimized', ...
    'Own channel - non optimized', 'Own channel - optimized', 'location', 'best')


figure;
plot(U,sr_gain_avg_1,'Marker','o','color','b'); hold on;
plot(U,sr_gain_avg_2,'Marker','square','color','r'); hold on;
plot(U,sr_gain_avg_5,'Marker','diamond','color','g'); hold on;
box on; grid on;
xlabel('BOR')
ylabel('Secrecy rate (bit/channel use)')
legend('Secrecy gain - Same decoder' , 'Secrecy gain - Matched filter',...
    'Secrecy gain - Own channel', 'location', 'best')


close(h);


