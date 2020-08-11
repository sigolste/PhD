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
nb_run = 100;              % number of experiments
alpha_global = 0.33; 
% Communication parameters
Q = 8;
U = [4];
N = Q./U;

M = 4;
k = log2(M);
nb_bit = k.*N;

% AWGN parameters
EbN0_b = 10; % energy per bit over noise psd @Bob - dB
EbN0_e = [10]; % energy per bit over noise psd @Eve - dB
snr_b  = EbN0_b + 10*log10(k);  % SNR @Bob
snr_e  = EbN0_e + 10*log10(k);  % SNR @Eve

% Channel parameters 
mu = 0;         % Channel mean
sigma = 1;      % Channel variance


alpha_opt           = zeros(Q,nb_run,length(U));
e_an_TX             = zeros(nb_run,length(U));
e_an_opt_TX         = zeros(nb_run,length(U));

e_sym_decod1_b      = zeros(nb_run,length(U));
e_noise_decod1_b    = zeros(nb_run,length(U));

e_sym_decod1_e      = zeros(nb_run,length(U));
e_noise_decod1_e    = zeros(nb_run,length(U));
e_an_decod1_e       = zeros(nb_run,length(U));
e_denom_decod1_e    = zeros(nb_run,length(U));

e_sym_decod2_e      = zeros(nb_run,length(U));
e_noise_decod2_e    = zeros(nb_run,length(U));
e_an_decod2_e       = zeros(nb_run,length(U));
e_denom_decod2_e    = zeros(nb_run,length(U));

e_sym_decod5_e      = zeros(nb_run,length(U));
e_noise_decod5_e    = zeros(nb_run,length(U));
e_an_decod5_e       = zeros(nb_run,length(U));
e_denom_decod5_e    = zeros(nb_run,length(U));

sinr1_b             = zeros(nb_run,length(U));
sinr1_e             = zeros(nb_run,length(U));
sinr2_e             = zeros(nb_run,length(U));
sinr5_e             = zeros(nb_run,length(U));

sr1             = zeros(nb_run,length(U));
sr2             = zeros(nb_run,length(U));
sr5             = zeros(nb_run,length(U));


%% Mainloop
for iter = 1:nb_run
for bb =1:length(U)
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

    
alpha_to_opt = alpha_global(bb)*ones(Q,1);                      % Coef alpha to optimize




%% Optimization problem

% Problems creations
prob = optimproblem('ObjectiveSense', 'max');                               % Objective is to maximize f --> minimization of -f

% Variables to optimize
x = optimvar('x', Q,'LowerBound',0,'UpperBound',1,'Type','continuous');


   
%  Objective function
f   = @(x) sum(abs(decod1*(diag(Hb_RX).*(sqrt(x).*sym_precoded))).^2);      % Maximization of Bob SINR numerator (only SINR term depending on alpha) - sum oer [Nx1]
fun = fcn2optimexpr(f,x,'OutputSize',[1,1]);
prob.Objective = fun;

% Nonlinear constraints
%sum_alpha = mean(x);                                                       % Constant total radiated energy
f2 = @(x) sum(abs(decod1*(diag(Hb_RX).*sqrt(ones(Q,1)-x).*an)).^2);         % W lies in null space of H_bob - sum over [Nx1]
f3 = @(x) sum(abs(sqrt(x).*sym_precoded + sqrt(ones(Q,1) - x).*an).^2)/Q;   % Total radiated energy constant - sum over [Qx1]
f4 = @(x) sum(abs(sqrt(ones(Q,1) - x).*an).^2)/Q;                           % AN energy constant - sum over [Qx1] 

cstr_ortho  = fcn2optimexpr(f2,x,'OutputSize',[1,1]);                       % Orthogonality constraints
cstr_e_tot  = fcn2optimexpr(f3,x,'OutputSize',[1,1]);                       % Total energy constraint
cstr_e_an   = fcn2optimexpr(f4,x,'OutputSize',[1,1]);                       % AN energy constraint

%prob.Constraints.cons_alpha = sum_alpha == alpha_global(nb_bor);
prob.Constraints.cons_energy_AN     = cstr_e_an     == sum(abs(sqrt(ones(Q,1) - alpha_to_opt).*an).^2)/Q;
prob.Constraints.cons_energy_total  = cstr_e_tot    == e_sym_transmitted;
prob.Constraints.cons_orthog        = cstr_ortho    <= 1e-15;

% Initial vector
x0.x = alpha_to_opt;


% Problem options
options = optimoptions('fmincon','Algorithm', 'interior-point','MaxIterations',1300,'MaxFunctionEvaluations',300000);
%options.Display = 'iter-detailed';                                         % Display the problem status at each iteration
options.StepTolerance       = 1e-7;
options.FunctionTolerance   = 1e-8;
options.ConstraintTolerance = 1e-8;
% show(prob)                                                                % Show the problem, i.e., function to optimization, constraints, initial point
% Solve the problem
[sol] = solve(prob,x0,'Options',options);
% disp(sol.x)                                                               % Display the solution
alpha_opt(:,iter,bb) = sol.x;                                               % Optimal coef alpha for a given data block, channel and BOR value
%% END of optimizaton 

% Non optimized
sym_precoded_TX     = sqrt(alpha_to_opt).*sym_precoded;     % weighted
an_TX               = sqrt(1-alpha_to_opt).*an;                       % weighted
e_an_TX(iter,bb)    = energy(an_TX);
sym_transmitted     = sym_precoded_TX + an_TX;

% Optimized 
sym_precoded_opt_TX     = sqrt(alpha_opt(:,iter,bb)).*sym_precoded;     % weighted
an_opt_TX               = sqrt(1-alpha_opt(:,iter,bb)).*an;                       % weighted
e_an_opt_TX(iter,bb)    = energy(an_opt_TX);
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

% AN symbol
an_e = He_RX*an_TX; % Only @Eve since no AN effect after decod1 @Bob

%% Decoder
decod1 = matrix_despread;                                   % despreading
decod2 = matrix_despread*Hb_RX*He_TX;                       % matched filter
decod5 = matrix_despread*He_TX;                             % Only He known by Eve
 
sym_decod1_b = decod1*sym_b;
sym_decod1_e = decod1*sym_e;
sym_decod2_e = decod2*sym_e;
sym_decod5_e = decod5*sym_e;

noise_decod1_b = decod1*noise_b;
noise_decod1_e = decod1*noise_e;
noise_decod2_e = decod2*noise_e;
noise_decod5_e = decod5*noise_e;

an_decod1_e = decod1*an_e;
an_decod2_e = decod2*an_e;
an_decod5_e = decod5*an_e;

%% Energy of the different RX components 
% @ Bob
e_sym_decod1_b(iter,bb)     = energy(sym_decod1_b);
e_noise_decod1_b(iter,bb)   = energy(noise_decod1_b);

% @ Eve : decod 1
e_sym_decod1_e(iter,bb)     = energy(sym_decod1_e);
e_noise_decod1_e(iter,bb)   = energy(noise_decod1_e);
e_an_decod1_e(iter,bb)      = energy(an_decod1_e);
e_denom_decod1_e(iter,bb)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve

% @ Eve : decod 2
e_sym_decod2_e(iter,bb)     = energy(sym_decod2_e);
e_noise_decod2_e(iter,bb)   = energy(noise_decod2_e);
e_an_decod2_e(iter,bb)      = energy(an_decod2_e);
e_denom_decod2_e(iter,bb)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve

% @ Eve : decod 5
e_sym_decod5_e(iter,bb)     = energy(sym_decod5_e);
e_noise_decod5_e(iter,bb)   = energy(noise_decod5_e);
e_an_decod5_e(iter,bb)      = energy(an_decod5_e);
e_denom_decod5_e(iter,bb)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 5 @Eve

% instantaneous SINRs
sinr1_b(iter,bb) = e_sym_decod1_b(iter,bb)/e_noise_decod1_b(iter,bb);
sinr1_e(iter,bb) = e_sym_decod1_e(iter,bb)/e_denom_decod1_e(iter,bb);
sinr2_e(iter,bb) = e_sym_decod2_e(iter,bb)/e_denom_decod2_e(iter,bb);
sinr5_e(iter,bb) = e_sym_decod5_e(iter,bb)/e_denom_decod5_e(iter,bb);


% instantaneous Secrecy capacity
sr1(iter,bb) = secrecyCapacity(sinr1_b(iter,bb),sinr1_e(iter,bb));
sr2(iter,bb) = secrecyCapacity(sinr1_b(iter,bb),sinr2_e(iter,bb));
sr5(iter,bb) = secrecyCapacity(sinr1_b(iter,bb),sinr5_e(iter,bb));


end
waitbar(iter / nb_run)
end



a = 1;



%% TRANSMISSION PART
for nb_bor = 1:length(BOR)
for cc = 1:nb_channels
msg_TX      = randi( [0 1] , nb_bit(nb_bor) , 1 ) ;                                 % Random bit data stream
sym_TX      = qammod(msg_TX,M,'gray','UnitAveragePower',true, 'InputType', 'bit');  % QAM modulation, can be changed to different modulation types
sym_TX      = reshape(sym_TX, [], nb_block) ;                                       % OFDM block per OFDM block
                                      


                                      
spread_matrix   = zeros( nb_subcar , nb_symb_per_block(nb_bor) );
spread_code     = reshape( 2.*randi([0 1], 1, nb_subcar).' - 1 , [] , BOR(nb_bor) )/sqrt(BOR(nb_bor));   % Diagonal elements of the spreading matrix. BPSK spreading code


for i = 1:size(spread_code,2)
    spread_matrix(i*nb_symb_per_block(nb_bor)-nb_symb_per_block(nb_bor)+1:i*end/BOR(nb_bor),:) = diag(spread_code(:,i));
end
despread_matrix = ctranspose(spread_matrix); 

%% 3. Symbol stream after spreading
sym_spread = spread_matrix*sym_TX;



H_bob = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,"decorrelated");
H_eve = TR_security_generate_rayleigh_channel(delta_t,nb_taps,nb_subcar,"decorrelated");
H_bob_TX    = ctranspose(H_bob).';
H_bob_RX    = H_bob;
H_eve_RX    = H_eve;



sym_useful_TX = H_bob_TX.*sym_spread;      % ATTENTION sqrt(alpha) needed  
W = TR_security_generate_AN(H_bob,nb_subcar,BOR(nb_bor),despread_matrix, "underdetermined"); % Resolution of the underdetermined system


tmp = sum(abs(W).^2)/nb_subcar;
energy_sym_useful_TX = mean(sum(abs(H_bob_TX.*sym_spread).^2))/nb_subcar;
W = W/sqrt(tmp)*sqrt(energy_sym_useful_TX);
energy_W = sum(abs(W).^2)/nb_subcar;




    
    
alpha_to_opt = alpha_global(nb_bor)*ones(nb_subcar,1);                      % Coef alpha to optimize



sym_raditated_noiseless = sqrt(alpha_to_opt).*sym_useful_TX + sqrt(ones(nb_subcar,1) - alpha_to_opt).*W;
test_noiseless(cc,nb_bor) = 1/nb_subcar*sum(abs(sym_raditated_noiseless).^2);

energy_sym_radiated_noiseless =  mean(sum(abs(sqrt(alpha_to_opt).*sym_useful_TX + sqrt(ones(nb_subcar,1) - alpha_to_opt).*W).^2))/nb_subcar;
%% RECEPTION BEFORE DESPREADING

% Total received symbol at B and E
sym_RX_bob = H_bob_RX.*sym_raditated_noiseless;
sym_RX_eve = H_eve_RX.*sym_raditated_noiseless;

test_energy_sym_RX_bob(cc,nb_bor) = 1/nb_subcar*sum(abs(sym_RX_bob).^2);
test_energy_sym_RX_eve(cc,nb_bor) = 1/nb_subcar*sum(abs(sym_RX_eve).^2);

% Noise addition on useful symbols + AN
[sym_RX_noisy_bob,energy_noise_bob] = TR_security_addNoise(sym_RX_bob, EbN0, M, var(sym_raditated_noiseless));   % Noisy received symbols (useful + AN) @ Bob
[sym_RX_noisy_eve,energy_noise_eve] = TR_security_addNoise(sym_RX_eve, EbN0, M, var(sym_raditated_noiseless));   % Noisy received symbols (useful + AN) @ Eve



% useful part of received symbol at B and E
sym_useful_RX_bob = H_bob_RX.*(sqrt(alpha_to_opt).*sym_useful_TX);
sym_useful_RX_eve = H_eve_RX.*(sqrt(alpha_to_opt).*sym_useful_TX);

% Artificial noise received symbol at E
sym_AN_RX_eve = H_eve_RX.*sqrt(ones(nb_subcar,1) - alpha_to_opt).*W;


tmp = sqrt(ones(nb_subcar,1) - alpha_to_opt).*W;
test_tmp(cc,nb_bor) = 1/nb_subcar*sum(abs(tmp).^2);


%% RECEPTION AFTER DESPREADING

% Useful symbols
sym_useful_RX_despread_bob = despread_matrix*(H_bob_RX.*(sqrt(alpha_to_opt).*sym_useful_TX)); % Despread useful symbols @ Bob
sym_useful_RX_despread_eve = despread_matrix*(H_eve_RX.*(sqrt(alpha_to_opt).*sym_useful_TX)); % Despread useful symbols @ Eve

% AN symbols (only @ Eve)
sym_AN_RX_despread_eve = despread_matrix*sym_AN_RX_eve; % Despread AN symbols @ Eve

test_energy_AN_eve(cc,nb_bor) = 1/nb_subcar*sum(abs(sym_AN_RX_eve).^2);
test_energy_AN_despread(cc,nb_bor) = 1/nb_subcar*sum(abs(sym_AN_RX_despread_eve).^2);
% useful symbols + AN
sym_RX_despread_bob = despread_matrix*sym_RX_noisy_bob; % Depsread symbols (useful + AN) @ Bob
sym_RX_despread_eve = despread_matrix*sym_RX_noisy_eve; % Depsread symbols (useful + AN) @ Eve

% noisy symbols
noise_despread_bob =  sym_RX_despread_bob - sym_useful_RX_despread_bob;
noise_despread_eve =  sym_RX_despread_eve - sym_AN_RX_despread_eve - sym_useful_RX_despread_eve;

%% Energy computation
energy_sym_useful_RX_despread_bob = mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_bob).^2)); 
energy_sym_useful_RX_despread_eve = mean(1/nb_subcar*sum(abs(sym_useful_RX_despread_eve).^2)); 
energy_noise_despread_bob         = mean(1/nb_subcar*sum(abs(noise_despread_bob).^2)); 
energy_denom_eve                  = mean(1/nb_subcar*sum(abs(sym_AN_RX_despread_eve + noise_despread_eve).^2));


%% SINR COMPUTATION
SINR_bob(cc,nb_bor) = energy_sym_useful_RX_despread_bob/energy_noise_despread_bob;%mean(noise_power_bob); %energy_noise_despread_bob;
SINR_eve(cc,nb_bor) = energy_sym_useful_RX_despread_eve/energy_denom_eve;

%% CAPA COMPUTATION

capa_bob(cc,nb_bor) = log2(1+SINR_bob(cc,nb_bor));
capa_eve(cc,nb_bor) = log2(1+SINR_eve(cc,nb_bor));


secrecy_capa(cc,nb_bor) = capa_bob(cc,nb_bor) - capa_eve(cc,nb_bor) ; 
% capa_bob_per_subcar(:,cc,nb_bor) = log2(1 +  1/nb_subcar*BOR(nb_bor)*sum(abs(despread_matrix*(H_bob_RX.*(sqrt(alpha_to_opt).*sym_useful_TX)))).^2 ./ sum(abs(noise_despread_bob).^2));    % Mean of capacity over all data blocks to have the mean capa per subcar
% capa_eve_per_subcar(:,cc,nb_bor) = log2(1 +  1/nb_subcar*BOR(nb_bor)*sum(abs(despread_matrix*(H_eve_RX.*(sqrt(alpha_to_opt).*sym_useful_TX)))).^2 ./ sum(abs(noise_despread_eve + despread_matrix*(H_eve_RX.*(sqrt(ones(nb_subcar,1) - alpha_to_opt).*W))).^2));
% 
% secrecy_capa_per_subcar(:,cc,nb_bor) = capa_bob_per_subcar(:,cc,nb_bor) - capa_eve_per_subcar(:,cc,nb_bor);

% 
% 

%% OPTIMIZATION PROBLEM

% Function to maximize

for ii=1:nb_block
    
% Problem creation
prob = optimproblem('ObjectiveSense', 'max');

% Variables to optimize
x = optimvar('x', nb_subcar,'LowerBound',0,'UpperBound',1,'Type','continuous');

   
%  Objective function
f = @(x) sum(abs(despread_matrix*(H_bob_RX.*(sqrt(x).*sym_useful_TX(:,ii)))).^2); % sum over subcarriers and for all block sent
%f = @(x) sum(abs(H_bob_RX.*sqrt(x).* H_bob_TX).^2); 
fun = fcn2optimexpr(f,x,'OutputSize',[1,1]);
prob.Objective = fun;

% Nonlinear constraints
%sum_alpha = mean(x);                                                       % Constant total radiated energy
f2 = @(x) sum(abs(despread_matrix*(H_bob_RX.*sqrt(ones(nb_subcar,1)-x).*W)).^2);            % W lies in null space of H_bob
f3 = @(x) sum(abs(sqrt(x).*sym_useful_TX + sqrt(ones(nb_subcar,1) - x).*W).^2)/nb_subcar;   % Total radiated energy constant   
f4 = @(x) sum(abs(sqrt(ones(nb_subcar,1) - x).*W).^2)/nb_subcar;                            % AN energy constant

orthogonality = fcn2optimexpr(f2,x,'OutputSize',[1,1]);
energy_total = fcn2optimexpr(f3,x,'OutputSize',[1,1]);
energy_AN = fcn2optimexpr(f4,x,'OutputSize',[1,1]);

%prob.Constraints.cons_alpha = sum_alpha == alpha_global(nb_bor);
prob.Constraints.cons_energy_AN = energy_AN == sum(abs(sqrt(ones(nb_subcar,1) - alpha_to_opt).*W).^2)/nb_subcar; 
prob.Constraints.cons_energy_total = energy_total == energy_sym_radiated_noiseless;
prob.Constraints.cons_orthog = orthogonality <= 1e-15;

% Initial vector
x0.x = alpha_to_opt;


% Problem options
options = optimoptions('fmincon','Algorithm', 'interior-point','MaxIterations',500,'MaxFunctionEvaluations',300000);
%options.Display = 'iter-detailed';
options.StepTolerance = 1e-6;
options.FunctionTolerance = 1e-9;
options.ConstraintTolerance = 1e-9;

% Solve the problem
[sol] = solve(prob,x0,'Options',options);
alpha_opt(:,ii,cc,nb_bor) = sol.x;                                          % Optimal coef alpha for a given data block, channel and BOR value

sym_raditated_noiseless_opt = sqrt(alpha_opt(:,ii,cc,nb_bor)).*sym_useful_TX + sqrt(ones(nb_subcar,1) - alpha_opt(:,ii,cc,nb_bor)).*W;

test_noiseless_opt(ii,cc,nb_bor) = 1/nb_subcar*sum(abs(sym_raditated_noiseless_opt).^2);
%% Eve signals when using optimized alpha coef
sym_RX_eve_opt                              = H_eve_RX.*( sym_raditated_noiseless_opt );

test_energy_sym_RX_eve_opt(ii,cc,nb_bor) = 1/nb_subcar*sum(abs(sym_RX_eve_opt).^2);


[sym_RX_noisy_eve_opt,energy_noise_eve_opt] = TR_security_addNoise( sym_RX_eve_opt , EbN0 , M , var(sym_raditated_noiseless_opt) );   

sym_AN_RX_eve_opt = H_eve_RX.*sqrt(ones(nb_subcar,1) - alpha_opt(:,ii,cc,nb_bor)).*W ;
sym_AN_RX_despread_eve_opt                  = despread_matrix*sym_AN_RX_eve_opt;

tmp_opt = sqrt(ones(nb_subcar,1) - alpha_opt(:,ii,cc,nb_bor)).*W;

test_tmp_opt(ii,cc,nb_bor) = 1/nb_subcar*sum(abs(tmp_opt).^2);
test_energy_AN_eve_opt(ii,cc,nb_bor) = 1/nb_subcar*sum(abs(sym_AN_RX_eve_opt).^2);
test_energy_AN_despread_opt(ii,cc,nb_bor) = 1/nb_subcar*sum(abs(sym_AN_RX_despread_eve_opt).^2);

sym_useful_RX_despread_eve_opt              = despread_matrix* (H_eve_RX.* sqrt( alpha_opt(:,ii,cc,nb_bor) )  .* sym_useful_TX(:,ii) ); 
sym_RX_despread_eve_opt                     = despread_matrix* sym_RX_noisy_eve_opt; 
noise_despread_eve_opt                      = sym_RX_despread_eve_opt - sym_AN_RX_despread_eve_opt - sym_useful_RX_despread_eve_opt;

energy_sym_useful_RX_despread_eve_opt       = sum( abs(sym_useful_RX_despread_eve_opt).^2)/nb_subcar ; 
energy_denom_eve_opt                        = sum( abs(sym_AN_RX_despread_eve_opt + noise_despread_eve_opt).^2)/nb_subcar ;


%% Bob signals when using optimized alpha coef
sym_RX_bob_opt                              = H_bob_RX.*( sym_raditated_noiseless_opt );

test_energy_sym_RX_bob_opt(ii,cc,nb_bor) = 1/nb_subcar*sum(abs(sym_RX_bob_opt).^2);



[sym_RX_noisy_bob_opt,energy_noise_bob_opt] = TR_security_addNoise( sym_RX_bob_opt , EbN0 , M , var(sym_raditated_noiseless_opt) );   
sym_useful_RX_despread_bob_opt              = despread_matrix * ( H_bob_RX .*  sqrt( alpha_opt(:,ii,cc,nb_bor)) .* sym_useful_TX(:,ii) ); 
sym_RX_despread_bob_opt                     = despread_matrix * sym_RX_noisy_bob_opt; 
noise_despread_bob_opt                      = sym_RX_despread_bob_opt - sym_useful_RX_despread_bob_opt;

energy_noise_despread_bob_opt               =sum( abs(noise_despread_bob_opt).^2)/nb_subcar; 
energy_sym_useful_RX_despread_bob_opt       =sum( abs(sym_useful_RX_despread_bob_opt).^2)/nb_subcar;


%% Capcity after optimization

capa_bob_opt(ii,cc,nb_bor) = log2(1+ energy_sym_useful_RX_despread_bob_opt/energy_noise_despread_bob_opt);
capa_eve_opt(ii,cc,nb_bor) = log2(1+ energy_sym_useful_RX_despread_eve_opt/energy_denom_eve_opt); 

% capa_bob_opt(ii,cc,nb_bor) = log2(1+ mean(1/nb_subcar*BOR(nb_bor)*sum(abs(despread_matrix*(H_bob_RX.*(sqrt(alpha_opt(:,ii,cc,nb_bor)).*sym_useful_TX(:,ii)))).^2)) ./ energy_noise_despread_bob);
% capa_eve_opt(ii,cc,nb_bor) = log2(1+ mean(1/nb_subcar*BOR(nb_bor)*sum(abs(despread_matrix*(H_eve_RX.*(sqrt(alpha_opt(:,ii,cc,nb_bor)).*sym_useful_TX(:,ii)))).^2)) ./  energy_denom_eve_opt); %mean(1/nb_subcar*BOR(nb_bor)*sum(abs( despread_matrix*(H_eve_RX.*(sqrt(ones(nb_subcar,1) - alpha_opt(:,ii,cc,nb_bor)).*W)) + noise_despread_eve).^2)));

secrecy_capa_opt(ii,cc,nb_bor) = capa_bob_opt(ii,cc,nb_bor)-capa_eve_opt(ii,cc,nb_bor);

clc;
fprintf(syyydlib.progress_bar(ii,nb_block,'Data block number:'))
fprintf(syyydlib.progress_bar(cc,nb_channels,'Channel Status'))
fprintf(syyydlib.progress_bar(nb_bor,length(BOR),'Back Of Rate'))

end
end
end


%% Mean over all blocks and all channels:
capa_bob_ergo_opt   = squeeze(mean(capa_bob_opt,2));
capa_eve_ergo_opt   = squeeze(mean(capa_eve_opt,2));
capa_bob_ergo       = squeeze(mean(capa_bob,1));
capa_eve_ergo       = squeeze(mean(capa_eve,1));

secrecy_capa_ergo       = capa_bob_ergo - capa_eve_ergo;
secrecy_capa_ergo_opt   = capa_bob_ergo_opt - capa_eve_ergo_opt;


%% TEST ENERGIES
test_noiseless_opt = squeeze(mean(test_noiseless_opt,2));
test_noiseless = squeeze(mean(test_noiseless,1))';

test_bob = squeeze(mean(test_energy_sym_RX_bob,1))';
test_eve = squeeze(mean(test_energy_sym_RX_eve,1))';
test_bob_opt = squeeze(mean(test_energy_sym_RX_bob_opt,2));
test_eve_opt = squeeze(mean(test_energy_sym_RX_eve_opt,2));

test_AN = squeeze(mean(test_energy_AN_despread,1))';
test_AN_opt = squeeze(mean(test_energy_AN_despread_opt,2));

test_energy_AN = mean(test_energy_AN_eve,1)';
test_energy_AN_opt = squeeze(mean(test_energy_AN_eve_opt,2));


test_energy_tmp_opt = squeeze(mean(test_tmp_opt,2));
test_energy_tmp = mean(test_tmp,1)';
%% Plot section
figure;
plot(BOR,capa_bob_ergo,'b'); hold on; plot(BOR,capa_bob_ergo_opt,'r')
hold on;
plot(BOR,capa_eve_ergo,'c'); hold on; plot(BOR,capa_eve_ergo_opt,'y')
legend('Bob no opt', 'Bob opt', 'Eve no opt', 'Eve opt')
ylabel('Capacity (bit/s/Hz)')
xlabel('BOR')


figure;
plot(BOR,secrecy_capa_ergo,'b'); hold on; plot(BOR,secrecy_capa_ergo_opt,'r')
legend('No opt', 'Opt')
ylabel('Secrecy Rate (bit/s/Hz)')
xlabel('BOR')

figure;
yyaxis left
bar(alpha_to_opt,'r'); hold on;
bar(squeeze(alpha_opt(:,1,end)),'m'); hold on;

yyaxis right
plot(abs(H_bob_TX.*H_bob_RX).^2,'b');
