

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







% Non optimized
% sym_decod1_b = decod1*sym_b;
% sym_decod1_e = decod1*sym_e;
% sym_decod2_e = decod2*sym_e;
% sym_decod5_e = decod5*sym_e;
% 
% noise_decod1_b = decod1*noise_b;
% noise_decod1_e = decod1*noise_e;
% noise_decod2_e = decod2*noise_e;
% noise_decod5_e = decod5*noise_e;
% 
% an_decod1_e = decod1*an_e;
% an_decod2_e = decod2*an_e;
% an_decod5_e = decod5*an_e;
% 
% % Optimized
% sym_decod1_opt_b = decod1*sym_opt_b;
% sym_decod1_opt_e = decod1*sym_opt_e;
% sym_decod2_opt_e = decod2*sym_opt_e;
% sym_decod5_opt_e = decod5*sym_opt_e;
% 
% noise_decod1_opt_b = decod1*noise_opt_b;
% noise_decod1_opt_e = decod1*noise_opt_e;
% noise_decod2_opt_e = decod2*noise_opt_e;
% noise_decod5_opt_e = decod5*noise_opt_e;
% 
% an_decod1_opt_e = decod1*an_opt_e;
% an_decod2_opt_e = decod2*an_opt_e;
% an_decod5_opt_e = decod5*an_opt_e;
% 
% 
% %% Energy of the different RX components 
% % @ Bob - non optimized
% e_sym_decod1_b(iter,bb)     = energy(sym_decod1_b);
% e_noise_decod1_b(iter,bb)   = energy(noise_decod1_b);
% 
% % @ Bob - optimized
% e_sym_decod1_opt_b(iter,bb)     = energy(sym_decod1_opt_b);
% e_noise_decod1_opt_b(iter,bb)   = energy(noise_decod1_opt_b);
% 
% 
% % @ Eve : decod 1 - non optimized
% e_sym_decod1_e(iter,bb)     = energy(sym_decod1_e);
% e_noise_decod1_e(iter,bb)   = energy(noise_decod1_e);
% e_an_decod1_e(iter,bb)      = energy(an_decod1_e);
% e_denom_decod1_e(iter,bb)   = energy(noise_decod1_e + an_decod1_e);        % energy of the sinr denominator for decoder 1 @Eve
% 
% % @ Eve : decod 1 - optimized
% e_sym_decod1_opt_e(iter,bb)     = energy(sym_decod1_opt_e);
% e_noise_decod1_opt_e(iter,bb)   = energy(noise_decod1_opt_e);
% e_an_decod1_opt_e(iter,bb)      = energy(an_decod1_opt_e);
% e_denom_decod1_opt_e(iter,bb)   = energy(noise_decod1_opt_e + an_decod1_opt_e);        % energy of the sinr denominator for decoder 1 @Eve
% 
% 
% 
% % @ Eve : decod 2 - non optimized
% e_sym_decod2_e(iter,bb)     = energy(sym_decod2_e);
% e_noise_decod2_e(iter,bb)   = energy(noise_decod2_e);
% e_an_decod2_e(iter,bb)      = energy(an_decod2_e);
% e_denom_decod2_e(iter,bb)   = energy(noise_decod2_e + an_decod2_e);        % energy of the sinr denominator for decoder 2 @Eve
% 
% % @ Eve : decod 2 - optimized
% e_sym_decod2_opt_e(iter,bb)     = energy(sym_decod2_opt_e);
% e_noise_decod2_opt_e(iter,bb)   = energy(noise_decod2_opt_e);
% e_an_decod2_opt_e(iter,bb)      = energy(an_decod2_opt_e);
% e_denom_decod2_opt_e(iter,bb)   = energy(noise_decod2_opt_e + an_decod2_opt_e);        % energy of the sinr denominator for decoder 2 @Eve
% 
% 
% 
% % @ Eve : decod 5 - non optimized
% e_sym_decod5_e(iter,bb)     = energy(sym_decod5_e);
% e_noise_decod5_e(iter,bb)   = energy(noise_decod5_e);
% e_an_decod5_e(iter,bb)      = energy(an_decod5_e);
% e_denom_decod5_e(iter,bb)   = energy(noise_decod5_e + an_decod5_e);        % energy of the sinr denominator for decoder 5 @Eve
% 
% % @ Eve : decod 5 - optimized
% e_sym_decod5_opt_e(iter,bb)     = energy(sym_decod5_opt_e);
% e_noise_decod5_opt_e(iter,bb)   = energy(noise_decod5_opt_e);
% e_an_decod5_opt_e(iter,bb)      = energy(an_decod5_opt_e);
% e_denom_decod5_opt_e(iter,bb)   = energy(noise_decod5_opt_e + an_decod5_opt_e);        % energy of the sinr denominator for decoder 5 @Eve
% 
% 
% 
% % instantaneous SINRs - non optimized 
% sinr1_b(iter,bb) = e_sym_decod1_b(iter,bb)/e_noise_decod1_b(iter,bb);
% sinr1_e(iter,bb) = e_sym_decod1_e(iter,bb)/e_denom_decod1_e(iter,bb);
% sinr2_e(iter,bb) = e_sym_decod2_e(iter,bb)/e_denom_decod2_e(iter,bb);
% sinr5_e(iter,bb) = e_sym_decod5_e(iter,bb)/e_denom_decod5_e(iter,bb);
% 
% % instantaneous SINRs - optimized 
% 
% % If optimization does not succeed, sinr bob optimized = sinr bob not optimized
% if e_sym_decod1_opt_b(iter,bb)/e_noise_decod1_opt_b(iter,bb) >= sinr1_b(iter,bb)
%     sinr1_opt_b(iter,bb) = e_sym_decod1_opt_b(iter,bb)/e_noise_decod1_opt_b(iter,bb) ;
% else
%     sinr1_opt_b(iter,bb) = sinr1_b(iter,bb);
% end
%     
% 
% sinr1_opt_e(iter,bb) = e_sym_decod1_opt_e(iter,bb)/e_denom_decod1_opt_e(iter,bb);
% sinr2_opt_e(iter,bb) = e_sym_decod2_opt_e(iter,bb)/e_denom_decod2_opt_e(iter,bb);
% sinr5_opt_e(iter,bb) = e_sym_decod5_opt_e(iter,bb)/e_denom_decod5_opt_e(iter,bb);


% % instantaneous capacity - non optimized 
% capa1_b(iter,bb) = capacity(sinr1_b(iter,bb));
% capa1_e(iter,bb) = capacity(sinr1_e(iter,bb));
% capa2_e(iter,bb) = capacity(sinr2_e(iter,bb));
% capa5_e(iter,bb) = capacity(sinr5_e(iter,bb));
% 
% % instantaneous capacity - optimized 
% capa1_opt_b(iter,bb) = capacity(sinr1_opt_b(iter,bb));
% capa1_opt_e(iter,bb) = capacity(sinr1_opt_e(iter,bb));
% capa2_opt_e(iter,bb) = capacity(sinr2_opt_e(iter,bb));
% capa5_opt_e(iter,bb) = capacity(sinr5_opt_e(iter,bb));
% 
% 
% % instantaneous Secrecy capacity - non optimized
% sr1(iter,bb) = secrecyCapacity(sinr1_b(iter,bb),sinr1_e(iter,bb));
% sr2(iter,bb) = secrecyCapacity(sinr1_b(iter,bb),sinr2_e(iter,bb));
% sr5(iter,bb) = secrecyCapacity(sinr1_b(iter,bb),sinr5_e(iter,bb));
% 
% % instantaneous Secrecy capacity - optimized
% sr_opt_1(iter,bb) = secrecyCapacity(sinr1_opt_b(iter,bb),sinr1_opt_e(iter,bb));
% sr_opt_2(iter,bb) = secrecyCapacity(sinr1_opt_b(iter,bb),sinr2_opt_e(iter,bb));
% sr_opt_5(iter,bb) = secrecyCapacity(sinr1_opt_b(iter,bb),sinr5_opt_e(iter,bb));

