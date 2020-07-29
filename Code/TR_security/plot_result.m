 % Plot of the FOM in function of the scenario investigated

%scenario_name of type "yparam_xparam_variableparam_constantparam"


% % Alwas plot a constellation of Bob and Eve whatever the scenario
% figure;
% subplot(1,2,1)
% stem(1:nb_taps, abs(h_bob)); hold on;
% stem(1:nb_taps, abs(h_eve),'r');
% legend('Bob PDP', 'Eve PDP')
% xlabel('Tap')
% ylabel('Amplitude')
% box on; grid on;
% subplot(1,2,2)
% plot(10*log10(abs(H_bob_TX))); hold on;
% plot(10*log10(abs(H_eve)),'r');
% legend('Bob FR', 'Eve FR')
% box on; grid on;
% xlabel('Bin')
% ylabel('Magnitude (dB)')


switch scenario
    % BER
    
    %1
    case "ber_ebno_bor_alpha"
        figure;
        semilogy(EbN0 , BER_bob,'Marker', 'o'); hold on;
        legendCell_bob = cellstr(num2str(BOR', 'Bob , BOR = %-d'));

        semilogy(EbN0 , BER_eve, 'Marker', 'x'); hold on;
        legendCell_eve = cellstr(num2str(BOR', 'Eve , BOR = %-d'));

        legendCell =[legendCell_bob;legendCell_eve]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
        legend(legendCell,'Location','best')

        box on; grid on;
        xlabel('$\frac{E_b}{N_0}$ (dB)')
        ylabel('Bit Error Rate (BER)')
        title(['Percentage of energy radiated dedicated for AN : ', num2str(100*(1-alpha)), '\%'])
        axis([min(EbN0) max(EbN0) 1e-6 1])
        
        
      
        
    %2    
    case "ber_ebno_alpha_bor"
        figure;
        semilogy(EbN0 , BER_bob,'Marker', 'o'); hold on;
        legendCell_bob = cellstr(num2str(100*(1-alpha'), 'Bob , AN energy = %-0.f %%'));

        semilogy(EbN0 , BER_eve, 'Marker', 'x'); hold on;
        legendCell_eve = cellstr(num2str(100*(1-alpha'), 'Eve , AN energy = %-0.f %%'));

        legendCell =[legendCell_bob;legendCell_eve]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
        legend(legendCell,'Location','SW')

        box on; grid on;
        xlabel('$\frac{E_b}{N_0}$ (dB)')
        ylabel('Bit Error Rate (BER)')
        title(['BOR = ', num2str(BOR)])
        axis([min(EbN0) max(EbN0) 1e-4 1])
        
        
    
    %3
    case "ber_alpha_bor_ebno"
        figure;
        semilogy(100*(1-alpha') , BER_bob,'Marker', 'o'); hold on;
        legendCell_bob = cellstr(num2str(BOR', 'Bob, BOR = %-d'));

        semilogy(100*(1-alpha') , BER_eve, 'Marker', 'x'); hold on;
        legendCell_eve = cellstr(num2str(BOR', 'Eve, BOR = %-d'));

        legendCell =[legendCell_bob;legendCell_eve];
        legend(legendCell,'Location','SW')

        box on; grid on;
        xlabel('Percentage of energy radiated dedicated for AN (\%)')
        ylabel('Bit Error Rate (BER)')
        title(['Noise level: $\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
        xlim([min(100*(1-alpha)) max(100*(1-alpha))])
        ylim([1e-6 1])
        
    
        
     %4   
     case "ber_alpha_ebno_bor" 
        figure;
        semilogy(100*(1-alpha') , BER_bob,'Marker', 'o'); hold on;
        legendCell_bob = cellstr(num2str(EbN0', 'Bob, $\\frac{Eb}{N_0}$ = %-d dB'));

        semilogy(100*(1-alpha') , BER_eve, 'Marker', 'x'); hold on;
        legendCell_eve = cellstr(num2str(EbN0', 'Eve, $\\frac{E_b}{N_0}$ = %-d dB'));

        legendCell =[legendCell_bob;legendCell_eve];
        legend(legendCell,'Location','best')

        box on; grid on;
        xlabel('Percentage of energy radiated dedicated for AN (\%)')
        ylabel('Bit Error Rate (BER)')
        title(['BOR = ', num2str(BOR)])
        xlim([min(100*(1-alpha')) max(100*(1-alpha))])
        ylim([1e-6 1])
        
        
       
        
    %5 
    case "ber_bor_ebno_alpha" 
        figure;
        semilogy(BOR , BER_bob,'Marker', 'o'); hold on;
        legendCell_bob = cellstr(num2str(EbN0', 'Bob, $\\frac{Eb}{N_0}$ = %-d dB'));

        semilogy(BOR , BER_eve, 'Marker', 'x'); hold on;
        legendCell_eve = cellstr(num2str(EbN0', 'Eve, $\\frac{E_b}{N_0}$ = %-d dB'));

        legendCell =[legendCell_bob;legendCell_eve]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
        legend(legendCell,'Location','SW')

        box on; grid on;
        xlabel('BOR')
        ylabel('Bit Error Rate (BER)')
        title(['Percentage of energy radiated dedicated for AN : ', num2str(100*(1-alpha)), '\%'])
        axis([min(BOR) max(BOR) 1e-7 1])
        
        
    %6
    case "ber_bor_alpha_ebno"
        figure;
        semilogy(BOR , BER_bob,'Marker', 'o'); hold on;
        legendCell_bob = cellstr(num2str(100*(1-alpha'), 'Bob , AN energy = %-0.f %%'));
        
        semilogy(BOR , BER_eve, 'Marker', 'x'); hold on;
        legendCell_eve = cellstr(num2str(100*(1-alpha'), 'Eve , AN energy = %-0.f %%'));

        legendCell =[legendCell_bob;legendCell_eve]; %legendCell_bob_no_AN;legendCell_eve_no_AN];
        legend(legendCell,'Location','SE')

        box on; grid on;
        xlabel('BOR')
        ylabel('Bit Error Rate (BER)')
        title(['Noise level: $\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
        axis([min(BOR) max(BOR) 1e-5 1])
        
        
        
        
        
        
    %% Secrecy rate
    
    %1
    case "secrecy_bor_ebno_alpha"
        figure;
        plot(BOR, secrecy_capa_despread ,'Marker', 'square'); hold on;
        legendCell_secrecy_despread = cellstr(num2str(EbN0', '$\\frac{E_b}{N_0}$ = %-d dB')); hold on;
        legend(legendCell_secrecy_despread,'Location','best')
        box on; grid on;
        title(['Secrecy before equalization, AN : ', num2str(100*(1-alpha)), '\%'])
        xlabel('BOR')
        ylabel('Secrecy rate (bit/s/Hz)')
        xlim([min(BOR) max(BOR)])        
        
        
        
    %2    
    case "secrecy_bor_alpha_ebno"
        figure;
        plot(BOR, secrecy_capa_despread ,'Marker', 'square'); hold on;
        legendCell_secrecy_despread = cellstr(num2str(100*(1-alpha'), 'AN energy = %-0.f %%')); hold on;
        legend(legendCell_secrecy_despread,'Location','best')
        box on; grid on;
        title(['Secrecy before equalization, noise level: $\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
        xlabel('BOR')
        ylabel('Secrecy rate (bit/s/Hz)')
        xlim([min(BOR) max(BOR)])

        
    %3    
    case "secrecy_alpha_bor_ebno"
        
        %%%%% Figure 2
        figure;
        plot(100*(1-alpha'), secrecy_capa_ergo ,'Marker', 'square'); hold on;
        for ii=1:length(energy_an_max_secrecy)
            max_value_theoretical(ii) = secrecy_capa_theoretical(ii,round(energy_an_max_secrecy(ii)));
        end
        plot(energy_an_max_secrecy,max_value_theoretical,'--','Marker','*')
         
        legendCell_secrecy_max=cellstr('Analytic max of SR');
        legendCell_secrecy_despread = cellstr(num2str(BOR', 'BOR = %-d')); hold on;
        legendCell = [legendCell_secrecy_despread;legendCell_secrecy_max];
        legend(legendCell,'Location','best')
        box on; grid on;
        title(['Secrecy before equalization, noise level: $\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
        xlabel('Percentage of energy radiated dedicated for AN (\%)')
        ylabel('Secrecy rate (bit/s/Hz)')
        xlim([min(100*(1-alpha')) max(100*(1-alpha'))])
%         filename = strcat('/graphs_2/scenarios/SISO/noEveKnowledge/Secrecy/secrecy_alpha_bor_ebno');
%         saveas(gcf,[pwd filename])
%         saveas(gcf,[pwd filename],'epsc')
        

        %%%%% Figure 2
        figure;
        plot(100*(1-alpha(2:end-1)') , secrecy_slope_theoretical(:,2:end-1)', '--')
        legendCell_secrecy_slope = cellstr(num2str(BOR', 'BOR = %-d')); hold on;
        legend(legendCell_secrecy_slope,'Location','best')
        box on; grid on;
        title(['Theoretical slope of secrecy rate function, noise level: $\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
        xlabel('Percentage of energy radiated dedicated for AN (\%)')
        ylabel('Slope of secrecy rate')
        xlim([min(100*(1-alpha')) max(100*(1-alpha'))])
        
%         filename = strcat('/graphs_2/scenarios/SISO/noEveKnowledge/Secrecy/secrecy_alpha_bor_ebno_slope');
%         saveas(gcf,[pwd filename])
%         saveas(gcf,[pwd filename],'epsc')
        
       
%         
        
        %%%%% Figure 3
        for ii = 1:size(secrecy_capa_ergo,1)
            [max_value(ii),max_idx(ii)] = max(secrecy_capa_ergo(ii,:));
        end
        figure;
        yyaxis left
        plot(BOR,energy_an_max_secrecy,'Marker','square') ; hold on;
        plot(BOR,max_idx,'Marker','o')
        ylabel('Optimal percentage of energy dedicated for AN (\%)')
        ylim([30 50])
    
        yyaxis right
        plot(BOR,max_value_theoretical,'Marker','square'); hold on
        plot(BOR, max_value,'Marker','o')
        ylabel('Max of secrecy rate (bit/s/Hz)')
        
        xlabel('BOR')
        title(['$\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
        box on; grid on;
        legend('Analytic', 'Simulation','Analytic', 'Simulation')
         
%         filename = strcat('/graphs_2/scenarios/SISO/noEveKnowledge/Secrecy/secrecy_alpha_bor_ebno_optimal_alpha_SR');
%         saveas(gcf,[pwd filename])
%         saveas(gcf,[pwd filename],'epsc')
%         

        
    %4    
    case "secrecy_alpha_ebno_bor"
        
        %%%%% Figure 1
        figure;
        plot(100*(1-alpha'), secrecy_capa_despread ,'Marker', 'square'); hold on;
        for ii=1:length(energy_an_max_secrecy)
            max_value_theoretical(ii) = secrecy_capa_despread(ii,round(energy_an_max_secrecy(ii)));
        end
        plot(energy_an_max_secrecy,max_value_theoretical,'--','Marker','*')
        
        legendCell_secrecy_despread = cellstr(num2str(EbN0', '$\\frac{E_b}{N_0}$ = %-d dB')); hold on;
        legendCell_secrecy_max=cellstr('Analytic max of SR');
        legendCell = [legendCell_secrecy_despread;legendCell_secrecy_max];
        legend(legendCell,'Location','best')
        
        box on; grid on;
        title(['Secrecy before equalization, BOR = ', num2str(BOR)])
        xlabel('Percentage of energy radiated dedicated for AN (\%)')
        ylabel('Secrecy rate (bit/s/Hz)')
        xlim([min(100*(1-alpha')) max(100*(1-alpha'))])
        filename = strcat('/graphs_2/scenarios/SISO/noEveKnowledge/Secrecy/secrecy_alpha_ebno_bor');
        saveas(gcf,[pwd filename])
        saveas(gcf,[pwd filename],'epsc')  
              
  
       
        %%%%% Figure 2
        figure;
        plot(100*(1-alpha(2:end-1)') , secrecy_slope_theoretical(:,2:end-1)', '--') ; hold on;
        legendCell_secrecy_slope = cellstr(num2str(EbN0', '$\\frac{E_b}{N_0}$ = %-d dB')); hold on;
        legend(legendCell_secrecy_slope,'Location','best') 
        box on; grid on;
        title(['Theoretical slope of secrecy rate function, BOR = ', num2str(BOR)])
        xlabel('Percentage of energy radiated dedicated for AN (\%)')
        ylabel('Slope of secrecy rate')
        xlim([min(100*(1-alpha')) max(100*(1-alpha'))])
        filename = strcat('/graphs_2/scenarios/SISO/noEveKnowledge/Secrecy/secrecy_alpha_ebno_bor_slope');
        saveas(gcf,[pwd filename])
        saveas(gcf,[pwd filename],'epsc')   
      
             
%   
        
        %%%%% Figure 3
        for ii = 1:size(secrecy_capa_despread,1)
            [max_value(ii),max_idx(ii)] = max(secrecy_capa_despread(ii,:));
        end
        figure;
        yyaxis left
        plot(EbN0,energy_an_max_secrecy,'Marker','square') ; hold on;
        plot(EbN0,max_idx,'Marker','o')
        ylabel('Optimal percentage of energy dedicated for AN (\%)')
    
        yyaxis right
        plot(EbN0,max_value_theoretical,'Marker','square'); hold on
        plot(EbN0, max_value,'Marker','o')
        ylabel('Max of secrecy rate (bit/s/Hz)')
        
        xlabel('$\frac{E_b}{N_0}$ (dB)')
        title(['BOR = ', num2str(BOR)])
        box on; grid on;
         legend('Analytic', 'Simulation','Analytic', 'Simulation')
        filename = strcat('/graphs_2/scenarios/SISO/noEveKnowledge/Secrecy/secrecy_alpha_ebno_bor_optimal_alpha_SR');
        saveas(gcf,[pwd filename])
        saveas(gcf,[pwd filename],'epsc')
        
        
        % figure;
%         plot(100*(1-alpha'), secrecy_capa_despread_eq ,'Marker', 'square'); hold on;
%         legendCell_secrecy_despread = cellstr(num2str(EbN0', '$\\frac{E_b}{N_0}$ = %-d dB')); hold on;
%         legend(legendCell_secrecy_despread,'Location','best')
%         box on; grid on;
%         title(['Secrecy after equalization, BOR = ', num2str(BOR)])
%         xlabel('Percentage of energy radiated dedicated for AN (\%)')
%         ylabel('Secrecy rate (bit/s/Hz)')
%         xlim([min(100*(1-alpha')) max(100*(1-alpha'))])
  
  
    %5
    case "secrecy_ebno_alpha_bor"
        figure;
        plot(EbN0, secrecy_capa_despread ,'Marker', 'square'); hold on;
        legendCell_secrecy_despread = cellstr(num2str(100*(1-alpha'), 'AN energy = %-0.f %%')); hold on;
        legend(legendCell_secrecy_despread,'Location','best')
        box on; grid on;
        title(['Secrecy before equalization, BOR = ', num2str(BOR)])
        xlabel('$E_b/N_0$ (dB)')
        ylabel('Secrecy rate (bit/s/Hz)')
        xlim([min(EbN0) max(EbN0)])
        
        
        
    %6
    case "secrecy_ebno_bor_alpha"    
        figure;
        plot(EbN0, secrecy_capa_despread ,'Marker', 'square'); hold on;
        legendCell_secrecy_despread = cellstr(num2str(BOR', 'BOR = %-d')); hold on;
        legend(legendCell_secrecy_despread,'Location','SE')
        box on; grid on;
        title(['Secrecy before equalization, AN : ', num2str(100*(1-alpha)), '\%'])
        xlabel('$E_b/N_0$ (dB)')
        ylabel('Secrecy rate (bit/s/Hz)')
        xlim([min(EbN0) max(EbN0)])
        
    %7.
    case "secrecy_alpha_bor_ebno_match_filt1"
        
        %%%%% FIGURE 1 : PERFORMANCE BW THE DIFFERENT DECODING STRUCTURES
        figure;
        plot(100*alpha', secrecy_capa_ergo_filt0 ,'Marker', 'square'); hold on;
        plot(100*alpha',secrecy_capa_ergo_filt1,'--','Marker','diamond') ; hold on;
        plot(100*alpha',secrecy_capa_ergo_filt2,'--','Marker','v') ; hold on;
        plot(100*alpha',secrecy_capa_ergo_filt3,'--','Marker','o') ; hold on;

      
        legendCell_secrecy_despread = cellstr(num2str(BOR', 'Despreading, BOR = %-d')); hold on;
        legendCell_secrecy_despread_filt1 = cellstr(num2str(BOR', 'Matched filtering, BOR = %-d')); hold on;
        legendCell_secrecy_despread_filt2 = cellstr(num2str(BOR', 'AN suppression, BOR = %-d')); hold on;
        legendCell_secrecy_despread_filt3 = cellstr(num2str(BOR', 'LMMSE, BOR = %-d')); hold on;


        legendCell = [legendCell_secrecy_despread;legendCell_secrecy_despread_filt1;legendCell_secrecy_despread_filt2;legendCell_secrecy_despread_filt3];
        legend(legendCell,'Location','best')
        box on; grid on;
        title(['Secrecy before equalization, noise level: $\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
        xlabel('Percentage of energy radiated dedicated for data signal (\%)')
        ylabel('Secrecy rate (bit/s/Hz)')
        xlim([min(100*alpha') max(100*alpha')])
        filename = strcat('/graphs_2/scenarios/SISO/noEveKnowledge/Secrecy/secrecy_alpha_bor_ebno');
%         saveas(gcf,[pwd filename])
%         saveas(gcf,[pwd filename],'epsc')
%         

        %%%%% FIGURE 2 : ANALYTIC vs SIMU - despreading
        figure;
        plot(100*alpha', secrecy_capa_ergo_filt0 ,'Marker', 'square'); hold on;
        plot(100*alpha',secrecy_capa_theoretical_filt0,'--','Marker','diamond') ; hold on;
        legendCell_secrecy_despread_filt0_simu = cellstr(num2str(BOR', 'Despreading simulation, BOR = %-d')); hold on;
        legendCell_secrecy_despread_filt0_anal = cellstr(num2str(BOR', 'Despreading analytic, BOR = %-d')); hold on;
        legendCell = [legendCell_secrecy_despread_filt0_simu;legendCell_secrecy_despread_filt0_anal];
        legend(legendCell,'Location','best')
        box on; grid on;
        title(['Secrecy before equalization, noise level: $\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
        xlabel('Percentage of energy radiated dedicated for data signal (\%)')
        ylabel('Secrecy rate (bit/s/Hz)')
        xlim([min(100*alpha') max(100*alpha')])
        
        
        %%%%% FIGURE 3 : ANALYTIC vs SIMU  - Matched filtering
        figure;
        plot(100*alpha', secrecy_capa_ergo_filt1 ,'Marker', 'square'); hold on;
        plot(100*alpha',secrecy_capa_theoretical_filt1,'--','Marker','diamond') ; hold on;
        legendCell_secrecy_despread_filt1_simu = cellstr(num2str(BOR', 'Matched filtering simulation, BOR = %-d')); hold on;
        legendCell_secrecy_despread_filt1_anal = cellstr(num2str(BOR', 'Matched filtering analytic, BOR = %-d')); hold on;
        legendCell = [legendCell_secrecy_despread_filt1_simu;legendCell_secrecy_despread_filt1_anal];
        legend(legendCell,'Location','best')
        box on; grid on;
        title(['Secrecy before equalization, noise level: $\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
        xlabel('Percentage of energy radiated dedicated for data signal (\%)')
        ylabel('Secrecy rate (bit/s/Hz)')
        xlim([min(100*alpha') max(100*alpha')])

%          %%%%% Figure 2
%         for ii = 1:size(secrecy_capa_ergo,1)
%             [max_value(ii),max_idx(ii)] = max(secrecy_capa_ergo(ii,:));
%             [max_value_filt1(ii),max_idx_filt1(ii)] = max(secrecy_capa_ergo_filt1(ii,:));
%             [max_value_filt2(ii),max_idx_filt2(ii)] = max(secrecy_capa_ergo_filt2(ii,:));
% 
%         end
%         figure;
%         yyaxis left
%         plot(BOR,100-max_idx,'Marker','square') ; hold on;
%         plot(BOR,100-max_idx_filt1,'Marker','o')
%         plot(BOR,100-max_idx_filt2,'Marker','x')
% 
%         ylabel('Optimal percentage of energy dedicated for data signal(\%)')
% 
%     
%         yyaxis right
%         plot(BOR,max_value,'Marker','v'); hold on
%         plot(BOR, max_value_filt1,'Marker','diamond')
%         plot(BOR,max_value_filt2,'Marker','x')
%         ylabel('Max of secrecy rate (bit/s/Hz)')
% 
%         
%         xlabel('BOR')
%         title(['$\frac{E_b}{N_0}$ = ', num2str(EbN0), ' dB'])
%         box on; grid on;
%         legend('No matched filtering','Matched filtering','AN killed','No matched filtering', 'Matched filtering','AN killed')
%          
%         filename = strcat('/graphs_2/scenarios/SISO/noEveKnowledge/Secrecy/secrecy_alpha_bor_ebno_optimal_alpha_SR');
% %         saveas(gcf,[pwd filename])
% %         saveas(gcf,[pwd filename],'epsc')
%         
        
        
end
