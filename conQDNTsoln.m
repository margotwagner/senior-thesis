function conQDNTsoln

% Plot estimates vs actual data
% Estimate
[tout,cout] = ode23s(@conQDnoTumor,[0 96],[40; 0; 0; 0; 0;0]);
[tout2,cout2] = ode23s(@conQDnoTumor,[0 2],[40; 0; 0; 0; 0;0]);


figure(1)
plot(tout,cout,'LineWidth',2)
axis([0 96 0 40])
title('Concentration of ConQDs in Various Organs over 4 days (in nM)')
xlabel('Time (hours)')
ylabel('Concentration (nM)')
legend('C_p_l_a_s_m_a','C_k_i_d_n_e_y','C_l_i_v_e_r','C_s_p_l_e_e_n','C_l_u_n_g','C_o_t_h_e_r','Location','northeast')

figure(2)
plot(tout2,cout2,'LineWidth',2)
title('Concentration of ConQDs in Various Organs over 2 hours (in nM)')
xlabel('Time (hours)')
ylabel('Concentration (nM)')
legend('C_p_l_a_s_m_a','C_k_i_d_n_e_y','C_l_i_v_e_r','C_s_p_l_e_e_n','C_l_u_n_g','C_o_t_h_e_r','Location','northeast')


%cout
cout(end,:)


% Cumulative equation
%dcdt = [cblp; ckp; clip; csp; clp; cop];  % Blood, Kidney, Liver, Spleen, Lung, Other