function TconQDsoln

% Plot estimates vs actual data
% Estimate
[tout,cout] = ode23s(@TconQD,[0 96],[40; 0; 0; 0; 0;0;0]);
[tout2,cout2] = ode23s(@TconQD,[0 10],[40; 0; 0; 0; 0;0;0]);

figure(1)
plot(tout,cout,'LineWidth',2)
axis([0 96 0 40])
title('ConQD Concentration in Tissues After 4 Days (Tumor)')
xlabel('Time (hours)')
ylabel('Concentration (nM)')
%legend('C_p_l_a_s_m_a','C_k_i_d_n_e_y','C_l_i_v_e_r','C_s_p_l_e_e_n','C_l_u_n_g','C_o_t_h_e_r','C_t_u_m_o_r','Location','northeast')

figure(2)
plot(tout2,cout2,'LineWidth',2)
title('ConQD Concentration in Tissues After 10 Hours (Tumor)')
xlabel('Time (hours)')
ylabel('Concentration (nM)')
%legend('C_p_l_a_s_m_a','C_k_i_d_n_e_y','C_l_i_v_e_r','C_s_p_l_e_e_n','C_l_u_n_g','C_o_t_h_e_r','C_t_u_m_o_r','Location','northeast')

%cout
cout(end,:)


% Cumulative equation
%dcdt = [cblp; ckp; clip; csp; clp; cop];  % Blood, Kidney, Liver, Spleen, Lung, Other