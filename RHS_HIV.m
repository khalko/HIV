function [UZ, UD, F] = RHS_HIV(coef)
syms A L I V E0 Q E 
syms A_w_E1 L_w_E1 I_w_E1 V_w_E1 E0_w_E1 Q_w_E1 E_w_E1 
syms A_w_E2 L_w_E2 I_w_E2 V_w_E2 E0_w_E2 Q_w_E2 E_w_E2
syms A_w_U L_w_U I_w_U V_w_U E0_w_U Q_w_U E_w_U
syms A_w_C L_w_C I_w_C V_w_C E0_w_C Q_w_C E_w_C
UZ = [A; L; I; V; E0; Q; E];

UD.w_U = [A_w_U; L_w_U; I_w_U; V_w_U; E0_w_U; Q_w_U; E_w_U];
UD.w_C = [A_w_C; L_w_C; I_w_C; V_w_C; E0_w_C; Q_w_C; E_w_C];
UD.w_E2 = [A_w_E2; L_w_E2; I_w_E2; V_w_E2; E0_w_E2; Q_w_E2; E_w_E2];
UD.w_E1 = [A_w_E1; L_w_E1; I_w_E1; V_w_E1; E0_w_E1; Q_w_E1; E_w_E1];






%((1-p_L)*g_AV*A*V+g_AI*A*I+g_LI*L*I+a_L*L)
%rho_C = (1 - coef.p_L) * coef.g_AV * A * V + coef.g_AI *  A * I + coef.g_LI * L * I + coef.a_L * L;
%rho_U = coef.nu_I * I;
rho_E1 = coef.g_E0AV * E0 * A * V;
rho_E2 = coef.g_QAV * Q * A * V;

% syms s;
% F_integral = exp(-coef.mu_C*(T_now - s))*rho_C;
% integral = int(F_integral,[(T_now - coef.w_C) T_now]);

rho_U_w_U = coef.nu_I * I_w_U;
rho_C_w_C = (1 - coef.p_L) * coef.g_AV * A_w_C * V_w_C + coef.g_AI *  A_w_C * I_w_C + coef.g_LI * L_w_C * I_w_C + coef.a_L * L_w_C;
rho_E1_w_E1 = coef.g_E0AV * E0_w_E1 * A_w_E1 * V_w_E1;
rho_E2_w_E2 = coef.g_QAV * Q_w_E2 * A_w_E2 * V_w_E2;

F(1) = coef.r_A-coef.mu_A*A-coef.g_AV*A*V-coef.g_AI*A*I;
F(2) = -(coef.mu_L+coef.a_L)*L-coef.g_LI*L*I+coef.p_L*coef.g_AV*A*V;
F(3) = -(coef.mu_I+coef.sigma_I*coef.nu_I)*I-coef.g_EI*E*I+exp(-coef.mu_C*coef.w_C)*rho_C_w_C;
F(4) = -coef.mu_V*V-coef.g_AV*A*V-coef.g_LV*L*V+exp(-coef.mu_U*coef.w_U)*rho_U_w_U ;
F(5) = coef.r_E0-coef.mu_E0*E0-rho_E1;
F(6) = -coef.mu_Q*Q-rho_E2+coef.n_Q1*rho_E1_w_E1+coef.n_Q2*rho_E2_w_E2;
F(7) = -coef.mu_E*E-coef.p_E*coef.g_EI*E*I+coef.n_E1*rho_E1_w_E1+coef.n_E2*rho_E2_w_E2;
%-coef.g_CV*V*integral
%
end

