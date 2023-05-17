function [coef, delays] = GetPar_HIV3
    coef.r_A = 2000;   %1
    coef.r_E0 = 120;   %2
    coef.mu_A = 0.01;    %3
    coef.mu_L = 0.01;    %4
    coef.mu_I = 0.01;    %5
    coef.mu_C = 0.01;    %6
    coef.mu_E0 = 0.01;   %7
    coef.mu_Q = 0.05;    %8
    coef.mu_E = 0.08;    %9
    coef.mu_U = 3.0;     %10
    coef.mu_V = 3.0;     %11
    coef.nu_I = 110.0;   %12
    coef.sigma_I = 0.0032;%13
    coef.p_L = 0.3;      %14
    coef.a_L = 0.02;      %15
    delays.w_C = 0.1;      %16
    delays.w_U = 0.02;     %17
    coef.n_Q1 = 12.0;    %18
    coef.n_Q2 = 10.0;    %19
    coef.n_E1 = 16.0;    %20
    coef.n_E2 = 12.0;    %21
    coef.p_E = 0.12;     %22
    coef.g_LV = 2e-6; %23
    coef.g_CV = 1.5e-6; %24
    coef.g_AV = 0.12e-6; %25
    coef.g_LI = 1e-6;    %26
    coef.g_AI = 0.15e-6;    %27
    coef.g_EI = 3.5e-5;  %28
    coef.g_E0AV = 7e-10; %29
    coef.g_QAV = 4e-10;  %30
    delays.w_E1 = 1.5;     %31
    delays.w_E2 = 1.2;     %32
    
    coef.w_C = 0.1;      
    coef.w_U = 0.02; 
  
end