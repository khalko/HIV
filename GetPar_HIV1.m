function [coef, delays] = GetPar_HIV1
    coef.r_A = 50000.0;   %1
    coef.r_E0 = 12.0;   %2
    coef.mu_A = 0.01;    %3
    coef.mu_L = 0.01;    %4
    coef.mu_I = 0.01;    %5
    coef.mu_C = 0.01;    %6
    coef.mu_E0 = 0.01;   %7
    coef.mu_Q = 0.01;    %8
    coef.mu_E = 0.08;    %9
    coef.mu_U = 3.0;     %10
    coef.mu_V = 3.0;     %11
    coef.nu_I = 110.0;   %12
    coef.sigma_I = 0.0032;%13
    coef.p_L = 0.3;      %14
    coef.a_L = 0.05;      %15
    
    coef.n_Q1 = 12.0;    %18
    coef.n_Q2 = 10.0;    %19
    coef.n_E1 = 16.0;    %20
    coef.n_E2 = 12.0;    %21
    coef.p_E = 0.12;     %22
    coef.g_LV = 0.1e-7; %23
    coef.g_CV = 0.1e-7; %24
    coef.g_AV = 0.1e-7; %25
    coef.g_LI = 0.2e-6;    %26
    coef.g_AI = 0.2e-6;    %27
    coef.g_EI = 1.5e-6;  %28
    coef.g_E0AV = 7e-14; %29
    coef.g_QAV = 4e-14;  %30
    coef.w_C = 0.2;      
    coef.w_U = 0.02; 
    %coef.delta = 0.5e-3;
    
    
    delays.w_C = 0.2;      %16
    delays.w_U = 0.02;     %17
    delays.w_E1 = 2.5;     %31
    delays.w_E2 = 2.3;     %32
    
end

