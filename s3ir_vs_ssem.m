%==========================================================================
%  S3IR CONTAGION MODEL vs SSCM and SCM
%--------------------------------------------------------------------------
%  Description:
%    This script compares the S3IR model using the non-standard Euler–
%    Maruyama discretisation and the SSCM using the (standard) Euler–
%    Maruyama and SCM.
%
%  Population Catagories: 
%    S: Susceptible (exposed) individuals
%    I: Infected (infectious) individuals
%    R: Recovered (immune) individuals
%
%    The total (constant) population density: S + I + R = 1.
%
%  Formulation:
%    For S31I and SSCM, the infected compartment is driven by a Brownian 
%    process. The diffusion (noise) coefficients sigma_k scale the 
%    stochastic perturbation of the transition rates beta_k for k=1,...,6 
%    (S3IR) and k=1,2 (SSCM).
%
%=========================================================================
clear all; close all; clc; 
rng('default')
%=========================================================================
%% ========================================================================
%  Simulation settings
%  ========================================================================
%  T  : simulation time 
%  dt : discretisation step
%  N  : number of time steps = T/dt
% =========================================================================
T = 100; 
dt= 1e-3; 
N = T/dt; 

%% ========================================================================
% Model Parameters
%%=========================================================================
% β1,...,β6 : infection transmission rates in pairwise, ternary, 
% quaternary, etc, contacts respectively
% γ         : recovery rate
% σ1,...,σ6 : noise intensities

beta = 0.05;    beta1 = 0.05;   beta2 = 0.05;   beta3 = 0.05;
                beta4 = 0.05;   beta5 = 0.05;

sgma = 0.13;    sgma1 = 0.13;   sgma2 = 0.13;   sgma3 = 0.13;
                sgma4 = 0.13;   sgma5 = 0.13;

bEpidemic = true;

if bEpidemic    fname = 'Epi';  else    fname = 'noEpi';    end

if bEpidemic                % outbreak will spread
    gamm = 0.001;
else                        % outbreak will die out
    gamm = 0.06;
end


%% ======================================================================= 
% Simulate S3IR model
% ========================================================================
% Alocate memory for arrays

NR = 1;                        % # realizations
S = zeros(N, NR);    I = zeros(N, NR);    R = zeros(N, NR);
S_s2 = zeros(N, NR); I_s2 = zeros(N, NR); R_s2 = zeros(N, NR);
S_s1 = zeros(N, NR); I_s1 = zeros(N, NR); R_s1 = zeros(N, NR);

sqrt_dt = sqrt(dt);
lam = 1;
phi = (1-exp(-lam*dt))/lam;    % denominator function (phi) of NSFD scheme

for r = 1:NR

    %Set initial population states
    I(1,r) = 0.3;    R(1,r) = 0.1;    S(1,r) = 1 - I(1,r) - R(1,r);
    I_s2(1,r) = 0.3; R_s2(1,r) = 0.1; S_s2(1,r) = 1 - I_s2(1,r) - R_s2(1,r);
    I_s1(1,r) = 0.3; R_s1(1,r) = 0.1; S_s1(1,r) = 1 - I_s1(1,r) - R_s1(1,r);


    for k = 1:(N-1)
    
        dB = sqrt_dt*randn;     % Brownian motion increment
        
        % S3IR model
        I(k+1,r) = I(k,r) + ...
                  (S(k,r)*(beta*I(k,r) + beta1*(I(k,r)^2) + beta2*(I(k,r)^3) + beta3*(I(k,r)^4) + beta4*(I(k,r)^5) + beta5*(I(k,r)^6)) - gamm*I(k,r))*phi ...
                 + S(k,r)*(sgma*I(k,r) + sgma1*(I(k,r)^2) + sgma2*(I(k,r)^3) + sgma3*(I(k,r)^4) + sgma4*(I(k,r)^5) + sgma5*(I(k,r)^6))*dB;

        R(k+1,r) = R(k,r) + gamm*I(k,r)*phi;
        S(k+1,r) = 1 - I(k+1,r) - R(k+1,r);  

        % SSCM model
        I_s2(k+1,r) = I_s2(k,r) + ...
                      (S_s2(k,r)*(beta*I_s2(k,r) + beta1*(I_s2(k,r)^2) ) - gamm*I_s2(k,r))*dt ...
                     + S_s2(k,r)*(sgma*I_s2(k,r) + sgma1*(I_s2(k,r)^2) )*dB;

        R_s2(k+1,r) = R_s2(k,r) + gamm*I_s2(k,r)*dt;
        S_s2(k+1,r) = 1 - I_s2(k+1,r) - R_s2(k+1,r);  
        
        % SCM model
        I_s1(k+1,r) = I_s1(k,r) + ...
                      (S_s1(k,r)*(beta*I_s1(k,r) + beta1*(I_s1(k,r)^2) ) - gamm*I_s1(k,r))*dt;

        R_s1(k+1,r) = R_s1(k,r) + gamm*I_s1(k,r)*dt;
        S_s1(k+1,r) = 1 - I_s1(k+1,r) - R_s1(k+1,r);  
    end
end


%% =======================================================================
% Plot sample paths
%=========================================================================
% susceptible 
f = figure;%('visible', 'off');
plot((1:N)*dt,S); axis([0 T 0 1]); 
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$S(t)$', 'Interpreter','latex', 'FontSize',15);
grid on;
hold on;
plot((1:N)*dt,S_s2); 
plot((1:N)*dt,S_s1); 
legend('S3IR','SSCM','SCM','Interpreter','latex','FontSize',12);
exportgraphics(f, strcat('Svs_',fname,'.png'));


% infected 
f = figure;%('visible', 'off');
plot((1:N)*dt,I); axis([0 T 0 1]);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$I(t)$', 'Interpreter','latex', 'FontSize',15);
grid on;
hold on;
plot((1:N)*dt,I_s2); 
plot((1:N)*dt,I_s1); 
legend('S3IR','SSCM','SCM','Interpreter','latex','FontSize',12,'Location','southeast');
exportgraphics(f, strcat('Ivs_',fname,'.png'));


% removed 
f = figure;%('visible', 'off');
plot((1:N)*dt,R); axis([0 T 0 1]);
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$R(t)$', 'Interpreter','latex', 'FontSize',15);
grid on;
hold on;
plot((1:N)*dt,R_s2); 
plot((1:N)*dt,R_s1); 
legend('S3IR','SSCM','SCM','Interpreter','latex','FontSize',12);
exportgraphics(f, strcat('Rvs_',fname,'.png'));
