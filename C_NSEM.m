%==========================================================================
%  Non-Standard Euler–Maruyama (NSEM) Scheme
%  Comparison Across Simplicial Dimensions  D = 1,...,7
%--------------------------------------------------------------------------
%  Stochastic Formulation:
%    Each compartment is driven by an independent Wiener process W_k(t).
%    The diffusion (noise) coefficients sigma_k scale the stochastic
%    perturbations applied to the corresponding transition rates beta_k for k=1,2,3.
%--------------------------------------------------------------------------
%  Non-Standard Scheme:
%    The denominator function phi = (exp(gamm*dt)-1)/gamm; replaces h in the
%    finite-difference approximation, ensuring dynamic consistency with
%    the continuous model's qualitative behaviour (positivity, boundedness)
%    at any step size.
%--------------------------------------------------------------------------
% DIMENSIONS COMPARED:
%  
%    D = 1  →  1-simplex only  (standard pairwise)
%    D = 2  →  2-simplex (triangular contagion)
%    D = 3  →  3-simplex (include tetrahedral connections too)
%    D = 4  →  4-simplex  
%    D = 5  →  5-simplex
%    D = 6  →  6-simplex
%    D = 7  →  7-simplex 
%==------------------------------------------------------------------------

clc;  
rng('default')
%% ========================================================================
%  Section 1: Simulation time settings
%  ========================================================================
%  T  : total simulation time (continuous-time units)
%  dt : discrete time step
%  N  : number of discrete time steps = T/dt
% =========================================================================

bEpidemic = true;

if bEpidemic    fname = 'Epi';  else    fname = 'noEpi';    end
T = 10;                     % simulation time
dt= 1e-3;                   % discretization step
N = int32(T/dt);            % # subintervals


%% ========================================================================
% Section 2:Model Parameters
%==========================================================================
% Parameters are fixed across all dimension runs so that differences in
%  dynamics are attributable solely to the simplicial structure.
if bEpidemic                % outbreak will spread
    gamm = 0.001;
    beta = 0.3;  beta1 = 0.3;  beta2 = 0.3;  beta3 = 0.3; beta4 = 0.3;  beta5 = 0.3;  beta6 = 0.3;
    sgma = 0.45; sgma1 = 0.35; sgma2 = 0.25; sgma3 = 0.4; sgma4 = 0.40; sgma5 = 0.40; sgma6 = 0.40;

else                        % outbreak will die out
    gamm = 0.3;
    beta = 0.3;  beta1 = 0.2;  beta2 = 0.1;  beta3 = 0.1;  beta4 = 0.1;  beta5 = 0.1;  beta6 = 0.1;
    sgma = 0.45; sgma1 = 0.35; sgma2 = 0.25; sgma3 = 0.40; sgma4 = 0.30; sgma5 = 0.25; sgma6 = 0.25; 
end
%% ======================================================================= 
% Section 3: Allocation and Initialization
% ========================================================================
% Preallocating arrays avoids dynamic memory reallocation inside the loop,
%  which is critical for performance in large-step simulations.


S = zeros(1, N+1);              % allocate arrays 
I = zeros(1, N+1);
R = zeros(1, N+1);
I_1 = zeros(1,N+1);
I_2 = I_1;  I_3 = I_1;  I_4 = I_1;  I_5 = I_1;  I_6 = I_1;  I_7 = I_1;  
                               % initialization
I_1(1)=0.3; I_2(1)=0.3; I_3(1)=0.3; I_4(1)=0.3; I_5(1) = 0.3; I_6(1)=0.3; I_7(1)=0.3;
R(1) = 0.2;
S(1) = 1-I_1(1)-R(1);
sqrt_dt = sqrt(dt);
phi = (exp(gamm*dt)-1)/gamm;    % denominator function (phi) of NSFD scheme


for k = 1:N

    dB = sqrt_dt*randn;     % Brownian motion increment
  
%% ======================================================================= 
% Section 4: Simulate S3IR model for D=1
% ========================================================================
    I_1(k+1) = I_1(k) + ( S(k)*(beta*I_1(k)) - gamm*I_1(k) )*phi ...
             + S(k)*( sgma*I_1(k) )*dB;
       
    R(k+1) = R(k) + gamm*I_1(k)*phi;
    S(k+1) = 1 - I_1(k+1) - R(k+1);
%% ======================================================================= 
% Section 5: Simulate S3IR model for D=2
% ========================================================================
    I_2(k+1) = I_2(k) +( S(k)*(beta*I_2(k)+ beta1*(I_2(k))^2) - gamm*I_2(k) )*phi ...
           + S(k)*( sgma*I_2(k)+ sgma1*(I_2(k)^2) )*dB;
    
    R(k+1) = R(k) + gamm*I_2(k)*phi;
    S(k+1) = 1 - I_2(k+1) - R(k+1);      
%% ======================================================================= 
% Section 6: Simulate S3IR model for D=3
% ========================================================================
    I_3(k+1) = I_3(k) + ( S(k)*(beta*I_3(k)+ beta1*(I_3(k))^2+ beta2*(I_3(k))^3) - gamm*I_3(k) )*phi ...
             + S(k)*( sgma*I_3(k)+ sgma1*(I_3(k)^2)+ sgma2*(I_3(k)^3) )*dB;
    
    R(k+1) = R(k) + gamm*I_3(k)*phi;
    S(k+1) = 1 - I_3(k+1) - R(k+1);
%% ======================================================================= 
% Section 7: Simulate S3IR model for D=4
% ========================================================================
    I_4(k+1) = I_4(k) + ( S(k)*(beta*I_4(k)+ beta1*(I_4(k))^2+ beta2*(I_4(k))^3+ beta3*(I_4(k))^4) - gamm*I_4(k) )*phi ...
             + S(k)*( sgma*I_4(k)+ sgma1*(I_4(k)^2)+ sgma2*(I_4(k)^3)+ sgma3*(I_4(k))^4 )*dB;
        
    R(k+1) = R(k) + gamm*I_4(k)*phi;
    S(k+1) = 1 - I_4(k+1) - R(k+1);
%% ======================================================================= 
% Section 8: Simulate S3IR model for D=5
% ========================================================================
    I_5(k+1) = I_5(k) + ( S(k)*(beta*I_5(k)+ beta1*(I_5(k))^2+ beta2*(I_5(k))^3+ beta3*(I_5(k))^4+ beta4*(I_5(k))^5) - gamm*I_5(k) )*phi ...
             + S(k)*( sgma*I_5(k)+ sgma1*(I_5(k)^2)+ sgma2*(I_5(k)^3)+ sgma3*(I_5(k))^4+ sgma4*(I_5(k))^5 )*dB;
    
    R(k+1) = R(k) + gamm*I_5(k)*phi;
    S(k+1) = 1 - I_5(k+1) - R(k+1);
%% ======================================================================= 
% Section 9: Simulate S3IR model for D=6
% ========================================================================
    I_6(k+1) = I_6(k) + ( S(k)*(beta*I_6(k)+ beta1*(I_6(k))^2+ beta2*(I_6(k))^3+ beta3*(I_6(k))^4+ beta4*(I_6(k))^5+ beta5*(I_6(k))^6) - gamm*I_6(k) )*phi ...
             + S(k)*( sgma*I_6(k)+ sgma1*(I_6(k)^2)+ sgma2*(I_6(k)^3)+ sgma3*(I_6(k))^4+ sgma4*(I_6(k))^5+ sgma5*(I_6(k))^6 )*dB;

    R(k+1) = R(k) + gamm*I_6(k)*phi;
    S(k+1) = 1 - I_6(k+1) - R(k+1);
 %% ======================================================================= 
% Section 10: Simulate S3IR model for D=7
% ========================================================================
    I_7(k+1) = I_7(k) + ( S(k)*(beta*I_7(k)+ beta1*(I_7(k))^2+ beta2*(I_7(k))^3+ beta3*(I_7(k))^4+ beta4*(I_7(k))^5+ beta5*(I_7(k))^6+ beta6*(I_7(k))^7) - gamm*I_7(k))*phi ...
             + S(k)*( sgma*I_7(k)+ sgma1*(I_7(k)^2)+ sgma2*(I_7(k)^3)+ sgma3*(I_7(k))^4+ sgma4*(I_7(k))^5+ sgma5*(I_7(k))^6+ sgma6*(I_7(k))^7 )*dB;

    R(k+1) = R(k) + gamm*I_7(k)*phi;
    S(k+1) = 1 - I_7(k+1) - R(k+1);
end

I = [I_1; I_2; I_3; I_4; I_5; I_6; I_7]';


%% =======================================================================
% Section 11: Construct plots
%=========================================================================

f = figure('visible', 'off');
plot(0:dt:T,I);
ylim([0 1])
ax = gca; ax.FontSize = 10; ax.PlotBoxAspectRatio = [1 1 1];
xlabel('$t$','Interpreter','latex', 'FontSize',15);
title('$I(t)$', 'Interpreter','latex', 'FontSize',15);
grid on;
legend('$D$=1','$D$=2','$D$=3','$D$=4','$D$=5','$D$=6','$D$=7','Interpreter','latex','FontSize',10,'location','southeast');
exportgraphics(f, strcat(fname,'_dim.png'));
