%% BESSELPACK demo
%
% This script illustrates the use of the main routines included in BESSELPACK.
% It provides examples for the evaluation of Bessel functions and their
% derivatives, and for the computation of zeros of general cylinder
% functions, derivative-related equations, and cross products.

clear; clc;
fprintf('---------------------------------------------\n');
fprintf('BESSELPACK demonstration script\n');
fprintf('---------------------------------------------\n\n');
%------------------------------------------------------
%% 1. Evaluation of Bessel functions and derivatives
%
% Compute J_nu(x), Y_nu(x) and their first derivatives
%-------------------------------------------------------
nu = 2.5;
x  = 10;
[J,Y,Jp,Yp,iflag] = bessel(nu,x);
fprintf('Example 1: Bessel functions\n');
fprintf('nu = %g, x = %g\n\n',nu,x);
fprintf('J_nu(x)  = % .16e\n',J);
fprintf('Y_nu(x)  = % .16e\n',Y);
fprintf('Jp_nu(x) = % .16e\n',Jp);
fprintf('Yp_nu(x) = % .16e\n',Yp);
fprintf('error flag = %d\n\n',iflag);
%----------------------------------------------------
%% 2. Zeros of a general cylinder function
%
% C_{nu,alpha}(x)=cos(alpha)J_nu(x)-sin(alpha)Y_nu(x)
%-----------------------------------------------------
nu = 0.5;
alpha = pi/3;
a = 0.1;
b = 50;
tol = 1e-12;
[z] = besselz(a,b,nu,alpha,tol);
fprintf('Example 2: zeros of a cylinder function\n');
fprintf('nu = %g, alpha = %g\n',nu,alpha);
disp(z);
% Plot the function and its zeros
xx = linspace(a,b,2000);
CC = zeros(size(xx));
for k=1:length(xx)
    [J,Y] = bessel(nu,xx(k));
    CC(k)=cos(alpha)*J-sin(alpha)*Y;
end
figure;
plot(xx,CC,'LineWidth',1);
hold on
plot(z,zeros(size(z)),'o');
grid on
xlabel('x');
ylabel('C_{\nu,\alpha}(x)');
title('Zeros of a general cylinder function');
%---------------------------------------------
%% 3. Zeros of 
%
%     x C'_{nu,alpha}(x)+gamma C_{nu,alpha}(x)
%----------------------------------------------
nu = 1;
alpha = pi/4;
gamma = 0.5;
a = 0.1;
b = 50;
tol = 1e-12;
[z] = besselzd(a,b,nu,gamma,alpha,tol);
fprintf('\nExample 3: derivative-related zeros\n');
disp(z);
%-------------------------------------------------
%% 4. Zeros of cross products of Bessel functions
%-------------------------------------------------
nu = 10;
lambda = 0.6;
icho = 1;      % selected cross product
a = 0.1;
b = 50;
tol = 1e-12;
[z] = crosspz(icho,nu,lambda,a,b,tol);
fprintf('\nExample 4: cross-product zeros\n');
fprintf('icho = %d, lambda = %g\n',icho,lambda);
disp(z);
fprintf('\nDemonstration completed successfully.\n');