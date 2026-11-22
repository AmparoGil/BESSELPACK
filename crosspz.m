function [xzer,xzbad,nit]=crosspz(icho,pnu,lambda,a,b,tol)
%--------------------------------------------------------------------
%  Computation of the zeros of the following cross products of
%  Bessel functions:
%
%  A) J_nu(x)*Y_nu(lambda*x)-J_nu(lambda*x)*Y_nu(x)
%
%  B) J'_nu(x)*Y'_nu(lambda*x)-J'_nu(lambda*x)*Y'_nu(x)
%
%  C) J'_nu(x)*Y_nu(lambda*x)-J_nu(lambda*x)*Y'_nu(x)
%
%  in the interval (a,b), a>0. 
%--------------------------------------------------------------------
% INPUTS:
%
%   icho     selects the cross product:
%              icho = 1  --> cross product A
%              icho = 2  --> cross product B
%              icho = 3  --> cross product C
%
%   pnu      order nu of the Bessel functions (real)
%   lambda   scaling parameter for the argument, lambda>0, lambda~=1
%   a        lower endpoint of the interval
%   b        upper endpoint of the interval
%   tol      requested relative accuracy for the computed zeros
%            (the routine uses max(tol,10*eps/|1-lambda|), with
%            lambda replaced by 1/lambda if lambda>1)
%
% OUTPUTS:
%
%   xzer      array containing the computed zeros
%   xzbad     zeros (also included in xzer) for which the requested
%             accuracy was not reached in itmax iterations
%   nit       (optional) number of iterations used for each zero
%--------------------------------------------------------------------
% Accompanying software of the paper 
% "Numerical Software for Bessel functions and Associated
%  Values" 
% Authors: Amparo Gil, Javier Segura and Nico M. Temme
%----------------------------------------------------------------
xzer=[];
xzbad=[];
nit=[];
iter=[];
if ~any(icho==[1 2 3])
  error('icho should be 1, 2 or 3');
end
if lambda==1
  return
end    
% The cross products are even functions of nu
pnu=abs(pnu);
% For lambda>1, the zeros are those of the cross product with 1/lambda
% (for icho=3, of the fourth cross product, icho=4 internally) in
% (lambda*a,lambda*b), divided by lambda
ilamg=0;
if lambda>1
  a=a*lambda;
  b=b*lambda;
  lambda=1/lambda;
  ilamg=1;
  if icho==3
    icho=4;
  end
end    
% The zeros are ill-conditioned with respect to lambda as lambda->1
tol=max(tol,10*eps/(1-lambda));
Pim=pi-eps;
itmax=50;
if b<=pnu
  return
end
AA=1.5;
BB=20;  
n0=4;
if lambda<0.7
   wp=2*BB*(1-lambda)^2;
else
   wp=BB*(1-1/lambda)^2;
end    
if icho==1
  if pnu>1/2
     x_0=max(pi/(1-lambda),pnu+1.855*pnu^(1/3));
  else
     x_0=(2+1.14*lambda)/(1-lambda);
  end    
  x_0=max(x_0,a);
elseif icho==3
  if pnu>1/2
     x_0=max(pi/2/(1-lambda)-0.42,pnu+0.8086*pnu^(1/3));
  else
     x_0=(0.5+lambda)/(1-lambda);
     if lambda<0.01
       x_0=min(x_0,0.9*max(sqrt(2/log(2/lambda)),sqrt(2*pnu)));
     end
  end    
  x_0=max(x_0,a);
elseif icho==4
  % first zero > pi/(2(1-lambda)) (numerical evidence) and > nu
  x_0=0.95*pi/2/(1-lambda);
  if pnu>1/2
     x_0=max(x_0,pnu+0.8086*pnu^(1/3));
  end    
  x_0=max(x_0,a);
else  
  x_0=max(a,pnu);
end
x_e=b;
x=x_0;
n=0;
rw=1;
while (x_e-x)>0 && rw>0 
  err=tol+1;
  i=0;
  while err>tol && (x_e-x)>0 && rw>0 && i<itmax
    i=i+1;
    h=fu(icho,lambda,pnu,x);
    rw=wp;
    if (rw > 0) 
      cortt=rw*h^2/3;
      if cortt < 0.001 
        errs=h*(1-cortt);
      else
        w=sqrt(abs(rw));
        angle=atan(w*h);
        if h>0 
          angle=angle-Pim;
        end       
        errs=angle/w;
      end           
      x=x-errs;
      err=abs(errs/x);
    end
  end
  if i==itmax
    if (x_e-x)>0 && rw>0 
      n=n+1;
      xzer(n)=x;    
      nb=numel(xzbad)+1;
      xzbad(nb)=x;
      if n>n0 
        wp=AA*(Pim/(xzer(n)-xzer(n-1)))^2;
      end
      iter(n)=i;
      rw=wp;
      w=sqrt(abs(rw));
      x=x+Pim/w;
    end 
  else
    if (x_e-x)>0 && rw>0 && err<tol
      n=n+1;
      xzer(n)=x;
      if n>n0 
        wp=AA*(Pim/(xzer(n)-xzer(n-1)))^2;
      end
      iter(n)=i;
      rw=wp;
      w=sqrt(abs(rw));
      x=x+Pim/w;
    end 
  end
end
nit=iter(xzer>=a);
xzer=xzer(xzer>=a);
xzbad=xzbad(xzbad>=a);
if ilamg
 xzer=xzer*lambda;
 xzbad=xzbad*lambda;
end    

function y=fu(icho,lambda,nu,z)
x=z;
lamz=lambda*z;
[J1,Y1,J1p,Y1p,ier1]=bessel(nu,z);
[J3,Y3,J3p,Y3p,ier3]=bessel(nu,lamz);
if icho==1
  if ier3==0
    M1=hypot(J1,Y1);
    N1=hypot(J1p,Y1p);
    Kx=N1/M1;
    M2=hypot(J3,Y3);
    N2=hypot(J3p,Y3p);
    Kl=N2/M2;
    thetax=atan2(Y1,J1);
    phix=atan2(Y1p,J1p);
    thetal=atan2(Y3,J3);
    phil=atan2(Y3p,J3p);
    Delta=thetax-thetal;
    Psi=phil-thetax;
    Psi2=thetal-phix;
    sdelta=sin(Delta);
    y=-x*sdelta/(-sdelta+x*Kx*sin(Psi2)+lambda*x*Kl*sin(Psi));
  else
    M1=hypot(J1,Y1);
    N1=hypot(J1p,Y1p);
    Kx=N1/M1;  
    thetax=atan2(Y1,J1);
    phix=atan2(Y1p,J1p);
    thetal=-pi/2;
    phil=pi/2;
    Kl=nu/lamz-lamz/(2*(nu-1));
    Delta=thetax-thetal;
    Psi=phil-thetax;
    Psi2=thetal-phix;
    sdelta=sin(Delta);
    y=-x*sdelta/(-sdelta+x*Kx*sin(Psi2)+lambda*x*Kl*sin(Psi));
  end  
elseif icho==2
   if ier3==0
    M1=hypot(J1,Y1);
    N1=hypot(J1p,Y1p);
    Kx=N1/M1;
    M2=hypot(J3,Y3);
    N2=hypot(J3p,Y3p);
    Kl=N2/M2;
    thetax=atan2(Y1,J1);
    phix=atan2(Y1p,J1p);
    thetal=atan2(Y3,J3);
    phil=atan2(Y3p,J3p);
    Phi=phix-phil;
    Psi=phil-thetax;
    Psi2=thetal-phix;
    sphi=sin(Phi);
    y=-x*sphi/(sphi-(x/Kx)*(1-nu^2/x^2)*sin(Psi)-...
       lambda*(x/Kl)*(1-nu^2/lamz^2)*sin(Psi2));
  else
    M1=hypot(J1,Y1);
    N1=hypot(J1p,Y1p);
    Kx=N1/M1;  
    thetax=atan2(Y1,J1);
    phix=atan2(Y1p,J1p);
    thetal=-pi/2;
    phil=pi/2;
    Kl=nu/lamz-lamz/(2*(nu-1));
    Phi=phix-phil;
    Psi=phil-thetax;
    Psi2=thetal-phix;
    sphi=sin(Phi);
    y=-x*sphi/(sphi-(x/Kx)*(1-nu^2/x^2)*sin(Psi)-...
       lambda*(x/Kl)*(1-nu^2/lamz^2)*sin(Psi2));
  end  
elseif icho==4
  % J(x)Y'(lambda x)-J'(lambda x)Y(x)=M(x)N(lambda x)sin(Psi),
  % Psi=phi(lambda x)-theta(x)
  M1=hypot(J1,Y1);
  N1=hypot(J1p,Y1p);
  Kx=N1/M1;
  thetax=atan2(Y1,J1);
  phix=atan2(Y1p,J1p);
  if ier3==0
    M2=hypot(J3,Y3);
    N2=hypot(J3p,Y3p);
    Kl=N2/M2;
    thetal=atan2(Y3,J3);
    phil=atan2(Y3p,J3p);
  else
    thetal=-pi/2;
    phil=pi/2;
    Kl=nu/lamz-lamz/(2*(nu-1));
  end
  Psi=phil-thetax;
  ca=cos(phil-thetal);
  sa=sin(phil-thetal);
  q=1-nu^2/lamz^2;
  A=Kx*cos(phix-thetax)+lambda*(-1/lamz-q*ca/Kl);
  dPsi=lambda*q*sa/Kl-Kx*sin(phix-thetax);
  sp=sin(Psi);
  y=sp/(A*sp+dPsi*cos(Psi));
else
  if ier3==0
    M1=hypot(J1,Y1);
    N1=hypot(J1p,Y1p);
    Kx=N1/M1;
    M2=hypot(J3,Y3);
    N2=hypot(J3p,Y3p);
    Kl=N2/M2;
    thetax=atan2(Y1,J1);
    phix=atan2(Y1p,J1p);
    thetal=atan2(Y3,J3);
    phil=atan2(Y3p,J3p);
    Delta=thetax-thetal;
    Phi=phix-phil;
    Psi=phil-thetax;
    Psi2=thetal-phix;
    sphi=sin(Phi);
    spsi2=sin(Psi2);
    y=spsi2/(-lambda*Kl*sphi+(1/Kx)*(1-nu^2/x^2)*sin(Delta));
  else
    M1=hypot(J1,Y1);
    N1=hypot(J1p,Y1p);
    Kx=N1/M1;  
    thetax=atan2(Y1,J1);
    phix=atan2(Y1p,J1p);
    thetal=-pi/2;
    phil=pi/2;
    Kl=nu/lamz-lamz/(2*(nu-1));
    Delta=thetax-thetal;
    Phi=phix-phil;
    Psi=phil-thetax;
    Psi2=thetal-phix;
    sphi=sin(Phi);
    spsi=sin(Psi);
    spsi2=sin(Psi2);
    y=spsi2/(-lambda*Kl*sphi+(1/Kx)*(1-nu^2/x^2)*sin(Delta));
  end  
end 
