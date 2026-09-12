function xzer=crosspz(icho,pnu,lambda,a,b,eps)
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
%   pnu      order nu of the Bessel functions
%   lambda   scaling parameter for the argument, 0<lambda
%   a        lower endpoint of the interval
%   b        upper endpoint of the interval
%   eps      requested accuracy for the computed zeros
%
% OUTPUTS:
%
%   xzer      array containing the computed zeros
%--------------------------------------------------------------------
% Accompanying software of the paper 
% "Numerical Software for Bessel functions and Associated
%  Values" 
% Authors: Amparo Gil, Javier Segura and Nico M. Temme
%----------------------------------------------------------------
if lambda==1
  xzer=[];
  return
end    
ilamg=0;
if lambda>1
  lambda=1/lambda;
  ilamg=1;
end    
Pim=pi-eps;
itmax=50;
if b<=pnu
  xzer=[];
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
  h=fu(icho,lambda,pnu,x);
  itp=0;
  itp2=0;
  err=eps+1;
  i=0;
  while err>eps && (x_e-x)>0 && rw>0 && i<itmax
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
  itp=i;
  if i==itmax
    if (x_e-x)>0 && rw>0 
      n=n+1;
      xzer(n)=x;    
      if n==10 
         AA=1;
      end   
      if n>n0 
        wp=AA*(Pim/(xzer(n)-xzer(n-1)))^2;
      end
      iter(n)=i;
      itp=0;
      h=fu(icho,lambda,pnu,x);              
      rw=wp;
      w=sqrt(abs(rw));
      x=x+Pim/w;
    end 
  else
    if (x_e-x)>0 && rw>0 && err<eps
      n=n+1;
      xzer(n)=x;
      if n>n0 
        wp=AA*(Pim/(xzer(n)-xzer(n-1)))^2;
      end
      iter(n)=i;
      itp=0;
      h=fu(icho,lambda,pnu,x);              
      rw=wp;
      w=sqrt(abs(rw));
      x=x+Pim/w;
    end 
  end
end
if n>0
  if(xzer(1)<a)
     if(n>1)
       xzer=xzer(2:end);  
     else
       xzer=[];
     end    
     n=n-1;
  end    
else
  xzer=[];
end
if ilamg
 xzer=xzer*lambda;
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