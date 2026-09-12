function xzer=besselz(a, b, pnu, alpha, eps)
%--------------------------------------------------------------------
% Computation of the zeros of the general cylinder function
%
%   C_{nu,alpha}(x) = cos(alpha)*J_nu(x) - sin(alpha)*Y_nu(x)
%
% in the interval (a,b), with a > 0.
%---------------------------------------------------------------------
% Inputs:
%
%   a       lower endpoint of the interval
%   b       upper endpoint of the interval
%   pnu     order nu of the Bessel functions
%   alpha   parameter defining the linear combination
%           C_{nu,alpha}(x) = cos(alpha)*J_nu(x) - sin(alpha)*Y_nu(x)
%   eps     required accuracy for the computed zeros
%
% Outputs:
%
%   xzer    vector containing all zeros of C_{nu,alpha}(x) in (a,b)
%           computed within the prescribed accuracy
%-----------------------------------------------------------------------
% Accompanying software of the paper 
% "Numerical Software for Bessel functions and Associated
%  Values" 
% Authors: Amparo Gil, Javier Segura and Nico M. Temme
%------------------------------------------------------------------------
dpi=pi;
tiny=10*realmin;
if a>=b
  error('a should be smaller than b');
end
if floor(pnu)==pnu
  pnu=abs(pnu);
end
calpha=cos(alpha);
salpha=sin(alpha);
pnuu=pnu^2;
pnu2=pnuu-0.25;
if pnu2~=0
  j=sign(pnu2);
else
  j=1;
end
if j==1
  xm=a;
  xc=b;
else
  xm=b;
  xc=a;
end
i=0;
omega2=1-(pnuu-0.25)/xc^2;
omegam=1-(pnuu-0.25)/xm^2;
if j > 0
  omegab=omega2;
else
  omegab=omegam;
end
ss=omega2*omegam;
if ss < 0
  k=0;
elseif omegab <= 0
  k=-1;
else
  k=1;
end
while (j*(xc-xm)>0)&&(k>=0)
  iter=0;
  dev=eps + 1;
  isal=0;
  h=1;
  while (dev>eps)&&(j*(xc-xm)>0)&&(k>=0) ...
          && (abs(h)/xc>1e-19)
    if iter > 1000
      error('convergence failure');
    end
    [J1,Y1,Jp1,Yp1,~]=bessel(pnu,xc);
    Cnu=calpha*J1-salpha*Y1; 
    if abs(Cnu)<tiny
      h=0;
    else
      Cnup=calpha*Jp1-salpha*Yp1;
      h=1/(1/(2*xc)+Cnup/Cnu);
    end
    omega=sqrt(omega2);
    atast=atan(omega*h);
    dest=atast/omega;
    dev=abs(dest)/xc;
    if j*dest < 0
      if dev > 1e-10
        dest=dest+j*dpi/omega;
      end
    end
    if isal==1
      dev=0.5*eps;
    end
    if log(dev)<0.25*log(eps)   
      isal=1;
    end
    xc=xc-dest;
    iter=iter+1;
    omega2=1-(pnuu-0.25)/xc^2;
    if omega2 <= 0
      k=-1;
    end
  end
  if (j*(xc-xm)>0)&&(k>=0)
    i=i+1;
    xzer(i)=xc;
    itc(i)=iter;
    xc=xc-j*dpi/omega;
    omega2=1-(pnuu - 0.25)/xc^2;
    if omega2 <= 0
      k = -1;
    end
  end
end
if k < 0
  [J1,Y1,Jp1,Yp1,~]=bessel(pnu,xm);
  Cnu=calpha*J1-salpha*Y1;
  if abs(Cnu)<tiny
    hm=0;
  else
    Cnup=calpha*Jp1-salpha*Yp1;
    hm=1/(1/(2*xm)+Cnup/Cnu);
  end  
  [J2,Y2,Jp2,Yp2,~]=bessel(pnu,xc);
  Cnub=calpha*J2-salpha*Y2;
  if abs(Cnub)<tiny
    h=0;
  else
    Cnupb=calpha*Jp2-salpha*Yp2;
    h=1/(1/(2*xc)+Cnupb/Cnub);
  end  
  if (j*h > 0) && (h*hm < 0)
    iter=0;
    dev=eps+1;
    omega=sqrt(-omega2);
    isal=0;
    while (dev > eps)&&(abs(h*omega) < 1) ...
              &&(j*(xc-xm) > 0)
      iter=iter+1;
      if iter>1000
        error('convergence failure');
      end
      if abs(h*omega)<1
        atast=0.5*log((1+omega*h)/(1-omega*h));
        dest=atast/omega;
        dev=abs(dest)/xc;
        if isal==1
          dev=0.5*eps;
        end
        if log(dev)<0.25*log(eps)
          isal=1;
        end
        xc=xc-dest;
        [J1,Y1,Jp1,Yp1,~]=bessel(pnu,xc);
        Cnu=calpha*J1-salpha*Y1; 
        if abs(Cnu)<tiny
          h=0;
        else
          Cnup=calpha*Jp1-salpha*Yp1;
          h=1/(1/(2*xc)+Cnup/Cnu);
        end
        omega2=1-(pnuu - 0.25)/xc^2;
        omega=sqrt(-omega2);
      end
    end
    if dev<eps
      i=i+1;
      xzer(i)=xc;
    end
  end
end
if j==1
  xzer=xzer(end:-1:1);
end
