function [xzer, xzbad] = besselzd(a, b, pnu, gamm, alpha, eps)
%---------------------------------------------------------------------
% Computation of the zeros of the equation
%
%   x C'_{nu,alpha}(x) + gamm C_{nu,alpha}(x) = 0
%
% in the interval (a,b), with a > 0.
%----------------------------------------------------------------------
% Inputs:
%
%   a       lower endpoint of the interval
%   b       upper endpoint of the interval
%   pnu     order nu of the Bessel functions
%   gamm    parameter in the equation
%           x C'_{nu,alpha}(x) + gamm C_{nu,alpha}(x) = 0
%   alpha   parameter defining the linear combination
%           C_{nu,alpha}(x) = cos(alpha)*J_nu(x) - sin(alpha)*Y_nu(x)
%   eps     required accuracy for the computed zeros
%
% Outputs:
%
%   xzer    vector containing all computed zeros of
%           x C'_{nu,alpha}(x) + gamm C_{nu,alpha}(x) = 0 in (a,b)
%           
%   xzbad   vector containing zeros of
%           x C'_{nu,alpha}(x) + gamm C_{nu,alpha}(x) = 0 in (a,b)
%           for which the required accuracy is not achieved (typically in
%           the case of nearly degenerate roots)
%-------------------------------------------------------------------------
% Accompanying software of the paper 
% "Numerical Software for Bessel functions and Associated
%  Values" 
% Authors: Amparo Gil, Javier Segura and Nico M. Temme
%--------------------------------------------------------------------------
talpha=tan(alpha);
xcert=zeros(1,1000);
itct=zeros(1,1000);
xzbad=zeros(1,3);
itbad=zeros(1,3);
pnu2=pnu*pnu;
gamm2=gamm*gamm;
[co, li, ui, signo] = intervals(a, b, pnu, gamm, talpha);
numt=0;
for i=1:co
  num=0;    
  if li(i)<ui(i)
    [xcer_temp, itc_temp, num]=...
            besjz(signo(i), li(i), ui(i), ...
                  pnu, gamm, talpha, eps);
    for j=numt+1:num+numt
      xcert(j)=xcer_temp(j-numt);
      itct(j)=itc_temp(j-numt);
    end
  end
  numt=numt+num;
end
if abs(pnu)>abs(gamm)  
  nuga=sqrt(pnu2-gamm2);
  if a < nuga      
    xt1=nuga*(1+2e-16); 
    xt2=nuga*(1-2e-16);    
 end
end
com=0;
i=0;
while i<numt
  i=i+1;  
  sip=xcert(i)^2-pnu2+gamm2;  
  if abs(sip*hfu(pnu, talpha, gamm, xcert(i))) > eps      
    com=com+1;
    xzbad(com)=xcert(i);
    itbad(com)=itct(i);
    for j = i:numt-1
      xcert(j)=xcert(j+1);
      itct(j)=itct(j+1);
    end     
    i=i-1;
    numt=numt-1;
  end
end
xcert=xcert(1:numt);
itct=itct(1:numt);
xzbad=xzbad(1:com);
itbad=itbad(1:com);
xzer=sort([xcert xzbad],'ascend');
end

function [co, li, ui, signo]=intervals(ai, bi, nu, g, talpha)
li=zeros(1,7);
ui=zeros(1,7);
signo=zeros(1,7);
xt=zeros(1,7);
nogo=0;
modv=0;     
nu2=nu^2;
lambda=nu2-g^2;
if lambda~=0
  si=lambda/abs(lambda);
  pi_=acos(-1);
  d=(3/4-2*g+nu2);
  a=si*(2*g-3*nu2+23/4)/d;
  b=3*(nu2-0.25)/d;
  c=si*(0.25-nu2)/d;
  q=(a^2-3*b)/9;
  r=(2*a^3-9*a*b+27*c)/54;
  delta=r^2-q^3;
  i=0;
  if delta > 0
    Ad=-r/abs(r)*(abs(r)+sqrt(delta))^(1/3);
    x=(Ad+q/Ad)-a/3;
    if x > 0
      i=i+1;
      xt(i)=x;
    end
  else
    the=acos(r/q^(3/2));
    sqq=sqrt(q);
    x1=-2*sqq*cos(the/3)-a/3;
    if x1>0
      i=i+1;
      xt(i)=x1;
    end
    x2=-2*sqq*cos((the+2*pi_)/3)-a/3;
    if x2>0
      i=i+1;
      xt(i)=x2;
    end
    x3=-2*sqq*cos((the-2*pi_)/3)-a/3;
    if x3 >0
      i=i+1;
      xt(i)=x3;
    end
  end
  alam=abs(lambda);
  al=ai^2/alam;
  bl=bi^2/alam;
  if lambda>0
    if modv==0
      i=i+1; xt(i)=1+1e-15;
      i=i+1; xt(i)=1-1e-15;
    elseif al<1
      if bl>1
        al=1+1e-15;
      else
        nogo=1;
      end
    end
  end
  i=i+1; xt(i)=al;
  i=i+1; xt(i)=bl;
  if nogo==0
    for ind=i:-1:2
      for k=1:ind-1
        if xt(k)>xt(k+1)
          xpp=xt(k);
          xt(k)=xt(k+1);
          xt(k+1)=xpp;
        end
      end
    end
    ja=1;
    if xt(1)<al
      while xt(ja)<al
        ja=ja+1;
      end
    end
    jb=ja;
    if xt(i)>bl
      while xt(jb)<bl
        jb=jb+1;
      end
    else
      jb=i;
    end
    co=0;
    for k=ja:jb-1
      if abs(1-xt(k)/xt(k+1))>1e-9
        co=co+1;
        li(co)=xt(k);
        ui(co)=xt(k+1);
      end
    end
    for j=1:co
      x=(li(j)+ui(j))/2;
      si2=(x-si)*d*(x^3+a*x^2+b*x+c);
      signo(j)=sign(si2);
      li(j)=sqrt(alam*li(j));
      ui(j)=sqrt(alam*ui(j));
    end
  else
    co=0;
  end
else
  co=1;
  li(1)=ai;
  ui(1)=bi;
  signo(1)=0;
  if nu==g
    signo(1)=(nu-1)^2-0.25;
  else
    signo(1)=(nu+1)^2-0.25;
  end
  if signo(1)~=0
    signo(1)=signo(1)/abs(signo(1));
  else
    signo(1)=1;
  end
end
if talpha~=0
  if li(1)<1e-6
    li(1)=1e-6;
  end
else
  li(1)=max(li(1), 1e-50);
end
li=li(1:co);
ui=ui(1:co);
signo=signo(1:co);
end

function [xcer, itc, num]=besjz(j, a, b, pnu, gamm, talpha, eps)
xcer=zeros(1,1000);
itc=zeros(1,1000);
dpi=pi;
if a>=b
  error('a should be smaller than b')
end
if floor(pnu)==pnu
  pnu=abs(pnu);
end
if j==1
  xm=a;
  xc=b;
else
  xm=b;
  xc=a;
end
i=0;
omega2=coef(pnu, gamm, xc);
omegam=coef(pnu, gamm, xm);
if j>0
  omegab=omega2;
else
  omegab=omegam;
end
ss=omega2*omegam;
if ss<0
  k=0;
elseif omegab <= 0
  k=-1;
else
  k=1;
end
while (j*(xc - xm) > 0)&&(k >= 0)
  iter=0;
  dev=eps+1;
  isal=0;
  hfut=1;
  sip=1;
  while (dev>eps) && ...
        (j*(xc - xm) > 0) && ...
        (k >= 0) && ...
        (abs(sip*hfut) > 1e-19)
    iter=iter+1;
    if iter>1000
      error('convergence failure')
    end
    sip=(xc^2-pnu^2+gamm^2);
    hfut=hfu(pnu, talpha, gamm, xc);
    si=sip/abs(sip);
    h=sqrt(omep(pnu, gamm, xc))*si*hfut;
    omega=sqrt(omega2);
    atast=atan(h);
    dest=atast/omega;
    dev=abs(dest)/xc;
    if isal==1
      dev=0.5*eps;
    end
    if log(dev) < 0.25*log(eps)
      isal=1;
    end
    if j*dest < 0
      if dev>1e-10
        dest=dest+j*dpi/omega;
      end
    end
    xc=xc-dest;
    omega2=coef(pnu, gamm, xc);
    if omega2 <= 0
      k=-1;
    end
  end
  if (j*(xc - xm) > 0)&&(k >= 0)
    i=i+1;
    xcer(i)=xc;
    itc(i)=iter;
    xc=xc-j*dpi/omega;
    omega2=coef(pnu, gamm, xc);
    if omega2 <= 0
      k=-1;
    end
  end
end
if k < 0
  h1t=hfu(pnu, talpha, gamm, xm);
  h2t=hfu(pnu, talpha, gamm, xc);
  si=(xc^2-pnu^2+gamm^2);
  si=si/abs(si);
  h=sqrt(-omep(pnu, gamm, xc))*si*h2t;
  if (j*h > 0) && (h1t*h2t < 0)
    iter=0;
    dev=eps + 1;
    isal=0;
    while (dev > eps) && ...
          (abs(h) < 1) && ...
          (j*(xc - xm) > 0)
      iter=iter+1;
      if iter > 1000
        error('convergence failure')
      end
      si=(xc^2 - pnu^2 + gamm^2);
      si=si/abs(si);
      h=sqrt(-omep(pnu, gamm, xc)) * ...
        si*hfu(pnu, talpha, gamm, xc);
      omega=sqrt(-omega2);
      if abs(h) < 1
        atast=atanh(h);
        dest=atast/omega;
        dev=abs(dest)/xc;
        if isal == 1
          dev=0.5*eps;
        end
        if log(dev) < 0.5*log(eps)
          isal=1;
        end
        xc=xc-dest;
        omega2=coef(pnu, gamm, xc);
      end
    end
    if dev<eps
      i=i+1;
      xcer(i)=xc;
      itc(i)=iter;
    end
  end
end
num=i;
xcer=xcer(1:num);
itc=itc(1:num);
end

function val=hfu(nu, talpha, gamm, x)
x2=x*x;
n2=nu*nu;
g2=gamm*gamm;
d=nu+gamm;
e=-x;
m=-x2^2 ...
    +d*(2*nu - gamm - 0.5)*x2 ...
    -d*(n2 - g2)*(nu + 0.5);
n=x*((0.5 - gamm)*x2 ...
     +(n2 - g2)*(gamm + 0.5));
r=fc(x, nu + 1.0, talpha);
val=(d + e*r) / (m + n*r);
end

function val=coef(nu, gamm, x)
x2=x*x;
g2=gamm*gamm;
n2=nu*nu;
val=omep(nu, gamm, x)/ ...
      (x2*(x2-n2+g2)^2);
end

function val=omep(nu, gamm, x)
n2=nu*nu;
g2=gamm*gamm;
x2=x*x;
P=(-3*n2+2*g2+2*gamm-0.75);
Q=(3*n2-g2-2*gamm-2.5)*(n2-g2);
R=(0.25-n2)*(n2-g2)^2;
val=x2^3+P*x2^2+Q*x2+R;
end

function val=fc(x, pnu, talpha)
% ---------------------------------------------------------
% fc = (cos(alpha)*J(pnu,x) - sin(alpha)*Y(pnu,x)) /
%      (cos(alpha)*J(pnu-1,x) - sin(alpha)*Y(pnu-1,x))
% ---------------------------------------------------------
if pnu<=0
  pnum=abs(pnu);
  pnumm=abs(pnu) + 1;
elseif pnu < 1
  pnum=pnu;
  pnumm=abs(pnu - 1);
else
  pnum=pnu;
  pnumm=pnu-1;
end
% --- Bessel J, Y ---
[fj1,fy1,~,~,~]=bessel(pnum,x);
[fj2,fy2,~,~,~]=bessel(pnumm,x);
if pnu <= 0
    fj1f=cos(pnum*pi)*fj1-sin(pnum*pi)*fy1;
    fy1f=sin(pnum*pi)*fj1+cos(pnum*pi)*fy1;
    fj2f=cos(pnumm*pi)*fj2-sin(pnumm*pi)*fy2;
    fy2f=sin(pnumm*pi)*fj2+cos(pnumm*pi)*fy2;
elseif pnu < 1
    fj1f=fj1;
    fy1f=fy1;
    fj2f=cos(pnum*pi)*fj2 - sin(pnum*pi)*fy2;
    fy2f=sin(pnum*pi)*fj2 + cos(pnum*pi)*fy2;
else
    fj1f=fj1;
    fy1f=fy1;
    fj2f=fj2;
    fy2f=fy2;
end
c=1.0;
talphai=talpha;
if abs(talpha)>1e250
  c=0.0;
  talphai=1.0;
end
val=(c*fj1f-talphai*fy1f)/ ...
    (c*fj2f-talphai*fy2f);
end



