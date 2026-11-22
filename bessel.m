function [J,Y,Jp,Yp,ier]=bessel(a,x)
%----------------------------------------------------
% Computation of the Bessel functions J_a(x), Y_a(x)
% and their derivatives for real orders (a) 
% and real positive arguments (x).
%----------------------------------------------------
% Inputs:
%   a ,    order of the Bessel functions
%   x,     argument of the Bessel functions.
% Outputs:
%   J,     function J_a(x)
%   Jp,    first order derivative of J_a(x)
%   Y,     function Y_a(x)
%   Yp,    first order derivative of Y_a(x)
%   ier, error flag:
%        ier=0, computation successful
%        ier=1, underflow/overflow values
%               if a>0 
%                 J, Jp are set to zero
%                 Y is set to -realmax, Yp to realmax
%               if a<0
%                 J=-s*realmax, Jp=s*realmax, s=sign(sin(a*pi))
%                 Y=-c*realmax, Yp=c*realmax, c=sign(cos(a*pi))
%        ier=2, the x value is out of range.
%               All the function values are set 
%               to zero.   
%-------------------------------------------------------
% Accompanying software of the paper 
% "Numerical Software for Bessel functions and Associated
%  Values" 
% Authors: Amparo Gil, Javier Segura and Nico M. Temme
%--------------------------------------------------------
ier=0;
tiny=realmin*10;
if x<=0 || isnan(a) || isnan(x)
  ier=2;
  J=0;
  Jp=0;
  Y=0;
  Yp=0;
else
  aa=abs(a);
  if aa == 0.5
  % Special case a = +-0.5
    s=sin(x);
    c=cos(x);
    r=sqrt(2/(pi*x));
    if a < 0
      jax=c*r;
      jap=a*jax/x-r*s;
      yax=s*r;
      yap=a*yax/x+r*c;
    else
      jax=r*s;
      jap=r*c-a*jax/x;
      yax=-c*r;
      yap=s*r-a*yax/x;
    end
    J=jax;
    Jp=jap;
    Y=yax;
    Yp=yap;
  elseif a < 0
    [s,c]=sincospi(a);
    [ia,yax,iad,yap,ier]=bessel(-a, x);
    if ier==0  
      J=c*ia+s*yax;
      Jp=c*iad+s*yap;
      Y=-s*ia+c*yax;
      Yp=-s*iad+c*yap;
    else
      ss=sign(s);
      cs=sign(c);
      J=-ss*realmax;
      Jp=ss*realmax;
      Y=-cs*realmax;
      Yp=cs*realmax;
    end
  else
    % power series for J if x<xmax (for a<=45, the recurrence region
    % below is more accurate for x>2*sqrt(1+a))
    if a<=45
      xmax=2*sqrt(1+a);
    else
      xmax=3*sqrt(1+a);
    end
    a1=1.5;
    x1=0.3*a;
    x2=a1*a;
    if a==0
      nunder=0;
    else
      z=x/a;
      if z<0.001
        logrmin = log(tiny);
        nunder = 2*exp((logrmin+gammaln(a+1))/a);
      elseif z<0.95
        eta = sqrt(1-z^2) - acosh(1/z);
        logJ = a*eta - 0.5*log(2*pi*a) - 0.25*log(1-z^2);
        logrmin = log(tiny);
        if logJ<logrmin
          nunder=x+1;
        else
          nunder=0;
        end
      else
        nunder=0;
      end    
    end
    if x>nunder
      if x<xmax 
        [Y,Yp]=bessyax(a, x);
        [J,Jp]=besseries(a,x);
%--------------------------------------------------------------
% Uniform Airy expansion
%--------------------------------------------------------------
      elseif (a>40) && (x1 <= x) && (x <= x2)
         [J,Y,Jp,Yp]=bessair(a,x);
%--------------------------------------------------------------
% Debye, x<a
%--------------------------------------------------------------
      elseif (a>45) && (x < 0.75*a)
        [J,Jp,Y,Yp,ier]=deb1(a,x);
%--------------------------------------------------------------
% Hankel (large x)
%--------------------------------------------------------------
      elseif (x > 1.9*(a/5.3)^2 + 20)&&((a<=20)||(x>a*a/6))
        [J,Y,Jp,Yp,ier]=expax(a,x);
%--------------------------------------------------------------
% Debye, x>a
%--------------------------------------------------------------
      elseif ((a>45)&&(x>1.35*a))||((a>20)&&(x>1.9*(a/5.3)^2+20))
        [J,Jp,Y,Yp]=deb2(a,x);
%--------------------------------------------------------------
% Recurrences for J
%--------------------------------------------------------------      
      else
        [Y,Yp]=bessyax(a,x);
        [J,Jp]=jrec(a,x);
        % Normalize J with the Wronskian J*Yp-Jp*Y=2/(pi*x)
        if J~=0
          rho=Jp/J;
          J=2/(pi*x)/(Yp-rho*Y);
          Jp=rho*J;
        end
      end  
    else  
      ier=1;
      J=0;
      Jp=0;
      Y=-realmax;
      Yp=realmax;   
    end	
  end
end

% Overflow of Y or Y' not detected before
if (ier==0)&&~(isfinite(J)&&isfinite(Jp)&&isfinite(Y)&&isfinite(Yp))
  ier=1;
  J=fixov(J); Jp=fixov(Jp); Y=fixov(Y); Yp=fixov(Yp);
end

function [jax, jap, yax, yap,ier]=deb1(a, x)
% Debye-type expansion a > x
y=x/a;
t2=1/(1-y^2);
t=sqrt(t2);
eta=log((1+t)/(t*y))-1/t;
if log(a*eta)<log(realmax)
  ier=0;  
  vv=exp(a*eta);
  c1=1;
  c2=0;
  c3=0;
  jj=0;
  jjj=0;
  u=t/(8*a);
  w=4*u;
  z=sqrt(t/(2*pi*a));
  zz=z/(y*t);
  k=1;
  m=1;
  n=-3;
  jax=1;
  jap=1;
  yax=1;
  yap=1;
  sig=1;
  while true
    v=m*c1-n*c3*t2;
    ck=m*u*v/k;
    if abs(c3) < abs(ck)
      jj=jj+10;
    else
      jj=0;
    end
    if abs(ck) < eps
      jjj=jjj+100;
    else
      jjj=0;
    end
    c3=c2;
    c2=c1;
    c1=ck;
    k=k+1;
    m=m+2;
    n=n+2;
    sig=-sig;
    v=ck-w*v;
    jax=jax+ck;
    jap=jap+v;
    yax=yax+sig*ck;
    yap=yap+sig*v;
    if (jj == 30) || (jjj == 400)
      break;
    end
  end
  jax=z*jax/vv;
  yax=-2*z*yax*vv;
  jap=zz*jap/vv;
  yap=2*zz*yap*vv;
else
  ier=1;
  jax=0;
  yax=-realmax;
  jap=0;
  yap=realmax;
end  

function [jax,jap,yax,yap]=deb2(a,x)
% Debye-type expansion x > a
a2=a*a;
x2=x*x;
w2= a2/(x2-a2);
ww=sqrt((x-a)*(x+a));
chi = ww - a*atan(ww/a) - pi/4;
c=sqrt(2.0/pi);
k2=sqrt(ww);
k1=c/k2;
k2=c*k2/x;
e1= 1.0;
e2=0.0;
e3=0.0;
jj=0;
jjj=0;
u=1.0/(8.0*ww);
w=4.0*u;
k=1;
m=1;
n=-3;
p=1.0;
q=0.0;
r=1;
s=0;
sig=1;
t=1.0;
while (jj~=30)&&(jjj~=400)
  v=n*e3*w2-m*e1;
  ek=m*u*v/k;
  if abs(e3)< abs(ek) 
    jj=jj+10;
  else
    jj=0;
  end   
  if abs(ek)<eps 
    jjj=jjj+100;
  else
    jjj=0;
  end   
  e3=e2;
  e2=e1;
  e1=ek;
  fk=ek-w*v;
  if (sig==1) 
    q=q+t*ek;
    s=s+t*fk;
    t=-t;
  else
    p=p+t*ek;
    r=r+t*fk;
  end
  k=k+1;
  m=m+2;
  n=n+2;
  sig=-sig;
end
c=cos(chi);
u=sin(chi);
jax=k1*(p*c-q*u);
yax=k1*(p*u+q*c);
jap=-k2*(r*u+s*c);
yap=k2*(r*c-s*u);

function [J,Jp]=besseries(a,x)
% Power series
eps=1.e-16;
if a==0
  % J0
  r=x*x/4.0;
  c=1.0;
  J=1.0;
  n=1;
  while abs(c)>eps
    c=-c*r/(n*n);
    J=J+c;
    n=n+1;
  end
  % J0p=-J1
  [j1,~]=besseries(1,x);
  Jp=-j1;
  return
else
  % (x/2)^a/gamma(a+1); the power is used when it does not
  % underflow (more accurate than the exponential form for small x)
  z=x/(2.0*a);
  if (a<700)&&(a*log(z)>-650)
    ia=z^a*exp(a)/(sqrt(2.0*pi*a)*gamstar(a));
  else
    ia=exp(a*(log(z)+1.0))/(sqrt(2.0*pi*a)*gamstar(a));
  end
  r=x*x/4.0;
  c=1.0;
  s=a;
  n=1;
  jax=0.0;
  jap=0.0;              
  % relative stopping criterion for both series
  while true
    jax=jax+c;
    jap=jap+s;
    c=-c*r/(n*(a+n));
    s=(a+2.0*n)*c;
    n=n+1;
    if (abs(c)<=eps*abs(jax))&&(abs(s)<=eps*abs(jap))
      break;
    end
  end 
  J=jax*ia;
  Jp=jap*ia/x;
end
  
function [jax,yax,jap,yap,ier]=expax(a,x)
% Asymptotic expansions, x large
ier=0;
piquart=0.78539816339744830962;
pihalf=pi/2;
giant=realmax;
m=0;
x2=x*x;
% phi=-(a/2+1/4)*pi; exact reduction of a/2 before adding 1/4
t=a/2-2*round(a/4)+0.25;
[sp,cp]=sincospi(t);
cphi=cp;
sphi=-sp;
cosw=cos(x)*cphi-sin(x)*sphi;
sinw=sin(x)*cphi+cos(x)*sphi;
fact=sqrt(1.0/(pihalf*x));
fnus=4.0*a*a;
err1=1.0; 
err1p=1.0;
y1=0.0; 
y1p=0.0;
k=0;
a2=1.0;
x2k=1.0;
x2kp=x;
while ((err1 > eps)||(err1p >eps))&&(m<500)
  acof=a2/x2k;
  k=k+1;
  l=2*k;
  a2=(fnus-(l-1.0)*(l-1.0))*a2/(8.0*k);   
  acofd=a2/x2kp;
  k=k+1;
  l=2*k;
  a2=(fnus-(l-1.0)*(l-1.0))*a2/(8.0*k);
  x2k=-x2*x2k;
  x2kp=-x2*x2kp;
  y1=y1+acof; 
  y1p=y1p+acofd;    
  err1=abs(acof)/(abs(y1)+eps);
  err1p=abs(acofd)/(abs(y1p)+eps);
  if (abs(a2)>giant) 
    m=500;
  end
end
if (m==500) 
  ier=1;
  jax=0.0; yax=-realmax;
  jap=0.0; yap=realmax;
else
  jax=fact*(cosw*y1-sinw*y1p);
  yax=fact*(sinw*y1+cosw*y1p);
end
if (ier==0) 
  err1=1.0; 
  err1p=1.0;
  y1=1.0; 
  a2=(fnus+3.0)/8.0;
  x2k=1.0;
  x2kp=x;
  y1p=a2/x2kp;
  a2=(fnus-1.0)*(fnus+15.0)/128.0;
  x2k=-x2*x2k;
  x2kp=-x2*x2kp;
  k=2;
  while ((err1 > eps)||(err1p >eps))&&(m<500)
    acof=a2/x2k;
    k=k+1;
    l=2*k;
    a2=(fnus-(l-3.0)*(l-3.0))*a2*(fnus+4*k*k-1.0)...
        /((8.0*k)*(fnus+4*(k-1)*(k-1)-1.0));
    acofd=a2/x2kp;
    k=k+1;
    l=2*k;
    a2=(fnus-(l-3.0)*(l-3.0))*a2*(fnus+4*k*k-1.0)...
        /((8.0*k)*(fnus+4*(k-1)*(k-1)-1.0));
    x2k=-x2*x2k;
    x2kp=-x2*x2kp;
    y1=y1+acof; 
    y1p=y1p+acofd;    
    err1=abs(acof/y1); 
    err1p=abs(acofd/y1p);
    if abs(a2)>giant 
      m=500;
    end
  end
  if m==500 
    ier=1;
    jax=0.0; yax=-realmax;
    jap=0.0; yap=realmax;
  else
    jap=-fact*(sinw*y1+cosw*y1p);
    yap=fact*(cosw*y1-sinw*y1p);
  end
end   

function  [yax, yap]=bessyax(a, x)
if (a <-0.5)     
  [f,g]=bessjax(-a, x);
  [p,q]=bessyax(-a, x);
  [s,c]=sincospi(a);
  yax=-s*f+c*p;
  yap=-s*g+c*q;
elseif (abs(a)==0.5) 
  s=sin(x);
  c=cos(x);
  t=sqrt(2.0/(pi*x));
  if (a < 0) 
    yax=s*t;
    yap=a*yax/x+t*c;   
  else
    yax=-c*t;
    yap=s*t-a*yax/x;
  end
else
  if (x<3)
    na=floor(a+0.5);
    r=a-na;  
    b=x/2.0;
    d=-log(b);
    e=r*d;
    c=r*pi;
    if (abs(c)<1.0e-5) 
      c=1.0 + c*c/6.0;
    else
      c=c/sin(c);
    end   
    c=c/pi;
    if (abs(e)<1.0e-5) 
      q=e*e;
      s=1.0+q/6.0;
      q=1.0+q/2.0;
    else
      [p,q]= hypfun(e);
      s=p/e;
    end
    e=exp(e);
    [gg,p,t]=recipgam(-r);
    g=e*gg;
    e=r*r;
    f=2.0*c*(-p*q+(1.0+e*t)*s*d);
    p=g*c;
    q=1.0/(g*pi);
    t=r*pi/2.0;
    if (abs(t)<1.0e-5) 
      s=1.0-t*t/6.0;
    else
      s=sin(t)/t;
    end   
    s=pi*t*s*s;
    c=1.0;
    d=-b*b;
    yax=f+s*q;
    ya1=p;
    n=1;
    h=ya1;
    g=yax;
    while ((abs(h/ya1)+ abs(g/yax))> eps)
      f=(f*n+p+q)/(n*n-e);
      c=c*d/n;
      p=p/(n-r);
      q=q/(n+r);
      g=c*(f+s*q);
      h=c*p-n*g;
      yax=yax+g;
      ya1=ya1+h;
      n=n+1;
    end   
    f=-yax;
    g=-ya1/b;
    b=2.0/x;
    fm=f;
    for n=1:na       
      h=b*(r+n)*g-f;
      fm=f;
      f=g;
      g=h;
    end
  else
    na=floor(a);
    r=a-na;
    b=x-pi*(r+0.5)/2.0;
    c=cos(b);
    s=sin(b);
    d=sqrt(2.0/(pi*x));
    [p, q, b, h]= besspqax(r, x);
    f=d*(p*s+q*c);
    g=d*(b*c-h*s);
    g=r*f/x-g;
    b=2.0/x;
    fm=f;
    for n=1:na       
      h=b*(r+n)*g-f;
      fm=f;
      f=g;
      g=h;
    end
  end  
  yax=f;
  if (na>=1)&&~isfinite(g)
    % Y_{a+1} overflows: use Y'_a=Y_{a-1}-a Y_a/x
    yap=fm-a*f/x;
  else
    yap=a*f/x-g;
  end
end 
   
function [pa, qa, ra, sa]= besspqax(a, x)  
r=abs(a);
na=floor(a);
r=r-na;
if (r==0.5) 
  pa=1.0;
  pa1=1.0;
  qa=0.0;
  qa1=1.0/x;
elseif (x >= 3) 
  c=0.25-r*r;
  b=x+x;
  f=1.0;
  g=1.0;
  p=1.0;
  q=0.0;
  % Number of terms of the backward recurrence (empirical, x>=3)
  n=floor(10+360/x);
  while (n > 0)
    t=(n+1.0)*(2.0-p)-2.0;
    s=b+(n+1)*q;
    d=(n-1.0+c/n)/(s*s+t*t);
    p=d*t;
    q=d*s;
    e=f;
    f=p*(e+1.0)-g*q;
    g=q*(e+1.0)+g*p;
    n=n-1.0;
  end   
  f=1.0+f;
  d=f*f+g*g;
  pa=f/d;
  qa=-g/d;
  d=r+0.5-p;
  q=q+x;
  pa1=(pa*q-qa*d)/x;
  qa1=(qa*q+pa*d)/x;
else 
  e=sqrt(pi*x/2.0);
  t=x-pi*(r/2.0+0.25);
  c=cos(t);
  s=sin(t);
  [p,q]= bessyax(r, x);
  [f,g]= bessjax(r, x);
  d=r/x;
  q=d*p-q;
  g=d*f-g;
  pa=e*(s*p+c*f);
  qa=e*(c*p-s*f);
  pa1=e*(s*g-c*q);
  qa1=e*(c*g+s*q);
end    
t=2.0/x;
b=(r+1.0)*t;
for n=1:na 
  c=pa-qa1*b;
  s=qa+pa1*b;
  pa=pa1;
  pa1=c;
  qa=qa1;
  qa1=s;
  b=b+t;
end
ra=abs(a)*qa/x+pa1;
sa=-abs(a)*pa/x+qa1;
  
function [sih,coh]=hypfun(x)   
ax=abs(x);
if ax<0.21 
  if ax<0.07 
    y=x*x;
  else
    y=x*x/9.0;
  end
  f=2.0+y*(y*28+2520.0)/(y*(y+420)+15120.0);
  f2=f*f;
  sih=2*x*f/(f2-y);
  coh=(f2+ y)/(f2-y);
  if (ax>=0.07) 
    ss=2.0*sih/3.0;
    f=ss*ss;
    sih=sih*(1.0+f/3.0);
    coh=coh*(1.0+f);
  end
else
  y=exp(x);
  f=1.0/y;
  coh=(y+f)/2.0;
  sih=(y-f)/2.0;
end
       
function [g] = gamstar(x)
giant=realmax/1000;
if (x>=3.0)
  g = exp(stirling(x));
elseif (x>0.0)
  g = gamma(x)/(exp(-x+(x-0.5)*log(x))*sqrt(2*pi));
else
  g = giant;
end
     
function [s]=stirling(x)
%Stirling series, function corresponding 
%with asymptotic series for log(gamma(x))
% that is:  1/(12x)-1/(360x**3)...; x>= 3
dwarf=realmin*1000.0;
giant=realmax/1000;
lnsqrttwopi=0.9189385332046727418;
if (x<dwarf)
  s =giant;
elseif (x<1.0)
  s = lngam1(x)-(x+0.5)*log(x)+x-lnsqrttwopi;
elseif (x<2.0)
  s =lngam1(x-1)-(x-0.5)*log(x)+x-lnsqrttwopi;
elseif (x<3.0)
  s =lngam1(x-2)-(x-0.5)*log(x)+x-lnsqrttwopi+log(x-1);
elseif (x<12.0)
  a=[1.996379051590076518221;
    -0.17971032528832887213e-2;
     0.131292857963846713e-4;
    -0.2340875228178749e-6;
     0.72291210671127e-8;
    -0.3280997607821e-9;
     0.198750709010e-10;
    -0.15092141830e-11;
     0.1375340084e-12;
    -0.145728923e-13;
     0.17532367e-14;
    -0.2351465e-15;
     0.346551e-16;
    -0.55471e-17;
     0.9548e-18;
    -0.1748e-18;
     0.332e-19;
    -0.58e-20];
     z=18.0/(x*x)-1.0;
  s=chepolsum(17,z,a)/(12.0*x);
else
  z=1.0/(x*x);
  if (x<1000.0)
    c=[0.25721014990011306473e-1;
          0.82475966166999631057e-1;
         -0.25328157302663562668e-2;
          0.60992926669463371e-3;
         -0.33543297638406e-3;
          0.250505279903e-3;
          0.30865217988013567769];
    s =((((((c(6)*z+c(5))*z+c(4))*z+c(3))*z+...
            c(2))*z+c(1))/(c(7)+z)/x);
  else
    s =(((-z*0.000595238095238095238095238095238+...
            0.000793650793650793650793650793651)*z...
            -0.00277777777777777777777777777778)*z+...
            0.0833333333333333333333333333333)/x;
  end
end

function y=lngam1(x)          
%ln(gamma(1+x)), -1<=x<=1
y=-logoneplusx(x*(x-1)*auxgam(x));

function y = auxgam(x)
% function g in 1/gamma(x+1)=1+x*(x-1)*g(x), -1<=x<=1
if x < 0
  y = -(1 + (1 + x)^2 * auxgam(1 + x)) / (1 - x);
else
  dr = zeros(1,18);
  dr=[-1.013609258009865776949;
     0.784903531024782283535e-1;
     0.67588668743258315530e-2;
     -0.12790434869623468120e-2;
     0.462939838642739585e-4;
     0.43381681744740352e-5;
     -0.5326872422618006e-6;
     0.172233457410539e-7;
     0.8300542107118e-9;
     -0.10553994239968e-9;
     0.39415842851e-11;
     0.362068537e-13;
     -0.107440229e-13;
     0.5000413e-15;
     -0.62452e-17;
     -0.5185e-18;
     0.347e-19;
     -0.9e-21];
  t=2*x-1;
  y=chepolsum(17, t, dr);
end

%% Function chepolsum
function [chep]=chepolsum(n,t,ak)
u0=0; u1=0; k=n; tt=t+t;
while k>=0
  u2=u1; 
  u1=u0; 
  u0=tt*u1-u2+ak(k+1); 
  k= k-1; 
end
s=(u0-u2)/2.0;
chep=s;
    
function y=logoneplusx(t)
% Computes log(1+t) with improved accuracy for small t
x=log(1+t);
if (t>-0.2928)&&(t<0.4142)
   [~,~,ex1]=hypfun1(x);
   s=ex1*x;
   x=x-2*(s-t)/(2+s+t);
end
y=x;

function [sinh_val, cosh_val, ex1] = hypfun1(x)
ax=abs(x);
if x==0
  sinh_val=0;
  cosh_val=1;
  ex1=1;
elseif ax<0.21
  if ax<0.07
    y=x^2;
  else
    y=x^2/9;
  end
  f=2+y*(y*28+2520)/(y*(y+420)+15120);
  f2=f^2;
  ex1=2/(f-x);
  sinh_val=2*x*f/(f2-y);
  cosh_val=(f2+y)/(f2-y);
  if ax >= 0.07
    ex1=2/(f-x/3);
    f=(2*sinh_val/3)^2;
    y=sinh_val/3 + cosh_val;
    sinh_val=sinh_val*(1+f/3);
    cosh_val=cosh_val*(1+f);
    ex1=ex1*(y^2+y+1)/3;
  end
else
  y=exp(x);
  f=1/y;
  cosh_val=(y+f)/2;
  ex1=(y-1)/x;
  sinh_val=(y-f)/2;
end

function [recip,q,r]=recipgam(x)  
%recipgam(x)=1/gamma(x+1)=1+x*(q + x * r); -0.5<=x<=0.5
if (x==0)
  q=0.5772156649015328606;
  r=-0.6558780715202538811;
else
  c=[+1.142022680371167841;
     -6.5165112670736881e-3;
     -3.087090173085368e-4;
     +3.4706269649043e-6;
     -6.9437664487e-9;
     -3.67795399e-11;
     +1.356395e-13;
     +3.68e-17;
     -5.5e-19];
  tx=2.0*x;
  t=2*tx*tx-1;
  q=chepolsum(8,t,c);
  c=[-1.270583625778727532;
     +2.05083241859700357e-2;
     -7.84761097993185e-5;
     -5.377798984020e-7;
     +3.8823289907e-9;
     -2.6758703e-12;
     -2.39860e-14;
     +3.80e-17;
     +4e-20];
  r=chepolsum(8,t,c);
end      
recip=1+x*(q+x*r);
   
function [jax, jap]= bessjax(a, x)
% Output: Ja(x) and derivative
aa=abs(a);
lneps=-log(eps);
if (a<-0.5) 
  [ia, iad]= bessjax(-a, x);
  [yax, yap]= bessyax(-a, x);
  [s,c]=sincospi(a);
  jax=c*ia+s*yax;
  jap=c*iad+s*yap;
elseif (aa==0.5) 
  s=sin(x);
  c=cos(x);
  r=sqrt(2.0/(pi*x));
  if (a<0) 
    jax=c*r;
    jap=a*jax/x-r*s;
  else
    jax=r*s;
    jap=r*c-a*jax/x;
  end  
else
  if ((x*x<4.0*(1.0 + aa))) 
    ia=exp(a*(log(x/(2.0*a))+1.0))...
       /(sqrt(2.0*pi*a)*gamstar(a)); 
    r=x*x/4.0;
    c=1.0;
    s=a;
    n=1;
    jax=0.0;
    jap=0.0;
    while true
      jax=jax+c;
      jap=jap+s;
      c=-c*r/(n*(a+n));
      s=(a+2.0*n)*c;
      n=n+1;
      if (abs(c)<=eps*abs(jax))&&(abs(s)<=eps*abs(jap))
        break;
      end
    end        
    jax=jax*ia;
    jap=jap*ia/x;    
  elseif (a>5.0+2.0*lneps*(0.4343+0.025*lneps)...
            /(2.0-3.6*x/a))&&(x<a/2.0) 
    [jax, jap, ~, ~, ~]= deb1(a, x);
  elseif (x>5.0+lneps*(0.4343+0.01*lneps)...
            /(2.0-3.4*a/x))&&(a<x/2.0) 
    [jax, jap, ~, ~]= deb2(a, x);
  else
    [jax, jap]= jrec(a, x);
  end
end

function [jax,jap]=jrec(a,x)
%Miller, with the standard recurrence relation
if x <= 0 
  if a == 0.0 
    jax=1.0;         
  else
    jax=0.0;
  end          
  jap=0.0; 
else
  % Function startijbes
  na=floor(a); 
  nu=startijbes(x,na+1,0,eps);
  nu2=nu+nu;
  x2=2.0/x;
  k=1;
  ap=a+1.0;
  am=a-1.0;
  % Logical variable
  notsimple=(a~=0.0)&&(a~=1.0);
  if a==0.0 
    la=2;
  elseif a==1.0 
    la=1;
  else
    g1=gamstar(nu+a);
    g2=gamstar(nu*1.0);
    g3=gamstar(a);
    la=g1*exp(a*log(1.0+nu/a)+nu*log(1.0+a/nu))...
           /(g2*g3*sqrt(2*pi*(nu+a)*nu*a));
  end
  j2=0.0;
  j1=1.0;
  n= nu2;
  m= nu; 
  s=0.0;
  while (n>0)
    j0=-j2+x2*(ap+n)*j1;
    if k==1 
      if a==0.0 
        s=s+2*j0;
      else
        s=s+la*(a+n)*j0;
      end
      if m>1 
        if (notsimple) 
          la=la*m/(m+am);
        end 
      end
      m=m-1;
    end
    j2=j1;
    j1=j0;
    k=-k;
    n=n-1;
  end
  j0=-j2+x2*ap*j1;
  s=s+j0;
  if a==0.0 
    s=1.0/s;
  elseif a==1.0 
    s=1.0/(s*x2);
  else   
    s=exp(a*(1.0-log(a*x2)))/(s*g3*sqrt(2.0*pi*a));
  end
  jax=s*j0;
  jap=s*(a/x*j0-j1);
end    
     
function s=startijbes(x,n,t,epss)
mactol=eps;
if x<=0.0 
  s=0;
else
  s=2*t-1;
  if (epss<mactol) 
    del=-log(mactol/2.0);
  else
    del=-log(epss/2.0);
  end
  p=del/x-t;
  r=n/x;
  if (r>1.0)||(t==1) 
    q=sqrt(r*r+s);
    r=r*log(q+r)-q;
  else
    r=0;
  end
  q=del/(2.0*x)+r;
  if (p>q)
    r=p;
  else
    r=q;
  end
  y=alfinv(t,r);
  [p,q]=hypfun(y);
  if (t==0) 
    s=floor(x*q)+1;
  else
    s=floor(x*p)+1;
  end
  if (mod(s,2)>0) 
    s=s+1;
  end
end
    
function alf=alfinv(t,r)
if ((t+r)<2.7) 
  if (t==0) 
    a=exp(log(3.0*r)/3.0);
    a2=a*a;
    b=a*(1.0+a2*(-1.0/30.0+0.004312*a2));
  else
    a=sqrt(2.0*(1.0+r));
    a2=a*a;
    b=a/(1.0+a2/8.0);
  end
else
  a=log(0.7357589*(r+t));
  lna=log(a)/a;
  b=1.0+a+log(a)*(1.0/a-1.0)+0.5*lna*lna;
end
while (abs(a/b-1.0)>1.0e-2)      
  a=b;
  b=fi(a,r,t); 
end
alf=b;
    
function [falf,df]=falfa(al,r,t) 
[sh,ch]=hypfun(al);
if (t==1) 
  falf=al*sh/ch-1.0-r/ch;
  df=(sh+(al+r*sh)/ch)/ch;
else
  falf=al-(sh+r)/ch;
  df=sh*(r+sh)/(ch*ch);
end
     
function fii=fi(al,r,t)
[p,q]=falfa(al,r,t);
fii=al-p/q;
     
function [J,Y,Jp,Yp]=bessair(nu,x)
%------------------------------------------
% Airy-type expansion for Bessel functions
%------------------------------------------
z= x/nu;
zeta= zetaz(z);
f= phizeta(zeta, z);
fhat= 2/(z*f);
nu2= nu*nu; nu13= nu^(1/3); nu23= nu13*nu13; 
w= zeta*nu23;
% number of terms (the k-th terms are O(nu^(-2k)))
if nu<100
  K=4;
elseif nu<400
  K=3;
else
  K=2;
end
sA=1; sB=Bkzeta(zeta,f,0); sC=Ckzeta(zeta,f,0); sD=1;
p=1;
for k=1:K
  p=p/nu2;
  sA=sA+Akzeta(zeta,f,k)*p;
  sD=sD+Dkzeta(zeta,f,k)*p;
  if k<K
    sB=sB+Bkzeta(zeta,f,k)*p;
    sC=sC+Ckzeta(zeta,f,k)*p;
  end
end
[Ai,Bi]=aibi(w);
[Aid,Bid]=aibip(w);
J= f*(Ai/nu13*sA+Aid/(nu23*nu)*sB);
Y=-f*(Bi/nu13*sA+Bid/(nu23*nu)*sB);
Jp=-fhat*(Ai/(nu*nu13)*sC+Aid/nu23*sD);
Yp=fhat*(Bi/(nu*nu13)*sC+Bid/nu23*sD);

function A=Akzeta(zeta, f, k)
two13=2^(1/3);
if abs(zeta) < 0.5
  eta= zeta/two13; 
  jmax= 18 - 2*k; 
  delta= eps*10^(2*k);
  if k==1
    ck=[-0.444444444444444444444444e-2;
       -0.184415584415584415584416e-2;
       0.112136752136752136752137e-2;
       0.134577521244187910854578e-2;
       0.388062656297950415597474e-3;
       -0.183068672378179913251358e-3;
       -0.199546088780673281140134e-3;
       -0.525619123404158729343881e-4;
       0.246061965245915785195825e-4;
       0.251924678092454136168976e-4;
       0.633315737653324223078736e-5;
       -0.295748573383020170629872e-5;
       -0.292525592056483802975188e-5;
       -0.715970261050200926824325e-6;
       0.333151072039094907664420e-6;
       0.322767047569230987018795e-6;
       0.776772938166419874312952e-7;
       -0.360095423792111979834763e-7;
       -0.344172444903422642801444e-7;
       -0.818819435639877179627180e-8;
       0.378314848515203785605629e-8;
       0.358157567915293191299802e-8];
  elseif k==2
    ck=[0.693735541354588973636593e-3;
       0.464483490365843307019778e-3;
       -0.428381301715351124588504e-3;
       -0.702670286877113311390726e-3;
       -0.263258004677881193972489e-3;
       0.166385366628870309272117e-3;
       0.221208768781858329671256e-3;
       0.702034561532966065008509e-4;
       -0.400042178254061367229691e-4;
       -0.478632496645396180982353e-4;
       -0.139460074147363092319545e-4;
       0.753618659127372677842368e-5;
       0.847850216106766137159346e-5;
       0.234535522845391156499770e-5;
       -0.122594329471088266135048e-5;
       -0.132508234340102697171799e-5;
       -0.353995477656999743787481e-6;
       0.180829171937667380735510e-6;
       0.190038351523365606575570e-6];
  elseif k==3
    ck=[-0.354211971457743840771126e-3;
       -0.312322527890318832782775e-3;
       0.371644223750229630164450e-3;
       0.753926915597773214092907e-3;
       0.340830005944473690438883e-3;
       -0.263496817206959254774719e-3;
       -0.408927572664843008395140e-3;
       -0.150110875956345964262796e-3;
       0.996401520553805536354607e-4;
       0.135249295575128253602523e-3;
       0.444311708727290144245992e-4;
       -0.271320507191411598209920e-4;
       -0.339679696977185933223319e-4;
       -0.104070886527304282710648e-4;
       0.602463906541441186398061e-5;
       0.714391960784688644874689e-5];
  elseif k==4
    ck=[0.378194199201772914026612e-3;
       0.404943905523632334772139e-3;
       -0.579130526946450830904535e-3;
       -0.138017901171009595794723e-2;
       -0.722520056780092714288715e-3;
       0.651265924036820506795172e-3;
       0.114674563328389541610941e-2;
       0.474423189340400087623894e-3;
       -0.356495172735465313237897e-3;
       -0.538157791035117996567133e-3;
       -0.195687390661222638146666e-3;
       0.132563525210291989094654e-3;
       0.181949256267290616373907e-3];
  elseif k==5
    ck=[-0.691141397288294167613745e-3;
       -0.859953266117743832936210e-3;
       0.142023355681435114892370e-2;
       0.385354269956030524409690e-2;
       0.227528116429013745967378e-2;
       -0.232195720345569883668239e-2;
       -0.454786433944346356182570e-2;
       -0.208244317582724499170221e-2;
       0.173704435731958087183915e-2;
       0.287424443211932650161305e-2];
  end
A= sumexpansion(ck, eta, jmax, delta);
else
  A=Akexact(k,zeta,f);
end

function A=Akexact(k,zeta,f)
% Exact coefficients (rational numbers rounded to double),
% polynomials in g=f^2 and zeta evaluated with Horner schemes
g=f*f;
if k==1
  s=((1.75781250000000000e-02)*g)*g;
  s=s*zeta+((((-2.50651041666666678e-02)*g)*g)*g-9.11458333333333391e-03)*g;
  s=s*zeta+((((((5.22189670138888899e-03)*g)*g)*g+3.79774305555555551e-03)*g)*g)*g-9.87413194444444475e-02;
  s=s/zeta^3;
elseif k==2
  s=((((7.00950622558593750e-03)*g)*g)*g)*g;
  s=s*zeta+((((((-3.69388580322265597e-02)*g)*g)*g-1.33514404296875000e-03)*g)*g)*g;
  s=s*zeta+((((((((3.43325138092041016e-02)*g)*g)*g+4.06150817871093767e-03)*g)*g)*g-1.73568725585937500e-03)*g)*g;
  s=s*zeta+((((((((((-1.09443384924052672e-02)*g)*g)*g-2.10371723881474246e-03)*g)*g)*g+2.47496145742910883e-03)*g)*g)*g-8.95700336974344057e-03)*g;
  s=s*zeta+((((((((((((1.14003525962554867e-03)*g)*g)*g+2.92182949835380909e-04)*g)*g)*g-5.15616970297731042e-04)*g)*g)*g+3.73208473739310067e-03)*g)*g)*g-3.17227202678413533e-01;
  s=s/zeta^6;
elseif k==3
  s=((((((8.94533470273017883e-03)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((-1.03482150339654508e-01)*g)*g)*g-1.03499740362167358e-03)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((2.13076671625206632e-01)*g)*g)*g+8.39543628195921557e-03)*g)*g)*g-6.92127893368403081e-04)*g)*g)*g)*g;
  s=s*zeta+((((((((((((-1.70795807464876104e-01)*g)*g)*g-1.21152747077522462e-02)*g)*g)*g+3.64739158087306562e-03)*g)*g)*g-1.31206104048976193e-03)*g)*g)*g;
  s=s*zeta+((((((((((((((6.46966828935546800e-02)*g)*g)*g+6.53816417239054450e-03)*g)*g)*g-3.39003771336542255e-03)*g)*g)*g+3.99128968516985604e-03)*g)*g)*g-5.57625942208148789e-03)*g)*g;
  s=s*zeta+((((((((((((((((-1.16768259909237927e-02)*g)*g)*g-1.50668722463532796e-03)*g)*g)*g+1.08065842318671806e-03)*g)*g)*g-2.06734654870502571e-03)*g)*g)*g+7.95133287963471422e-03)*g)*g)*g-5.89018217473200184e-02)*g;
  s=s*zeta+((((((((((((((((((8.10890693814152273e-04)*g)*g)*g+1.25557268719610681e-04)*g)*g)*g-1.12568585748616464e-04)*g)*g)*g+2.87131465097920232e-04)*g)*g)*g-1.65652768325723213e-03)*g)*g)*g+2.45424257280500054e-02)*g)*g)*g-3.51120304082635437e+00;
  s=s/zeta^9;
elseif k==4
  s=((((((((2.37267265674745431e-02)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((-4.82339164817468746e-01)*g)*g)*g-1.96843953744973987e-03)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((1.73572126525619241e+00)*g)*g)*g+3.07876154996241276e-02)*g)*g)*g-8.83274151419755071e-04)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((-2.51420013237900086e+00)*g)*g)*g-8.55135131074469673e-02)*g)*g)*g+1.02179640634858515e-02)*g)*g)*g-1.01710356830153614e-03)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((1.86463111851528107e+00)*g)*g)*g+9.44507061012110904e-02)*g)*g)*g-2.10394716991035187e-02)*g)*g)*g+8.25028948860069035e-03)*g)*g)*g-2.22360605209955160e-03)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((-7.75910100099241440e-01)*g)*g)*g-5.18720966072411457e-02)*g)*g)*g+1.68646033846611594e-02)*g)*g)*g-1.19058164717024198e-02)*g)*g)*g+1.17180106036982798e-02)*g)*g)*g-8.62819654501758085e-03)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((1.83627129776507886e-01)*g)*g)*g+1.51166677546563082e-02)*g)*g)*g-6.38823583258840699e-03)*g)*g)*g+6.42512733520873500e-03)*g)*g)*g-1.08912073166118218e-02)*g)*g)*g+2.62469738899434799e-02)*g)*g)*g-6.17203659520257561e-02)*g)*g;
  s=s*zeta+((((((((((((((((((((((-2.31219764682382388e-02)*g)*g)*g-2.24240840562825059e-03)*g)*g)*g+1.15298520526699766e-03)*g)*g)*g-1.48063845100339484e-03)*g)*g)*g+3.47184188511140890e-03)*g)*g)*g-1.35950019080152001e-02)*g)*g)*g+8.80086699686293272e-02)*g)*g)*g-9.82953976273002850e-01)*g;
  s=s*zeta+((((((((((((((((((((((((1.20426960772074157e-03)*g)*g)*g+1.33476690811205407e-04)*g)*g)*g-8.00684170324303916e-05)*g)*g)*g+1.23386537583616254e-04)*g)*g)*g-3.61650196365771769e-04)*g)*g)*g+1.88819470944655570e-03)*g)*g)*g-1.83351395767977765e-02)*g)*g)*g+4.09564156780417854e-01)*g)*g)*g-8.22814390971859382e+01;
  s=s/zeta^12;
elseif k==5
  s=((((((((((1.07438613544186268e-01)*g)*g)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((-3.39015863127857431e+00)*g)*g)*g-6.94432405179412626e-03)*g)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((1.88102053596595695e+01)*g)*g)*g+1.78007134731836630e-01)*g)*g)*g-2.34280828736999069e-03)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((-4.25051594251473190e+01)*g)*g)*g-8.04980933811287103e-01)*g)*g)*g+4.76268055538082227e-02)*g)*g)*g-1.93440763283092565e-03)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((5.06926237732765870e+01)*g)*g)*g+1.47598025474036221e+00)*g)*g)*g-1.71387407919176976e-01)*g)*g)*g+3.02553353994787187e-02)*g)*g)*g-2.83770350476923250e-03)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((-3.58268515212663203e+01)*g)*g)*g-1.41120861581983137e+00)*g)*g)*g+2.48255438418499452e-01)*g)*g)*g-8.40350893782311342e-02)*g)*g)*g+3.28273530793956464e-02)*g)*g)*g-6.68853751702632327e-03)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((((1.58177643477270671e+01)*g)*g)*g+7.82606704477697201e-01)*g)*g)*g-1.84116136919369122e-01)*g)*g)*g+9.28177692697448459e-02)*g)*g)*g-6.75937164956911946e-02)*g)*g)*g+5.42544264818408048e-02)*g)*g)*g-2.46117995739686055e-02)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((((((-4.41701103517607763e+00)*g)*g)*g-2.61693190371938833e-01)*g)*g)*g+7.66143870540700600e-02)*g)*g)*g-5.09752916962797273e-02)*g)*g)*g+5.41810762312835495e-02)*g)*g)*g-7.82934035663547245e-02)*g)*g)*g+1.29699830647406894e-01)*g)*g)*g-1.43987398868115651e-01)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((((((((7.59252942989239732e-01)*g)*g)*g+5.21195715411671701e-02)*g)*g)*g-1.81315850799286214e-02)*g)*g)*g+1.48553191150901306e-02)*g)*g)*g-2.05235477368947219e-02)*g)*g)*g+4.22520445041585715e-02)*g)*g)*g-1.20548426886090243e-01)*g)*g)*g+4.38009667356807797e-01)*g)*g)*g-1.44635342163022163e+00)*g)*g;
  s=s*zeta+((((((((((((((((((((((((((((-7.34136220754393470e-02)*g)*g)*g-5.70115031217587487e-03)*g)*g)*g+2.28309446463723911e-03)*g)*g)*g-2.20363991539784332e-03)*g)*g)*g+3.70420684526334872e-03)*g)*g)*g-9.73677227897794609e-03)*g)*g)*g+3.84277945943662935e-02)*g)*g)*g-2.26873478383215171e-01)*g)*g)*g+2.06239284195420502e+00)*g)*g)*g-3.07722106577294063e+01)*g;
  s=s*zeta+((((((((((((((((((((((((((((((3.05890091980997279e-03)*g)*g)*g+2.63942144082216455e-04)*g)*g)*g-1.18911170033189541e-04)*g)*g)*g+1.31169042583204965e-04)*g)*g)*g-2.57236586476621433e-04)*g)*g)*g+8.11397689914828768e-04)*g)*g)*g-4.00289527024648862e-03)*g)*g)*g+3.15102053310021071e-02)*g)*g)*g-4.29665175407126065e-01)*g)*g)*g+1.28217544407205857e+01)*g)*g)*g-3.31621856854797261e+03;
  s=s/zeta^15;
else
  s=0;
end
A=s;

function B=Bkexact(k,zeta,f)
% Exact coefficients (rational numbers rounded to double),
% polynomials in g=f^2 and zeta evaluated with Horner schemes
g=f*f;
if k==0
  s=(-6.25000000000000000e-02)*g;
  s=s*zeta+(((2.60416666666666678e-02)*g)*g)*g-1.04166666666666671e-01;
  s=s/zeta^2;
elseif k==1
  s=(((-9.15527343750000000e-03)*g)*g)*g;
  s=s*zeta+(((((2.78503417968750007e-02)*g)*g)*g-1.83105468750000000e-03)*g)*g;
  s=s*zeta+(((((((-1.44254896375868059e-02)*g)*g)*g+2.61094835069444449e-03)*g)*g)*g-5.22189670138888899e-03)*g;
  s=s*zeta+(((((((((2.00354022744261189e-03)*g)*g)*g-5.43947573061342548e-04)*g)*g)*g+2.17579029224537019e-03)*g)*g)*g-1.28226574556327161e-01;
  s=s/zeta^5;
elseif k==2
  s=(((((-7.09712505340576172e-03)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((5.75687059334346260e-02)*g)*g)*g-7.30156898498535156e-04)*g)*g)*g)*g;
  s=s*zeta+(((((((((-8.30761694245868293e-02)*g)*g)*g+3.84779771169026707e-03)*g)*g)*g-7.64926274617513057e-04)*g)*g)*g;
  s=s*zeta+(((((((((((4.48331257535351663e-02)*g)*g)*g-3.57630352179209377e-03)*g)*g)*g+2.32690572738647461e-03)*g)*g)*g-2.25398275587293832e-03)*g)*g;
  s=s*zeta+(((((((((((((-1.03315695403565343e-02)*g)*g)*g+1.14003525962554867e-03)*g)*g)*g-1.20525466807094622e-03)*g)*g)*g+3.21401244818918979e-03)*g)*g)*g-1.82405641540087787e-02)*g;
  s=s*zeta+(((((((((((((((8.60964128363044564e-04)*g)*g)*g-1.18753672877661324e-04)*g)*g)*g+1.67396481676520306e-04)*g)*g)*g-6.69585926706081224e-04)*g)*g)*g+7.60023506417032476e-03)*g)*g)*g-8.81627267443757634e-01;
  s=s/zeta^8;
elseif k==3
  s=(((((((-1.34978711139410734e-02)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((2.11115077711708315e-01)*g)*g)*g-9.31805698201060295e-04)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((-5.86378375593922141e-01)*g)*g)*g+1.07793906603806788e-02)*g)*g)*g-5.92967262491583824e-04)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((6.47661984694018922e-01)*g)*g)*g-2.21954866276256914e-02)*g)*g)*g+4.80988536987246758e-03)*g)*g)*g-8.98804972638134647e-04)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((-3.55694376735367856e-01)*g)*g)*g+1.77912299442579275e-02)*g)*g)*g-6.94104280131639054e-03)*g)*g)*g+4.73654323349488333e-03)*g)*g)*g-2.67195763974737962e-03)*g)*g)*g;
  s=s*zeta+(((((((((((((((((1.03657150317643262e-01)*g)*g)*g-6.73923780141194584e-03)*g)*g)*g+3.74582322376541647e-03)*g)*g)*g-4.40234064166204146e-03)*g)*g)*g+8.12809514011152839e-03)*g)*g)*g-1.54973543105348024e-02)*g)*g;
  s=s*zeta+(((((((((((((((((((-1.53765147814508615e-02)*g)*g)*g+1.21633604072122841e-03)*g)*g)*g-8.63206222447323300e-04)*g)*g)*g+1.40335503566608511e-03)*g)*g)*g-4.21006510699825556e-03)*g)*g)*g+2.20980792946514765e-02)*g)*g)*g-2.07588017616422982e-01)*g;
  s=s*zeta+(((((((((((((((((((((9.15268736991122758e-04)*g)*g)*g-8.44677806056408527e-05)*g)*g)*g+7.19338518706102840e-05)*g)*g)*g-1.46182816215217204e-04)*g)*g)*g+5.84731264860868815e-04)*g)*g)*g-4.60376651971905818e-03)*g)*g)*g+8.64950073401762332e-02)*g)*g)*g-1.49957629868625553e+01;
  s=s/zeta^11;
elseif k==4
  s=(((((((((-4.76182220694454372e-02)*g)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((1.22062035244687972e+00)*g)*g)*g-2.47153401744526491e-03)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((-5.51986926042025505e+00)*g)*g)*g+5.02436630018196634e-02)*g)*g)*g-1.12775181833058014e-03)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((1.01210074610767702e+01)*g)*g)*g-1.80804298464186691e-01)*g)*g)*g+1.76387380466596563e-02)*g)*g)*g-1.14702962719093193e-03)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((((-9.67685907990741612e+00)*g)*g)*g+2.61895847122812608e-01)*g)*g)*g-4.89921168844748289e-02)*g)*g)*g+1.32691616657767651e-02)*g)*g)*g-2.07128903753073246e-03)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((((((5.36644597356135211e+00)*g)*g)*g-1.94232408178675098e-01)*g)*g)*g+5.41123837038188574e-02)*g)*g)*g-2.73220917203635984e-02)*g)*g)*g+1.68013707814732814e-02)*g)*g)*g-6.17977181979333795e-03)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((((((((-1.79446759112186616e+00)*g)*g)*g+8.08239687603376500e-02)*g)*g)*g-2.97183886812319061e-02)*g)*g)*g+2.19005613398030353e-02)*g)*g)*g-2.42456991689356562e-02)*g)*g)*g+3.25663044694447998e-02)*g)*g)*g-3.04084010180307104e-02)*g)*g)*g;
  s=s*zeta+(((((((((((((((((((((((3.57391347710860563e-01)*g)*g)*g-1.91278260183862375e-02)*g)*g)*g+8.66059090110517758e-03)*g)*g)*g-8.29583403259744577e-03)*g)*g)*g+1.30845041045136220e-02)*g)*g)*g-3.02684803340836882e-02)*g)*g)*g+9.25023558968494225e-02)*g)*g)*g-2.63597396253443361e-01)*g)*g;
  s=s*zeta+(((((((((((((((((((((((((-3.90936021406345716e-02)*g)*g)*g+2.40853921544148314e-03)*g)*g)*g-1.28471314905785204e-03)*g)*g)*g+1.49727939850644834e-03)*g)*g)*g-3.01525851219962195e-03)*g)*g)*g+9.64882723903879094e-03)*g)*g)*g-4.79129407522063469e-02)*g)*g)*g+3.75870361324354418e-01)*g)*g)*g-4.93268831322415746e+00)*g;
  s=s*zeta+(((((((((((((((((((((((((((1.80988898799234132e-03)*g)*g)*g-1.25444750804243918e-04)*g)*g)*g+7.64710207772530939e-05)*g)*g)*g-1.03977736007392250e-04)*g)*g)*g+2.51271542683301847e-04)*g)*g)*g-1.00508617073320739e-03)*g)*g)*g+6.65457510447310403e-03)*g)*g)*g-7.83063252759071682e-02)*g)*g)*g+2.05528679717673235e+00)*g)*g)*g-4.74451538868264322e+02;
  s=s/zeta^14;
else
  s=0;
end
B=s;

function B=Bkzeta(zeta, f, k) 
two13=2^(1/3);
if abs(zeta) < 0.5
  eta= zeta/two13; 
  jmax= 18 - 2*k; 
  delta= eps*10^(2*k);
  if k==0
    ck=[0.142857142857142857142857e-1;
       0.888888888888888888888889e-2;
       0.204823747680890538033395e-2;
       -0.578266178266178266178266e-3;
       -0.604120897998449018857182e-3;
       -0.147268574562692209751033e-3;
       0.532410214800978351784867e-4;
       0.520656100658341554700098e-4;
       0.123311505089493920144250e-4;
       -0.490593272853136642707052e-5;
       -0.463223098713634979781493e-5;
       -0.107717452345523495409932e-5;
       0.447596397893282157770878e-6;
       0.415258618846462448367911e-6;
       0.955581929358923426517202e-7;
       -0.406059920840305938391248e-7;
       -0.373136718798848249861612e-7;
       -0.853267064555377755230178e-8;
       0.367301724557362439212528e-8;
       0.335596046078453629690337e-8;
       0.764310709511047504468419e-9;
       -0.301966912392930035742926e-9];
  elseif k==1
    ck=[-0.118485958485958485958486e-2;
       -0.139406307977736549165121e-2;
       -0.481410055863837376442418e-3;
       0.268417053660161429584988e-3;
       0.341970698270990239629547e-3;
       0.103454823490207815639056e-3;
       -0.541819198209550495372548e-4;
       -0.620218483069016491051219e-4;
       -0.172488588605608661593333e-4;
       0.874467599288705317859166e-5;
       0.942068421618092915552474e-5;
       0.249492211208585079667788e-5;
       -0.123845860883635664944695e-5;
       -0.128546171380976866640804e-5;
       -0.329971086253750597036374e-6;
       0.161344110578831478602255e-6;
       0.163362319440237459346168e-6;
       0.410425294960577916294255e-7;
       -0.198431704232698899355524e-7;
       -0.197394814276970718379204e-7];
  elseif k==2
    ck=[0.438291809448988109926169e-3;
       0.711048651167086689443780e-3;
       0.318583839453875805781720e-3;
       -0.240480942680445795531986e-3;
       -0.372296603862153505390841e-3;
       -0.135275205959561728461416e-3;
       0.869169437270413932413888e-4;
       0.115875075359137640445758e-3;
       0.372496592784644574529580e-4;
       -0.219833494960693482739066e-4;
       -0.268644963387045129596591e-4;
       -0.802306161203252181331827e-5;
       0.449475659218012475517714e-5;
       0.519350476385601375304097e-5;
       0.147715619152961744841117e-5;
       -0.798879382609678364636004e-6;
       -0.887808168959802029947558e-6];
  elseif k==3
    ck=[-0.376704394771054542175466e-3;
       -0.758562716587986423627298e-3;
       -0.410325396877504542517502e-3;
       0.379126331042900768687048e-3;
       0.685098167390344516518718e-3;
       0.287831057193221563639176e-3;
       -0.215701063611570529397680e-3;
       -0.326086399137349916645268e-3;
       -0.118131700874867756897152e-3;
       0.788752684158257871094102e-4;
       0.107208183342068435757316e-3;
       0.354459525128873402372981e-4;
       -0.220144792073382420302000e-4;
       -0.278933635962081525723724e-4];
  elseif k==4
    ck=[0.584533301220761872533486e-3;
       0.138546904223724012512638e-2;
       0.868303741849469002543831e-3;
       -0.935029048013459516962523e-3;
       -0.191754860055254920576692e-2;
       -0.907950471133081379467134e-3;
       0.770504298063922351074396e-3;
       0.129531001282559834871519e-2;
       0.519338694718995507588524e-3;
       -0.384826319487598346504950e-3;
       -0.573353930990124765012774e-3];
  end
B= two13*sumexpansion(ck, eta, jmax, delta);
else
  B=Bkexact(k,zeta,f);
end

function C=Ckexact(k,zeta,f)
% Exact coefficients (rational numbers rounded to double),
% polynomials in g=f^2 and zeta evaluated with Horner schemes
g=f*f;
if k==0
  s=(1.87500000000000000e-01)*g;
  s=s*zeta+(((-3.64583333333333356e-02)*g)*g)*g+1.45833333333333343e-01;
  s=s/zeta^1;
elseif k==1
  s=(((1.28173828125000000e-02)*g)*g)*g;
  s=s*zeta+(((((-3.40393066406250028e-02)*g)*g)*g-4.27246093750000000e-03)*g)*g;
  s=s*zeta+(((((((1.66447957356770822e-02)*g)*g)*g+4.69970703125000000e-03)*g)*g)*g-1.85139973958333322e-02)*g;
  s=s*zeta+(((((((((-2.23925084243586014e-03)*g)*g)*g-8.99985984519675918e-04)*g)*g)*g+3.59994393807870367e-03)*g)*g)*g+1.43312053915895049e-01;
  s=s/zeta^4;
elseif k==2
  s=(((((8.67426395416259766e-03)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((-6.64254299231937934e-02)*g)*g)*g-1.31428241729736328e-03)*g)*g)*g)*g;
  s=s*zeta+(((((((((9.28498364157146899e-02)*g)*g)*g+6.36635621388753278e-03)*g)*g)*g-1.26560529073079420e-03)*g)*g)*g;
  s=s*zeta+(((((((((((-4.91029472538718481e-02)*g)*g)*g-5.67440158791012242e-03)*g)*g)*g+3.36108605066935207e-03)*g)*g)*g-4.19859532956723836e-03)*g)*g;
  s=s*zeta+(((((((((((((1.11580951035850570e-02)*g)*g)*g+1.76405455963111221e-03)*g)*g)*g-1.64352909282401748e-03)*g)*g)*g+4.61845486252396211e-03)*g)*g)*g-5.94801005022025375e-02)*g;
  s=s*zeta+(((((((((((((((-9.20340964801875288e-04)*g)*g)*g-1.80712110900788959e-04)*g)*g)*g+2.21106582749200614e-04)*g)*g)*g-8.84426330996802455e-04)*g)*g)*g+1.15655750976504934e-02)*g)*g)*g+9.42429147957120295e-01;
  s=s/zeta^7;
elseif k==3
  s=(((((((1.55744666699320078e-02)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((-2.35952145677791630e-01)*g)*g)*g-1.54171488247811794e-03)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((6.42223935174295657e-01)*g)*g)*g+1.71032998478040099e-02)*g)*g)*g-8.56508268043398857e-04)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((-6.99474943469540422e-01)*g)*g)*g-3.43445950974839639e-02)*g)*g)*g+6.55893459528063705e-03)*g)*g)*g-1.29156008673210932e-03)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((3.80225023406772500e-01)*g)*g)*g+2.70736107847403262e-02)*g)*g)*g-9.16811535788849413e-03)*g)*g)*g+6.25628972552532726e-03)*g)*g)*g-4.06602249526775202e-03)*g)*g)*g;
  s=s*zeta+(((((((((((((((((-1.09939401852045884e-01)*g)*g)*g-1.01338168421231476e-02)*g)*g)*g+4.84848980045826574e-03)*g)*g)*g-5.57629814610525239e-03)*g)*g)*g+1.07981940267182150e-02)*g)*g)*g-2.76102289440562573e-02)*g)*g;
  s=s*zeta+(((((((((((((((((((1.62076777426103692e-02)*g)*g)*g+1.81273306713937909e-03)*g)*g)*g-1.10176503301458368e-03)*g)*g)*g+1.73355622052869343e-03)*g)*g)*g-5.28018199038242722e-03)*g)*g)*g+3.03712518384618817e-02)*g)*g)*g-6.58350570154941472e-01)*g;
  s=s*zeta+(((((((((((((((((((((-9.59915992454104345e-04)*g)*g)*g-1.25012315296348466e-04)*g)*g)*g+9.08756812033101633e-05)*g)*g)*g-1.77587820210302230e-04)*g)*g)*g+7.10351280841208920e-04)*g)*g)*g-5.81604359701185045e-03)*g)*g)*g+1.28012610863460830e-01)*g)*g)*g+1.57272636203680456e+01;
  s=s/zeta^10;
elseif k==4
  s=(((((((((5.32203658423213710e-02)*g)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((-1.33686990982277298e+00)*g)*g)*g-3.92150064101315365e-03)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((5.96145880125387517e+00)*g)*g)*g+7.77454574870262149e-02)*g)*g)*g-1.53784338863260928e-03)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((-1.08190079756337880e+01)*g)*g)*g-2.75136975923762350e-01)*g)*g)*g+2.32982261899729176e-02)*g)*g)*g-1.51506052361582988e-03)*g)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((((1.02633353877805931e+01)*g)*g)*g+3.93813755303192292e-01)*g)*g)*g-6.34140387379133058e-02)*g)*g)*g+1.68076047766505705e-02)*g)*g)*g-2.75171248947319539e-03)*g)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((((((-5.65652413429439793e+00)*g)*g)*g-2.89468943801767409e-01)*g)*g)*g+6.90670788365106147e-02)*g)*g)*g-3.37508191839785629e-02)*g)*g)*g+2.10719533212457531e-02)*g)*g)*g-8.49338097400168142e-03)*g)*g)*g)*g;
  s=s*zeta+(((((((((((((((((((((1.88200259556683536e+00)*g)*g)*g+1.19619473765299716e-01)*g)*g)*g-3.75439204969794937e-02)*g)*g)*g+2.66055412695816937e-02)*g)*g)*g-2.94544938753054653e-02)*g)*g)*g+4.11417576839694515e-02)*g)*g)*g-4.50044335066854470e-02)*g)*g)*g;
  s=s*zeta+(((((((((((((((((((((((-3.73275407609121035e-01)*g)*g)*g-2.81522362424453873e-02)*g)*g)*g+1.08555615978040091e-02)*g)*g)*g-9.95861558161479851e-03)*g)*g)*g+1.55767906006114547e-02)*g)*g)*g-3.66700899679588549e-02)*g)*g)*g+1.19518916984183210e-01)*g)*g)*g-4.60759676377970051e-01)*g)*g;
  s=s*zeta+(((((((((((((((((((((((((4.06892593708645528e-02)*g)*g)*g+3.52879001332124251e-03)*g)*g)*g-1.60036748543570259e-03)*g)*g)*g+1.78139313636345941e-03)*g)*g)*g-3.53965129692999101e-03)*g)*g)*g+1.13999755582835788e-02)*g)*g)*g-5.84432574010429126e-02)*g)*g)*g+5.06835644015767084e-01)*g)*g)*g-1.54277698307223652e+01)*g;
  s=s*zeta+(((((((((((((((((((((((((((-1.87818668565242959e-03)*g)*g)*g-1.83095955429173041e-04)*g)*g)*g+9.47833716507416390e-05)*g)*g)*g-1.22851005736465143e-04)*g)*g)*g+2.91957189774451169e-04)*g)*g)*g-1.16782875909780468e-03)*g)*g)*g+7.86246436713376913e-03)*g)*g)*g-9.70581725703594383e-02)*g)*g)*g+2.99984413375157111e+00)*g)*g)*g+4.92355370523670501e+02;
  s=s/zeta^13;
else
  s=0;
end
C=s;

function D=Dkexact(k,zeta,f)
% Exact coefficients (rational numbers rounded to double),
% polynomials in g=f^2 and zeta evaluated with Horner schemes
g=f*f;
if k==1
  s=((-2.92968750000000000e-02)*g)*g;
  s=s*zeta+((((3.22265625000000000e-02)*g)*g)*g-1.95312500000000000e-02)*g;
  s=s*zeta+((((((-6.17133246527777797e-03)*g)*g)*g+3.79774305555555551e-03)*g)*g)*g+8.35503472222222238e-02;
  s=s/zeta^3;
elseif k==2
  s=((((-9.01222229003906250e-03)*g)*g)*g)*g;
  s=s*zeta+((((((4.36550140380859361e-02)*g)*g)*g-1.33514404296875000e-03)*g)*g)*g;
  s=s*zeta+((((((((-3.89101823170979813e-02)*g)*g)*g+3.54576110839843750e-03)*g)*g)*g-2.44776407877604152e-03)*g)*g;
  s=s*zeta+((((((((((1.20963741231847697e-02)*g)*g)*g-1.73383288913302961e-03)*g)*g)*g+2.69254048665364598e-03)*g)*g)*g-2.40424827293113409e-02)*g;
  s=s*zeta+((((((((((((-1.23916876046255287e-03)*g)*g)*g+2.33255296087068792e-04)*g)*g)*g-5.15616970297731042e-04)*g)*g)*g+4.67492719736609455e-03)*g)*g)*g+2.91849026464140460e-01;
  s=s/zeta^6;
elseif k==3
  s=((((((-1.05717591941356659e-02)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((1.17279770384941778e-01)*g)*g)*g-9.03569161891937256e-04)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((-2.35505794954175757e-01)*g)*g)*g+6.91931561699935405e-03)*g)*g)*g-7.52974301576614380e-04)*g)*g)*g)*g;
  s=s*zeta+((((((((((((1.85647616809647947e-01)*g)*g)*g-9.67185795997027933e-03)*g)*g)*g+3.64739158087306562e-03)*g)*g)*g-1.64352909282401748e-03)*g)*g)*g;
  s=s*zeta+((((((((((((((-6.94890297745587304e-02)*g)*g)*g+5.11489033894498447e-03)*g)*g)*g-3.25095924307350738e-03)*g)*g)*g+4.36474369079978362e-03)*g)*g)*g-8.55026444719161546e-03)*g)*g;
  s=s*zeta+((((((((((((((((1.24301696032414564e-02)*g)*g)*g-1.16230157329011018e-03)*g)*g)*g+1.01065625812199141e-03)*g)*g)*g-2.13430514137563383e-03)*g)*g)*g+9.40529089191077632e-03)*g)*g)*g-1.65305112645704549e-01)*g;
  s=s*zeta+((((((((((((((((((-8.57227304889246673e-04)*g)*g)*g+9.58688505001953334e-05)*g)*g)*g-1.03532980203577008e-04)*g)*g)*g+2.87131465097920232e-04)*g)*g)*g-1.80109737197786343e-03)*g)*g)*g+3.21426607922203345e-02)*g)*g)*g+3.32140828186276771e+00;
  s=s/zeta^9;
elseif k==4
  s=((((((((-2.68902901098044822e-02)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((5.33111708482465474e-01)*g)*g)*g-1.62234027811791748e-03)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((-1.88665354919151329e+00)*g)*g)*g+2.45783485081032954e-02)*g)*g)*g-8.83274151419755071e-04)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((2.70043717922188975e+00)*g)*g)*g-6.68983265806557953e-02)*g)*g)*g+9.79876553780438010e-03)*g)*g)*g-1.11227115363969164e-03)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((-1.98492990035497652e+00)*g)*g)*g+7.28619732780771273e-02)*g)*g)*g-1.96765909412668534e-02)*g)*g)*g+8.51750534248249362e-03)*g)*g)*g-2.63020830162632698e-03)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((8.20247820104912373e-01)*g)*g)*g-3.96067732715388032e-02)*g)*g)*g+1.55109228454241438e-02)*g)*g)*g-1.19058164717024198e-02)*g)*g)*g+1.27406733472937659e-02)*g)*g)*g-1.13001541847649600e-02)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((-1.93043905662482640e-01)*g)*g)*g+1.14520210262547793e-02)*g)*g)*g-5.80583256579971993e-03)*g)*g)*g+6.29630272698399832e-03)*g)*g)*g-1.13558988287872590e-02)*g)*g)*g+3.00099808992543744e-02)*g)*g)*g-9.73068832576982623e-02)*g)*g;
  s=s*zeta+((((((((((((((((((((((2.41974172342028053e-02)*g)*g)*g-1.68829976485524672e-03)*g)*g)*g+1.03854498638193595e-03)*g)*g)*g-1.43076431370643838e-03)*g)*g)*g+3.53031501159749555e-03)*g)*g)*g-1.46745057816044973e-02)*g)*g)*g+1.07037571583468091e-01)*g)*g)*g-2.81170556003672889e+00)*g;
  s=s*zeta+((((((((((((((((((((((((-1.25551512294290067e-03)*g)*g)*g+9.99912492139692026e-05)*g)*g)*g-7.16216389718663145e-05)*g)*g)*g+1.18012169340409731e-04)*g)*g)*g-3.61650196365771769e-04)*g)*g)*g+1.97418460133786007e-03)*g)*g)*g-2.04975147603021803e-02)*g)*g)*g+5.46720525562697346e-01)*g)*g)*g+7.89230130115865194e+01;
  s=s/zeta^12;
elseif k==5
  s=((((((((((-1.18747941285679559e-01)*g)*g)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((3.68495503399845026e+00)*g)*g)*g-5.54378810857514281e-03)*g)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((-2.02035539048195396e+01)*g)*g)*g+1.39257282273205524e-01)*g)*g)*g-2.24669307558045261e-03)*g)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((4.52474277751568223e+01)*g)*g)*g-6.20985291797278682e-01)*g)*g)*g+4.45416683519421000e-02)*g)*g)*g-1.99706051162706899e-03)*g)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((-5.35893451317495391e+01)*g)*g)*g+1.12697999746185284e+00)*g)*g)*g-1.57630559122988845e-01)*g)*g)*g+3.02553353994787187e-02)*g)*g)*g-3.08535762882182008e-03)*g)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((3.76641259582543384e+01)*g)*g)*g-1.06909743622714504e+00)*g)*g)*g+2.25622463975787252e-01)*g)*g)*g-8.23501753054846353e-02)*g)*g)*g+3.42279868107831925e-02)*g)*g)*g-7.64746762699425493e-03)*g)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((((-1.65534743173887904e+01)*g)*g)*g+5.89221263988999766e-01)*g)*g)*g-1.65841582386429265e-01)*g)*g)*g+8.96912759890797495e-02)*g)*g)*g-6.87321369840396823e-02)*g)*g)*g+5.85624702719621595e-02)*g)*g)*g-2.99332697521239796e-02)*g)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((((((4.60496895156654951e+00)*g)*g)*g-1.96041937038212016e-01)*g)*g)*g+6.85319901780362950e-02)*g)*g)*g-4.87549523120497583e-02)*g)*g)*g+5.41810762312835495e-02)*g)*g)*g-8.18589475617864382e-02)*g)*g)*g+1.44996125170934020e-01)*g)*g)*g-1.92206434768135759e-01)*g)*g)*g;
  s=s*zeta+((((((((((((((((((((((((((-7.89027568204504037e-01)*g)*g)*g+3.88828549592834458e-02)*g)*g)*g-1.61288853472343376e-02)*g)*g)*g+1.40971529082593739e-02)*g)*g)*g-2.02803056896426366e-02)*g)*g)*g+4.32904972108660024e-02)*g)*g)*g-1.29236601796799438e-01)*g)*g)*g+5.10445374619949099e-01)*g)*g)*g-2.31219764682382367e+00)*g)*g;
  s=s*zeta+((((((((((((((((((((((((((((7.60832083327280528e-02)*g)*g)*g-4.23846451779839092e-03)*g)*g)*g+2.02170261179862847e-03)*g)*g)*g-2.07825499844775273e-03)*g)*g)*g+3.62773289749017012e-03)*g)*g)*g-9.83728089605126650e-03)*g)*g)*g+4.01769971932563683e-02)*g)*g)*g-2.49601411816954083e-01)*g)*g)*g+2.54341741150620626e+00)*g)*g)*g-8.89596635377995568e+01)*g;
  s=s*zeta+((((((((((((((((((((((((((((((-3.16259247641370068e-03)*g)*g)*g+1.95644446422128077e-04)*g)*g)*g-1.04898724464630387e-04)*g)*g)*g+1.23086739574226981e-04)*g)*g)*g-2.50180954390405539e-04)*g)*g)*g+8.11397689914828768e-04)*g)*g)*g-4.11578538362594292e-03)*g)*g)*g+3.35792749013004710e-02)*g)*g)*g-4.87060152455944362e-01)*g)*g)*g+1.72977123545721376e+01)*g)*g)*g+3.20749009089066203e+03;
  s=s/zeta^15;
else
  s=0;
end
D=s;

function D=Dkzeta(zeta, f, k)
two13=2^(1/3);
if abs(zeta) < 0.5
  eta= zeta/two13; 
  jmax= 18 - 2*k; 
  delta= eps*10^(2*k);
  if k==1
    ck=[0.730158730158730158730159e-2;
       0.419336219336219336219336e-2;
       -0.450495536209821924107638e-3;
       -0.152302582778773254963731e-2;
       -0.602260874305692232863101e-3;
       0.134516202087850843046701e-3;
       0.227957335910397990094302e-3;
       0.759662380303310963431066e-4;
       -0.208929293145702503498759e-4;
       -0.288989984337042437996113e-4;
       -0.882085288915952350246854e-5;
       0.268699754140326485858770e-5;
       0.336422008038391207907941e-5;
       0.973596713460319029167637e-6;
       -0.315173324724484568747324e-6;
       -0.371902658726380167612929e-6;
       -0.103842224770769608750491e-6;
       0.350150687274032091141572e-7;
       0.397171235257758884576377e-7;
       0.108068143474318368422919e-7;
       -0.375354806995033684725386e-8;
       -0.413840555021600212423505e-8];
  elseif k==2
    ck=[-0.937299455394693489931585e-3;
       -0.790690489715699799733413e-3;
       0.293543703725772342692795e-3;
       0.759272400940109832633704e-3;
       0.350160108830627551277944e-3;
       -0.139795674313775324610736e-3;
       -0.238208912033003463924301e-3;
       -0.885682782902038147880529e-4;
       0.358182517107331492611466e-4;
       0.514900346471815453486348e-4;
       0.171249859424247232793150e-4;
       -0.697628383205670347703814e-5;
       -0.911796118302585196984065e-5;
       -0.283069240173736478532390e-5;
       0.115921686495136845230928e-5;
       0.142494761926653767630103e-5;
       0.422052804623112329027247e-6;
       -0.173575341276739650538416e-6;
       -0.204376944779354029041699e-6];
  elseif k==3
    ck=[0.444495041599140470657888e-3;
       0.472063930501695695068578e-3;
       -0.288641537836417409885435e-3;
       -0.801448906698328116457619e-3;
       -0.427270049238916578621937e-3;
       0.231003885183789413048089e-3;
       0.431796240309483025358877e-3;
       0.180109473648234133944025e-3;
       -0.910271534839798903976819e-4;
       -0.142460091487209173987685e-3;
       -0.521208962926222839685499e-4;
       0.253553321898728396521221e-4;
       0.357360245119576212047157e-4;
       0.120336367723061296449446e-4;
       -0.571467368238561493453327e-5;
       -0.751083430778512389626646e-5];
  elseif k==4
    ck=[-0.455709396340424418035779e-3;
       -0.571724624632071801184617e-3;
       0.476698510527782146839932e-3;
       0.145362409347314380791410e-2;
       0.875271889787213278878581e-3;
       -0.584212244697429717425129e-3;
       -0.119863142824806225281440e-2;
       -0.553361926292205772487962e-3;
       0.329748357174476289919883e-3;
       0.560846441567589440030798e-3;
       0.223805653823594922878522e-3;
       -0.124764472641943425453482e-3;
       -0.189331882520153576169346e-3];
  elseif k==5
    ck=[0.811234305193098332019330e-3;
       0.116042840661402233286282e-2;
       -0.120961038552699514663161e-2;
       -0.403337872044996403015012e-2;
       -0.269322476027294353847792e-2;
       0.211422330764814095297867e-2;
       0.472359714331690567866988e-2;
       0.238370161988457179934623e-2;
       -0.162089313723573923967785e-2;
       -0.297606926592315334080523e-2];
  end
D= sumexpansion(ck, eta, jmax, delta);
else
  D=Dkexact(k,zeta,f);
end

function C=Ckzeta(zeta, f, k) 
two13=2^(1/3);
if abs(zeta) < 0.5
  eta= zeta/two13; 
  jmax= 18 - 2*k; 
  delta= eps*10^(2*k);
  if k==0
    ck=[0.1;
       0.200000000000000000000000e-1;
       -0.330158730158730158730159e-2;
       -0.296911976911976911976912e-2;
       -0.137695954838811981669124e-3;
       0.555048929906072763215620e-3;
       0.242456146115529869031270e-3;
       -0.154566037237636020689912e-4;
       -0.560494560611812641811244e-4;
       -0.205938726629285257807890e-4;
       0.260832410695472146322360e-5;
       0.532729210222188857946138e-5;
       0.178881182990467430870839e-5;
       -0.300617408480663364532565e-6;
       -0.495357296368007386699085e-6;
       -0.157154550860651784075191e-6;
       0.310582008908507027210460e-7;
       0.455703940740490910695069e-7;
       0.138952925679589253258334e-7;
       -0.305767396855684609868392e-8;
       -0.416704885308313929164001e-8;
       -0.123351938981635545635728e-8;
       0.293335110754148084704709e-9;
       0.379626674718002001131840e-9];
  elseif k==1
    ck=[-0.136652236652236652236652e-2;
       -0.273304473304473304473304e-3;
       0.780378351806923235494664e-3;
       0.480481095831235887258296e-3;
       -0.149133522849003299556727e-3;
       -0.300753222613497123442870e-3;
       -0.116371247406915680355006e-3;
       0.396104305305290465475782e-4;
       0.596307815250464462047999e-4;
       0.198810673627150872495934e-4;
       -0.711571245181156693105293e-5;
       -0.941398869171534444964037e-5;
       -0.289613648251761565880795e-5;
       0.106874541035013210932373e-5;
       0.131276867676816411452532e-5;
       0.383761251199166055079476e-6;
       -0.144594124257249176458698e-6;
       -0.169190789079476125912616e-6;
       -0.477325260434205248333610e-7;
       0.182632080752898752616217e-7;
       0.206465177093998576833282e-7];
  elseif k==2
    ck=[0.301615299318380550873548e-3;
       0.603230598636761101747096e-4;
       -0.391597684545969512977759e-3;
       -0.289790085717604851031592e-3;
       0.148338746317245900757649e-3;
       0.318187930024330786834948e-3;
       0.140863378363957417331528e-3;
       -0.661957688156824128532334e-4;
       -0.107844623154343122498834e-3;
       -0.400103296947872344389066e-4;
       0.182052411474350468881140e-4;
       0.259353860787712690660584e-4;
       0.871634320891265678415106e-5;
       -0.389911063231578420704130e-5;
       -0.511759079066164065972108e-5;
       -0.161273941555264126239528e-5;
       0.714125189263839431850448e-6;
       0.886393151169713542004539e-6];
  elseif k==3
    ck=[-0.191582461090933800468500e-3;
       -0.383164922181867600937000e-4;
       0.412025391628760368304863e-3;
       0.354435556191377360085550e-3;
       -0.244344099186617866886315e-3;
       -0.577690224079697003990733e-3;
       -0.288229407397989651107383e-3;
       0.167679384035986552008787e-3;
       0.299177428518370852808697e-3;
       0.122592827219627459365658e-3;
       -0.660336620798892415061247e-4;
       -0.101973509296252262006631e-3;
       -0.373065702856469195365027e-4;
       0.192016880205339806075821e-4;
       0.270693197224160526115489e-4];
  elseif k==4
    ck=[0.240291372681993458788730e-3;
       0.480582745363986917577461e-4;
       -0.745008930228873560876256e-3;
       -0.724897584623060371465682e-3;
       0.618191768954996746066377e-3;
       0.160403374732237116145661e-2;
       0.887021131887307152191308e-3;
       -0.606942883053824288153770e-3;
       -0.117879817494336191988832e-2;
       -0.527593162723709314899751e-3;
       0.324720140692314086985285e-3;
       0.540845643878935522856430e-3];
  elseif k==5
    ck=[-0.499090772787701333229479e-3;
       -0.998181545575402666458959e-4;
       0.206206753490600447118091e-2;
       0.222570153585730233716748e-2;
       -0.223914974046656858808991e-2;
       -0.632084572564640421312756e-2;
       -0.382265214430849401331132e-2;
       0.298287351772147095269266e-2;
       0.625707189095209801465069e-2];
  end
C= two13^2*sumexpansion(ck, eta, jmax, delta);
else
  C=Ckexact(k,zeta,f);
end

function f=phizeta(zeta, z)
two13=2^(1/3);
if abs(zeta) < 0.2
  eta= zeta/two13;
  phik=[1.0;
     0.200000000000000000000000;
     0.257142857142857142857143e-1;
     -0.565079365079365079365079e-2;
     -0.393679653679653679653680e-2;
     -0.520936206650492364778079e-3;
     0.370854126473174092221711e-3;
     0.212382784029362660815242e-3;
     0.215062978875314527278662e-4;
     -0.263690406298179894036302e-4;
     -0.140546982649312867820944e-4;
     -0.114932889902944139647916e-5;
     0.197264119393862359423994e-5;
     0.101432430559332922363478e-5;
     0.703433110019219963923458e-7;
     -0.152504477739267649835292e-6;
     -0.767786625690057235680721e-7;
     -0.466984263869301799437010e-8;
     0.120667364596532865221556e-7;
     0.599087766809234353800582e-8;
     0.326910215007771473164964e-9;
     -0.971383502444283350951986e-9;
     -0.477459342952322338344740e-9;
     -0.237503186428391557787830e-10;
     0.792445981091066555670358e-10;
     0.386535842308178655279643e-10;
     0.177341058464268078733321e-11];
  f= two13*sumexpansion(phik,eta, 18, eps);
else
  f=(4*zeta/(1-z)/(1+z))^(1/4);
end

function y=zetaz(z) 
iz=0;
if z == 1 
  y=0;
  iz=1;
else
  if z < 1 
    theta=acosh(1/z); 
    theta2= theta^2;
  else
    theta= acos(1/z); 
    theta2= -theta^2;
  end
end
if iz<1
  if z < 0.85  
    y= (3*(theta-tanh(theta))/2)^(2/3); 
  elseif z > 1.15 
    y= -(3*(tan(theta)-theta)/2)^(2/3); 
  else 
    s= 1+theta2/10+theta2^2/280+theta2^3/15120+...
       theta2^4/1330560 + theta2^5/172972800 +...   
       theta2^6/31135104000 + theta2^7/7410154752000 +...
     theta2^8/2252687044608000;
    y= theta2*(s*z/2)^(2/3);
  end
end
function s=sumexpansion(ck, z, kmax, eps)
%summation ck(k)*z^k, k=1,kmax+1; 
%stop when abs(ck(k)/s) < eps
d= 1; k= 3; zk= z; 
s= ck(1)+ck(2)*z;
while (d > eps)&&(k <= kmax+1) 
  zk= z*zk; term= ck(k)*zk;
  s= s + term; 
  d= abs(term/s); 
  k= k+1; 
end
 
function [ai,bi]=aibi(x)
piquart=0.78539816339744830962;
sqrt3=1.732050807568877;
if x < -5
  ar=[+1.1280094441222097012e+0;
     -1.836709914558369e-4;
     +1.1677286575962e-6;
     -2.19524008170e-8;
     +7.698908177e-10;
     -4.09915521e-11;
     +2.9533495e-12;
     -2.679714e-13;
     +2.91808e-14;
     -3.6863e-15;
     +5.269e-16;
     -8.36e-17;
     +1.45e-17;
     -2.7e-18;
     +5.0e-19;
     -1.0e-19];
  br=[+7.7988242533135881502e-2;
     -1.83421439318478417e-4;
     +2.226279540000655e-6;
     -6.0307865841937e-8;
     +2.708298845975e-9;
     -1.73155804053e-10;
     +1.4396760940e-11;
     -1.467799419e-12;
     +1.76234030e-13;
     -2.4205568e-14;
     +3.721795e-15;
     -6.30106e-16;
     +1.15949e-16;
     -2.2951e-17;
     +4.846e-18;
     -1.084e-18;
     +2.55e-19;
     -6.3e-20;
     +1.6e-20;
     -4.0e-21;
     +1.0e-21];
  z=(2.0/3.0)*(-x)^(1.5); la= 5.0; 
  t=-2*(la/x)^3-1.0;
  s=sin(z+piquart); c=cos(z+piquart);
  f1=chepolsumC(15,t,ar);
  f2=chepolsumC(20,t,br);
  y=(-x)^(-0.25);
  ai=y*(s*f1-c*f2/z);
  bi=y*(c*f1+s*f2/z);
elseif x <0  
  ar=[+1.5217399882834532859e-1;
     +9.463304395658582345e-2;
     +8.835933286664339033e-2;
     -1.8809032021891572589e-1;
     +7.873696950590187481e-2;
     -1.616162609419079569e-2;
     +2.03792657403144947e-3;
     -1.7619321508091265e-4;
     +1.116557634724685e-5;
     -5.4278037289400e-7;
     +2.092500233007e-8;
     -6.5627359901e-10;
     +1.708749758e-11;
     -3.7550458e-13;
     +7.06091e-15;
     -1.1494e-16;
     +1.64e-18;
     -2.1e-20];
  br=[+7.5365743540308812655e-2;
     -4.3238169422348489365e-2;
     +9.5247283336721394897e-2;
     -6.1705567716412224149e-2;
     +1.7697290707409225048e-2;
     -2.909533805902076483e-3;
     +3.13306256975299299e-4;
     -2.3946459531381245e-5;
     +1.370954782252896e-6;
     -6.1112783503340e-8;
     +2.184145572027e-9;
     -6.4038637539e-11;
     +1.569118251e-12;
     -3.2625307e-14;
     +5.83052e-16;
     -9.054e-18;
     +1.23e-19;
     -1.5e-21];
  la= 5.0; 
  t= -2*(x/la)^(3)-1.0;
  f1=chepolsumC(17,t,ar);
  f2=chepolsumC(17,t,br);
  ai=f1-x*f2;
  bi=sqrt3*(f1+x*f2);
elseif x==0.0                                                  
  ai=0.35502805388781723926;
  bi=0.61492662744600073515;
elseif x<4.5                           
  ar=[+9.7934675132923263719e-1;
     -5.817935703328346930e-2;
     -1.8573106989330210396e-1;
     +2.466671790187645891e-2;
     +1.523182915947318079e-2;
     -3.18268812403824637e-3;
     -6.4421595859367550e-4;
     +2.1680741331029091e-4;
     +1.036058369682457e-5;
     -9.16356774348920e-6;
     +3.6678407013177e-7;
     +2.5198203640692e-7;
     -2.903896001941e-8;
     -4.22648003774e-9;
     +9.7569578821e-10;
     +2.168534257e-11;
     -2.083233561e-11;
     +1.00889184e-12;
     +2.9268193e-13;
     -3.543870e-14;
     -2.19178e-15;
     +6.4842e-16;
     -1.036e-17;
     -7.75e-18;
     +5.9e-19;
     +6.0e-20];
  br=[+7.2340133270078824597e-1;
     -6.186433877265966711e-2;
     +1.6859154747460945662e-1;
     -1.516452843795166719e-2;
     +1.268194063516318154e-2;
     +2.37665412524311142e-3;
     -1.35486356206889915e-3;
     +9.4811517986800841e-4;
     -3.0184629675683576e-4;
     +8.26176092954571e-5;
     -1.29237026410167e-5;
     +3.86277876708896e-7;
     +6.7349791262055e-7;
     -2.6650662327003e-7;
     +6.313398211272e-8;
     -9.28482785276e-9;
     +4.2243068981e-10;
     +2.2149462027e-10;
     -7.885556616e-11;
     +1.489486494e-11;
     -1.61561338e-12;
     -1.99887e-15;
     +4.310890e-14;
     -1.026039e-14;
     +1.36524e-15;
     -7.870e-17;
     -1.152e-17;
     +3.99e-18;
     -6.0e-19;
     +5.0e-20];
  la= 4.5; t=2*x/la-1.0;
  f1=chepolsumC(25,t,ar);
  f2=chepolsumC(29,t, br);
  ai=f1*exp(-1.5*x);
  bi=f2*exp(1.375*x);
elseif x<9
  ar=[+5.0491255733788020503e+1;
     -6.257510019590865924e+0;
     -5.510390933031668771e+0;
     +9.23785747211947637e-1;
     +2.69701933599952531e-1;
     -6.0503610861932189e-2;
     -7.272383379486571e-3;
     +2.432739101814126e-3;
     +8.9407383924725e-5;
     -6.7626485765665e-5;
     +1.143189712477e-6;
     +1.369896864000e-6;
     -8.37850105221e-8;
     -2.0445153671e-8;
     +2.247360374e-9;
     +2.14042502e-10;
     -4.0617972e-11;
     -1.169040e-12;
     +5.51909e-13;
     -8.169e-15;
     -5.801e-15;
     +3.09e-16;
     +4.6e-17;
     -5.0e-18];
  br=[+5.1865644423401596677e-3;
     +3.190078113025389266e-4;
     +6.589260266019625254e-4;
     +1.77323651950275976e-5;
     +4.24724564662038450e-5;
     -3.095767914834698e-7;
     +1.8968968954113726e-6;
     -8.26612616574416e-8;
     +6.82206567197979e-8;
     -5.5232442059920e-9;
     +2.1841788459370e-9;
     -2.532110361946e-10;
     +6.69024065663e-11;
     -9.6266259764e-12;
     +2.0421501469e-12;
     -3.287463020e-13;
     +6.17455626e-14;
     -1.01838450e-14;
     +1.7333508e-15;
     -2.689068e-16;
     +3.97206e-17;
     -5.1324e-18;
     +5.434e-19;
     -3.01e-20;
     -5.0e-21;
     +2.3e-21;
     -5.0e-22;
     +1.0e-22];
  la=4.5; t=2*(x-la)/(9.0-la)-1.0;
  f1=chepolsumC(23,t,ar);
  f2=chepolsumC(27,t,br);
  ai=f1*exp(-2.5*x);
  bi=f2*exp(2.5*x);
else
  ar=[+5.6312443149444830912e-1;
     -5.2879891960803332e-4;
     +3.72709673104629e-6;
     -4.905124032702e-8;
     +9.3521089895e-10;
     -2.313916112e-11;
     +6.9870927e-13;
     -2.476246e-14;
     +1.00261e-15;
     -4.547e-17;
     +2.28e-18;
     -1.2e-19];
  br=[+1.1306068068318665650e+0;
     +1.1225343035826004e-3;
     +8.8491541977453e-6;
     +1.379157291154e-7;
     +3.3014426995e-9;
     +1.089423404e-10;
     +4.6726552e-12;
     +2.514325e-13;
     +1.66086e-14;
     +1.3312e-15;
     +1.292e-16;
     +1.53e-17;
     +2.2e-18;
     +4.0e-19];
  la= 9.0; y=(x)^(-0.25); 
  z= 2.0/3.0*(x)^1.5;
  t=2*(la/x)^1.5-1;
  f1=chepolsumC(11,t,ar);
  f2=chepolsumC(13, t, br);
  la=exp(z);
  ai=y/la*f1;
  bi=y*la*f2;
end

function  [aip,bip]=aibip(x)
piquart=0.78539816339744830962;
sqrt3=1.732050807568877;
if x<-5 
  cr=[+1.1288167982669480212e+0;
     +2.175180808271174e-4;
     -1.2733503718850e-6;
     +2.33020578931e-8;
     -8.064911115e-10;
     +4.26054712e-11;
     -3.0537799e-12;
     +2.760646e-13;
     -2.99790e-14;
     +3.7789e-15;
     -5.392e-16;
     +8.55e-17;
     -1.48e-17;
     +2.8e-18;
     -6.0e-19;
     +1.0e-19];
  dr=[-1.0928773748447830587e-1;
     +2.0544241243325013e-4;
     -2.38739034179860e-6;
     +6.350265424372e-8;
     -2.82350803428e-9;
     +1.7939539824e-10;
     -1.485151455e-11;
     +1.50940584e-12;
     -1.8079342e-13;
     +2.478440e-14;
     -3.80487e-15;
     +6.4333e-16;
     -1.1825e-16;
     +2.339e-17;
     -4.93e-18;
     +1.10e-18;
     -2.60e-19;
     +6.0e-20;
     -2.0e-20];
  z= 2.0/3.0*(-x)^1.5; 
  la=5.0; t= -2*(la/x)^3.0-1.0;
  s=sin(z+piquart); c=cos(z+piquart);
  f3=chepolsumC(15,t,cr);
  f4=chepolsumC(18,t,dr);
  y=(-x)^0.25;
  aip=-y*(c*f3+s*f4/z);
  bip=y*(s*f3-c*f4/z);
elseif x < 0
  cr=[+5.2961692024697341466e-2;
     -4.8382429177635177789e-2;
     +6.2046464244529580518e-2;
     -3.1417437267239646847e-2;
     +7.876452021481851462e-3;
     -1.182440976973326921e-3;
     +1.18871496270269531e-4;
     -8.595270331212026e-6;
     +4.69655735896232e-7;
     -2.0107696526448e-8;
     +6.93493715818e-10;
     -1.9694289584e-11;
     +4.68795260e-13;
     -9.492371e-15;
     +1.65543e-16;
     -2.513e-18;
     +3.4e-20];
  dr=[+1.9598877757490956561e-1;
     +2.0709937287929773169e-1;
     -9.618436611491528534e-2;
     -2.6927080774066725642e-1;
     +1.5311467164195350955e-1;
     -3.621663427173734531e-2;
     +5.00970025411579034e-3;
     -4.6418943778727143e-4;
     +3.110071191129936e-5;
     -1.58422527393619e-6;
     +6.359244256854e-8;
     -2.06683376305e-9;
     +5.556077451e-11;
     -1.25683910e-12;
     +2.426805e-14;
     -4.0481e-16;
     +5.89e-18;
     -7.6e-20];
  la= 5.0; t= -2.0*(x/la)^3.0-1.0; z= x*x;
  f3=chepolsumC(16,t,cr);
  f4=chepolsumC(17,t,dr);
  aip=z*f3-f4;
  bip=sqrt3*(z*f3+f4);
elseif x==0.0
  aip=-0.25881940379280679841;
  bip=0.44828835735382635791;
elseif x < 4.5
  cr=[-1.0838288082415277302e+0;
     -2.75303325521307795e-2;
     +2.630996799464136716e-1;
     -2.20070838039707907e-2;
     -2.65686371901896143e-2;
     +4.6607302595359725e-3;
     +1.3331627500919505e-3;
     -3.993412832236136e-4;
     -2.57477942743769e-5;
     +1.92965112669086e-5;
     -8.179057497221e-7;
     -5.675046562446e-7;
     +7.38399111711e-8;
     +9.1592819462e-9;
     -2.5909300717e-9;
     -5.596039e-13;
     +5.48262416e-11;
     -4.2085028e-12;
     -6.882068e-13;
     +1.230089e-13;
     +2.2115e-15;
     -2.0268e-15;
     +1.092e-16;
     +2.01e-17;
     -2.8e-18];
      dr=[+6.5492503637358356971e-1;
         +7.398273602101549954e-2;
         +1.3916014396498561737e-1;
         -1.318037729688479868e-2;
         +2.721562752659802946e-2;
         -8.95735391942756359e-3;
         +4.64318555084023870e-3;
         -1.31982848083974091e-3;
         +3.2861337599865999e-4;
         -2.823577298366577e-5;
         -1.263661132378531e-5;
         +8.58819079977971e-6;
         -2.83605845640361e-6;
         +6.4357170147486e-7;
         -9.054345148828e-8;
         +1.02361813059e-9;
         +3.75568593189e-9;
         -1.26430133998e-9;
         +2.5000891998e-10;
         -2.967290357e-11;
         +2.6790190e-13;
         +8.3103543e-13;
         -2.2570250e-13;
         +3.487651e-14;
         -2.83929e-15;
         -1.5678e-16;
         +9.941e-17;
         -1.867e-17;
         +2.01e-18;
         -6.0e-20;
         -2.0e-20];
      la= 4.5; t= 2*x/la-1.0;
      f3=chepolsumC(24, t, cr);
      f4=chepolsumC(30, t, dr);
      aip=f3*exp(-1.375*x);
      bip=f4*exp(1.5*x);
elseif x < 9
      cr=[-1.2958123754838538881e+2;
         +6.76832374624909949e+0;
         +1.598511024718904852e+1;
         -1.38866512314607882e+0;
         -9.2855057862186510e-1;
         +1.1311806368323397e-1;
         +3.279015099087572e-2;
         -5.43676653553677e-3;
         -7.5136921782930e-4;
         +1.7836159217021e-4;
         +1.030315382664e-5;
         -4.29105073708e-6;
         -2.395648053e-8;
         +7.851108600e-8;
         -2.78229868e-9;
         -1.10405576e-9;
         +8.374716e-11;
         +1.165092e-11;
         -1.51205e-12;
         -8.180e-14;
         +2.018e-14;
         +1.2e-16;
         -2.1e-16];
      dr=[+5.5860183624034104454e-3;
         +2.87080283197283049e-5;
         +6.566960781790676321e-4;
         -1.29004618056797811e-5;
         +3.94734989962719105e-5;
         -1.7468162471848648e-6;
         +1.6312939322265420e-6;
         -1.095384680289487e-7;
         +5.23022431650618e-8;
         -4.4609927148093e-9;
         +1.3623610617155e-9;
         -1.277857258067e-10;
         +2.77519636083e-11;
         -2.2855189985e-12;
         +3.227684345e-13;
         +2.5537222e-15;
         -4.8715113e-15;
         +1.8406926e-15;
         -3.499972e-16;
         +5.20075e-17;
         -3.2250e-18;
         -9.750e-19;
         +5.043e-19;
         -1.464e-19;
         +3.39e-20;
         -6.7e-21;
         +1.2e-21;
         -2.0e-22];
      la= 4.5; t= 2*(x-la)/(9.0-la)-1.0;
      f3= chepolsumC(22, t, cr);
      f4= chepolsumC(27, t, dr);
      aip= f3*exp(-2.5*x);
      bip= f4*exp(2.625*x);
else
      cr=[+5.6568578662442917280e-1;
         +7.4362138971729004e-4;
         -4.42404425936355e-6;
         +5.505791762929e-8;
         -1.02083536642e-9;
         +2.483795569e-11;
         -7.4166337e-13;
         +2.607525e-14;
         -1.04943e-15;
         +4.738e-17;
         -2.36e-18;
         +1.3e-19];
      dr=[+1.1252717606343474943e+0;
         -1.5639596736331532e-3;
         -1.04063425452027e-5;
         -1.533581036012e-7;
         -3.5697709625e-9;
         -1.158278459e-10;
         -4.9121069e-12;
         -2.621683e-13;
         -1.72103e-14;
         -1.3725e-15;
         -1.326e-16;
         -1.56e-17;
         -2.2e-18;
         -4.0e-19];
      la= 9.0; y=(x)^0.25; 
      z=2.0/3.0*(x)^1.5;
      t=2*(la/x)^1.5-1.0;
      f3=chepolsumC(11,t,cr);
      f4=chepolsumC(13,t,dr);
      la=exp(z);
      aip=-y/la*f3;
      bip= y*la*f4;
end
  
function [chep]=chepolsumC(n,t,ak)
d2=0.0;
half=0.5;
tval=0.6;
u1=0.0;
u2=0.0;
if abs(t) < tval %Clenshaw method
  u0=0;
  tt=t + t;
  for k=n:-1:1
    u2=u1;
    u1=u0;
    u0=tt*u1+ak(k)-u2;
  end
  chep=(u0-u2)/2;
else
  d1=0;
  if t>0
    tt=(t-half)-half;
    tt=tt + tt;
    for k=n:-1:1
      d2=d1;
      u2=u1;
      d1=tt*u2+ak(k)+d2;
      u1=d1+u2;
    end
    chep=(d1+d2)/2;
  else
    tt=(t+half)+half;
    tt=tt+tt;
    for k=n:-1:1
      d2=d1;
      u2=u1;
      d1=tt*u2+ak(k)-d2;
      u1=d1-u2;
    end
    chep=(d1-d2)/2;
  end
end

     

function [s,c]=sincospi(a)
% sin(pi*a) and cos(pi*a) with exact argument reduction:
% a=n/2+r, |r|<=1/4, so that both values have small relative error
n=round(2*a);
r=a-n/2;
sr=sin(pi*r);
cr=cos(pi*r);
m=mod(n,4);
if m==0
  s=sr; c=cr;
elseif m==1
  s=cr; c=-sr;
elseif m==2
  s=-sr; c=-cr;
else
  s=-cr; c=sr;
end

function w=fixov(v)
% Inf -> sign*realmax, NaN -> realmax
w=v;
if isnan(v)
  w=realmax;
elseif isinf(v)
  w=sign(v)*realmax;
end
