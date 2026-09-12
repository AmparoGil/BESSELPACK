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
%        ier=0, computation succesful
%        ier=1, underflow/overflow values
%               if a>0 
%                 J, Jp are set to zero
%                 Y, Yp are set to realmax
%               if a<0
%                 J, Jp are set to sign(sin(a*pi))*realmax
%                 Y, Yp are set to sign(cos(a*pi))*realmax
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
  % Special case a = ±0.5
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
    n=round(a);
    epsil=a-n;  
    if abs(epsil)<0.1
      s=(-1)^n*sin(pi*epsil);
    else
      s=sin(a*pi);
    end
    c=cos(a*pi);
    [ia,yax,iad,yap,ier]=bessel(-a, x);
    if ier==0  
      J=c*ia+s*yax;
      Jp=c*iad+s*yap;
      Y=-s*ia+c*yax;
      Yp=-s*iad+c*yap;
    else
      ss=sign(s);
      cs=sign(c);
      J=ss*realmax;
      Jp=ss*realmax;
      Y=cs*realmax;
      Yp=cs*realmax;
    end
  else
    xmax=3*sqrt((1+a));
    a1=1.41;
    x1=0.38*a;
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
      elseif x > 1.9*(a/5.3)^2 + 20
        [J,Y,Jp,Yp,ier]=expax(a,x);
%--------------------------------------------------------------
% Debye, x>a
%--------------------------------------------------------------
      elseif (a>45) && (x > 1.35*a)
        [J,Jp,Y,Yp]=deb2(a,x);
%--------------------------------------------------------------
% Recurrences for J
%--------------------------------------------------------------      
      else
        [J,Jp]=jrec(a,x);
        [Y,Yp]=bessyax(a,x);
      end  
    else  
      ier=1;
      J=0;
      Jp=0;
      Y=realmax;
      Yp=realmax;   
    end	
  end
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
  yax=realmax;
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
eps=5.e-16;
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
  ia=exp(a*(log(x/(2.0*a))+1.0))...
    /(sqrt(2.0*pi*a)*gamstar(a));     
  r=x*x/4.0;
  c=1.0;
  s=a;
  n=1;
  jax=0.0;
  jap=0.0;              
  while abs(s)>eps 
    jax=jax+c;
    jap=jap+s;
    c=-c*r/(n*(a+n));
    s=(a+2.0*n)*c;
    n=n+1;
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
phi=-a*pihalf-piquart;
cosw=cos(x)*cos(phi)-sin(x)*sin(phi);
sinw=sin(x)*cos(phi)+cos(x)*sin(phi);
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
while ((err1 > eps)&&(err1p >eps))&&(m<500)
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
  jax=0.0; yax=realmax;
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
    jax=0.0; yax=realmax;
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
  c=cos(a*pi);
  s=sin(a*pi);
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
  if (x<8)
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
    for n=1:na       
      h=b*(r+n)*g-f;
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
    for n=1:na       
      h=b*(r+n)*g-f;
      f=g;
      g=h;
    end
  end  
  yax=f;
  yap=a*f/x-g;
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
elseif (x >= 5) 
  c=0.25-r*r;
  b=x+x;
  f=1.0;
  g=1.0;
  p=1.0;
  q=0.0;
  n=startpqbes(x, eps);
  n=n+20;
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
  
function st=startpqbes(x,epss)
mactol=eps;
if (epss<mactol) 
  del=-log(mactol/2.0);
else
  del=-log(epss/2.0);
end
t=del/(2.0*x);
if (t<0.5) 
  s=1.0/(2.0*t);
  c=sqrt(1.0+s*s);
  a=log(s+c)+s*c/(1.0+c*c);
  a=a*(1.0+7.16*t)/(1.0+5.33*t);
else
  q=1.0/t*t;
  a=(1.0+q*q*(-7.0/360.0+4.1e-3*q))/t;
end
b=0.0;
while (abs(b/a-1.0)>1.0e-2)
  b=a;
  [s,c]=hypfun(a);
  a=a+(a*c*s+s*s*(1.0-2.0*t*s))/(a*(1.0+c*c));
end  
[s,c]=hypfun(a);
st=1.0+floor(x*c/s*s);

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
  dr(1)=-1.013609258009865776949;
  dr(2)=0.784903531024782283535e-1;
  dr(3)=0.67588668743258315530e-2;
  dr(4)=-0.12790434869623468120e-2;
  dr(5)=0.462939838642739585e-4;
  dr(6)=0.43381681744740352e-5;
  dr(7)=-0.5326872422618006e-6;
  dr(8)=0.172233457410539e-7;
  dr(9)=0.8300542107118e-9;
  dr(10)=-0.10553994239968e-9;
  dr(11)=0.39415842851e-11;
  dr(12)=0.362068537e-13;
  dr(13)=-0.107440229e-13;
  dr(14)=0.5000413e-15;
  dr(15)=-0.62452e-17;
  dr(16)=-0.5185e-18;
  dr(17)=0.347e-19;
  dr(18)=-0.9e-21;
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
  c(1)=+1.142022680371167841;
  c(2)=-6.5165112670736881e-3;
  c(3)=-3.087090173085368e-4;
  c(4)=+3.4706269649043e-6;
  c(5)=-6.9437664487e-9;
  c(6)=-3.67795399e-11;
  c(7)=+1.356395e-13;
  c(8)=+3.68e-17;
  c(9)=-5.5e-19;
  tx=2.0*x;
  t=2*tx*tx-1;
  q=chepolsum(8,t,c);
  c(1)=-1.270583625778727532;
  c(2)=+2.05083241859700357e-2;
  c(3)=-7.84761097993185e-5;
  c(4)=-5.377798984020e-7;
  c(5)=+3.8823289907e-9;
  c(6)=-2.6758703e-12;
  c(7)=-2.39860e-14;
  c(8)=+3.80e-17;
  c(9)=+4e-20;
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
  s=sin(a*pi);
  c=cos(a*pi);
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
    while (abs(s)>eps)
      jax=jax+c;
      jap=jap+s;
      c=-c*r/(n*(a+n));
      s=(a+2.0*n)*c;
      n=n+1;
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
nu2p2=nu2*nu2;
nu2p3=nu2p2*nu2;
nu2p4=nu2p3*nu2;
nu2p5=nu2p4*nu2;
w= zeta*nu23;
sA= 1+Akzeta(zeta,f,1)/nu2+Akzeta(zeta,f,2)/nu2p2+...
    Akzeta(zeta,f,3)/nu2p3+Akzeta(zeta,f,4)/nu2p4+...
    Akzeta(zeta,f,5)/nu2p5;
sB= Bkzeta(zeta,f,0)+Bkzeta(zeta,f,1)/nu2+...
    Bkzeta(zeta,f,2)/nu2p2+Bkzeta(zeta,f,3)/nu2p3+...
    Bkzeta(zeta,f,4)/nu2p4;
sC= Ckzeta(zeta,f,0)+Ckzeta(zeta,f,1)/nu2+...
    Ckzeta(zeta,f,2)/nu2p2+Ckzeta(zeta,f,3)/nu2p3+...
    Ckzeta(zeta,f,4)/nu2p4;
sD= 1+Dkzeta(zeta,f,1)/nu2+Dkzeta(zeta,f,2)/nu2p2+...
    Dkzeta(zeta,f,3)/nu2p3+Dkzeta(zeta,f,4)/nu2p4+...
    Dkzeta(zeta,f,5)/nu2p5;
[Ai,Bi]=aibi(w);
[Aid,Bid]=aibip(w);
J= f*(Ai/nu13*sA+Aid/(nu23*nu)*sB);
Y=-f*(Bi/nu13*sA+Bid/(nu23*nu)*sB);
Jp=-fhat*(Ai/(nu*nu13)*sC+Aid/nu23*sD);
Yp=fhat*(Bi/(nu*nu13)*sC+Bid/nu23*sD);

function A=Akzeta(zeta, f, k)
two13=2^(1/3);
%ki=0, first component. Changed to ki=1
Ajk(1,1)= -0.444444444444444444444444e-2;
Ajk(2,1)= -0.184415584415584415584416e-2;
Ajk(3,1)= 0.112136752136752136752137e-2;
Ajk(4,1)= 0.134577521244187910854578e-2;
Ajk(5,1)= 0.388062656297950415597474e-3;
Ajk(6,1)= -0.183068672378179913251358e-3;
Ajk(7,1)= -0.199546088780673281140134e-3;
Ajk(8,1)= -0.525619123404158729343881e-4;
Ajk(9,1)= 0.246061965245915785195825e-4;
Ajk(10,1)= 0.251924678092454136168976e-4;
Ajk(11,1)= 0.633315737653324223078736e-5;
Ajk(12,1)= -0.295748573383020170629872e-5;
Ajk(13,1)= -0.292525592056483802975188e-5;
Ajk(14,1)= -0.715970261050200926824325e-6;
Ajk(15,1)= 0.333151072039094907664420e-6;
Ajk(16,1)= 0.322767047569230987018795e-6;
Ajk(17,1)= 0.776772938166419874312952e-7;
Ajk(18,1)= -0.360095423792111979834763e-7;
Ajk(19,1)= -0.344172444903422642801444e-7;
Ajk(20,1)= -0.818819435639877179627180e-8;
Ajk(21,1)= 0.378314848515203785605629e-8;
Ajk(22,1)= 0.358157567915293191299802e-8;
%----------------------
Ajk(1,2)= 0.693735541354588973636593e-3;
Ajk(2,2)= 0.464483490365843307019778e-3;
Ajk(3,2)= -0.428381301715351124588504e-3;
Ajk(4,2)= -0.702670286877113311390726e-3;
Ajk(5,2)= -0.263258004677881193972489e-3;
Ajk(6,2)= 0.166385366628870309272117e-3;
Ajk(7,2)= 0.221208768781858329671256e-3;
Ajk(8,2)= 0.702034561532966065008509e-4;
Ajk(9,2)= -0.400042178254061367229691e-4;
Ajk(10,2)= -0.478632496645396180982353e-4;
Ajk(11,2)= -0.139460074147363092319545e-4;
Ajk(12,2)= 0.753618659127372677842368e-5;
Ajk(13,2)= 0.847850216106766137159346e-5;
Ajk(14,2)= 0.234535522845391156499770e-5;
Ajk(15,2)= -0.122594329471088266135048e-5;
Ajk(16,2)= -0.132508234340102697171799e-5;
Ajk(17,2)= -0.353995477656999743787481e-6;
Ajk(18,2)= 0.180829171937667380735510e-6;
Ajk(19,2)= 0.190038351523365606575570e-6;
%----------------------
Ajk(1,3)= -0.354211971457743840771126e-3;
Ajk(2,3)= -0.312322527890318832782775e-3;
Ajk(3,3)= 0.371644223750229630164450e-3;
Ajk(4,3)= 0.753926915597773214092907e-3;
Ajk(5,3)= 0.340830005944473690438883e-3;
Ajk(6,3)= -0.263496817206959254774719e-3;
Ajk(7,3)= -0.408927572664843008395140e-3;
Ajk(8,3)= -0.150110875956345964262796e-3;
Ajk(9,3)= 0.996401520553805536354607e-4;
Ajk(10,3)= 0.135249295575128253602523e-3;
Ajk(11,3)= 0.444311708727290144245992e-4;
Ajk(12,3)= -0.271320507191411598209920e-4;
Ajk(13,3)= -0.339679696977185933223319e-4;
Ajk(14,3)= -0.104070886527304282710648e-4;
Ajk(15,3)= 0.602463906541441186398061e-5;
Ajk(16,3)= 0.714391960784688644874689e-5;
%-------------------------
Ajk(1,4)= 0.378194199201772914026612e-3;
Ajk(2,4)= 0.404943905523632334772139e-3;
Ajk(3,4)= -0.579130526946450830904535e-3;
Ajk(4,4)= -0.138017901171009595794723e-2;
Ajk(5,4)= -0.722520056780092714288715e-3;
Ajk(6,4)= 0.651265924036820506795172e-3;
Ajk(7,4)= 0.114674563328389541610941e-2;
Ajk(8,4)= 0.474423189340400087623894e-3;
Ajk(9,4)= -0.356495172735465313237897e-3;
Ajk(10,4)= -0.538157791035117996567133e-3;
Ajk(11,4)= -0.195687390661222638146666e-3;
Ajk(12,4)= 0.132563525210291989094654e-3;
Ajk(13,4)= 0.181949256267290616373907e-3;
%----------------------------
Ajk(1,5)= -0.691141397288294167613745e-3;
Ajk(2,5)= -0.859953266117743832936210e-3;
Ajk(3,5)= 0.142023355681435114892370e-2;
Ajk(4,5)= 0.385354269956030524409690e-2;
Ajk(5,5)= 0.227528116429013745967378e-2;
Ajk(6,5)= -0.232195720345569883668239e-2;
Ajk(7,5)= -0.454786433944346356182570e-2;
Ajk(8,5)= -0.208244317582724499170221e-2;
Ajk(9,5)= 0.173704435731958087183915e-2;
Ajk(10,5)= 0.287424443211932650161305e-2;
if abs(zeta) < 0.3
  eta= zeta/two13; 
  jmax= 18 - 2*k; 
  delta= eps*10^(2*k);
  for j=0:jmax 
    ck(j+1)= Ajk(j+1,k); 
  end 
  A= sumexpansion(ck, eta, jmax, delta);
else
  A=Akexact(k,zeta,f);
end

function A=Akexact(k,zeta, f) 
if k==1
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  A= 1/73728*(385*f12-1848*f8*zeta+280*f6+...
     1296*f4*zeta2-672*f2*zeta-7280)/zeta3;
elseif (k==2) 
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  A= 0.6132165100e-11*(185910725.*f24-1784742960.*f20*zeta+...
     47647600.*f18+5598758880.*f16*zeta2-...
     343062720.*f14*zeta-6023787264.*f12*zeta3-...
     84084000.*f12+662328576.*f10*zeta2+...
     1143072000.*f8*zeta4+403603200.*f8*zeta-...
     217728000.*f6*zeta3+608608000.*f6-...
     283046400.*f4*zeta2-1460659200.*f2*zeta-...
     0.5173168000e11)/zeta6;
elseif (k==3)
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  zeta9=zeta8*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f28=f26*f2;
  f30=f26*f4;
  f32=f30*f2;
  f34=f32*f2;
  f36=f34*f2;
  A=0.7921220213e-18*(-0.1474119602e17*f32*zeta+...
    0.8167514746e17*f28*zeta2-0.1902089810e16*f26*zeta-...
    0.2156180524e18*f24*zeta3+0.8253986124e16*f22*zeta2+...
    0.2689947583e18*f20*zeta4+0.1364257519e16*f20*zeta-...
    0.1529470761e17*f18*zeta3-0.1306391535e18*f16*zeta5-...
    0.4279691288e16*f16*zeta2+0.1059866543e17*f14*zeta4+...
    0.1129287466e17*f12*zeta6-0.2609883949e16*f14*zeta+...
    0.4604582985e16*f12*zeta3-0.1306613597e16*f10*zeta5+...
    0.5038730875e16*f10*zeta2-0.8737642368e15*f8*zeta4+...
    0.1003801519e17*f8*zeta-0.1656387533e16*f6*zeta3-...
    0.7039647014e16*f4*zeta2-0.7435953069e17*f2*zeta+...
    0.1023694168e16*f36+0.1585074841e15*f30-...
    0.1421101582e15*f24+0.3624838818e15*f18-...
    0.2091253164e16*f12+0.3098313779e17*f6-...
    0.4432654246e19)/zeta9;
elseif (k==4)
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  zeta9=zeta8*zeta;
  zeta10=zeta9*zeta;
  zeta11=zeta10*zeta;
  zeta12=zeta11*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f28=f26*f2;
  f30=f26*f4;
  f32=f30*f2;
  f34=f32*f2;
  f36=f34*f2;
  f38=f36*f2;
  f40=f38*f2;
  f42=f40*f2;
  f44=f42*f2;
  f48=f44*f4;
  A=0.2685960630e-23*( -0.8608456956e22*f44*zeta + ...
    0.6836553289e23*f40*zeta2 - 0.8348627231e21*f38*zeta - ...
    0.2888762000e24*f36*zeta3 + 0.5628030280e22*f34*zeta2 + ...
    0.6942138681e24*f32*zeta4 + 0.4292636282e21*f32*zeta - ...
    0.1931230712e23*f30*zeta3 - 0.9360524888e24*f28*zeta5 - ...
    0.2378380294e22*f28*zeta2 + 0.3516459067e23*f26*zeta4 + ...
    0.6462199207e24*f24*zeta6 - 0.5512509880e21*f26*zeta + ...
    0.6278797685e22*f24*zeta3 - 0.3183721763e23*f22*zeta5 - ...
    0.1795778983e24*f20*zeta7 + 0.2392115232e22*f22*zeta2 - ...
    0.7833127361e22*f20*zeta4 + 0.1146242248e23*f18*zeta6 + ...
    0.8833609214e22*f16*zeta8 + 0.1292588524e22*f20*zeta - ...
    0.4432610195e22*f18*zeta3 + 0.3804212150e22*f16*zeta5 - ...
    0.7328623940e21*f14*zeta7 - 0.4054864839e22*f16*zeta2 + ...
    0.3071634557e22*f14*zeta4 - 0.3288485101e21*f12*zeta6 - ...
    0.5061504535e22*f14*zeta + 0.4362688892e22*f12*zeta3 - ...
    0.3786740420e21*f10*zeta5 + 0.9771913109e22*f10*zeta2 - ...
    0.8278624889e21*f8*zeta4 + 0.3276618019e23*f8*zeta - ...
    0.3212331726e22*f6*zeta3 - 0.2297887961e23*f4*zeta2 - ...
    0.3659599345e24*f2*zeta + 0.4483571331e21*f48 + ...
    0.4969420971e20*f42 - 0.2980997418e20*f36 + 0.4593758234e20*f30 - ...
    0.1346446379e21*f24 + 0.7029867409e21*f18 - 0.6826287539e22*f12 + ...
    0.1524833061e24*f6 - 0.3063389619e26 )/zeta12;
else
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  zeta9=zeta8*zeta;
  zeta10=zeta9*zeta;
  zeta11=zeta10*zeta;
  zeta12=zeta11*zeta;
  zeta13=zeta12*zeta;
  zeta15=zeta13*zeta2;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f28=f26*f2;
  f30=f26*f4;
  f32=f30*f2;
  f34=f32*f2;
  f36=f34*f2;
  f38=f36*f2;
  f40=f38*f2;
  f42=f40*f2;
  f44=f42*f2;
  f46=f44*f2;
  f48=f44*f4;
  f50=f48*f2;
  f52=f50*f2;
  f54=f52*f2;
  f56=f54*f2;
  f58=f56*f2;
  f60=f58*f2;
  A = 0.1051390263e-31*( ...
    0.2171500484e30*f44*zeta - ...
    0.1724534240e31*f40*zeta2 - ...
    0.2095929545e30*f38*zeta + ...
    0.7286959920e31*f36*zeta3 + ...
    0.1412921504e31*f34*zeta2 - ...
    0.1751168367e32*f32*zeta4 + ...
    0.3523151180e30*f32*zeta - ...
    0.4848370153e31*f30*zeta3 + ...
    0.2361211124e32*f28*zeta5 - ...
    0.1952038978e31*f28*zeta2 + ...
    0.8828098619e31*f26*zeta4 - ...
    0.1630102674e32*f24*zeta6 - ...
    0.9260854531e30*f26*zeta + ...
    0.5153279248e31*f24*zeta3 - ...
    0.7992758955e31*f22*zeta5 + ...
    0.4529888400e31*f20*zeta7 + ...
    0.4018683262e31*f22*zeta2 - ...
    0.6428984449e31*f20*zeta4 + ...
    0.2877650333e31*f18*zeta6 - ...
    0.2228295591e30*f16*zeta8 + ...
    0.3654950588e31*f20*zeta - ...
    0.7446654809e31*f18*zeta3 + ...
    0.3122280493e31*f16*zeta5 - ...
    0.1839856903e30*f14*zeta7 - ...
    0.1146562139e32*f16*zeta2 + ...
    0.5160255749e31*f14*zeta4 - ...
    0.2699001128e30*f12*zeta6 - ...
    0.2157842681e32*f14*zeta + ...
    0.1233603118e32*f12*zeta3 - ...
    0.6361612575e30*f10*zeta5 + ...
    0.4166004602e32*f10*zeta2 - ...
    0.2340881445e31*f8*zeta4 + ...
    0.1961586401e33*f8*zeta - ...
    0.1369495267e32*f6*zeta3 - ...
    0.1375657996e33*f4*zeta2 - ...
    0.2926811456e34*f2*zeta - ...
    0.5422487267e30*f50*zeta + ...
    0.4957205082e31*f46*zeta2 - ...
    0.2489020487e32*f42*zeta3 + ...
    0.7443541491e32*f38*zeta4 - ...
    0.1342231011e33*f34*zeta5 + ...
    0.1403836717e33*f30*zeta6 - ...
    0.7656347622e32*f26*zeta7 + ...
    0.1693064327e32*f22*zeta8 - ...
    0.6604896676e30*f18*zeta9 - ...
    0.6982528246e31*f56*zeta + ...
    0.7221418819e32*f52*zeta2 - ...
    0.4201114649e33*f48*zeta3 + ...
    0.1504461750e34*f44*zeta4 - ...
    0.3407569272e34*f40*zeta5 + ...
    0.4821484996e34*f36*zeta6 - ...
    0.4042757569e34*f32*zeta7 + ...
    0.1789079282e34*f28*zeta8 - ...
    0.3224453137e33*f24*zeta9 + ...
    0.1021871871e32*f20*zeta10 - ...
    0.1130989835e29*f48 + ...
    0.1247577110e29*f42 - ...
    0.2446632764e29*f36 + ...
    0.7717378776e29*f30 - ...
    0.3807240196e30*f24 + ...
    0.2997003724e31*f18 - ...
    0.4086638336e32*f12 + ...
    0.1219504773e34*f6 + ...
    0.2510410772e29*f54 + ...
    0.2909386769e30*f60 - ...
    0.3154127146e36)/zeta15;  
end

function B=Bkexact(k,zeta,f) 
if k==0   
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2; 
  B= 1/192*(5*f6-12*f2*zeta-20)/zeta^2;
elseif k==1
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f18=f14*f4;
  B= 1/212336640*(425425*f18-3063060*f14*zeta-...
    115500*f12+5913648*f10*zeta2+554400*f8*zeta-...
    1944000*f6*zeta3+462000*f6-388800*f4*zeta2-...
    1108800*f2*zeta-27227200)/zeta5;
elseif k==2
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f30=f26*f4;
  B=0.4562622843e-14*(0.1886993859e12*f30-...
    0.2264392630e13*f26*zeta-0.2602750150e11*f24+...
    0.9826173958e13*f22*zeta2+0.2498640144e12*f20*zeta-...
    0.1820798525e14*f18*zeta3+0.3668865200e11*f18-...
    0.7838262432e12*f16*zeta2+0.1261745884e14*f14*zeta4-...
    0.2641582944e12*f14*zeta+0.8433302170e12*f12*zeta3-...
    0.1555492378e13*f10*zeta5-0.1467546080e12*f12+...
    0.5099930035e12*f10*zeta2-0.1600300800e12*f8*zeta4+...
    0.7044221184e12*f8*zeta-0.1676505600e12*f6*zeta3+...
    0.1665760096e13*f6-0.4940103168e12*f4*zeta2-...
    0.3997824230e13*f2*zeta-0.1932281711e15)/zeta8;
elseif k==3  
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  zeta9=zeta8*zeta;
  zeta10=zeta9*zeta;
  zeta11=zeta10*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f28=f26*f2;
  f30=f26*f4;
  f32=f30*f2;
  f34=f32*f2;
  f36=f34*f2;
  f38=f34*f4;
  f42=f38*f4;
  B=0.4125635527e-20*(-0.3727065728e19*f38*zeta+...
    0.2512513518e20*f34*zeta2+0.2948239205e18*f32*zeta-...
    0.8621565680e20*f30*zeta3-0.1633502949e19*f28*zeta2+...
    0.1569847798e21*f26*zeta4-0.2092298791e18*f26*zeta+...
    0.4312361047e19*f24*zeta3-0.1421304358e21*f22*zeta5+...
    0.9079384737e18*f22*zeta2-0.5379895165e19*f20*zeta4+...
    0.5117152892e20*f18*zeta6+0.3401548746e18*f20*zeta-...
    0.1682417837e19*f18*zeta3+0.2612783070e19*f16*zeta5-...
    0.3271707116e19*f14*zeta7-0.1067069694e19*f16*zeta2+...
    0.1165853197e19*f14*zeta4-0.2258574932e18*f12*zeta6-...
    0.1020464624e19*f14*zeta+0.1148076024e19*f12*zeta3-...
    0.1437274957e18*f10*zeta5+0.1970143772e19*f10*zeta2-...
    0.2178585497e18*f8*zeta4+0.5356284904e19*f8*zeta-...
    0.6476475253e18*f6*zeta3-0.3756355647e19*f4*zeta2-...
    0.5031661576e20*f2*zeta+0.2218491505e18*f42-...
    0.2047388337e17*f36+0.1743582325e17*f30-...
    0.3543279944e17*f24+0.1417311978e18*f18-...
    0.1115892688e19*f12+0.2096525657e20*f6-...
    0.3634776482e22)/zeta11;
elseif k==4
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  zeta9=zeta8*zeta;
  zeta10=zeta9*zeta;
  zeta11=zeta10*zeta;
  zeta12=zeta11*zeta;
  zeta13=zeta12*zeta;
  zeta14=zeta13*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f28=f26*f2;
  f30=f26*f4;
  f32=f30*f2;
  f34=f32*f2;
  f36=f34*f2;
  f38=f36*f2;
  f40=f38*f2;
  f42=f40*f2;
  f44=f42*f2;
  f46=f44*f2;
  f48=f46*f2;
  f50=f48*f2;
  f52=f50*f2;
  f54=f52*f2;
  B = 0.2826137026e-28*(0.8522372387e26*f44*zeta - ...
    0.6768187756e27*f40*zeta2 - 0.4545827527e26*f38*zeta + ...
    0.2859874380e28*f36*zeta3 + 0.3064462487e27*f34*zeta2 - ...
    0.6872717295e28*f32*zeta4 + 0.5297971700e26*f32*zeta - ...
    0.1051555123e28*f30*zeta3 + 0.9266919639e28*f28*zeta5 - ...
    0.2935396959e27*f28*zeta2 + 0.1914711962e28*f26*zeta4 - ...
    0.6397577215e28*f24*zeta6 - 0.1066918725e27*f26*zeta + ...
    0.7749292103e27*f24*zeta3 - 0.1733536500e28*f22*zeta5 + ...
    0.1777821193e28*f20*zeta7 + 0.4629819426e27*f22*zeta2 - ...
    0.9667645788e27*f20*zeta4 + 0.6241289039e27*f18*zeta6 - ...
    0.8745273122e26*f16*zeta8 + 0.3414139919e27*f20*zeta - ...
    0.8579095402e27*f18*zeta3 + 0.4695158636e27*f16*zeta5 - ...
    0.3990435735e26*f14*zeta7 - 0.1071019560e28*f16*zeta2 + ...
    0.5944995103e27*f14*zeta4 - 0.4058648312e26*f12*zeta6 - ...
    0.1695350944e28*f14*zeta + 0.1152325742e28*f12*zeta3 - ...
    0.7329046745e26*f10*zeta5 + 0.3273102296e28*f10*zeta2 - ...
    0.2186649749e27*f8*zeta4 + 0.1329979254e29*f8*zeta - ...
    0.1075970511e28*f6*zeta3 - 0.9327127234e28*f4*zeta2 - ...
    0.1745381865e30*f2*zeta - 0.1383287568e28*f50*zeta + ...
    0.1264593133e29*f46*zeta2 - 0.6349542058e29*f42*zeta3 + ...
    0.1898862625e30*f38*zeta4 - 0.3424058703e30*f34*zeta5 + ...
    0.3581216115e30*f30*zeta6 - 0.1953149904e30*f26*zeta7 + ...
    0.4319041650e29*f22*zeta8 - 0.1684922621e28*f18*zeta9 - ...
    0.4438735618e25*f48 + 0.2705849719e25*f42 - 0.3679147014e25*f36 + ...
    0.8890989373e25*f30 - 0.3556395749e26*f24 + 0.2354654089e27*f18 - ...
    0.2770790112e28*f12 + 0.7272424437e29*f6 + 0.6404109111e26*f54 - ...
    0.1678798779e32) / zeta14;
else
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  zeta9=zeta8*zeta;
  zeta10=zeta9*zeta;
  zeta11=zeta10*zeta;
  zeta12=zeta11*zeta;
  zeta13=zeta12*zeta;
  zeta15=zeta13*zeta2;
  zeta17=zeta15*zeta2;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f28=f26*f2;
  f30=f26*f4;
  f32=f30*f2;
  f34=f32*f2;
  f36=f34*f2;
  f38=f36*f2;
  f40=f38*f2;
  f42=f40*f2;
  f44=f42*f2;
  f46=f44*f2;
  f48=f44*f4;
  f50=f48*f2;
  f52=f50*f2;
  f54=f52*f2;
  f56=f54*f2;
  f58=f56*f2;
  f60=f58*f2;
  f62=f60*f2;
  f64=f62*f2;
  f66=f64*f2; 
  B = 0.4212300732e-35*( ...
    0.7038556903e33*f44*zeta - ...
    0.5589790317e34*f40*zeta2 - ...
    0.1065360988e34*f38*zeta + ...
    0.2361946609e35*f36*zeta3 + ...
    0.7181880005e34*f34*zeta2 - ...
    0.5676120399e35*f32*zeta4 + ...
    0.2443939511e34*f32*zeta - ...
    0.2464426549e35*f30*zeta3 + ...
    0.7653472323e35*f28*zeta5 - ...
    0.1354090398e35*f28*zeta2 + ...
    0.4487322528e35*f26*zeta4 - ...
    0.5283706135e35*f24*zeta6 - ...
    0.8146465036e34*f26*zeta + ...
    0.3574726749e35*f24*zeta3 - ...
    0.4062719377e35*f22*zeta5 + ...
    0.1468287827e35*f20*zeta7 + ...
    0.3535101709e35*f22*zeta2 - ...
    0.4459657933e35*f20*zeta4 + ...
    0.1462709664e35*f18*zeta6 - ...
    0.7222648777e33*f16*zeta8 + ...
    0.3896177327e35*f20*zeta - ...
    0.6550574014e35*f18*zeta3 + ...
    0.2165863532e35*f16*zeta5 - ...
    0.9351992638e33*f14*zeta7 - ...
    0.1222235240e36*f16*zeta2 + ...
    0.4539304974e35*f14*zeta4 - ...
    0.1872243102e34*f12*zeta6 - ...
    0.2702805851e36*f14*zeta + ...
    0.1315020924e36*f12*zeta3 - ...
    0.5596098529e34*f10*zeta5 + ...
    0.5218129065e36*f10*zeta2 - ...
    0.2495379620e35*f8*zeta4 + ...
    0.2823202330e37*f8*zeta - ...
    0.1715361297e36*f6*zeta3 - ...
    0.1979908128e37*f4*zeta2 - ...
    0.4759112500e38*f2*zeta - ...
    0.7754156791e33*f50*zeta + ...
    0.7088803268e34*f46*zeta2 - ...
    0.3559299296e35*f42*zeta3 + ...
    0.1064426433e36*f38*zeta4 - ...
    0.1919390346e36*f34*zeta5 + ...
    0.2007486505e36*f30*zeta6 - ...
    0.1094857710e36*f26*zeta7 + ...
    0.2421081987e35*f22*zeta8 - ...
    0.9445002247e33*f18*zeta9 + ...
    0.1815457344e34*f56*zeta - ...
    0.1877568893e35*f52*zeta2 + ...
    0.1092289809e36*f48*zeta3 - ...
    0.3911600550e36*f44*zeta4 + ...
    0.8859680107e36*f40*zeta5 - ...
    0.1253586099e37*f36*zeta6 + ...
    0.1051116968e37*f32*zeta7 - ...
    0.4651606132e36*f28*zeta8 + ...
    0.8383578157e35*f24*zeta9 - ...
    0.2656866866e34*f20*zeta10 - ...
    0.3599144184e35*f62*zeta + ...
    0.4154226202e36*f58*zeta2 - ...
    0.2747951690e37*f54*zeta3 + ...
    0.1146925832e38*f50*zeta4 - ...
    0.3132396147e38*f46*zeta5 + ...
    0.5613527869e38*f42*zeta6 - ...
    0.6432672076e38*f38*zeta7 + ...
    0.4432786075e38*f34*zeta8 - ...
    0.1625581244e38*f30*zeta9 + ...
    0.2434434335e37*f26*zeta10 - ...
    0.6390972479e35*f22*zeta11 - ...
    0.3665915053e32*f48 + ...
    0.6341434452e32*f42 - ...
    0.1697180216e33*f36 + ...
    0.6788720863e33*f30 - ...
    0.4058518049e34*f24 + ...
    0.3753897015e35*f18 - ...
    0.5881671522e36*f12 + ...
    0.1982963542e38*f6 + ...
    0.3589887403e32*f54 - ...
    0.7564405600e32*f60 + ...
    0.1363312191e34*f66 - ...
    0.5718145776e40)/zeta17; 
end

function B=Bkzeta(zeta, f, k) 
two13=2^(1/3);
%ki=0, first and second component. 
%Changed to ki=1
Bjk(1,1)= 0.142857142857142857142857e-1;
Bjk(2,1)= 0.888888888888888888888889e-2;
Bjk(3,1)= 0.204823747680890538033395e-2;
Bjk(4,1)= -0.578266178266178266178266e-3;
Bjk(5,1)= -0.604120897998449018857182e-3;
Bjk(6,1)= -0.147268574562692209751033e-3;
Bjk(7,1)= 0.532410214800978351784867e-4;
Bjk(8,1)= 0.520656100658341554700098e-4;
Bjk(9,1)= 0.123311505089493920144250e-4;
Bjk(10,1)= -0.490593272853136642707052e-5;
Bjk(11,1)= -0.463223098713634979781493e-5;
Bjk(12,1)= -0.107717452345523495409932e-5;
Bjk(13,1)= 0.447596397893282157770878e-6;
Bjk(14,1)= 0.415258618846462448367911e-6;
Bjk(15,1)= 0.955581929358923426517202e-7;
Bjk(16,1)= -0.406059920840305938391248e-7;
Bjk(17,1)= -0.373136718798848249861612e-7;
Bjk(18,1)= -0.853267064555377755230178e-8;
Bjk(19,1)= 0.367301724557362439212528e-8;
Bjk(20,1)= 0.335596046078453629690337e-8;
Bjk(21,1)= 0.764310709511047504468419e-9;
Bjk(22,1)= -0.331727355222099784736441e-9;
Bjk(22,1)= -0.301966912392930035742926e-9;
%------------------------
Bjk(1,2)= -0.118485958485958485958486e-2;
Bjk(2,2)= -0.139406307977736549165121e-2;
Bjk(3,2)= -0.481410055863837376442418e-3;
Bjk(4,2)= 0.268417053660161429584988e-3;
Bjk(5,2)= 0.341970698270990239629547e-3;
Bjk(6,2)= 0.103454823490207815639056e-3;
Bjk(7,2)= -0.541819198209550495372548e-4;
Bjk(8,2)= -0.620218483069016491051219e-4;
Bjk(9,2)= -0.172488588605608661593333e-4;
Bjk(10,2)= 0.874467599288705317859166e-5;
Bjk(11,2)= 0.942068421618092915552474e-5;
Bjk(12,2)= 0.249492211208585079667788e-5;
Bjk(13,2)= -0.123845860883635664944695e-5;
Bjk(14,2)= -0.128546171380976866640804e-5;
Bjk(15,2)= -0.329971086253750597036374e-6;
Bjk(16,2)= 0.161344110578831478602255e-6;
Bjk(17,2)= 0.163362319440237459346168e-6;
Bjk(18,2)= 0.410425294960577916294255e-7;
Bjk(19,2)= -0.198431704232698899355524e-7;
Bjk(20,2)= -0.197394814276970718379204e-7;
%------------------------
Bjk(1,3)= 0.438291809448988109926169e-3;
Bjk(2,3)= 0.711048651167086689443780e-3;
Bjk(3,3)= 0.318583839453875805781720e-3;
Bjk(4,3)= -0.240480942680445795531986e-3;
Bjk(5,3)= -0.372296603862153505390841e-3;
Bjk(6,3)= -0.135275205959561728461416e-3;
Bjk(7,3)= 0.869169437270413932413888e-4;
Bjk(8,3)= 0.115875075359137640445758e-3;
Bjk(9,3)= 0.372496592784644574529580e-4;
Bjk(10,3)= -0.219833494960693482739066e-4;
Bjk(11,3)= -0.268644963387045129596591e-4;
Bjk(12,3)= -0.802306161203252181331827e-5;
Bjk(13,3)= 0.449475659218012475517714e-5;
Bjk(14,3)= 0.519350476385601375304097e-5;
Bjk(15,3)= 0.147715619152961744841117e-5;
Bjk(16,3)= -0.798879382609678364636004e-6;
Bjk(17,3)= -0.887808168959802029947558e-6;
%---------------------
Bjk(1,4)= -0.376704394771054542175466e-3;
Bjk(2,4)= -0.758562716587986423627298e-3;
Bjk(3,4)= -0.410325396877504542517502e-3;
Bjk(4,4)= 0.379126331042900768687048e-3;
Bjk(5,4)= 0.685098167390344516518718e-3;
Bjk(6,4)= 0.287831057193221563639176e-3;
Bjk(7,4)= -0.215701063611570529397680e-3;
Bjk(8,4)= -0.326086399137349916645268e-3;
Bjk(9,4)= -0.118131700874867756897152e-3;
Bjk(10,4)= 0.788752684158257871094102e-4;
Bjk(11,4)= 0.107208183342068435757316e-3;
Bjk(12,4)= 0.354459525128873402372981e-4;
Bjk(13,4)= -0.220144792073382420302000e-4;
Bjk(14,4)= -0.278933635962081525723724e-4;
%-----------------------
Bjk(1,5)= 0.584533301220761872533486e-3;
Bjk(2,5)= 0.138546904223724012512638e-2;
Bjk(3,5)= 0.868303741849469002543831e-3;
Bjk(4,5)= -0.935029048013459516962523e-3;
Bjk(5,5)= -0.191754860055254920576692e-2;
Bjk(6,5)= -0.907950471133081379467134e-3;
Bjk(7,5)= 0.770504298063922351074396e-3;
Bjk(8,5)= 0.129531001282559834871519e-2;
Bjk(9,5)= 0.519338694718995507588524e-3;
Bjk(10,5)= -0.384826319487598346504950e-3;
Bjk(11,5)= -0.573353930990124765012774e-3;
if abs(zeta) < 0.3
  eta= zeta/two13; 
  jmax= 18 - 2*k; 
  delta= eps*10^(2*k);
  for j=0:jmax 
    ck(j+1)= Bjk(j+1,k+1); 
  end
  B= two13*sumexpansion(ck, eta, jmax, delta);
else
  B=Bkexact(k,zeta,f);
end

function C=Ckexact(k,zeta,f) 
if k==0  
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  C=-0.52083333333333333333e-2*(7.*f6-36.*f2*zeta-28.)/zeta;
elseif k==1
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  C=-0.47095027970679012346e-8*(475475.*f18-3534300.*f14*zeta+...
     191100.*f12+7227792.*f10*zeta2-997920.*f8*zeta-...
     2721600.*f6*zeta3-764400.*f6+907200.*f4*zeta2+...
     3931200.*f2*zeta-30430400.)/zeta4;
elseif k==2
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f30=f26*f4;
  C=-0.45626228425214506938e-14*(201713136625.*f30-...
     2445544040940.*f26*zeta+39607067500.*f24+...
     10762000048800.*f22*zeta2-386631685440.*f20*zeta-...
     20350101163392.*f18*zeta3-48460412000.*f18+...
     1243670972544.*f16*zeta2+14558606357760.*f14*zeta4+...
     360215856000.*f14*zeta-1395328177152.*f12*zeta3-...
     1901157350400.*f10*zeta5+193841648000.*f12-...
     736656560640.*f10*zeta2+288054144000.*f8*zeta4-...
     1012236825600.*f8*zeta+277385472000.*f6*zeta3-...
     2534852320000.*f6+920215296000.*f4*zeta2+...
     13036383360000.*f2*zeta-206554251904000.)/zeta7;
elseif k==3
   zeta2=zeta*zeta;
   zeta3=zeta2*zeta;
   zeta4=zeta3*zeta;
   zeta5=zeta4*zeta;
   zeta6=zeta5*zeta;
   zeta7=zeta6*zeta;
   zeta8=zeta7*zeta;
   zeta9=zeta8*zeta;
   zeta10=zeta9*zeta;
   f2=f*f;
   f4=f2*f2;
   f6=f4*f2;
   f8=f6*f2;
   f10=f8*f2;
   f12=f10*f2;
   f14=f12*f2;
   f16=f14*f2;
   f18=f16*f2;
   f20=f18*f2;
   f22=f20*f2;
   f24=f22*f2;
   f26=f24*f2;
   f28=f26*f2;
   f30=f26*f4;
   f32=f30*f2;
   f34=f32*f2;
   f36=f34*f2;
   f38=f36*f2;
   f40=f38*f2;
   f42=f40*f2;
   C=-0.41256355274535686974e-20*...
       (0.15957555285095424000e21*f2*zeta+...
        1279846936368000000.*f14*zeta-...
        2617340759953920000.*f10*zeta2-...
        7361593537858560000.*f8*zeta+...
        985550582016000000.*f6*zeta3+...
        6692357761689600000.*f4*zeta2+...
        267053409270648000.*f26*zeta-...
        1175210405328960000.*f22*zeta2-...
        420191315736192000.*f20*zeta+...
        2222231047042406400.*f18*zeta3+...
        1351621612960819200.*f16*zeta2-...
        1589799814267392000.*f14*zeta4-...
        1516442662928793600.*f12*zeta3+...
        207606382663680000.*f10*zeta5+...
        313057243699200000.*f8*zeta4-...
        3928528740543907500.*f38*zeta+...
        26647870642102710000.*f34*zeta2-...
        439382746022220000.*f32*zeta-...
        92161564170322046400.*f30*zeta3+...
        2456304434720136000.*f28*zeta2+...
        0.16954356215302215091e21*f26*zeta4-...
        6562288550353536000.*f24*zeta3-...
        0.15566666781413192602e21*f22*zeta5+...
        8324679887242048512.*f20*zeta4+...
        57191708794360320000.*f18*zeta6-...
        4145615804884377600.*f16*zeta5-...
        3775046672517120000.*f14*zeta7+...
        373691488794624000.*f12*zeta6-...
        31028579721018880000.*f6-...
        172179843836000000.*f18+...
        1409732769244800000.*f12-...
        22027074519450000.*f30+43044960959000000.*f24+...
        232671060268521875.*f42+...
        30301347383807500.*f36-...
        0.38120826514394624000e22)/zeta10;
elseif k==4  
   zeta2=zeta*zeta;
   zeta3=zeta2*zeta;
   zeta4=zeta3*zeta;
   zeta5=zeta4*zeta;
   zeta6=zeta5*zeta;
   zeta7=zeta6*zeta;
   zeta8=zeta7*zeta;
   zeta9=zeta8*zeta;
   zeta10=zeta9*zeta;
   zeta11=zeta10*zeta;
   zeta12=zeta11*zeta;
   zeta13=zeta12*zeta;
   f2=f*f;
   f4=f2*f2;
   f6=f4*f2;
   f8=f6*f2;
   f10=f8*f2;
   f12=f10*f2;
   f14=f12*f2;
   f16=f14*f2;
   f18=f16*f2;
   f20=f18*f2;
   f22=f20*f2;
   f24=f22*f2;
   f26=f24*f2;
   f28=f26*f2;
   f30=f26*f4;
   f32=f30*f2;
   f34=f32*f2;
   f36=f34*f2;
   f38=f36*f2;
   f40=f38*f2;
   f42=f40*f2;
   f44=f42*f2;
   f46=f44*f2;
   f48=f44*f4;
   f50=f48*f2;
   f52=f50*f2;
   f54=f52*f2;   
   C=-0.28261370263600418638e-28*...
       (0.54589603005176120366e30*f2*zeta+...
        0.20679555469508010916e28*f14*zeta-...
        0.42290559824029155769e28*f10*zeta2-...
        0.17933866592043922072e29*f8*zeta+...
        0.15924363570102425518e28*f6*zeta3+...
        0.16303515083676292792e29*f4*zeta2+...
        0.12524698073429682941e27*f26*zeta-...
        0.55116897925766092416e27*f22*zeta2-...
        0.40337660389263995197e27*f20*zeta+...
        0.10422174721387004320e28*f18*zeta3+...
        0.12975340411992886902e28*f16*zeta2-...
        0.74560975369214977843e27*f14*zeta4-...
        0.14557594801749046523e28*f12*zeta3+...
        0.97366563043735265280e26*f10*zeta5+...
        0.30052969458953785344e27*f8*zeta4+...
        0.56627384677696100268e26*f38*zeta-...
        0.38411306658352530302e27*f34*zeta2-...
        0.63032794225757227872e26*f32*zeta+...
        0.13284536505766901056e28*f30*zeta3+...
        0.35237553907430738223e27*f28*zeta2-...
        0.24438687222985224921e28*f26*zeta4-...
        0.94141016594119742607e27*f24*zeta3+...
        0.22438416165400232344e28*f22*zeta5+...
        0.11942385973920148986e28*f20*zeta4-...
        0.82438416724542739661e27*f18*zeta6-...
        0.59472009389078108799e27*f16*zeta5+...
        0.54415032756330774528e26*f14*zeta7+...
        0.53608884122903651942e26*f12*zeta6-...
        0.14397482850741596886e28*f50*zeta+...
        0.13207972724871225663e29*f46*zeta2-...
        0.12486266520014393707e27*f44*zeta-...
        0.66592758171771445484e29*f42*zeta3+...
        0.99613840305204192856e27*f40*zeta2+...
        0.20015038483748922317e30*f38*zeta4-...
        0.42326140823881103376e28*f36*zeta3-...
        0.36315774118707124699e30*f34*zeta5+...
        0.10242565774476707856e29*f32*zeta4+...
        0.38281965363754011064e30*f30*zeta6-...
        0.13934701383195476995e29*f28*zeta5-...
        0.21094018958210280550e30*f26*zeta7+...
        0.97354435881026066072e28*f24*zeta6+...
        0.47303789496173550313e29*f22*zeta8-...
        0.27509443725437275934e28*f20*zeta7-...
        0.18831488121744471613e28*f18*zeta9+...
        0.13875833352864347505e27*f16*zeta8-...
        0.10614645028784245627e30*f6-...
        0.27820534863662737997e27*f18+...
        0.34343052606818116762e28*f12-...
        0.10330609841323972200e26*f30+...
        0.41322439365295888800e26*f24-...
        0.33538137311345817150e25*f42+...
        0.43469585724473028120e25*f36+...
        0.66457736059297289034e26*f54+...
        0.64786651787013217937e25*f48-...
        0.17421496761528428537e32)/zeta13;    
end  

function D=Dkexact(k,zeta,f) 
if k==1
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  D=-0.13563368055555555556e-4*(455.*f12-...
     2376.*f8*zeta-280.*f6+2160.*f4*zeta2+...
     1440.*f2*zeta-6160.)/zeta3;
elseif k==2
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  D=-0.61321651003488297325e-11*(202076875.*f24-...
     1972610640.*f20*zeta-38038000.*f18+...
     6345260064.*f16*zeta2+282744000.*f14*zeta-...
     7119021312.*f12*zeta3+84084000.*f12-...
     578223360.*f10*zeta2+1469664000.*f8*zeta4-...
     439084800.*f8*zeta+217728000.*f6*zeta3-...
     762361600.*f6+399168000.*f4*zeta2+...
     3920716800.*f2*zeta-47593145600.)/zeta6;
elseif k==3
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  zeta9=zeta8*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f28=f26*f2;
  f30=f26*f4;
  f32=f30*f2;
  f34=f32*f2;
  f36=f34*f2;
  D= -0.79212202127108518990e-18*( ...
    208686424826880000*f2*zeta + ...
    2694414602880000*f14*zeta - ...
    5510191073587200*f10*zeta2 - ...
    11873537964288000*f8*zeta + ...
    2074843330560000*f6*zeta3 + ...
    10794125422080000*f4*zeta2 + ...
    1467326424564000*f26*zeta - ...
    6457200029280000*f22*zeta2 - ...
    1275884561952000*f20*zeta + ...
    12210060698035200*f18*zeta3 + ...
    4104114209395200*f16*zeta2 - ...
    8735163814656000*f14*zeta4 - ...
    4604582984601600*f12*zeta3 + ...
    1140694410240000*f10*zeta5 + ...
    950578675200000*f8*zeta4 - ...
    15692240929365000*f32*zeta + ...
    87725158382862000*f28*zeta2 - ...
    234367448226912000*f24*zeta3 + ...
    297309995972930304*f20*zeta4 - ...
    148057707317299200*f16*zeta5 + ...
    13346124599808000*f12*zeta6 - ...
    40577915938560000*f6 - ...
    362483881760000*f18 + ...
    2273762531040000*f12 - ...
    121027881975000*f30 + ...
    130703322750000*f24 + ...
    1082190977993125*f36 - ...
    4193051313651200000 ) / zeta9;
elseif k==4
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  zeta9=zeta8*zeta;
  zeta10=zeta9*zeta;
  zeta11=zeta10*zeta;
  zeta12=zeta11*zeta;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f28=f26*f2;
  f30=f26*f4;
  f32=f30*f2;
  f34=f32*f2;
  f36=f34*f2;
  f38=f36*f2;
  f40=f38*f2;
  f42=f40*f2;
  f44=f42*f2;
  f48=f44*f4;
  D= -0.26859606298525837874e-23*( ...
    0.10468156267022598144e25*f2*zeta + ...
    0.54634106019677184000e22*f14*zeta - ...
    0.11172904236091293696e23*f10*zeta2 - ...
    0.39850759684941004800e23*f8*zeta + ...
    0.42071183245099008000e22*f6*zeta3 + ...
    0.36227963349946368000e23*f4*zeta2 + ...
    0.53268253369185254400e21*f26*zeta - ...
    0.23441530218294988800e22*f22*zeta2 - ...
    0.13143584356228085760e22*f20*zeta + ...
    0.44326101951672532992e22*f18*zeta3 + ...
    0.42278724053414424576e22*f16*zeta2 - ...
    0.31711206961920245760e22*f14*zeta4 - ...
    0.47434326496412663808e22*f12*zeta3 + ...
    0.41410553128648704000e21*f10*zeta5 + ...
    0.97924305829109760000e21*f8*zeta4 + ...
    0.62856459848702520000e21*f38*zeta - ...
    0.42636593027364336000e22*f34*zeta2 - ...
    0.38665681649955360000e21*f32*zeta + ...
    0.14745850267251527424e23*f30*zeta3 + ...
    0.21615479025537196800e22*f28*zeta2 - ...
    0.27126969944483544146e23*f26*zeta4 - ...
    0.57748139243111116800e22*f24*zeta3 + ...
    0.24906666850261108163e23*f22*zeta5 + ...
    0.73257183007730026906e22*f20*zeta4 - ...
    0.91506734070976512000e22*f18*zeta6 - ...
    0.36481419082982522880e22*f16*zeta5 + ...
    0.60400746760273920000e21*f14*zeta7 + ...
    0.32884851013926912000e21*f12*zeta6 - ...
    0.90088503030406881000e22*f44*zeta + ...
    0.71871457651662476808e23*f40*zeta2 - ...
    0.30538341142771358857e24*f36*zeta3 + ...
    0.73900185963035410218e24*f32*zeta4 - ...
    0.10053897101872638524e25*f28*zeta5 + ...
    0.70241295729456036127e24*f24*zeta6 - ...
    0.19848083495986490573e24*f20*zeta7 + ...
    0.10011423775515402240e23*f16*zeta8 - ...
    0.20354748296988385280e24*f6 - ...
    0.73500131736711680000e21*f18 + ...
    0.76313533908451840000e22*f12 - ...
    43936671308129600000*f30 + ...
    0.13464463787975200000e21*f24 - ...
    37227369642963500000*f42 + ...
    26665185697750600000*f36 + ...
    0.46743616007946044688e21*f48 - ...
    0.29383533077295376179e26 ) / zeta12;
elseif k==5
  zeta2=zeta*zeta;
  zeta3=zeta2*zeta;
  zeta4=zeta3*zeta;
  zeta5=zeta4*zeta;
  zeta6=zeta5*zeta;
  zeta7=zeta6*zeta;
  zeta8=zeta7*zeta;
  zeta9=zeta8*zeta;
  zeta10=zeta9*zeta;
  zeta11=zeta10*zeta;
  zeta12=zeta11*zeta;
  zeta13=zeta12*zeta;
  zeta15=zeta13*zeta2;
  f2=f*f;
  f4=f2*f2;
  f6=f4*f2;
  f8=f6*f2;
  f10=f8*f2;
  f12=f10*f2;
  f14=f12*f2;
  f16=f14*f2;
  f18=f16*f2;
  f20=f18*f2;
  f22=f20*f2;
  f24=f22*f2;
  f26=f24*f2;
  f28=f26*f2;
  f30=f26*f4;
  f32=f30*f2;
  f34=f32*f2;
  f36=f34*f2;
  f38=f36*f2;
  f40=f38*f2;
  f42=f40*f2;
  f44=f42*f2;
  f46=f44*f2;
  f48=f44*f4;
  f50=f48*f2;
  f52=f50*f2;
  f54=f52*f2;
  f56=f54*f2;
  f58=f56*f2;
  f60=f58*f2; 
  D= -0.10513902627827536696e-31*( ...
    0.84611458453444978740e34*f2*zeta + ...
    0.23740129678995196531e32*f14*zeta - ...
    0.48549562677985470823e32*f10*zeta2 - ...
    0.24190992646008046482e33*f8*zeta + ...
    0.18281169378477584495e32*f6*zeta3 + ...
    0.21991811496370951348e33*f4*zeta2 + ...
    0.93564504487749103441e30*f26*zeta - ...
    0.41174527426464301678e31*f22*zeta2 - ...
    0.38213210275429424784e31*f20*zeta + ...
    0.77857814038649477070e31*f18*zeta3 + ...
    0.12291972483627928192e32*f16*zeta2 - ...
    0.55700031039818357048e31*f14*zeta4 - ...
    0.13790894808856930072e32*f12*zeta3 + ...
    0.72736717256191992575e30*f10*zeta5 + ...
    0.28470179734115552649e31*f8*zeta4 + ...
    0.19766732411494452067e30*f38*zeta - ...
    0.13408106777542256578e31*f34*zeta2 - ...
    0.34504151559179506537e30*f32*zeta + ...
    0.46371888762796995954e31*f30*zeta3 + ...
    0.19289037008927586103e31*f28*zeta2 - ...
    0.85307310866367091791e31*f26*zeta4 - ...
    0.51532792483621147103e31*f24*zeta3 + ...
    0.78325031361357077701e31*f22*zeta5 + ...
    0.65372620821238895550e31*f20*zeta4 - ...
    0.28776503331313718991e31*f18*zeta6 - ...
    0.32554977939581356757e31*f16*zeta5 + ...
    0.18994474100809862362e30*f14*zeta7 + ...
    0.29345503168877459073e30*f12*zeta6 + ...
    0.40312951982076471281e30*f50*zeta - ...
    0.36982323629639431856e31*f46*zeta2 - ...
    0.19228850440822166308e30*f44*zeta + ...
    0.18645972288096004735e32*f42*zeta3 + ...
    0.15340531407001445700e31*f40*zeta2 - ...
    0.56042107754496982486e32*f38*zeta4 - ...
    0.65182256868776899199e31*f36*zeta3 + ...
    0.10168416753237994916e33*f34*zeta5 + ...
    0.15773551292694130099e32*f32*zeta4 - ...
    0.10718950301851123098e33*f30*zeta6 - ...
    0.21459440130121034572e32*f28*zeta5 + ...
    0.59063253082988785541e32*f26*zeta7 + ...
    0.14992583125678014175e32*f24*zeta6 - ...
    0.13245061058928594088e32*f22*zeta8 - ...
    0.42364543337173404938e31*f20*zeta7 + ...
    0.52728166740884520518e30*f18*zeta9 + ...
    0.21368783363411095157e30*f16*zeta8 - ...
    0.72364383641290146765e31*f56*zeta + ...
    0.75046117139810240390e32*f52*zeta2 - ...
    0.43798854855079282757e33*f48*zeta3 + ...
    0.15744367152094499372e34*f44*zeta4 - ...
    0.35823164139419832070e34*f40*zeta5 + ...
    0.50969984247250518544e34*f36*zeta6 - ...
    0.43035806376405631899e34*f32*zeta7 + ...
    0.19216036727738034843e34*f28*zeta8 - ...
    0.35048403665498508495e33*f24*zeta9 + ...
    0.11294373315897464295e32*f20*zeta10 - ...
    0.16452228032614301422e34*f6 - ...
    0.31937974023484823220e31*f18 + ...
    0.46325343661336957700e32*f12 - ...
    0.77173787758626601923e29*f30 + ...
    0.39146124225390305323e30*f24 - ...
    0.11707045797480446573e29*f42 + ...
    0.23795251225576535593e29*f36 - ...
    0.18608166096603240930e29*f54 + ...
    0.99771443752000355623e28*f48 + ...
    0.30080100495159138963e30*f60 - ...
    0.30507131409047661779e36 ) / zeta15;
end

function D=Dkzeta(zeta, f, k)
two13=2^(1/3);
%ki=0, first component. Changed to ki=1
Djk(1,1)= 0.730158730158730158730159e-2;
Djk(2,1)= 0.419336219336219336219336e-2;
Djk(3,1)= -0.450495536209821924107638e-3;
Djk(4,1)= -0.152302582778773254963731e-2;
Djk(5,1)= -0.602260874305692232863101e-3;
Djk(6,1)= 0.134516202087850843046701e-3;
Djk(7,1)= 0.227957335910397990094302e-3;
Djk(8,1)= 0.759662380303310963431066e-4;
Djk(9,1)= -0.208929293145702503498759e-4;
Djk(10,1)= -0.288989984337042437996113e-4;
Djk(11,1)= -0.882085288915952350246854e-5;
Djk(12,1)= 0.268699754140326485858770e-5;
Djk(13,1)= 0.336422008038391207907941e-5;
Djk(14,1)= 0.973596713460319029167637e-6;
Djk(15,1)= -0.315173324724484568747324e-6;
Djk(16,1)= -0.371902658726380167612929e-6;
Djk(17,1)= -0.103842224770769608750491e-6;
Djk(18,1)= 0.350150687274032091141572e-7;
Djk(19,1)= 0.397171235257758884576377e-7;
Djk(20,1)= 0.108068143474318368422919e-7;
Djk(21,1)= -0.375354806995033684725386e-8;
Djk(22,1)= -0.413840555021600212423505e-8;
%---------------------------------
Djk(1,2)= -0.937299455394693489931585e-3;
Djk(2,2)= -0.790690489715699799733413e-3;
Djk(3,2)= 0.293543703725772342692795e-3;
Djk(4,2)= 0.759272400940109832633704e-3;
Djk(5,2)= 0.350160108830627551277944e-3;
Djk(6,2)= -0.139795674313775324610736e-3;
Djk(7,2)= -0.238208912033003463924301e-3;
Djk(8,2)= -0.885682782902038147880529e-4;
Djk(9,2)= 0.358182517107331492611466e-4;
Djk(10,2)= 0.514900346471815453486348e-4;
Djk(11,2)= 0.171249859424247232793150e-4;
Djk(12,2)= -0.697628383205670347703814e-5;
Djk(13,2)= -0.911796118302585196984065e-5;
Djk(14,2)= -0.283069240173736478532390e-5;
Djk(15,2)= 0.115921686495136845230928e-5;
Djk(16,2)= 0.142494761926653767630103e-5;
Djk(17,2)= 0.422052804623112329027247e-6;
Djk(18,2)= -0.173575341276739650538416e-6;
Djk(19,2)= -0.204376944779354029041699e-6;
%----------------------
Djk(1,3)= 0.444495041599140470657888e-3;
Djk(2,3)= 0.472063930501695695068578e-3;
Djk(3,3)= -0.288641537836417409885435e-3;
Djk(4,3)= -0.801448906698328116457619e-3;
Djk(5,3)= -0.427270049238916578621937e-3;
Djk(6,3)= 0.231003885183789413048089e-3;
Djk(7,3)= 0.431796240309483025358877e-3;
Djk(8,3)= 0.180109473648234133944025e-3;
Djk(9,3)= -0.910271534839798903976819e-4;
Djk(10,3)= -0.142460091487209173987685e-3;
Djk(11,3)= -0.521208962926222839685499e-4;
Djk(12,3)= 0.253553321898728396521221e-4;
Djk(13,3)= 0.357360245119576212047157e-4;
Djk(14,3)= 0.120336367723061296449446e-4;
Djk(15,3)= -0.571467368238561493453327e-5;
Djk(16,3)= -0.751083430778512389626646e-5;
%------------------
Djk(1,4)= -0.455709396340424418035779e-3;
Djk(2,4)= -0.571724624632071801184617e-3;
Djk(3,4)= 0.476698510527782146839932e-3;
Djk(4,4)= 0.145362409347314380791410e-2;
Djk(5,4)= 0.875271889787213278878581e-3;
Djk(6,4)= -0.584212244697429717425129e-3;
Djk(7,4)= -0.119863142824806225281440e-2;
Djk(8,4)= -0.553361926292205772487962e-3;
Djk(9,4)= 0.329748357174476289919883e-3;
Djk(10,4)= 0.560846441567589440030798e-3;
Djk(11,4)= 0.223805653823594922878522e-3;
Djk(12,4)= -0.124764472641943425453482e-3;
Djk(13,4)= -0.189331882520153576169346e-3;
%------------------------------------------
Djk(1,5)= 0.811234305193098332019330e-3;
Djk(2,5)= 0.116042840661402233286282e-2;
Djk(3,5)= -0.120961038552699514663161e-2;
Djk(4,5)= -0.403337872044996403015012e-2;
Djk(5,5)= -0.269322476027294353847792e-2;
Djk(6,5)= 0.211422330764814095297867e-2;
Djk(7,5)= 0.472359714331690567866988e-2;
Djk(8,5)= 0.238370161988457179934623e-2;
Djk(9,5)= -0.162089313723573923967785e-2;
Djk(10,5)= -0.297606926592315334080523e-2;
if abs(zeta) < 0.3
  eta= zeta/two13; 
  jmax= 18 - 2*k; 
  delta= eps*10^(2*k);
  for j=0:jmax 
    ck(j+1)= Djk(j+1,k); 
  end 
  D= sumexpansion(ck, eta, jmax, delta);
else
  D=Dkexact(k,zeta,f);
end

function C=Ckzeta(zeta, f, k) 
two13=2^(1/3);
%ki=0, first and second component. 
%Changed to ki=1
Cjk(1,1)= 0.1;
Cjk(2,1)= 0.200000000000000000000000e-1;
Cjk(3,1)= -0.330158730158730158730159e-2;
Cjk(4,1)= -0.296911976911976911976912e-2;
Cjk(5,1)= -0.137695954838811981669124e-3;
Cjk(6,1)= 0.555048929906072763215620e-3;
Cjk(7,1)= 0.242456146115529869031270e-3;
Cjk(8,1)= -0.154566037237636020689912e-4;
Cjk(9,1)= -0.560494560611812641811244e-4;
Cjk(10,1)= -0.205938726629285257807890e-4;
Cjk(11,1)= 0.260832410695472146322360e-5;
Cjk(12,1)= 0.532729210222188857946138e-5;
Cjk(13,1)= 0.178881182990467430870839e-5;
Cjk(14,1)= -0.300617408480663364532565e-6;
Cjk(15,1)= -0.495357296368007386699085e-6;
Cjk(16,1)= -0.157154550860651784075191e-6;
Cjk(17,1)= 0.310582008908507027210460e-7;
Cjk(18,1)= 0.455703940740490910695069e-7;
Cjk(19,1)= 0.138952925679589253258334e-7;
Cjk(20,1)= -0.305767396855684609868392e-8;
Cjk(21,1)= -0.416704885308313929164001e-8;
Cjk(22,1)= -0.123351938981635545635728e-8;
Cjk(23,1)= 0.293335110754148084704709e-9;
Cjk(24,1)= 0.379626674718002001131840e-9;
%---------
Cjk(1,2)= -0.136652236652236652236652e-2;
Cjk(2,2)= -0.273304473304473304473304e-3;
Cjk(3,2)= 0.780378351806923235494664e-3;
Cjk(4,2)= 0.480481095831235887258296e-3;
Cjk(5,2)= -0.149133522849003299556727e-3;
Cjk(6,2)= -0.300753222613497123442870e-3;
Cjk(7,2)= -0.116371247406915680355006e-3;
Cjk(8,2)= 0.396104305305290465475782e-4;
Cjk(9,2)= 0.596307815250464462047999e-4;
Cjk(10,2)= 0.198810673627150872495934e-4;
Cjk(11,2)= -0.711571245181156693105293e-5;
Cjk(12,2)= -0.941398869171534444964037e-5;
Cjk(13,2)= -0.289613648251761565880795e-5;
Cjk(14,2)= 0.106874541035013210932373e-5;
Cjk(15,2)= 0.131276867676816411452532e-5;
Cjk(16,2)= 0.383761251199166055079476e-6;
Cjk(17,2)= -0.144594124257249176458698e-6;
Cjk(18,2)= -0.169190789079476125912616e-6;
Cjk(19,2)= -0.477325260434205248333610e-7;
Cjk(20,2)= 0.182632080752898752616217e-7;
Cjk(21,2)= 0.206465177093998576833282e-7;
%------------------------
Cjk(1,3)= 0.301615299318380550873548e-3;
Cjk(2,3)= 0.603230598636761101747096e-4;
Cjk(3,3)= -0.391597684545969512977759e-3;
Cjk(4,3)= -0.289790085717604851031592e-3;
Cjk(5,3)= 0.148338746317245900757649e-3;
Cjk(6,3)= 0.318187930024330786834948e-3;
Cjk(7,3)= 0.140863378363957417331528e-3;
Cjk(8,3)= -0.661957688156824128532334e-4;
Cjk(9,3)= -0.107844623154343122498834e-3;
Cjk(10,3)= -0.400103296947872344389066e-4;
Cjk(11,3)= 0.182052411474350468881140e-4;
Cjk(12,3)= 0.259353860787712690660584e-4;
Cjk(13,3)= 0.871634320891265678415106e-5;
Cjk(14,3)= -0.389911063231578420704130e-5;
Cjk(15,3)= -0.511759079066164065972108e-5;
Cjk(16,3)= -0.161273941555264126239528e-5;
Cjk(17,3)= 0.714125189263839431850448e-6;
Cjk(18,3)= 0.886393151169713542004539e-6;
%----------------------
Cjk(1,4)= -0.191582461090933800468500e-3;
Cjk(2,4)= -0.383164922181867600937000e-4;
Cjk(3,4)= 0.412025391628760368304863e-3;
Cjk(4,4)= 0.354435556191377360085550e-3;
Cjk(5,4)= -0.244344099186617866886315e-3;
Cjk(6,4)= -0.577690224079697003990733e-3;
Cjk(7,4)= -0.288229407397989651107383e-3;
Cjk(8,4)= 0.167679384035986552008787e-3;
Cjk(9,4)= 0.299177428518370852808697e-3;
Cjk(10,4)= 0.122592827219627459365658e-3;
Cjk(11,4)= -0.660336620798892415061247e-4;
Cjk(12,4)= -0.101973509296252262006631e-3;
Cjk(13,4)= -0.373065702856469195365027e-4;
Cjk(14,4)= 0.192016880205339806075821e-4;
Cjk(15,4)= 0.270693197224160526115489e-4;
%-----------------
Cjk(1,5)= 0.240291372681993458788730e-3;
Cjk(2,5)= 0.480582745363986917577461e-4;
Cjk(3,5)= -0.745008930228873560876256e-3;
Cjk(4,5)= -0.724897584623060371465682e-3;
Cjk(5,5)= 0.618191768954996746066377e-3;
Cjk(6,5)= 0.160403374732237116145661e-2;
Cjk(7,5)= 0.887021131887307152191308e-3;
Cjk(8,5)= -0.606942883053824288153770e-3;
Cjk(9,5)= -0.117879817494336191988832e-2;
Cjk(10,5)= -0.527593162723709314899751e-3;
Cjk(11,5)= 0.324720140692314086985285e-3;
Cjk(12,5)= 0.540845643878935522856430e-3;
%-----------------------
Cjk(1,6)= -0.499090772787701333229479e-3;
Cjk(2,6)= -0.998181545575402666458959e-4;
Cjk(3,6)= 0.206206753490600447118091e-2;
Cjk(4,6)= 0.222570153585730233716748e-2;
Cjk(5,6)= -0.223914974046656858808991e-2;
Cjk(6,6)= -0.632084572564640421312756e-2;
Cjk(7,6)= -0.382265214430849401331132e-2;
Cjk(8,6)= 0.298287351772147095269266e-2;
Cjk(9,6)= 0.625707189095209801465069e-2;
if abs(zeta) < 0.3
  eta= zeta/two13; 
  jmax= 18 - 2*k; 
  delta= eps*10^(2*k);
  for j=0:jmax 
    ck(j+1)= Cjk(j+1,k+1); 
  end
  C= two13^2*sumexpansion(ck, eta, jmax, delta);
else
  C=Ckexact(k,zeta,f);
end

function f=phizeta(zeta, z)
two13=2^(1/3);
%ki=0. Change to ki=1
phik(1)= 1.0;
phik(2)= 0.200000000000000000000000;
phik(3)= 0.257142857142857142857143e-1;
phik(4)= -0.565079365079365079365079e-2;
phik(5)= -0.393679653679653679653680e-2;
phik(6)= -0.520936206650492364778079e-3;
phik(7)= 0.370854126473174092221711e-3;
phik(8)= 0.212382784029362660815242e-3;
phik(9)= 0.215062978875314527278662e-4;
phik(10)= -0.263690406298179894036302e-4;
phik(11)= -0.140546982649312867820944e-4;
phik(12)= -0.114932889902944139647916e-5;
phik(13)= 0.197264119393862359423994e-5;
phik(14)= 0.101432430559332922363478e-5;
phik(15)= 0.703433110019219963923458e-7;
phik(16)= -0.152504477739267649835292e-6;
phik(17)= -0.767786625690057235680721e-7;
phik(18)= -0.466984263869301799437010e-8;
phik(19)= 0.120667364596532865221556e-7;
phik(20)= 0.599087766809234353800582e-8;
phik(21)= 0.326910215007771473164964e-9;
phik(22)= -0.971383502444283350951986e-9;
phik(23)= -0.477459342952322338344740e-9;
phik(24)= -0.237503186428391557787830e-10;
phik(25)= 0.792445981091066555670358e-10;
phik(26)= 0.386535842308178655279643e-10;
phik(27)= 0.177341058464268078733321e-11;
%if abs(zeta) < 0.15
if abs(zeta) < 0.2
  eta= zeta/two13;
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
  ar(1)=+1.1280094441222097012e+0;
  ar(2)=-1.836709914558369e-4;
  ar(3)=+1.1677286575962e-6;
  ar(4)=-2.19524008170e-8;
  ar(5)=+7.698908177e-10;
  ar(6)=-4.09915521e-11;
  ar(7)=+2.9533495e-12;
  ar(8)=-2.679714e-13;
  ar(9)=+2.91808e-14;
  ar(10)=-3.6863e-15;
  ar(11)=+5.269e-16;
  ar(12)=-8.36e-17;
  ar(13)=+1.45e-17;
  ar(14)=-2.7e-18;
  ar(15)=+5.0e-19;
  ar(16)=-1.0e-19;
  br(1)=+7.7988242533135881502e-2;
  br(2)=-1.83421439318478417e-4;
  br(3)=+2.226279540000655e-6;
  br(4)=-6.0307865841937e-8;
  br(5)=+2.708298845975e-9;
  br(6)=-1.73155804053e-10;
  br(7)=+1.4396760940e-11;
  br(8)=-1.467799419e-12;
  br(9)=+1.76234030e-13;
  br(10)=-2.4205568e-14;
  br(11)=+3.721795e-15;
  br(12)=-6.30106e-16;
  br(13)=+1.15949e-16;
  br(14)=-2.2951e-17;
  br(15)=+4.846e-18;
  br(16)=-1.084e-18;
  br(17)=+2.55e-19;
  br(18)=-6.3e-20;
  br(19)=+1.6e-20;
  br(20)=-4.0e-21;
  br(21)=+1.0e-21;
  z=(2.0/3.0)*(-x)^(1.5); la= 5.0; 
  t=-2*(la/x)^3-1.0;
  s=sin(z+piquart); c=cos(z+piquart);
  f1=chepolsumC(15,t,ar);
  f2=chepolsumC(20,t,br);
  y=(-x)^(-0.25);
  ai=y*(s*f1-c*f2/z);
  bi=y*(c*f1+s*f2/z);
elseif x <0  
  ar(1)=+1.5217399882834532859e-1;
  ar(2)=+9.463304395658582345e-2;
  ar(3)=+8.835933286664339033e-2;
  ar(4)=-1.8809032021891572589e-1;
  ar(5)=+7.873696950590187481e-2;
  ar(6)=-1.616162609419079569e-2;
  ar(7)=+2.03792657403144947e-3;
  ar(8)=-1.7619321508091265e-4;
  ar(9)=+1.116557634724685e-5;
  ar(10)=-5.4278037289400e-7;
  ar(11)=+2.092500233007e-8;
  ar(12)=-6.5627359901e-10;
  ar(13)=+1.708749758e-11;
  ar(14)=-3.7550458e-13;
  ar(15)=+7.06091e-15;
  ar(16)=-1.1494e-16;
  ar(17)=+1.64e-18;
  ar(18)=-2.1e-20;    
  br(1)=+7.5365743540308812655e-2;
  br(2)=-4.3238169422348489365e-2;
  br(3)=+9.5247283336721394897e-2;
  br(4)=-6.1705567716412224149e-2;
  br(5)=+1.7697290707409225048e-2;
  br(6)=-2.909533805902076483e-3;
  br(7)=+3.13306256975299299e-4;
  br(8)=-2.3946459531381245e-5;
  br(9)=+1.370954782252896e-6;
  br(10)=-6.1112783503340e-8;
  br(11)=+2.184145572027e-9;
  br(12)=-6.4038637539e-11;
  br(13)=+1.569118251e-12;
  br(14)=-3.2625307e-14;
  br(15)=+5.83052e-16;
  br(16)=-9.054e-18;
  br(17)=+1.23e-19;
  br(18)=-1.5e-21;
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
  ar(1)=+9.7934675132923263719e-1;
  ar(2)=-5.817935703328346930e-2;
  ar(3)=-1.8573106989330210396e-1;
  ar(4)=+2.466671790187645891e-2;
  ar(5)=+1.523182915947318079e-2;
  ar(6)=-3.18268812403824637e-3;
  ar(7)=-6.4421595859367550e-4;
  ar(8)=+2.1680741331029091e-4;
  ar(9)=+1.036058369682457e-5;
  ar(10)=-9.16356774348920e-6;
  ar(11)=+3.6678407013177e-7;
  ar(12)=+2.5198203640692e-7;
  ar(13)=-2.903896001941e-8;
  ar(14)=-4.22648003774e-9;
  ar(15)=+9.7569578821e-10;
  ar(16)=+2.168534257e-11;
  ar(17)=-2.083233561e-11;
  ar(18)=+1.00889184e-12;
  ar(19)=+2.9268193e-13;
  ar(20)=-3.543870e-14;
  ar(21)=-2.19178e-15;
  ar(22)=+6.4842e-16;
  ar(23)=-1.036e-17;
  ar(24)=-7.75e-18;
  ar(25)=+5.9e-19;
  ar(26)=+6.0e-20;
  br(1)=+7.2340133270078824597e-1;
  br(2)=-6.186433877265966711e-2;
  br(3)=+1.6859154747460945662e-1;
  br(4)=-1.516452843795166719e-2;
  br(5)=+1.268194063516318154e-2;
  br(6)=+2.37665412524311142e-3;
  br(7)=-1.35486356206889915e-3;
  br(8)=+9.4811517986800841e-4;
  br(9)=-3.0184629675683576e-4;
  br(10)=+8.26176092954571e-5;
  br(11)=-1.29237026410167e-5;
  br(12)=+3.86277876708896e-7;
  br(13)=+6.7349791262055e-7;
  br(14)=-2.6650662327003e-7;
  br(15)=+6.313398211272e-8;
  br(16)=-9.28482785276e-9;
  br(17)=+4.2243068981e-10;
  br(18)=+2.2149462027e-10;
  br(19)=-7.885556616e-11;
  br(20)=+1.489486494e-11;
  br(21)=-1.61561338e-12;
  br(22)=-1.99887e-15;
  br(23)=+4.310890e-14;
  br(24)=-1.026039e-14;
  br(25)=+1.36524e-15;
  br(26)=-7.870e-17;
  br(27)=-1.152e-17;
  br(28)=+3.99e-18;
  br(29)=-6.0e-19;
  br(30)=+5.0e-20; 
  la= 4.5; t=2*x/la-1.0;
  f1=chepolsumC(25,t,ar);
  f2=chepolsumC(29,t, br);
  ai=f1*exp(-1.5*x);
  bi=f2*exp(1.375*x);
elseif x<9
  ar(1)=+5.0491255733788020503e+1;
  ar(2)=-6.257510019590865924e+0;
  ar(3)=-5.510390933031668771e+0;
  ar(4)=+9.23785747211947637e-1;
  ar(5)=+2.69701933599952531e-1;
  ar(6)=-6.0503610861932189e-2;
  ar(7)=-7.272383379486571e-3;
  ar(8)=+2.432739101814126e-3;
  ar(9)=+8.9407383924725e-5;
  ar(10)=-6.7626485765665e-5;
  ar(11)=+1.143189712477e-6;
  ar(12)=+1.369896864000e-6;
  ar(13)=-8.37850105221e-8;
  ar(14)=-2.0445153671e-8;
  ar(15)=+2.247360374e-9;
  ar(16)=+2.14042502e-10;
  ar(17)=-4.0617972e-11;
  ar(18)=-1.169040e-12;
  ar(19)=+5.51909e-13;
  ar(20)=-8.169e-15;
  ar(21)=-5.801e-15;
  ar(22)=+3.09e-16;
  ar(23)=+4.6e-17;
  ar(24)=-5.0e-18;
  br(1)=+5.1865644423401596677e-3;
  br(2)=+3.190078113025389266e-4;
  br(3)=+6.589260266019625254e-4;
  br(4)=+1.77323651950275976e-5;
  br(5)=+4.24724564662038450e-5;
  br(6)=-3.095767914834698e-7;
  br(7)=+1.8968968954113726e-6;
  br(8)=-8.26612616574416e-8;
  br(9)=+6.82206567197979e-8;
  br(10)=-5.5232442059920e-9;
  br(11)=+2.1841788459370e-9;
  br(12)=-2.532110361946e-10;
  br(13)=+6.69024065663e-11;
  br(14)=-9.6266259764e-12;
  br(15)=+2.0421501469e-12;
  br(16)=-3.287463020e-13;
  br(17)=+6.17455626e-14;
  br(18)=-1.01838450e-14;
  br(19)=+1.7333508e-15;
  br(20)=-2.689068e-16;
  br(21)=+3.97206e-17;
  br(22)=-5.1324e-18;
  br(23)=+5.434e-19;
  br(24)=-3.01e-20;
  br(25)=-5.0e-21;
  br(26)=+2.3e-21;
  br(27)=-5.0e-22;
  br(28)=+1.0e-22;
  la=4.5; t=2*(x-la)/(9.0-la)-1.0;
  f1=chepolsumC(23,t,ar);
  f2=chepolsumC(27,t,br);
  ai=f1*exp(-2.5*x);
  bi=f2*exp(2.5*x);
else
  ar(1)=+5.6312443149444830912e-1;
  ar(2)=-5.2879891960803332e-4;
  ar(3)=+3.72709673104629e-6;
  ar(4)=-4.905124032702e-8;
  ar(5)=+9.3521089895e-10;
  ar(6)=-2.313916112e-11;
  ar(7)=+6.9870927e-13;
  ar(8)=-2.476246e-14;
  ar(9)=+1.00261e-15;
  ar(10)=-4.547e-17;
  ar(11)=+2.28e-18;
  ar(12)=-1.2e-19;
  br(1)=+1.1306068068318665650e+0;
  br(2)=+1.1225343035826004e-3;
  br(3)=+8.8491541977453e-6;
  br(4)=+1.379157291154e-7;
  br(5)=+3.3014426995e-9;
  br(6)=+1.089423404e-10;
  br(7)=+4.6726552e-12;
  br(8)=+2.514325e-13;
  br(9)=+1.66086e-14;
  br(10)=+1.3312e-15;
  br(11)=+1.292e-16;
  br(12)=+1.53e-17;
  br(13)=+2.2e-18;
  br(14)=+4.0e-19;
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
  cr(1)=+1.1288167982669480212e+0;
  cr(2)=+2.175180808271174e-4;
  cr(3)=-1.2733503718850e-6;
  cr(4)=+2.33020578931e-8;
  cr(5)=-8.064911115e-10;
  cr(6)=+4.26054712e-11;
  cr(7)=-3.0537799e-12;
  cr(8)=+2.760646e-13;
  cr(9)=-2.99790e-14;
  cr(10)=+3.7789e-15;
  cr(11)=-5.392e-16;
  cr(12)=+8.55e-17;
  cr(13)=-1.48e-17;
  cr(14)=+2.8e-18;
  cr(15)=-6.0e-19;
  cr(16)=+1.0e-19;
  dr(1)=-1.0928773748447830587e-1;
  dr(2)=+2.0544241243325013e-4;
  dr(3)=-2.38739034179860e-6;
  dr(4)=+6.350265424372e-8;
  dr(5)=-2.82350803428e-9;
  dr(6)=+1.7939539824e-10;
  dr(7)=-1.485151455e-11;
  dr(8)=+1.50940584e-12;
  dr(9)=-1.8079342e-13;
  dr(10)=+2.478440e-14;
  dr(11)=-3.80487e-15;
  dr(12)=+6.4333e-16;
  dr(13)=-1.1825e-16;
  dr(14)=+2.339e-17;
  dr(15)=-4.93e-18;
  dr(16)=+1.10e-18;
  dr(17)=-2.60e-19;
  dr(18)=+6.0e-20;
  dr(19)=-2.0e-20;
  z= 2.0/3.0*(-x)^1.5; 
  la=5.0; t= -2*(la/x)^3.0-1.0;
  s=sin(z+piquart); c=cos(z+piquart);
  f3=chepolsumC(15,t,cr);
  f4=chepolsumC(18,t,dr);
  y=(-x)^0.25;
  aip=-y*(c*f3+s*f4/z);
  bip=y*(s*f3-c*f4/z);
elseif x < 0
  cr(1)=+5.2961692024697341466e-2;
  cr(2)=-4.8382429177635177789e-2;
  cr(3)=+6.2046464244529580518e-2;
  cr(4)=-3.1417437267239646847e-2;
  cr(5)=+7.876452021481851462e-3;
  cr(6)=-1.182440976973326921e-3;
  cr(7)=+1.18871496270269531e-4;
  cr(8)=-8.595270331212026e-6;
  cr(9)=+4.69655735896232e-7;
  cr(10)=-2.0107696526448e-8;
  cr(11)=+6.93493715818e-10;
  cr(12)=-1.9694289584e-11;
  cr(13)=+4.68795260e-13;
  cr(14)=-9.492371e-15;
  cr(15)=+1.65543e-16;
  cr(16)=-2.513e-18;
  cr(17)=+3.4e-20;
  dr(1)=+1.9598877757490956561e-1;
  dr(2)=+2.0709937287929773169e-1;
  dr(3)=-9.618436611491528534e-2;
  dr(4)=-2.6927080774066725642e-1;
  dr(5)=+1.5311467164195350955e-1;
  dr(6)=-3.621663427173734531e-2;
  dr(7)=+5.00970025411579034e-3;
  dr(8)=-4.6418943778727143e-4;
  dr(9)=+3.110071191129936e-5;
  dr(10)=-1.58422527393619e-6;
  dr(11)=+6.359244256854e-8;
  dr(12)=-2.06683376305e-9;
  dr(13)=+5.556077451e-11;
  dr(14)=-1.25683910e-12;
  dr(15)=+2.426805e-14;
  dr(16)=-4.0481e-16;
  dr(17)=+5.89e-18;
  dr(18)=-7.6e-20;
  la= 5.0; t= -2.0*(x/la)^3.0-1.0; z= x*x;
  f3=chepolsumC(16,t,cr);
  f4=chepolsumC(17,t,dr);
  aip=z*f3-f4;
  bip=sqrt3*(z*f3+f4);
elseif x==0.0
  aip=-0.25881940379280679841;
  bip=0.44828835735382635791;
elseif x < 4.5
  cr(1)=-1.0838288082415277302e+0;
      cr(2)=-2.75303325521307795e-2;
      cr(3)=+2.630996799464136716e-1;
      cr(4)=-2.20070838039707907e-2;
      cr(5)=-2.65686371901896143e-2;
      cr(6)=+4.6607302595359725e-3;
      cr(7)=+1.3331627500919505e-3;
      cr(8)=-3.993412832236136e-4;
      cr(9)=-2.57477942743769e-5;
      cr(10)=+1.92965112669086e-5;
      cr(11)=-8.179057497221e-7;
      cr(12)=-5.675046562446e-7;
      cr(13)=+7.38399111711e-8;
      cr(14)=+9.1592819462e-9;
      cr(15)=-2.5909300717e-9;
      cr(16)=-5.596039e-13;
      cr(17)=+5.48262416e-11;
      cr(18)=-4.2085028e-12;
      cr(19)=-6.882068e-13;
      cr(20)=+1.230089e-13;
      cr(21)=+2.2115e-15;
      cr(22)=-2.0268e-15;
      cr(23)=+1.092e-16;
      cr(24)=+2.01e-17;
      cr(25)=-2.8e-18;
      dr(1)=+6.5492503637358356971e-1;
      dr(2)=+7.398273602101549954e-2;
      dr(3)=+1.3916014396498561737e-1;
      dr(4)=-1.318037729688479868e-2;
      dr(5)=+2.721562752659802946e-2;
      dr(6)=-8.95735391942756359e-3;
      dr(7)=+4.64318555084023870e-3;
      dr(8)=-1.31982848083974091e-3;
      dr(9)=+3.2861337599865999e-4;
      dr(10)=-2.823577298366577e-5;
      dr(11)=-1.263661132378531e-5;
      dr(12)=+8.58819079977971e-6;
      dr(13)=-2.83605845640361e-6;
      dr(14)=+6.4357170147486e-7;
      dr(15)=-9.054345148828e-8;
      dr(16)=+1.02361813059e-9;
      dr(17)=+3.75568593189e-9;
      dr(18)=-1.26430133998e-9;
      dr(19)=+2.5000891998e-10;
      dr(20)=-2.967290357e-11;
      dr(21)=+2.6790190e-13;
      dr(22)=+8.3103543e-13;
      dr(23)=-2.2570250e-13;
      dr(24)=+3.487651e-14;
      dr(25)=-2.83929e-15;
      dr(26)=-1.5678e-16;
      dr(27)=+9.941e-17;
      dr(28)=-1.867e-17;
      dr(29)=+2.01e-18;
      dr(30)=-6.0e-20;
      dr(31)=-2.0e-20;
      la= 4.5; t= 2*x/la-1.0;
      f3=chepolsumC(24, t, cr);
      f4=chepolsumC(30, t, dr);
      aip=f3*exp(-1.375*x);
      bip=f4*exp(1.5*x);
elseif x < 9
      cr(1)=-1.2958123754838538881e+2;
      cr(2)=+6.76832374624909949e+0;
      cr(3)=+1.598511024718904852e+1;
      cr(4)=-1.38866512314607882e+0;
      cr(5)=-9.2855057862186510e-1;
      cr(6)=+1.1311806368323397e-1;
      cr(7)=+3.279015099087572e-2;
      cr(8)=-5.43676653553677e-3;
      cr(9)=-7.5136921782930e-4;
      cr(10)=+1.7836159217021e-4;
      cr(11)=+1.030315382664e-5;
      cr(12)=-4.29105073708e-6;
      cr(13)=-2.395648053e-8;
      cr(14)=+7.851108600e-8;
      cr(15)=-2.78229868e-9;
      cr(16)=-1.10405576e-9;
      cr(17)=+8.374716e-11;
      cr(18)=+1.165092e-11;
      cr(19)=-1.51205e-12;
      cr(20)=-8.180e-14;
      cr(21)=+2.018e-14;
      cr(22)=+1.2e-16;
      cr(23)=-2.1e-16;
      dr(1)=+5.5860183624034104454e-3;
      dr(2)=+2.87080283197283049e-5;
      dr(3)=+6.566960781790676321e-4;
      dr(4)=-1.29004618056797811e-5;
      dr(5)=+3.94734989962719105e-5;
      dr(6)=-1.7468162471848648e-6;
      dr(7)=+1.6312939322265420e-6;
      dr(8)=-1.095384680289487e-7;
      dr(9)=+5.23022431650618e-8;
      dr(10)=-4.4609927148093e-9;
      dr(11)=+1.3623610617155e-9;
      dr(12)=-1.277857258067e-10;
      dr(13)=+2.77519636083e-11;
      dr(14)=-2.2855189985e-12;
      dr(15)=+3.227684345e-13;
      dr(16)=+2.5537222e-15;
      dr(17)=-4.8715113e-15;
      dr(18)=+1.8406926e-15;
      dr(19)=-3.499972e-16;
      dr(20)=+5.20075e-17;
      dr(21)=-3.2250e-18;
      dr(22)=-9.750e-19;
      dr(23)=+5.043e-19;
      dr(24)=-1.464e-19;
      dr(25)=+3.39e-20;
      dr(26)=-6.7e-21;
      dr(27)=+1.2e-21;
      dr(28)=-2.0e-22;
      la= 4.5; t= 2*(x-la)/(9.0-la)-1.0;
      f3= chepolsumC(22, t, cr);
      f4= chepolsumC(27, t, dr);
      aip= f3*exp(-2.5*x);
      bip= f4*exp(2.625*x);
else
      cr(1)=+5.6568578662442917280e-1;
      cr(2)=+7.4362138971729004e-4;
      cr(3)=-4.42404425936355e-6;
      cr(4)=+5.505791762929e-8;
      cr(5)=-1.02083536642e-9;
      cr(6)=+2.483795569e-11;
      cr(7)=-7.4166337e-13;
      cr(8)=+2.607525e-14;
      cr(9)=-1.04943e-15;
      cr(10)=+4.738e-17;
      cr(11)=-2.36e-18;
      cr(12)=+1.3e-19;
      dr(1)=+1.1252717606343474943e+0;
      dr(2)=-1.5639596736331532e-3;
      dr(3)=-1.04063425452027e-5;
      dr(4)=-1.533581036012e-7;
      dr(5)=-3.5697709625e-9;
      dr(6)=-1.158278459e-10;
      dr(7)=-4.9121069e-12;
      dr(8)=-2.621683e-13;
      dr(9)=-1.72103e-14;
      dr(10)=-1.3725e-15;
      dr(11)=-1.326e-16;
      dr(12)=-1.56e-17;
      dr(13)=-2.2e-18;
      dr(14)=-4.0e-19;
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

     