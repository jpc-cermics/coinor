// -*- Mode: scilab -*-
// Copyright (C) 2017-2017 Aurelien Alfonsi, Jean-Philippe Chancelier Enpc/Cermics
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//
// a test for a quadratic problem for which we have
// an other way to compute the solution

function [rp,rw]=convex_hull_psi(p_m,w_m,p_n,w_n)
  //#calcule la conction psi enveloppe convexe de int_0^q F^{-1}_m(u)-F^{-1}(u)du
  L=[ ]
  VAL=[]
  L.concatr[0.]
  VAL.concatr[0.]
  s_m=cumsum(w_m)
  s_n=cumsum(w_n)
  l_m=length(w_m)
  l_n=length(w_n)
  i_m=1
  i_n=1
  ta=0.
  val=0.
  while (i_m <= l_m)&(i_n <= l_n) do
    t=min(s_m(i_m),s_n(i_n))
    val=val+(t-ta)*(p_m(i_m)-p_n(i_n))
    while length(L)>1 && (VAL($)-VAL($-1))/(L($)-L($-1))>=(val-VAL($))/(t-L($)) do        
      L($)=[]
      VAL($)=[]
    end
    L.concatr[t];
    VAL.concatr[val];
    ta=t
    if t==s_m(i_m) then
      i_m=i_m+1
    end
    if t==s_n(i_n) then
      i_n=i_n+1
    end
  end
  [rp,rw]=(L,VAL);
endfunction;

function [rp,rw]=curlywedge(p_m,w_m,p_n,w_n)
  [L,VAL]=convex_hull_psi(p_m,w_m,p_n,w_n)   
  DL= L ; 
  DVAL= VAL;
  // DL=DL[1:]-DL[:-1]
  DL=DL(2:$) -DL(1:$-1);
  DVAL=(DVAL(2:$) -DVAL(1:$-1))./ DL
  i_m=1
  l_m=length(w_m)
  i_L=1
  l_L=length(DL)
  s_m=cumsum(w_m)
  s_L=cumsum(DL)
  ta=0.
  L=[]
  VAL=[]
  while (i_m<=l_m)&(i_L<=l_L) do
    t=min(s_m(i_m),s_L(i_L))
    VAL.concatr[p_m(i_m)-DVAL(i_L)]
    L.concatr[t-ta]
    if t==s_m(i_m) then
      i_m=i_m+1
    end
    if t==s_L(i_L) then
      i_L=i_L+1
    end
    ta=t
  end
  [rp,rw]=(VAL,L);
endfunction

function res=wasserstein(p_m,w_m,p_n,w_n,rho)    
  i_n=1; l_n=length(w_n);
  i_m=1; l_m=length(w_m)
  s_n=cumsum(w_n);
  s_m=cumsum(w_m);
  ta=0.;
  wass=0.;
  while (i_m < l_m)&&(i_n < l_n) do
    t=min(s_m(i_m),s_n(i_n))
    wass=wass+(t-ta)*abs(p_m(i_m)-p_n(i_n))^rho
    if t==s_n(i_n) then
      i_n=i_n+1
    end
    if t==s_m(i_m) then 
      i_m=i_m+1
    end
    ta=t;
  end
  res= wass^(1.0/rho)   
endfunction

function [p,m]=filter(p_m,w_m)
  // remove points with very low proba.
  // we assume that p_m and w_m must have the same size
  I=find(w_m> 1e-12)
  p=p_m(I);
  m=w_m(I);
endfunction

// Create a clp quadratic problem
N=10
N2=N*N;

X=2*rand(1,N)-1;
Y=4*rand(1,N)-2;

p=zeros(1,N2);
M=zeros(N2,N2);
C=zeros(2*N,N2);
b=ones(2*N,1); 
me=2*N;

for i=1:N,
  for j=1:N,
    p((i-1)*N+j)=-X(i)*Y(j);
    C(i,(i-1)*N+j)=1;
    C(N+j,(i-1)*N+j)=1;
    for k=1:N,
      M((i-1)*N+j,(i-1)*N+k)=Y(j)*Y(k);
    end;
  end;
end;

// solve the  quadratic problem
[xopt,fopt,flag,lambda] = linprog_clp(p,zeros(0,N2),[],C,b,sense="min",Q=M);

Qopt=zeros(N,N);
for i=1:N,
  for j=1:N,
    Qopt(i,j)=xopt((i-1)*N+j)
  end;
end;

copt=norm(Qopt*Y'-X')^2

// Same problem using sparse representation

if %f then 
  T=timer();
  N=200
  N2=N*N;
  X=grand(N,"mn",0,1);
  Y=grand(N,"mn",0,1.1);
  p=zeros(1,N2);
  M=sparse([],[],[N2,N2]);
  C=sparse([],[],[2*N,N2]);
  b=ones(2*N,1); 
  for i=1:N,
    for j=1:N,
      p((i-1)*N+j)=-X(i)*Y(j);
      C(i,(i-1)*N+j)=1;
      C(N+j,(i-1)*N+j)=1;
      for k=1:N,
	M((i-1)*N+j,(i-1)*N+k)=Y(j)*Y(k);
      end;
    end;
  end;
  p=p';
  [xopt,fopt,flag,lambda] = linprog_clp(p',sparse(zeros(0,N2)),[],C,b,sense="min",Q=M);
  T=timer()
  
  Qopt=zeros(N,N);
  for i=1:N,
    for j=1:N,
      Qopt(i,j)=xopt((i-1)*N+j)
    end;
  end;
  copt=norm(Qopt*Y'-X')^2
  fprintfMat("newX_1d.txt",Qopt*Y')
  fprintfMat("Y_1d.txt",Y')
  fprintfMat("X_1d.txt",X')
end

// Compare to known solution

newX= Y*Qopt';
newX=sort(newX,dir='i');
X=sort(X,dir='i');
Y=sort(Y,dir='i' );

wX=ones(1,length(X))*(1./length(X))
wY=ones(1,length(Y))*(1./length(Y))

p_m=X;w_m=wX;p_n=Y;w_n=wY;
[L,VAL]=convex_hull_psi(p_m,w_m,p_n,w_n)

[pt_cw,w_cw]=curlywedge(X,wX,Y,wY)
[pt_cw,w_cw]=filter(pt_cw,w_cw)

if %f then
  function y=phi(P,W,t)
    y=sum(W.*max(t-P,0.))
  endfunction

  // graphics 
  Npt=100
  a=-3.0
  b=3.0
  pas=(b-a)/Npt
  T= linspace(a,b,Npt);
  phi_nu= 1:Npt
  phi_mu=1:Npt
  phi_cw=1:Npt
  phi_CO=1:Npt
  for i = 1:Npt do
    phi_nu(i)=phi(Y,wY,T(i))
    phi_mu(i)=phi(X,wX,T(i))
    phi_cw(i)=phi(pt_cw,w_cw,T(i))
    phi_CO(i)=phi(newX,wX,T(i))    
  end
  hold('on');
  plot(T,phi_nu,'r');
  plot(T,phi_mu,'g')   
  plot(T,phi_CO,'b')
  plot(T,phi_cw,'m')
end

w1=wasserstein(pt_cw,w_cw,X,wX,2);
w2=wasserstein(newX,wX,X,wX,2);
w3=wasserstein(pt_cw,w_cw,newX,wX,2);
if abs(w1-w2) >= 1.e-5 then pause;end

//printf("W2 curlywedge %f\n",wasserstein(pt_cw,w_cw,X,wX,2))
//printf("W2 avec Coin-OR %f\n",wasserstein(newX,wX,X,wX,2))
//printf("W2 curlywedge/Coin-OR %f\n",wasserstein(pt_cw,w_cw,newX,wX,2))

