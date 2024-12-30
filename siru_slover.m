
%
function [T,S,I,R,U]=siru_slover(t0,tau,S0,I0,R0,U0,f,v,n,h)

% here we are working with RUNGE-KUTTA method

%S0=19458310;
%S0=11.081*10^6; %wuham

% calculo de parameters
v2=(1-f)*v;
v1=f*v;
%R0=(tau0*S0/v)*(1+(v2/n));
%R0 = 0;

% % define functions handles
% fS=@(t,S,I,R,U) -tau*S*(I+U);
% fI=@(t,S,I,R,U) tau*S*(I+U) - v*I;
% fR=@(t,S,I,R,U) v1*I-n*R ;
% fU=@(t,S,I,R,U) v2*I-n*U;

% initial conditions
t(1)=t0;
S(1)=S0;
I(1)=I0;
R(1)=R0;
U(1)=U0;

N=length(tau)-1;

% update loop
T=[t0];
for i=1:N
    % update time
    t(i+1)=t(i)+h;
    T=[T t(i+1)];
    % define functions handles

    fS=@(t,S,I,R,U) -tau(i+1)*S*(I+U);
    fI=@(t,S,I,R,U) tau(i+1)*S*(I+U) - v*I;
    fR=@(t,S,I,R,U) v1*I-n*R ;
    fU=@(t,S,I,R,U) v2*I-n*U;

    % update functions
    k1S=h*fS(t(i),S(i),I(i),R(i),U(i));
    k1I=h*fI(t(i),S(i),I(i),R(i),U(i));
    k1R=h*fR(t(i),S(i),I(i),R(i),U(i));
    k1U=h*fU(t(i),S(i),I(i),R(i),U(i));

    %
    k2S=h*fS(t(i)+h/2,S(i)  + 1/2 *k1S    ,I(i) + 1/2 *k1I  ,R(i) + 1/2*k1R  ,U(i) + 1/2*k1U);
    k2I=h*fI(t(i)+h/2,S(i)  + 1/2 *k1S    ,I(i) + 1/2 *k1I  ,R(i) + 1/2*k1R  ,U(i) + 1/2*k1U);
    k2R=h*fR(t(i)+h/2,S(i)  + 1/2 *k1S    ,I(i) + 1/2 *k1I  ,R(i) + 1/2*k1R  ,U(i) + 1/2*k1U);
    k2U=h*fU(t(i)+h/2,S(i)  + 1/2 *k1S    ,I(i) + 1/2 *k1I  ,R(i) + 1/2*k1R  ,U(i) + 1/2*k1U);

    %
    k3S=h*fS(t(i)+h/2,S(i)  + 1/2 *k2S    ,I(i) + 1/2*k2I  ,R(i) + 1/2*k2R  ,U(i) + 1/2*k2U);
    k3I=h*fI(t(i)+h/2,S(i)  + 1/2 *k2S    ,I(i) + 1/2*k2I  ,R(i) + 1/2*k2R  ,U(i) + 1/2*k2U);
    k3R=h*fR(t(i)+h/2,S(i)  + 1/2 *k2S    ,I(i) + 1/2*k2I  ,R(i) + 1/2*k2R  ,U(i) + 1/2*k2U);
    k3U=h*fU(t(i)+h/2,S(i)  + 1/2 *k2S    ,I(i) + 1/2*k2I  ,R(i) + 1/2*k2R  ,U(i) + 1/2*k2U);
    %
    k4S=h*fS(t(i)+h ,S(i)  + h*k3S    ,I(i) + 1*k3I  ,R(i) + 1*k3R  ,U(i) + 1*k3U);
    k4I=h*fI(t(i)+h ,S(i)  + h*k3S    ,I(i) + 1*k3I  ,R(i) + 1*k3R  ,U(i) + 1*k3U);
    k4R=h*fR(t(i)+h ,S(i)  + h*k3S    ,I(i) + 1*k3I  ,R(i) + 1*k3R  ,U(i) + 1*k3U);
    k4U=h*fU(t(i)+h ,S(i)  + h*k3S    ,I(i) + 1*k3I  ,R(i) + 1*k3R  ,U(i) + 1*k3U );
    % already the functions
    S(i+1)=S(i)+1/6*(k1S + 2*k2S + 2*k3S + k4S);
    I(i+1)=I(i)+1/6*(k1I + 2*k2I + 2*k3I + k4I);
    R(i+1)=R(i)+1/6*(k1R + 2*k2R + 2*k3R + k4R);
    U(i+1)=U(i)+1/6*(k1U + 2*k2U + 2*k3U + k4U);
end
end