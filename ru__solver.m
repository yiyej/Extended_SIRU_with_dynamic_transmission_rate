
function [T,R,U]=ru__solver(t_0,I_hat,h,U_0,R_0,f,v,n,appx_md)

% here we are working with RUNGE-KUTTA method

% calculo de parameters
v2=(1-f)*v;
v1=f*v;

% initial conditions
t(1)=t_0;
R(1)=R_0;
U(1)=U_0;

N=length(I_hat)-1;

% update loop
T=[t_0];
for i=1:N 
    % update time
    t(i+1)=t(i)+h;
    T=[T t(i+1)];
    % define functions handles

    if appx_md == "poly"
        I = polyval(I_hat,t(i+1));
    elseif  appx_md == "spline"
        I = ppval(I_hat,t(i+1))./v1;
    else
        I = I_hat(i+1);
    end

    fR=@(t,R) v1*I-n*R ;
    fU=@(t,U) v2*I-n*U;

    % update functions
    k1R=h*fR(t(i),R(i));
    k1U=h*fU(t(i),U(i));
    %
    k2R=h*fR(t(i)+h/2, R(i) + 1/2*k1R);
    k2U=h*fU(t(i)+h/2, U(i) + 1/2*k1U);
    %
    k3R=h*fR(t(i)+h/2, R(i) + 1/2*k2R);
    k3U=h*fU(t(i)+h/2, U(i) + 1/2*k2U);
    %
    k4R=h*fR(t(i)+h , R(i) + 1*k3R );
    k4U=h*fU(t(i)+h ,U(i) + 1*k3U);
    % already the functions
    R(i+1)=R(i)+1/6*(k1R + 2*k2R + 2*k3R + k4R);
    U(i+1)=U(i)+1/6*(k1U + 2*k2U + 2*k3U + k4U);
end
end