%% load data
Q = xlsread('poblacioon_cuarentena.xlsx', 1, 'B2:B281');
Q=100/19458310.*Q;
%% analyze the effect of governement measurement wrt transmission rate
t1=13:(13+length(Q));
tau1 = interp1q(T(t_span)', tau_hat', t1');

%period
fig=figure();
plot(diff(log(tau1)), '-*b')
xlabel('t')
ylabel('log\tau(t) - log\tau(t-1)')
fig=figure();
plot(diff(log(tau1),5), '-*b') %the most regular with periodi 7 (a week)
fig=figure();
plot(diff(log(tau1),9), '-*b') %cloest to white noise 

% ARCH model
d5 = diff(log(tau1),2)';
nb_week = floor(length(d5)/7);
ind = repmat(1:7, 1, nb_week);
m = [];
for i=1:7
    m = [m sum((ind == i).*d5(1:(7*nb_week)))/nb_week];
end
d5 = d5(1:7*nb_week) - repmat(m, 1,nb_week);
fig=figure();
plot(d5,'*-')





%period
Q = Q(101: length(Q));
i = 5; %1: saturday
dlog_tau = diff(log(tau1));
dQ = [0 diff(Q)'];
fig=figure();
scatter(Q(find(ind == i)), dlog_tau(find(ind == i)))
b = regress(dlog_tau(find(ind == i)), [ones(size(find(ind == i)))' Q(find(ind == i))]);
hold on 
plot(Q(find(ind == i)), [ones(size(find(ind == i)))' Q(find(ind == i))]*b)
xlabel('quarantine percentage')
ylabel('log\tau(t) - log\tau(t-1)')
title('Friday')



fig=figure();
plot(Q)
scatter(Q'.*t1, log(tau1)) %we can see periodic trend
%we want to extend pierre's method log(tau(t)) = -mu(t)*t + f(mu(t)') + C

dQ = [0 diff(Q)']; 
figure()
plot(dQ)
scatter(cumsum(tau1), cumsum(Q))
figure()
scatter(cumsum(Q)'./(1:numel(Q)), tau1)
figure()
scatter(diff(Q(20:250), 2), diff(tau1(20:249)))
scatter(Q(1:100), tau1)

%% linear regression
b = regress(diff(log(tau1(1: 51))), [ones(size(Q(1: 50))), Q(1: 50), dQ(1: 50)']);
figure()
plot(diff(log(tau1(1: 51))))
hold on
plot([ones(size(Q(1: 50))), Q(1: 50), dQ(1: 50)']*b)

figure()
plot(diff(log(tau1(51: 71))))
hold on
plot([ones(size(Q(51: 70))), Q(51: 70), dQ(51: 70)']*b)


%% find the functional relationship between Q and tau
n_Q = length(Q);
d = 10;
Y = ones(n_Q-d,1);
for i=0:d
   Y = [Y Q((i+1):(n_Q-d+i),1)];
end
for i=0:d
    for j=i:d
        Y = [Y Q((i+1):(n_Q-d+i),1).*Q((j+1):(n_Q-d+j),1)];
    end
end
n_pred = 50;
b = regress(tau1((d+1):(n_Q-n_pred)),Y(1:(length(Y)-n_pred),:));

norm(Y*b- tau1((d+1):n_Q),2)

figure()
plot((d+1):n_Q, Y*b, '-b',(d+1):n_Q, tau1((d+1):n_Q), 'r')
hold on
plot((n_Q-n_pred+1):n_Q, Y((length(Y)-n_pred+1):length(Y),:)*b, ...
    '.b',(n_Q-n_pred+1):n_Q, tau1((n_Q-n_pred+1):n_Q,:), '.r')

%% slove SIRU with tau_hat
f=0.3; %affect the final size of epidemic 
v=1/7; 
n=1/7;
v2=(1-f)*v;
v1=f*v;

X = 1:20;
[x1,x2,S0,I_0,U_0,t_0]=conin(X,CR_data(1, X));

t11 = [t_0 1:(13+length(Q)-1)];
tau11 = interp1q(T(t_span)', tau_hat', t11');
tau11(24:length(tau11),1) = Y*b;

[T,S,I,R,U]=siru_slover(t_0,t11(length(t11)),tau11,I_0,U_0,f,v,n,h);
CR=[];
CU=[];
for i=1:length(T)
    CRi =v1*trapz(T(1,1:i),I(1,1:i),2);
    CUi =v2*trapz(T(1,1:i),I(1,1:i),2);
    CR=[CR CRi];
    CU=[CU CUi];
end

fig=figure();
plot(T, I, 'r')

fig=figure();
plot(T, CR, 'r')

fig=figure();
plot(T, R, 'r')
hold on
plot(T, U, 'k')

figure()
plot(tau11)
figure()
plot(tau_hat)

fig=figure();
plot(R_hat, 'r')














