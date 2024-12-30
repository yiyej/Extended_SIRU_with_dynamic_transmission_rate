%% validate the estimation of tau_hat
% test to slove SIRU with the estimated tau_hat 
% it should recover the true CR

[T_pred,S_pred,I_pred,R_pred,U_pred]=...
    siru_slover(t_0,ppval(pp, T(1:300/h))./10^10,S0,I_hat(1),R0,U_0,f,v,n,h);

CR_pred=[];
CU_pred=[];
for i=1:length(I_pred)
    CRi =v1*trapz(T(1:i),I_pred(1:i),2);
    CUi =v2*trapz(T(1:i),I_pred(1:i),2);
    CR_pred=[CR_pred CRi];
    CU_pred=[CU_pred CUi];
end
fig=figure();
plot([t_0 [X x]], [0 CR_data(1,[X x])],'b*',T_pred,CR_pred,'.r')
hold on 
plot(T_pred,CU_pred,'.k')
xlabel('t')
legend({'CR data','CR hat', 'CU hat'},'Location','northwest')
figure()
plot(T(t_span), I_hat(t_span), 'r')
hold on
plot(T_pred, I_pred, 'k')

%% focus on the second phase 
CR2_0 = interp1q(T_pred', CR_pred', t_int(1));
I2_0 = interp1q(T(t_span)', I_hat(t_span)', t_int(1));
S2_0 = interp1q(T(t_span)', S_hat(t_span)', t_int(1));
U2_0 = interp1q(T(t_span)', U_hat(t_span)', t_int(1));
R2_0 = interp1q(T(t_span)', R_hat(t_span)', t_int(1));
tau2 = interp1q(t_int', tau_int, (t_int(1):h:t_int(length(t_int)))');
tau2_pred = interp1q(t_int', tau_pred', (t_int(1):h:t_int(length(t_int)))');
[T_pred,S_pred,I_pred,R_pred,U_pred]=...
    siru_slover(t_int(1),tau2_pred,S2_0,I2_0,R2_0,U2_0,f,v,n,h);
figure()
plot(T(t_span), I_hat(t_span), 'r')
hold on
plot(T_pred, I_pred, 'k.')

CR_pred=[];
CU_pred=[];
for i=1:length(I_pred)
    CRi =v1*trapz(T_pred(1:i),I_pred(1:i),2);
    CUi =v2*trapz(T_pred(1:i),I_pred(1:i),2);
    CR_pred=[CR_pred CRi];
    CU_pred=[CU_pred CUi];
end
CR_pred = CR2_0+CR_pred;
fig=figure();
plot([t_0 [X x]], [0 CR_data(1,[X x])],'b*',T_pred,CR_pred,'.r')
hold on 
plot(T_pred,CU_pred,'.k')
xlabel('t')
legend({'CR data','CR hat', 'CU hat'},'Location','northwest')
%% training test split 
% t=13 <-> 14 Mars (Saturday)
% training set: the data of following 20 weeks since t=13 
% test set: the data of the following 10 weeks since t=13+20wk

% load observations of the covariate variable: quarantine percentage  
Q = xlsread('poblacioon_cuarentena.xlsx', 1, 'B2:B281');
Q=100/19458310.*Q;

% interpolate tau_hat to the common recording time stamps as Q
nb_week = 30;
D = -90; %or 15
t_int=(110+D):(110+D+nb_week*7);
tau_int = interp1q(T(t_span)', tau_hat', t_int');
dlogtau = diff(log(tau_int));
Q = Q((98+D):(98+D+nb_week*7-1));
dQ = [0 diff(Q)'];

figure();
scatter(Q, tau_int(1:210))
figure();
plot(t_int(1:nb_week*7), dlogtau,'b*-')

nb_week_tn = 20;
t_tn = t_int(1:nb_week_tn*7);
dlogtau_tn = dlogtau(1:nb_week_tn*7);
Q_tn = Q(1:nb_week_tn*7);
t_test = t_int((nb_week_tn+1)*7:nb_week*7);
dlogtau_test = dlogtau((nb_week_tn*7+1):nb_week*7);
Q_test = Q((nb_week_tn*7+1):nb_week*7);

% we observe that the impact of covariate variable on transmission rate
% varies across the specific days on a week, thus we propose to set up 
% 7 regression models, each corresponding to a day. When forcasting the 
% transmission rate of a date, we need not only to know the covariate value 
% of that day but also which day it is in a week

% label all dates in the data set with their day information
ind = repmat([6 7 1:5], 1, nb_week);
%dlogtau_pred = Q_test; %initialization

%set up the regression models and forcast 
dlogtau_pred = dlogtau;
blist = ones([6,7]);
for i=1:7
    %fitting/parameter estimation
    ind_tn_i = find(ind(1:nb_week_tn*7) == i);
    b =regress(dlogtau_tn(ind_tn_i),[ones(size(ind_tn_i))' ...
        Q_tn(ind_tn_i) Q(ind_tn_i).^2 Q(ind_tn_i).^3 Q(ind_tn_i).^4 dQ(ind_tn_i)']);
    %forcasting 
%     ind_test_i = find(ind((nb_week_tn*7+1):nb_week*7) == i);
%     dlogtau_pred(ind_test_i)=[ones(size(ind_test_i))' Q_test(ind_test_i)]*b;
%     scatter(Q_test(ind_test_i), dlogtau_test(ind_test_i), 'g')
%     hold on
%     scatter(Q_test(ind_test_i), dlogtau_pred(ind_test_i),'d')
%     title(i)
    %predict
    ind_i = find(ind == i);
    dlogtau_pred(ind_i)=[ones(size(ind_i))' ...
        Q(ind_i) Q(ind_i).^2 Q(ind_i).^3 Q(ind_i).^4 dQ(ind_i)']*b;
    blist(:,i) = b;
end 


figure();
plot(dlogtau, '-*')
hold on
plot(dlogtau_pred, '-r*')
hold on
plot(zeros(size(exp(dlogtau_pred))), '-k')

%recover tau(t)
logtau = 1:(length(t_test)+1);
logtau(1) = log(tau_int((nb_week_tn*7+1)));
for j=1:length(t_test)
     logtau(j+1) = logtau(j) + dlogtau_pred(j);
end
figure();
plot(tau_int((nb_week_tn*7+1):(nb_week*7+1)), 'b*-')
hold on
figure();
plot(exp(logtau), '-*r')
figure();
plot(Q((nb_week_tn*7+1):(nb_week*7+1)))

tau_test = log(tau_int(1)) + cumsum(dlogtau);
tau_test = [tau_int(1) exp(tau_test)'];
tau_pred = log(tau_int(1)) + cumsum(dlogtau_pred);
tau_pred = [tau_int(1) exp(tau_pred)'];
figure();
plot(tau_test)
hold on
plot(tau_pred)

figure()
plot(13:(13+nb_week*7-1),Q(1:nb_week*7),'*-')







