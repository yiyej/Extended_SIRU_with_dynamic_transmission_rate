xx = 0:350;
y= 0.1*tau0*exp(-0.014*xx) + tau0*0.001.*xx + 0.025*tau0;
figure()
plot(xx, y)
hold on
plot(T(t_span), tau_hat, '.')
tau_hat2 = tau0*exp(-0.014*T(t_span)) + tau0*0.001.*T(t_span) + 0.025*tau0;

Q = xlsread('poblacioon_cuarentena.xlsx', 1, 'B2:B281');
Q=100/19458310.*Q;
t_int=1:(length(Q)+13-1);
Q = [zeros(1,12) Q'];
p = ksr(t_int, Q, 3, 1:25:(length(Q)+13-1)); %29
fig = figure();
plot(t_int, Q, '*-')
hold on
plot(t_int, ppval(pp, t_int))
xticks([1 60 120 180 240 300])
xticklabels({'3/2/2020','4/30/2020','6/29/2020','8/28/2020','10/27/2020',...
    '12/26/2020'})
saveas(fig,"Q",'pdf')

plot(p.x, p.f)
pp=spline(p.x, p.f);
figure()
plot(t_int, ppval(pp, t_int))
[breaks,coefs,l,k,d] = unmkpp(pp);
pp2_Q = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
figure()
plot(t_int, ppval(pp2_Q, t_int))
[breaks,coefs,l,k,d] = unmkpp(pp2_Q);
pp3 = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
[breaks,coefs,l,k,d] = unmkpp(pp3);
pp4 = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
figure()
plot(t_int, ppval(pp3, t_int), '*-')
hold on
plot(t_int, zeros(size(t_int)))
figure()
plot(t_int, ppval(pp3, t_int).*sign(ppval(pp4, t_int)), '*-')
hold on
plot(t_int, zeros(size(t_int)))
inf_pt_Q = [];
for i=1:30
    root = bisectionMethod(pp3, (i-1)*10, i*10, 0.0001);
    if ~isnan(root)
        inf_pt_Q = [inf_pt_Q floor(root)];
    end
end

figure()
plot(t_int, ppval(pp3, t_int), '*-')
hold on
plot(t_int, zeros(size(t_int)))
hold on
plot(inf_pt_tau, ppval(pp3, inf_pt_tau), '*r')
hold on
plot(noinf_pt_tau, ppval(pp3, noinf_pt_tau), '*k')

tau_int = interp1q(T(t_span)', tau_hat', t_int');
%p = ksr(t_int, tau_int.*10^10, 6, t_int(1:10:length(t_int)));
%figure()
%plot(t_int, tau_int.*10^10)
%hold on
%plot(p.x, p.f)
%lambda  = 10 for ksr initial and here
pp=spline(t_int(1:6:length(t_int)), tau_int(1:6:length(t_int)).*10^10);
figure()
plot(t_int, tau_int.*10^10)
hold on
plot(t_int, ppval(pp, t_int), 'r')
[breaks,coefs,l,k,d] = unmkpp(pp);
pp2 = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
figure()
plot(t_int, ppval(pp2, t_int))
hold on
plot(t_int, zeros(size(t_int)))
inf_pt_tau  = [];
for i=1:30
    root = bisectionMethod(pp2, (i-1)*10, i*10, 0.001);
    if ~isnan(root)
        inf_pt_tau = [inf_pt_tau floor(root)];
    end
end
inf_pt_tau
figure()
plot(t_int, ppval(pp3, t_int), '*-')
hold on
plot(t_int, zeros(size(t_int)))
hold on
plot(noinf_pt_tau, ppval(pp3, noinf_pt_tau), '*k')
hold on
plot(inf_pt_tau_tn, ppval(pp3, inf_pt_tau_tn), '*r')

%% logistic regression
noinf_pt_tau = [floor((1 + inf_pt_tau(1))/2)-3: ...
    floor((1 + inf_pt_tau(1))/2)+3];
for i=1:3
    tmp = (inf_pt_tau(i) + inf_pt_tau(i+1))/2;
    noinf_pt_tau = [noinf_pt_tau floor(tmp)-3:floor(tmp)+3];
end 
inf_pt_tau_tn = [];
for i=1:4
    inf_pt_tau_tn = [inf_pt_tau_tn (inf_pt_tau(i)-3):(inf_pt_tau(i)+3)];
end 

X = ppval(pp3, [inf_pt_tau_tn noinf_pt_tau]);
X = [ones(size(X')) X' X'.^2 ...
    ppval(pp2_Q, [inf_pt_tau_tn noinf_pt_tau])' ...
    ppval(pp2_Q, [inf_pt_tau_tn noinf_pt_tau])'.^2];
y = [ones(size(inf_pt_tau_tn)),zeros(size(noinf_pt_tau))]';

theta = ones(5, 1);
iter_max = 200000;
eta = 0.005;
J = 1:iter_max;
%gradient descent 
for n=1:iter_max
    [J(n), grad] = costFunction(theta, X, y);
    theta = theta - eta*grad;
end
figure()
plot(J)

prob = sigmoid(X*theta);
figure()
plot(prob, '.-')
d2 = -1:0.01:1;
prob = sigmoid([ones(size(d2')),d2', d2'.^2]*theta(1:3));
figure()
plot(d2,prob)

X = ppval(pp3, t_int);
X = [ones(size(X')) X' X'.^2 ppval(pp2_Q, t_int)' ...
    ppval(pp2_Q, t_int)'.^2];
prob = sigmoid(X*theta);
figure()
plot(t_int, prob, '*-')
title([eta 25 8]')














