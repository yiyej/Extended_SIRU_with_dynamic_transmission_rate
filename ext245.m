figure()
plot(CR_data)

%% predict to the next extreme point of tau(t)
% assume we have CR data until t_0
% we now want to predict CR into the future 
% method
% step 1: from t, we predict the next extreme point of tau using 
%         logistic regression
% step 2: we use the LS estimation with step 1 constraint to get 
%         to prediction of CR 

t_inv_eq =120+10:185;
%chebpoints
N = length(t_inv_eq);
t_inv = floor((-cos((0:N).*pi/N) + 1)*N/2+t_inv_eq(1));

CR = @(x, t) [t'.^5 t'.^4 t'.^3 t'.^2 t' ones(size(t'))]*x;
data_term = @(x) mean((CR(x,t_inv) - CR_data(t_inv)').^2);

I = @(x, t) v1*[5*t'.^4 4*t'.^3 3*t'.^2 2*t' ones(size(t')) zeros(size(t'))]*x;
dI = @(x, t) v1*[20*t'.^3 12*t'.^2 6*t' 2*ones(size(t')) zeros(size(t')) zeros(size(t'))]*x;
d2I = @(x, t) v1*[60*t'.^2 24*t' 6*ones(size(t')) zeros(size(t')) ...
    zeros(size(t')) zeros(size(t'))]*x;
d3I = @(x, t) v1*[120*t' 24*ones(size(t')) 6*ones(size(t')) zeros(size(t')) ...
    zeros(size(t')) zeros(size(t'))]*x;

S_ = @(x, t) -I(x, t)-CR(x, t)/f;
S = @(x, t) S_(x, t) ...
    + mean(interp1q(T(t_span)', S_hat(t_span)', t_inv(1:end-10)')...
      - S_(x,t_inv(1:end-10)));
dS = @(x, t) -dI(x, t)-v*I(x, t);

U_ = @(x, t) v1*I(x, t)/n -v1*dI(x, t)/n/n + ...
    v1*d2I(x, t)/n/n/n - v1*d3I(x, t)/n/n/n/n;
U = @(x, t) U_(x, t) ...
    + mean(interp1q(T(t_span)', U_hat(t_span)', t_inv(1:end-10)')...
  - U_(x,t_inv(1:end-10)));
dU = @(x, t) v1*dI(x, t)/n -v1*d2I(x, t)/n/n;

dtau=@(x, t) (d2I(x,t)+v*dI(x,t))./S(x,t)./(I(x,t)+U(x,t))...
    -(dI(x,t)+v*I(x,t)).*dS(x,t)./(S(x,t).^2)./(I(x,t)+U(x,t))...
    -(dI(x,t)+v*I(x,t)).*(dI(x,t)+dU(x,t))./S(x,t)./((I(x,t)+U(x,t)).^2);

p = polyfit(t_inv, CR_data(t_inv), 5);
x0 = p';
options = optimset('MaxFunEvals',50000);
x_star = fminsearch(data_term,x0,options);
figure()
plot(CR_data)
hold on
plot(t_inv(1):t_inv(end)+50, CR(x_star, t_inv(1):t_inv(end)+50))
hold on
plot(t_inv(1):t_inv(end)+50, polyval(p,(t_inv(1):t_inv(end)+50)'))

% the "next" extreme points
t_E = 245;

lambda = 1e+31;
L = @(x) data_term(x) + lambda*dtau(x,t_E).^2;
x0 = p';
options = optimset('MaxFunEvals',50000);
x_star = fminsearch(L,x0,options);
figure()
plot(CR_data)
hold on
plot(t_inv(1):t_E(end), CR(x_star, t_inv(1):t_E(end)),'r-')
hold on
plot(t_inv(end):t_E(end), CR(x_star, t_inv(end):t_E(end)),'r.-')
title([t_E(end) lambda t_inv(1)])

%tune the hyperparameter: lambda 
lg_vld = 10;
grid = [1e+26 1e+27 1e+28]; 
err_vd = grid;
for n=1:length(grid)
    L = @(x) data_term(x) + grid(n)*dtau(x,t_E).^2;
    x0 = p';
    options = optimset('MaxFunEvals',50000);
    x_star = fminsearch(L,x0,options);
    err_vd(n) = mean((CR_data(t_inv(end):t_inv(end)+lg_vld)' ...
                - CR(x_star, t_inv(end):t_inv(end)+lg_vld)).^2);
end 

lambda = 1e+27;
L = @(x) data_term(x) + lambda*dtau(x,t_E).^2;
x0 = p';
options = optimset('MaxFunEvals',50000);
x_star = fminsearch(L,x0,options);
figure()
plot(CR_data(1:310))
hold on
plot(t_inv(1):t_inv(end)+40, CR(x_star, t_inv(1):t_inv(end)+40),'r-')
hold on
plot(t_inv(end)+10:t_inv(end)+40, CR(x_star, t_inv(end)+10:t_inv(end)+40),'r.-')
title([t_E(end) lambda t_inv(1)])

%% S I R U
dCR = diff(CR(x_star, t_inv(end)+10:h:t_inv(end)+40))/h; 
I = dCR./v1;
fig = figure();
plot(t_inv(end)+10:h:t_inv(end)+40-h, I)
hold on 
plot(t_inv(end)+10:t_inv(end)+40, ...
    newcases(t_inv(end)+10:t_inv(end)+40)./v1, 'k.')
legend({'I\_pred','daily new cases'},'Location','northwest')
saveas(fig,"I_pred",'epsc')

%obtain initial values
ind = find(T <= t_inv(end)+10);
t2_0 = t_inv(end)+10;
S2_0 = S_hat(ind(end));
U2_0 = U_hat(ind(end));
R2_0 = R_hat(ind(end));
dI = diff(I)/h;

[T2,R,U]=ru__solver(t2_0,I,h,U2_0,R2_0,f,v,n,"kernel");
fig = figure();
plot(T(23000:30000), U_hat(23000:30000), 'k')
hold on
plot(T(23000:30000), R_hat(23000:30000), 'r')
hold on
plot(T2, R, 'r-.')
hold on
plot(T2, U, 'k-.')
legend({'R\_data','U\_data','R\_pred','U\_pred'},'Location','northwest')
saveas(fig,"RU_pred",'epsc')





