function [Y,Nbins]=GenerateThalamicSpikeTrains(bins,pfail,pspont,Nrep)

load('spike_trains.mat');

% read spike trains
A       = cell2mat(spike_trains(1,1));
a       = length(A);
B       = cell2mat(spike_trains(1,2));
b       = length(B);
C       = cell2mat(spike_trains(1,3));
c       = length(C);
D       = cell2mat(spike_trains(1,4));
dur       = length(D);
E       = cell2mat(spike_trains(1,5));
e       = length(E);
V       = zeros(a+b+c+dur+e,1);

% Creating a vector containg all the time between the spikes of spike
% trains A,B,C,D and E
V(1,1)  = A(1,1);
for i= 1:a-1
    V(i+1,1) = A(i+1,1)-A(i,1);
end

V(a+1,1) = B(1,1);
for i= 1:b-1
    V(a+i+1,1) = B(i+1,1)-B(i,1);
end

V(a+b+1,1) = C(1,1);
for i= 1:c-1
    V(a+b+i+1,1) = C(i+1,1)-C(i,1);
end

V(a+b+c+1,1) = D(1,1);
for i= 1:dur-1
    V(a+b+c+i+1,1) = D(i+1,1)-D(i,1);
end

V(a+b+c+dur+1,1) = E(1,1);
for i= 1:e-1
    V(a+b+c+dur+i+1,1) = E(i+1,1)-E(i,1);
end

% Times in ms
V = V*1000;

% Plot histogram and cdf
figure('Name','Histogram and cdf');
fig.a = axes; 
hold(fig.a,'all');
fig.hist = histogram(V,round(max(V)/3),'Normalization','probability');
fig.cdf = cdfplot(V);

% fit cdf
x   = transpose(get(fig.cdf,'XData'));
y   = transpose(get(fig.cdf,'YData'));
t   = ~isinf(x) & ~isinf(y);
% coeffNames      = {'p1','p2','p3','p4','p5','q1','q2','q3','q4'};
coeffNames = {'p1','p2','q1'};
% myfun   = fittype('(p1*x^4 + p2*x^3 + p3*x^2 + p4*x + p5)/(x^4 + q1*x^3 + q2*x^2 + q3*x + q4)','independent','x','coefficients',coeffNames);
myfun   = fittype('(p1*x + p2)/(x + q1)','independent','x','coefficients',coeffNames);
options = fitoptions(...
   'method','NonLinearLeastSquares',...
   'StartPoint',[1 -3 15],...
   'MaxFunEvals',5000,...
   'TolFun',1e-07,...
   'TolX',1e-07,...
   'Lower',[-Inf -Inf],...
   'Upper',[+Inf +Inf]);
[fitcdf,gof] = fit(x(t),y(t),myfun,options);
fig.fit = plot(x(t),fitcdf(x(t)));
legend('Histogram','cdf','fit of cdf');
xlabel('time (ms)');
ylabel('probability');

% plot inverse of cdf, which is the distribution of the time intervals
syms u
%f(u)    = (fit.p1*u^4 + fit.p2*u^3 + fit.p3*u^2 + fit.p4*u + fit.p5)/(u^4 + fit.q1*u^3 + fit.q2*u^2 + fit.q3*u + fit.q4);
f(u)     = (fitcdf.p1*u + fitcdf.p2)/(u + fitcdf.q1);
distr    = finverse(f);
% 
% figure('Name','Inverse of cdf');
% z = 0:0.01:1;
% fig.invcdf = plot(z,distr(z));
% xlabel('probability');
% ylabel('time interval (ms)');

% create spike trains
durtot     = 128000;           % ms
Nbins   = floor(durtot/bins);  % #bins
X       = zeros(1,Nbins);
Y       = zeros(Nrep,Nbins) > 1;
dur     = 0;

while dur <= 128000
    var     = rand;
    dur     = dur + distr(var);
    dbin    = floor(dur/3);
    if dbin < Nbins
        X(1,dbin)  = 1;
    end
end

Ydel    = rand(Nrep,Nbins) < (1-pfail);
Yadd    = rand(Nrep,Nbins) < pspont;

for i = 1:1:Nrep
    Y(i,:) = X.*Ydel(i,:)+Yadd(i,:);
end
end