function [YMI,XTE,YTE]=GenerateThalamicSpikeTrains(durtot,Nbins,pfail,pspont,Nmax)

load('spike_trains.mat');

now = tic();

% read spike trains
A       = cell2mat(spike_trains(1,1));
a       = length(A);
B       = cell2mat(spike_trains(1,2));
b       = length(B);
C       = cell2mat(spike_trains(1,3));
c       = length(C);
D       = cell2mat(spike_trains(1,4));
durMI       = length(D);
E       = cell2mat(spike_trains(1,5));
e       = length(E);
V       = zeros(a+b+c+durMI+e,1);

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
for i= 1:durMI-1
    V(a+b+c+i+1,1) = D(i+1,1)-D(i,1);
end

V(a+b+c+durMI+1,1) = E(1,1);
for i= 1:e-1
    V(a+b+c+durMI+i+1,1) = E(i+1,1)-E(i,1);
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
coeffNames = {'a','b','n'};
myfun1   = fittype('a*x/sqrt(b+x^n)','independent','x','coefficients',coeffNames);
options = fitoptions(...
   'method','NonLinearLeastSquares',...
   'StartPoint',[1,560,2],...
   'MaxFunEvals',5000,...
   'TolFun',1e-07,...
   'TolX',1e-07,...
   'Lower',[-Inf -Inf],...
   'Upper',[+Inf +Inf]);
[fitcdf,gof] = fit(x(t),y(t),myfun1,options);
fig.fit = plot(x(t),fitcdf(x(t)));
legend('Histogram','cdf','fit of cdf');
xlabel('time (ms)');
ylabel('probability');

% plot inverse of cdf, which is the distribution of the time intervals
figure('Name','Inverse of cdf');
fig.a = axes; 
hold(fig.a,'all');
fig.invcdf = plot(fitcdf(x(t)),x(t));
xlabel('probability');
ylabel('time interval (ms)');


% fit inverse of cdf
invx    = transpose(get(fig.invcdf,'XData'));
invy    = transpose(get(fig.invcdf,'YData'));
coeffNames = {'c','d','e','f'};
myfun2   = fittype('c*exp(d*x)+e*exp(f*x)','independent','x','coefficients',coeffNames);
options = fitoptions(...
   'method','NonLinearLeastSquares',...
   'StartPoint',[22,1,0,6],...
   'MaxFunEvals',5000,...
   'TolFun',1e-07,...
   'TolX',1e-07,...
   'Lower',[-Inf -Inf],...
   'Upper',[+Inf +Inf]);
[fitinvcdf,gof] = fit(invx,invy,myfun2,options);
fig.invfit = plot(invx,fitinvcdf(invx),'r');
legend('inverse of cdf','fit of inverse of cdf');
xlabel('probability');
ylabel('time interval (ms)');


% create spike trains for MI calculation

XMI      = zeros(1,Nbins);
YMI      = zeros(Nmax,Nbins) > 1;
durMI     = 0;

while durMI <= durtot
    var     = rand;
    durMI     = durMI + fitinvcdf(var);
    dbin    = floor(durMI/3);
    if dbin < Nbins
        XMI(1,dbin)  = 1;
    end
end

XMI       = logical(XMI);
YdelMI    = rand(Nmax,Nbins) < (1-pfail);
YaddMI    = rand(Nmax,Nbins) < pspont;

for i = 1:1:Nmax
    YMI(i,:) = XMI.*YdelMI(i,:)+YaddMI(i,:);
end

% create spike trains for TE calculation

XTE      = zeros(1,Nbins*Nmax);
YTE      = zeros(1,Nbins*Nmax) > 1;
durTE     = 0;

while durTE <= durtot*Nmax
    var     = rand;
    durTE     = durTE + fitinvcdf(var);
    dbin    = floor(durTE/3);
    if dbin < Nbins*Nmax
        XTE(1,dbin)  = 1;
    end
end

XTE     = logical(XTE);
YdelTE    = rand(1,Nbins*Nmax) < (1-pfail);
YaddTE    = rand(1,Nbins*Nmax) < pspont;
YTE       = XTE.*YdelTE+YaddTE;



totaltime=toc(now);
end