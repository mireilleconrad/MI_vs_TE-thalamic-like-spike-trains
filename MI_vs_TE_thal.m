function MI_vs_TE_thal(MIcalc,TEcalc)

% Last Massive change: 26.05.2017

% Begin global timer and timer for MI calculation
now1 = tic();

% Parameters
rep         = [2 5 10 15 20 30 40 50 75 100 150 200 300 400 500]; % Vector containing the number of repetitions
bins        = 3;                    % ms
dur         = 128000;               % ms
Nbins       = floor(dur/bins);      % #bins
pfail       = 0.9;                  % Probability that a spike fails being transmitted
pspont      = 0.8*3/1000;      % Probability that a spike is spontaneously generated (From non-spike or spike not transmitted)
%Nmean       = 10;                   % Number of calculation for MI with the Strong method (MI is then the mean over all the calculations)
Nmean = 2;
words       = 7;               % word lengths to compute

saveName = [datestr(datetime('now'),'dd-mm-yyyy-HH:MM') ')_pfail(' num2str(pfail)... 
    ')_pspont(' num2str(pspont) ')_rep(' num2str(min(rep)) '-' num2str(max(rep)) ')_dur(' num2str(dur) ')'];

% Preallocating variables
MItot           = zeros(Nmean,length(rep));
MItheo1tot      = zeros(Nmean,length(rep));
%MItheo2tot      = zeros(Nmean,length(rep));
TEin_out_mean   = zeros(length(words),length(rep));
TEout_in_mean   = zeros(length(words),length(rep));
STDin_out       = zeros(length(words),length(rep));
STDout_in       = zeros(length(words),length(rep));

%% Create Thalamic like Input spike train

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


%% Mutual Information Calculation

if MIcalc == 1

%     % Calculation of Theoretical value of MI
%     [~, ~, MItheo1, MItheo2] = MITheory(bins,freq,pfail,pspont);
%     MItheo1tot(:,:)          = MItheo1;
%     %MItheo2tot(:,:)          = MItheo2;

% Creat input spike train for MI calculation
X       = zeros(1,Nbins);
dur     = 0;

while dur <= 128000
    var     = rand;
    dur     = dur + distr(var);
    dbin    = floor(dur/3);
    if dbin < Nbins
        X(1,dbin)  = 1;
    end
end

    % Calculation of MI with the Strong method
    for j = [1:Nmean]
    
        fprintf(['\n calculations for the ' num2str(j) ' time ...']); 
        MI      = zeros(1,length(rep));
        x       = 0;

        for i = rep
  
            fprintf(['\n calculations for Nrep = ' num2str(i) ' ...']);
            x = x+1;
            [MI(x),~,~] = MutualInformation_thal(0, i, bins, Nbins, pfail, pspont,X);

        end

    MItot(j,:)        = MI;

    end

    MImean      = mean(MItot);
    MItheo1mean = mean(MItheo1tot);
    %MItheo2mean = mean(MItheo2tot);

    % End timer for MI calculation
    MITime = toc(now1);

end

%% TE calculation
if TEcalc == 1

    % Begin timer for TE calculation
    now2 = tic();

    % Calculation of Tranfer Entropy
    fprintf(['\n calculations of Transfer Entropy ... \n']);
    y = 0;
    
    for k = rep
        
    % Creat input spike train for TE calculation
    X       = zeros(1,Nbins*k);
    dur     = 0;

    while dur <= 128000*k
        var     = rand;
        dur     = dur + distr(var);
        dbin    = floor(dur/3);
        if dbin < Nbins
            X(1,dbin)  = 1;
        end
    end    
    fprintf(['\n calculations for Nrep = ' num2str(k) ' ...']);
    y = y+1;
        for l = words
        fprintf(['\n \t calculations for words of length ' num2str(l) '...']);
            [TEin_out_mean(l,y), TEout_in_mean(l,y), STDin_out(l,y), STDout_in(l,y)] = TransferEntropy_thal(k,Nbins,bins,pfail,pspont,l,X);
        end
    end

    % End timer for TE calculation
    TETime = toc(now2);
    
end

%% Code end

%end global timer
TotalTime = toc(now1);

% plot
figure('Name','Information');
h.a = axes;
hold(h.a,'all');
if MIcalc == 1
    if TEcalc == 1
        h.mi        = errorbar(rep, MImean, std(MItot));
        h.mitheo1   = plot(rep, MItheo1mean);
        %h.mitheo2   = plot(rep, MItheo2mean);
        for l = words
            h.teio      = errorbar(rep, TEin_out_mean(l,:), STDin_out(l,:));
            h.teoi      = errorbar(rep, TEout_in_mean(l,:), STDout_in(l,:));
        end
        legend('MI', 'MItheo1','TEin->out','TEout->in');
    else
        h.mi        = errorbar(rep, MImean, std(MItot));
        h.mitheo1   = plot(rep, MItheo1mean);
        %h.mitheo2   = plot(rep, MItheo2mean);
        legend('MI', 'MItheo1');
    end
else
    if TEcalc == 1
         for l = words
            h.teio      = errorbar(rep, TEin_out_mean(l,:), STDin_out(l,:));
            h.teoi      = errorbar(rep, TEout_in_mean(l,:), STDout_in(l,:));
        end
        legend('TEin->out','TEout->in');
    end
end
xlabel('Number of repetitions');
ylabel('[bit/sec]');

%Save figure
saveas(gcf, ['Results/' saveName '.fig']);

%Save values
if MIcalc == 1
    if TEcalc ==1
        save(['Results/' saveName '.mat'],'rep','MImean','MItheo1mean','TEin_out_mean','TEout_in_mean',...
        'MITime','TETime','TotalTime')
    else
        save(['Results/' saveName '.mat'],'rep','MImean','MItheo1mean','MITime','TotalTime')
    end
else
    if TEcalc ==1
        save(['Results/' saveName '.mat'],'rep','TEin_out_mean','TEout_in_mean','TETime','TotalTime')
    end
end

fprintf(['\n Done! \n \n']);

end