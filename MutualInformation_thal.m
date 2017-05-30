function [MI,Entropy,NoiseEntropy] = MutualInformation_thal(FirstCorr, Nrep, bins, Nbins, pfail, pspont,X)
%% Parameters

% Last massive change: 26.05.2017

% Sim         = 0;               % Determine if the results of the simulation are used (1)
% RandomData  = 1;               % Determine if the random generate data are used (1)
% FirstCorr   = 0;               % Determine if the First Correction is applied (1) or not (0)

%bins    = 3;                % ms
%freq    = 20;               % Hz
fracs   = 1:1:5;            % fractions of the dataset to test in the first correction
Nsmpl   = 30;               % #samples for averaging
Words   = [1:10];           % word lengths to compute
%pfail   = 0.7;              % Probability that a spike fails being transmitted
%pspont  = 4*0.07/(1000*3);  % Probability that a spike is spontaneously generated (From non-spike or spike not transmitted)




%% Calculation with Thalamic like data

Y       = zeros(Nrep,Nbins) > 1;

Ydel    = rand(Nrep,Nbins) < (1-pfail);
Yadd    = rand(Nrep,Nbins) < pspont;

for i = 1:1:Nrep
    Y(i,:) = X.*Ydel(i,:)+Yadd(i,:);
end


%% Corrections for the Entropy

if FirstCorr == 1
    %First Correction
    fprintf(['\n calculations of the 1st correction ...']);
    [HTotCorr,HNoiseCorr]=FirstCorrection(fracs,Nbins,Y,Words,Nsmpl);
else
    %Entropy Calculation without the first correction
    frequencyCount          = Ptable(Y,Words);
    [HtotalRaw, HtotalRate] = Htot(Y,Words,frequencyCount,bins);
    [HnoiseRaw, HnoiseRate] = Hnoise(Words, Y,bins);
    HTotCorr(:,1)                = HtotalRate;
    HNoiseCorr(:,1)              = HnoiseRate;
end

%Second Correction
% figure('Name','2nd correction');
% fig.a = axes; 
% hold(fig.a,'all');
% fig.htc = plot(1./Words,HTotCorr,'-','LineWidth',2,'Color','r','Parent',fig.a);
% fig.hnc = plot(1./Words,HNoiseCorr,'-','LineWidth',2,'Color','b','Parent',fig.a);
% legend('Htot','Hnoise');
%set(fig.htc,'MarkerFaceColor','r');
%set(fig.hnc,'MarkerFaceColor','b');
% axis([0 max(1./Words)+1 0 1.2*max(HTotCorr)]);
% xticks([0 sort(1./Words)]);
WordsStr = strsplit(rats(sort(1./Words))); WordsStr(length(WordsStr)) = []; WordsStr(1) = [];
% xticklabels([0 WordsStr]);
% xlabel('Inverse of Word lengths');
% ylabel('Entropy [bits/sec]');

% Fit line through total entropy
coeffNames      = {'a','b'};
myfun   = fittype(...
   'a*x+b',...
   'independent','x',...
   'coefficients',coeffNames);
options = fitoptions(...
   'method','NonLinearLeastSquares',...
   'StartPoint',[0.1 22],...
   'MaxFunEvals',5000,...
   'TolFun',1e-07,...
   'TolX',1e-07,...
   'Lower',[-Inf -Inf],...
   'Upper',[+Inf +Inf]);
[ht.fun,gof] = fit(1./Words',HTotCorr,myfun,options);
% fig.htp = plot([0 1./Words],ht.fun([0 1./Words]),'-k','LineWidth',1);

% Fit line through noise entropy
options = fitoptions(...
    'method','NonLinearLeastSquares',...
    'StartPoint',[0.1 1],...
    'MaxFunEvals',5000,...
    'TolFun',1e-07,...
    'TolX',1e-07,...
    'Lower',[-Inf -Inf],...
    'Upper',[+Inf +Inf]);
[hn.fun,gof] = fit(1./Words',HNoiseCorr,myfun,options);
% fig.hnp = plot([0 1./Words],hn.fun([0 1./Words]),'-k','LineWidth',1);

fprintf('\n');

% Mutual Information Calculation
Entropy         = ht.fun(0);
NoiseEntropy    = hn.fun(0);
MI              = Entropy - NoiseEntropy;

end

%% Subfonctions
% calculate entropy tables
function frequencyCount = Ptable(Data,Words)
%
[Nrep, Nbins] = size(Data);
for i = Words
    for k =1:1:Nrep
        frequencyCount{i,k}      = zeros(2^i,1);
        for j = 1:Nbins-i+1
            word                         = myfasterbin2dec(Data(k,j:(j+i-1)));
            frequencyCount{i,k}(word+1)  = frequencyCount{i,k}(word+1)+1;
        end
        frequencyCount{i,k}      = frequencyCount{i,k}/sum(frequencyCount{i,k});
    end
end
end

% calculate total entropy
function [HtotalRaw, HtotalRate]= Htot(Data,Words,frequencyCount,bins)
%
[Nrep, Nbins] = size(Data);
fprintf(['\n calculations for total entropy ...']);
for i = Words
    for k = 1:1:Nrep
        prob            = frequencyCount{i,k};
        prob(prob == 0) = [];
        prob(prob == 1) = [];
        if isempty(prob)
            HRaw(k)   = 0;
            HRate(k)  = 0;
        else
            HRaw(k)   = -(prob'*log2(prob));
            HRate(k)  = HRaw(k)/(i*bins/1000);
        end
    end
HtotalRaw(i)	= mean(HRaw);
HtotalRate(i)	= mean(HRate);
end
end

% calculate the entropy tables for noise entropy
function [HnoiseRaw, HnoiseRate] = Hnoise(Words, Data, bins)
[Nrep, Nbins] = size(Data);
 fprintf(['\n calculations for noise entropy ...']);
for i = Words
    fprintf(['\n \t calculations for words of length ' num2str(i) '...']);
    for j = 1:Nbins-i+1
        frequencyCount = zeros(2^i,1);
        for k = 1:1:Nrep
                word                        = myfasterbin2dec(Data(k,j:(j+i-1)));
                frequencyCount(word+1)      = frequencyCount(word+1)+1;
        end
        prob = frequencyCount/sum(frequencyCount);            
        prob(prob == 0) = [];
        prob(prob == 1) = [];
        if isempty(prob)
            HRaw(j)   = 0;
            HRate(j)  = 0;
        else
            HRaw(j)   = -(prob'*log2(prob));
            HRate(j)  = HRaw(j)/(i*bins/1000);
        end
    end
    HnoiseRaw(i)   = mean(HRaw);
    HnoiseRate(i)  = mean(HRate);                
end
end