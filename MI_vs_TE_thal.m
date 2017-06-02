function MI_vs_TE_thal(MIcalc,TEcalc)

% Last Massive change: 2.06.2017

% Begin global timer and timer for MI calculation
now1 = tic();

% Parameters
rep         = [2 5 10 15 20 30 40 50 75 100 150 200 300 400 500]; % Vector containing the number of repetitions
%rep = [2 3 10];
Nmax        = max(rep);
bins        = 3;                    % ms
dur         = 128000;               % ms
Nbins       = floor(dur/bins);      % #bins
pfail       = 0.8;                  % Probability that a spike fails being transmitted
pspont      = 4*0.07*3/1000;      % Probability that a spike is spontaneously generated (From non-spike or spike not transmitted)
Nmean       = 10;                   % Number of calculation for MI with the Strong method (MI is then the mean over all the calculations)
%Nmean = 2;
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

%% Create Thalamic like Input and output spike train for size of maximum repetitions number

[YMI,XTE,YTE]=GenerateThalamicSpikeTrains(dur,Nbins,pfail,pspont,Nmax);


%% Mutual Information Calculation

if MIcalc == 1

%     % Calculation of Theoretical value of MI
%     [~, ~, MItheo1, MItheo2] = MITheory(bins,freq,pfail,pspont);
%     MItheo1tot(:,:)          = MItheo1;
%     %MItheo2tot(:,:)          = MItheo2;

% Creat input spike train for MI calculation

    % Calculation of MI with the Strong method
    for j = [1:Nmean]
    
        fprintf(['\n calculations for the ' num2str(j) ' time ...']); 
        MI      = zeros(1,length(rep));
        x       = 0;

        for i = rep
            
            YMIshaped     = YMI(1:i,:);
            fprintf(['\n calculations for Nrep = ' num2str(i) ' ...']);
            x           = x+1;
            [MI(x),~,~] = MutualInformation_thal(0,bins, Nbins,YMIshaped);

        end

    MItot(j,:)        = MI;

    end

    MImean      = mean(MItot);
    %MItheo1mean = mean(MItheo1tot);
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
    fprintf(['\n calculations for Nrep = ' num2str(k) ' ...']);
    
    XTEshaped       = XTE(:,1:Nbins*k);
    YTEshaped       = YTE(:,1:Nbins*k);
    y = y+1;
        for l = words
        fprintf(['\n \t calculations for words of length ' num2str(l) '...']);
            [TEin_out_mean(l,y), TEout_in_mean(l,y), STDin_out(l,y), STDout_in(l,y)] = TransferEntropy_thal(bins,l,XTEshaped,YTEshaped);
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
        %h.mitheo1   = plot(rep, MItheo1mean);
        %h.mitheo2   = plot(rep, MItheo2mean);
        for l = words
            h.teio      = errorbar(rep, TEin_out_mean(l,:), STDin_out(l,:));
            h.teoi      = errorbar(rep, TEout_in_mean(l,:), STDout_in(l,:));
        end
        legend('MI','TEin->out','TEout->in');
    else
        h.mi        = errorbar(rep, MImean, std(MItot));
        %h.mitheo1   = plot(rep, MItheo1mean);
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
        save(['Results/' saveName '.mat'],'rep','MImean','TEin_out_mean','TEout_in_mean',...
        'MITime','TETime','TotalTime')
    else
        save(['Results/' saveName '.mat'],'rep','MImean','MITime','TotalTime')
    end
else
    if TEcalc ==1
        save(['Results/' saveName '.mat'],'rep','TEin_out_mean','TEout_in_mean','TETime','TotalTime')
    end
end

fprintf(['\n Done! \n \n']);

end