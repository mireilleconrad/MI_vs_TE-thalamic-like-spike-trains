function [TEin_out_mean, TEout_in_mean, STDin_out, STDout_in] = TransferEntropy_thal(bins,wordlength,X,Y)

%parameter
nMean   = 50;                   % Number of calculation of TE before taking the mean

%preallocating variables
TEin_out     = zeros(1,nMean);         % Transfer entropy from input to output
TEout_in     = zeros(1,nMean);         % Transfer entropy from output to input

    
for j = 1:nMean

    %Calculate tranfer entropy
    trains                  = [X;Y];
    asdf                    = SparseToASDF(trains,1);
    [peakTE, CI, TEdelay]   = ASDFTE(asdf,0:1,wordlength,wordlength);
    %peakTE                  = peakTE/(3e-03);
    TEin_out(j)             = TEdelay(2,1,1)/(bins/1000);
    TEout_in(j)             = TEdelay(1,2,1)/(bins/1000);
end
    
TEin_out_mean   = mean(TEin_out);
STDin_out       = std(TEin_out);
TEout_in_mean   = mean(TEout_in);
STDout_in       = std(TEout_in);

end