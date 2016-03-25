%% Low pass Filter

scanRate = GRABinfo.scanFrameRate;

Nqst = scanRate/2;


% bands is a vector of pairs of normalized frequency points, specified in the
% range between 0 and 1, where 1 corresponds to the Nyquist frequency. The 
%frequencies must be in increasing order.

%
%cutoff_freq = 8*scan

bands = [0,0.6,0.65]%/scanRate;
bands(end+1) = 1;
ampl = [1,1,0,0];


low_pass = firpm(512,bands,ampl);

%%

trx = traces(689,:);


dwn_spld = upfirdn(trx,low_pass,1,1);

plot(trx,'bo-')
hold on;
plot(dwn_spld(257:end-256),'r','LineWidth',3);

legend('Original','Low Pass filtered')

%%
freqz(low_pass);
