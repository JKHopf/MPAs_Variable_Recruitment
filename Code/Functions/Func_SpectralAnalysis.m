function [period, AS, PS, tsDTstd] = Func_SpectralAnalysis(tsdata, tstime, tsselect, selects)
figure
hold on

tsN = NaN(1,length(selects));
lag1 = NaN(1,length(selects));
    
for s = 1:length(selects)
    
    % create time series object
    if length(selects) > 1
    ts = timeseries(tsdata(tsselect == selects(s)),...
               tstime(tsselect == selects(s)));
    else
        ts = timeseries(tsdata, tstime);
    end
    
    % scale to unit variance (i.e. normalise)
    tsStd = ts/std(ts);  
      
    % detrend
        % fit linear model
        Model = fitglm(tsStd.Time, tsStd.Data);
        
        % mean plus the redsiduals 
        % (difference between the observed and values predicted by the fitted model)
        tsDT = mean(tsStd.Data) + tsStd.Data - predict(Model,tsStd.Time);
        
        % scale to unit variance 
        tsDTstd = tsDT/std(tsDT);
    
    % find how much variance occurs at specific frequencies  
        % frequencies (see matlab help doc fft)
        freqi = linspace(0,0.5,length(tsStd.Data)/2+1);
        % get length of fft output
        tsN(s) = length(freqi);
        
        % get discrete fourier transform (DFT) of the standardised-detrended data
        % using the Fast Fourier Transformation (fft)
        fhat = fft(tsDTstd);
        
        % normalise the DFT output (as calculated by the fft)
            % this standarises the fft output for different lengthed samples
            % (when analysis amplitudes)
            % see https://dsp.stackexchange.com/questions/31984/matlab-tt-fft-and-tt-ifft-scaling
            % and Ahrens et al (2020) (https://appliedacousticschalmers.github.io/scaling-of-the-dft/AES2020_eBrief/)
%         fhat = fhat/length(tsStd.Data);
        
        % the amplitude spectrum (Ahrens et al 2020)
        % this give the amplitude of each component, scaled by the length
        % of the time series
        amp_spec = abs(fhat)/length(tsStd.Data);
        
        % multiple by two to get single sided spectrum
        AS = 2.*amp_spec(1:length(freqi));
                
    % get periods (and inc only those up to half the length of time series)
    period = 1./freqi; 
    inc = true(1,length(period));% period <= tsN(s);
    period = period(inc);
    
    % figure out single-sided power spectrum
    PS = 2.*amp_spec(inc).^2;
    
    % plot variance (single-sided power spectrum) vs. period (1/freq index)
        % have only up to periods half the length of time series
    line(period, PS, 'Color','k')
    
    % Calcluate the lag-1 autocorrelation (used to calc red noise sig
    % value)
    ac = autocorr(tsStd.Data);
    lag1(s) = ac(2);
    
end

% Red noise sig ----
% Plot red noise significance (based on Torrence & Compo 1998). General
% priciple is to build a background spectrum (for which different 
% realisations of the process will be randomly distributed around) and then 
% compare the observed specturm to the background one; if a peak in the 
% observed spectrum is significantly above the background spectrum, then it
% can be assumed to be a true feature, with a set percentage confidence 
% (i.e. 95%). I.e. we are testing against the null hyposthesis that the 
% observed spectrum stems from a red-noise process, with a 5% significance level.

% first determine the lag value that corresponds to the assumed AR(1)
% (lag-1 autocorrelation value)
a = mean(abs(lag1)); % 0.43; % 0.72; % 

% get the mean number of fft output point as the sample size
tsNmean = mean(tsN);

% define background spectrum (in this case we are using the normalised 
% discrete Fourier power spectrum from the AR(1) [Markov] process, 
% if a = 0 this is white noise)

BGspec = (1-a^2)./(1-2.*a.*cos(2.*pi.*freqi)+a^2);
BGspec = BGspec./(2*tsNmean);

% get chi2 value, with 2 dregrees of freedom
% i.e. If you generate random numbers from this chi-square distribution,
% you would observe numbers greater than chi22 only 5% of the time.
chi22 = chi2inv(0.95,2);

% get significance level for specturm
sigif = chi22*BGspec/2;

% add to figure (point above the red line are significant)
sigifx = 1./linspace(0,0.5,max(tsN));
line(sigifx(inc), sigif(inc),'Color','r')

end

