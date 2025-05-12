function cohr = fn_getCrossFieldCoherence(x,y,t)

K = size(x,1); %Define the number of trials.

N = size(x,2); %Define the number of indices per trial.

dt = t(2)-t(1); %Define the sampling interval.

T = t(end); %Define the duration of data.

Sxx = zeros(K,N); %Create variables to save the spectra.
Syy = zeros(K,N);
Sxy = zeros(K,N);

for k=1:K %Compute the spectra for each trial.
  
  Sxx(k,:) = 2*dt^2/T * fft(x(k,:)) .* conj(fft(x(k,:)));
  Syy(k,:) = 2*dt^2/T * fft(y(k,:)) .* conj(fft(y(k,:)));
  Sxy(k,:) = 2*dt^2/T * fft(x(k,:)) .* conj(fft(y(k,:)));
  
end

Sxx = Sxx(:,1:N/2+1); %Ignore negative frequencies.
Syy = Syy(:,1:N/2+1);
Sxy = Sxy(:,1:N/2+1);

Sxx = mean(Sxx,1); %Average the spectra across trials.
Syy = mean(Syy,1);
Sxy = mean(Sxy,1);

cohr = abs(Sxy) ./ (sqrt(Sxx) .* sqrt(Syy));

end

% %Compute the coherence.
% df = 1/max(T); %Determine the frequency resolution.
% fNQ = 1/ dt / 2; %Determine the Nyquist frequency.
%
% faxis = (0:df:fNQ); %Construct frequency axis.
% plot(faxis, real(cohr)); %Plot the results
% xlim([0 50]); ylim([0 1]) %Set the axes limits
% xlabel('Frequency [Hz]') %Label axes.
% ylabel('Coherence [ ]')