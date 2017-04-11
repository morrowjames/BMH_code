            

        cut_up_freq=3.37;  %in ppm
        cut_down_freq=3.19;     %in ppm
        
        diff_freq=freq(2)-freq(1);
        cut_low = (cut_up_freq-freq(min(freqbounds)))/diff_freq+1;
        cut_high = (cut_down_freq-freq(min(freqbounds)))/diff_freq+1;
%         cut_low = 350;
%         cut_high = 580;
        freqbounds_chop=[freqbounds(1:cut_low) freqbounds(cut_high:end)];
        freqbounds_orig=freqbounds;
        freqbounds=freqbounds_chop;
        plotbounds=(lowerbound-150):(upperbound+150);
        plotbounds_chop=[plotbounds(1:(cut_low+150)) plotbounds((cut_high+150):end)];
        plotbounds_orig=plotbounds;
        plotbounds=plotbounds_chop;

        plot(freq(freqbounds),GaussModel_area(GaussModelParam(ii,:),freq(freqbounds)),'r',...
                freq(plotbounds),real(GABAData(ii,plotbounds)), 'b', ...
                freq(freqbounds),residg,'k',...
                cut_up_freq,GABAData(ii,freqbounds_orig(round(cut_low))),'co',...
                cut_down_freq,GABAData(ii,freqbounds_orig(round(cut_high))),'co');
            hold on;
            line ([cut_down_freq cut_down_freq],[YL(1) real(GABAData(ii,freqbounds_orig(round(cut_high))))],...
                'Color', 'c', 'LineStyle',':')
            line ([cut_up_freq cut_up_freq],[YL(1) real(GABAData(ii,freqbounds_orig(round(cut_low))))],...
                'Color', 'c', 'LineStyle',':')