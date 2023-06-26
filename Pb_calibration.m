function Pb_calibration
clear all
close all

% Use this single wavelength to plot signal variation as spectra are acquired.
wavelength_variance_display = 515;
% Samples will be collected in the order shown here.  "blank" is a special term for computing absorbance.
% All "blank" acquisitions will be averaged for calculating absorbance of other samples, and
% variance across all "blank" acqusitions will be plotted to show spectrometer drift over the series of acqusitions.
sample_names = {"blank", "2000 ug/l", "750 ug/l", "500 ug/l", "250 ug/l", "100 ug/l", "25 ug/l", "0 ug/l", "blank"};


%%%%%%%%%%%%%%%%% Internal state variables
files_read = 1;
files_total = 0;
current_acquisition = 0;
scans_collected = 0;
currently_collecting = 0;
done_collecting = 0;
data_lengths = 0;

fig1 = figure(1);
hax_live_data = axes('Parent',fig1);

fig2 = figure(2);
hax_live_variance = axes('Parent',fig2);

disp(["Press n to start acquiring " ,num2str( sample_names{current_acquisition+1})]);

%%%%%%%%%%%%%%%%%%%%%    Realtime data collection loop
while(done_collecting == 0)
  S = dir("sample_*.txt");
  [~,idx] = sort([S.datenum]);
  S = S(idx);
  files_total = length(S);
  if( files_read < (files_total - 1 ))
    fid = fopen(S(files_total - 1).name,'r','l');
    newdata = cell2mat(textscan(fid,"%f %f"));
    fclose(fid);
    files_read = files_total;
    plot(hax_live_data,newdata(:,1),newdata(:,2));
    title(hax_live_data,"Live data");
    ylabel(hax_live_data,"Counts");
    xlabel(hax_live_data,"Wavelength (nm)");
    axis(hax_live_data, [400 700]);

    switch(kbhit(1))
      case "n"
        if (currently_collecting == 0)
          scans_collected = 0;
          current_acquisition++;
          disp("Acquiring now... press x to stop");
          currently_collecting = 1;
        endif

      case "x"
        if (currently_collecting == 1)
          currently_collecting = 0;
          data_lengths(1,current_acquisition) = scans_collected;
          disp("Stopped.");
          if( current_acquisition == length(sample_names))
            done_collecting = 1;
          else
            disp(["Press n to start acquiring " ,num2str( sample_names{current_acquisition+1})]);
          endif
        endif
    endswitch


    if (currently_collecting == 1)

      raw_spectra(current_acquisition,++scans_collected,:) = newdata(:,2);
      raw_wavelengths(:,1) = newdata(:,1);
    endif

    if(scans_collected >= 2)
      [~,wl_index] = min(abs(raw_wavelengths(:,1)-wavelength_variance_display));
      clear tempdata
      for (i = 1:scans_collected)
        tempdata(1,i) = raw_spectra(current_acquisition,i,wl_index);
      end
     plot(hax_live_variance,1:scans_collected,tempdata(1,:) , [1 scans_collected],[mean(tempdata(1,:)) mean(tempdata(1,:))]    );
     title(hax_live_variance,["Variance at ",num2str(wavelength_variance_display), " nm"] );
     ylabel(hax_live_variance,"Counts");
     xlabel(hax_live_variance,"Samples");

    endif

    pause(0.05);
   endif

   pause(0.05);
  endwhile



%%%%%%%%%%%%%%%%%%%%%    Average collected spectra for each sample
for (i = 1 : length(sample_names))
  %TODO if there was only one spectrum acquired for a sample, the mean function collapses the dimension we want to preserve.
  # Need to add a test if there is only one sample, or figure out a parameter to add to the mean function to make sure it preserves the important dimension
  avg_spectra(i,:) = mean(squeeze(raw_spectra(i,1:data_lengths(i),:)));
endfor

%%%%%%%%%%%%%%%%%%%%%    Compute blank average and plot variance for blank acqusitions
fig3 = figure(3);
hax_baseline_drifts = axes('Parent',fig3);
title(hax_baseline_drifts,"Blank drift");
ylabel(hax_baseline_drifts,"Counts");
xlabel(hax_baseline_drifts,"Wavelength (nm)");
axis(hax_baseline_drifts, [400 700]);
hold on

blank_sum = 0;
blank_index = find(strcmp(sample_names, "blank"));
for (i = 1 : (length(blank_index)))
  plot(hax_baseline_drifts,raw_wavelengths(:,1),(avg_spectra(blank_index(1),:) - avg_spectra(blank_index(i),:))');
  blank_sum = blank_sum + avg_spectra(blank_index(i),:)';
endfor
blankavg = blank_sum ./ length(blank_index);
legend
hold off



%%%%%%%%%%%%%%%%%%%%%%    Compute and plot absorbance
fig4 = figure(4);
hax_abs = axes('Parent', fig4);
title(hax_abs, "Samples");
ylabel(hax_abs,"Absorbance");
xlabel(hax_abs,"Wavelength (nm)");
hold on
for (i = 1 : (length(sample_names)))
  %temp_startindex = find(
  absorbances(:,i) = log10(10 ./ (avg_spectra(i,:)' ./ blankavg(:,1)));
  plot(raw_wavelengths(:,1), absorbances(:,i))
endfor
legend(sample_names{1:end})
hold off


%%%%%%%%%%%%%%%%%%%%%%%  Compute and plot absorbance normalized to absorbance measured 750 to 800
fig5 = figure(5);
hax_normabs = axes('Parent', fig5);
title(hax_normabs, "Normalized Absorbance");
ylabel(hax_normabs, "Absorbance");
xlabel(hax_normabs, "Wavelength (nm)");
hold on
[~,start_norm_index] = min(abs(raw_wavelengths(:,1)-750));
[~,end_norm_index] = min(abs(raw_wavelengths(:,1)-800));

for (i = 1 : (length(sample_names)))

  norm_average = mean(absorbances([start_norm_index end_norm_index],i));
  norm_absorbances(:,i) = absorbances(:,i) - norm_average + 1;
  plot(raw_wavelengths(:,1), norm_absorbances(:,i))
endfor
legend(sample_names{1:end})
hold off



%%%%%%%%%%%  Compute  calibration curve for single measurement at 515nm
fig6 = figure(6);
hax_calcurve = axes('Parent', fig6);
title(hax_calcurve, "Calibration curve");
ylabel(hax_calcurve, "Absorbance");
xlabel(hax_calcurve, "Concentration (ug/l)");

[~, pbtarget_index] = min(abs(raw_wavelengths(:,1)-515));

sample_index_no_blanks = find(~strcmp(sample_names, "blank"));
for (i = 1 : (length(sample_index_no_blanks)))
  calcurve(i,1) = sscanf(cell2mat(sample_names(sample_index_no_blanks(i))),'%f');
  calcurve(i,2) = norm_absorbances(pbtarget_index,sample_index_no_blanks(i));
endfor
hold on
plot(calcurve(:,1), calcurve(:,2),'*')
P =  polyfit(calcurve(:,1),calcurve(:,2),3);
concentrations = 0:1000;
absfit = P(1)*concentrations.^3 + P(2)*concentrations.^2 + P(3)*concentrations + P(4);
plot(concentrations, absfit)
hold off

save calsession
save caldata concentrations absfit P calcurve


%%%%%%%%%%%%%%%  Compute calibration curve for whole-spectrum fit

[~,startfit] = min(abs(raw_wavelengths(:,1)-400));
[~,endfit] = min(abs(raw_wavelengths(:,1)-800));

function [a, b] = optimize_ab(X, Y, S)
res = fmincon(@MSE, [0,0], [], [], [], [], [0,0], []);  %Second argument is the starting point, second to the last argument is the lower bound to ensure the solutions are all positive
a = res(1);
b = res(2);
    function mse = MSE(x_)
        S_ = x_(1)*X + x_(2)*Y;
        mse = norm(S_-S);
    end
end


sample_index_no_blanks = find(~strcmp(sample_names, "blank"));
minindex = sample_index_no_blanks(1);
maxindex = sample_index_no_blanks(length(sample_index_no_blanks));

for (i = 1 : (length(sample_index_no_blanks)))
  calfitzero(i,1) = sscanf(cell2mat(sample_names(sample_index_no_blanks(i))),'%f');
  calfitmax(i,1) = sscanf(cell2mat(sample_names(sample_index_no_blanks(i))),'%f');
  [calfitzero(i,2), calfitmax(i,2)]  = optimize_ab(norm_absorbances(startfit:endfit,minindex), norm_absorbances(startfit:endfit,maxindex), norm_absorbances(startfit:endfit,sample_index_no_blanks(i)));
endfor

fig7 = figure(7);
hax_spectfit = axes('Parent', fig7);
title(hax_spectfit, "Calibration curve");
ylabel(hax_spectfit, "Concentration (ug/l)");
xlabel(hax_spectfit, "Fit");
hold on
plot(calfitzero(:,2), calfitzero(:,1),'*')

L =  polyfit(calfitzero(2:7,2),calfitzero(2:7,1),1)
rangefits = 0:.01:1;
fittedconc = L(1)*rangefits+L(2);
plot(rangefits, fittedconc)
axis([0 1 -100 750] )

save calfit norm_absorbances L rangefits fittedconc minindex maxindex startfit endfit

end
