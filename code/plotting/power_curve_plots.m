opt_arclen_purcell_previous=9.648*r;
opt_wave_length_purcell_previous=3.687*r;
opt_hrad_purcell_previous=0.5695*r;

wave_length_array=[linspace(0.5*r,12*r,11) opt_wave_length_purcell_previous];
wave_length_array=sort(wave_length_array,"descend");

hrad_array=[linspace(0.2*r,0.8*r,11) opt_hrad_purcell_previous];
hrad_array=sort(hrad_array);

arclen_array=[linspace(4*r,20*r,11) opt_arclen_purcell_previous];
arclen_array=sort(arclen_array,"descend");

energy_per_dist_wl=readmatrix('/Users/nathan/Downloads/power_curve_wavelength07-21-2025_14-48/energy_per_dist.txt');
energy_per_dist_hrad=readmatrix('/Users/nathan/Downloads/power_curve_hrad07-22-2025_09-24/energy_per_dist.txt');
energy_per_dist_arclen=readmatrix('/Users/nathan/Downloads/power_curve_arclen07-21-2025_14-49/energy_per_dist.txt');

bigdir='./power_curve_plots';
mkdir(bigdir)

figure('Theme','light')
plot(wave_length_array/r,diag(energy_per_dist_wl))
xline(opt_wave_length_purcell_previous/r,"Color","k","Label","Purcell-inefficiency-minimizing wavelength","LabelVerticalAlignment","middle")
xlabel("\lambda/r")
ylabel("Energy per distance (10^{-12} J/m)")
saveas(gcf,[bigdir,'/energy_per_dist_vs_wl.png'])

figure('Theme','light')
plot(hrad_array/r,diag(energy_per_dist_hrad))
xline(opt_hrad_purcell_previous/r,"Color","k","Label","Purcell-inefficiency-minimizing radius","LabelVerticalAlignment","middle")
xlabel("R/r")
ylabel("Energy per distance (10^{-12} J/m)")
saveas(gcf,[bigdir,'/energy_per_dist_vs_hrad.png'])

figure('Theme','light')
plot(arclen_array/r,diag(energy_per_dist_arclen))
xline(opt_arclen_purcell_previous/r,"Color","k","Label","Purcell-inefficiency-minimizing arc length","LabelVerticalAlignment","middle")
xlabel("S/r")
ylabel("Energy per distance (10^{-12} J/m)")
saveas(gcf,[bigdir,'/energy_per_dist_vs_arclen.png'])

