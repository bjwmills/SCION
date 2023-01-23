%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SCION - Spatial Continuous Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Earth Evolution Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Coded by BJW Mills
%%%% b.mills@leeds.ac.uk
%%%%
%%%% plot sensitivity analysis

%%%%%% define colours
c_mean = [255 132 34]./255 ;
c_std = [255 225 192]./255 ;
c_range = [255 225 192]./255 ;

%%%% Proxy color chart
pc1 = [65 195 199]./255 ;
pc2 = [73 167 187]./255 ;
pc3 = [82 144 170]./255 ;
pc4 = [88 119 149]./255 ;
pc5 = [89 96 125]./255 ;
pc6 = [82 56 100]./255 ;

%%%% output to screen
fprintf('running sens plotting script... \t')
tic

%%%%%%% make figure
figure('Color',[0.80 0.80 0.70])

%%%% load geochem data
load('data/geochem_data_2020.mat')

%%%% make column vector
sens.time_myr = sens.time_myr(:,1) ;

%%%% Forcings
subplot(5,2,1)
hold on
box on
xlabel('Time (Ma)')
ylabel('\delta^{13}C_{carb}')
%%%% plot DEGASS
plot((sens.time_myr),mean(sens.DEGASS,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.DEGASS,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.DEGASS,[],2),'linewidth',0.5,'color',c_range)
%%%% plot GRAN_AREA
plot((sens.time_myr),mean(sens.GRAN_AREA,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.GRAN_AREA,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.GRAN_AREA,[],2),'linewidth',0.5,'color',c_range)
%%%% plot BAS_AREA
plot((sens.time_myr),mean(sens.BAS_AREA,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.BAS_AREA,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.BAS_AREA,[],2),'linewidth',0.5,'color',c_range)

%%%% d13C record
subplot(5,2,2)
hold on
box on
xlabel('Time (Ma)')
ylabel('\delta^{13}C_{carb}')
%%%% plot data comparison
plot(d13c_x,d13c_y,'.','color',pc2)
%%%% plot this model
plot((sens.time_myr),mean(sens.delta_mccb,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.delta_mccb,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.delta_mccb,[],2),'linewidth',0.5,'color',c_range)

%%%% d34S record
subplot(5,2,3)
hold on
box on
xlabel('Time (Ma)')
ylabel('\delta^{34}S_{sw}')
%%%% plot data comparison
plot(d34s_x,d34s_y,'.','color',pc2)
%%%% plot this model
plot((sens.time_myr),mean(sens.d34s_S,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.d34s_S,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.d34s_S,[],2),'linewidth',0.5,'color',c_range)

%%%% Ocean 87Sr/86Sr 
subplot(5,2,4)
hold on
box on
ylim([0.706 0.71])
xlabel('Time (Ma)')
ylabel('^{87}Sr/^{86}Sr seawater')
%%%% plot data comparison
plot(sr_x,sr_y,'color',pc2)
%%%% plot this model
plot((sens.time_myr),mean(sens.delta_OSr,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.delta_OSr,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.delta_OSr,[],2),'linewidth',0.5,'color',c_range)

%%%% SO4
subplot(5,2,5)
hold on
box on
xlabel('Time (Ma)')
ylabel('Marine SO_{4} (mM)')
%%%% plot algeo data window comparison
plot(sconc_max_x,sconc_max_y,'color',pc1)
plot(sconc_min_x,sconc_min_y,'color',pc1)
plot(sconc_mid_x,sconc_mid_y,'color',pc2)
%%%% plot fluid inclusion data comparison
for u = 1:2:length(SO4_x-1)
   plot( [SO4_x(u) SO4_x(u)] , [SO4_y(u) SO4_y(u+1)], 'color' , pc3 ) ;     
end
%%%% plot this model
plot((sens.time_myr),mean(sens.SmM,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.SmM,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.SmM,[],2),'linewidth',0.5,'color',c_range)

%%%% O2 (%) 
subplot(5,2,6)
hold on
box on
xlabel('Time (Ma)')
ylabel('Atmospheric O_{2} (%)')
%%%% plot data comparison
for u = 1:2:length(O2_x) - 1
   plot( [O2_x(u) O2_x(u)] , [O2_y(u) O2_y(u+1)] , 'color' , pc2  ) ;     
end
%%%% plot this model
plot((sens.time_myr),mean(sens.mrO2.*100,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.mrO2.*100,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.mrO2.*100,[],2),'linewidth',0.5,'color',c_range)


%%%% CO2ppm
subplot(5,2,7)
set(gca, 'YScale', 'log')
hold on
box on
ylim([100 10000])
xlabel('Time (Ma)')
ylabel('Atmospheric CO_{2} (ppm)')
%%%% plot data comparison
%%%% paleosol
% errorbar(paleosol_age,paleosol_co2,paleosol_low,paleosol_high,'color',[0.4 0.7 0.7],'linestyle','none')
plot(paleosol_age, paleosol_co2,'.','markerfacecolor',pc1,'markeredgecolor',pc1)
%%%% alkenone
% errorbar(alkenone_age,alkenone_co2,alkenone_low,alkenone_high,'color',[0.4 0.7 0.4],'linestyle','none')
plot(alkenone_age, alkenone_co2,'.','markerfacecolor',pc2,'markeredgecolor',pc2)
%%%% boron
% errorbar(boron_age,boron_co2,boron_low,boron_high,'color',[0.7 0.4 0.4],'linestyle','none')
plot(boron_age, boron_co2,'.','markerfacecolor',pc3,'markeredgecolor',pc3)
%%%% stomata
% errorbar(stomata_age,stomata_co2,stomata_low,stomata_high,'color',[0.7 0.7 0.4],'linestyle','none')
plot(stomata_age, stomata_co2,'.','markerfacecolor',pc4,'markeredgecolor',pc4)
%%%% liverwort
% errorbar(liverwort_age,liverwort_co2,liverwort_low,liverwort_high,'color',[0.7 0.7 0.4],'linestyle','none')
plot(liverwort_age, liverwort_co2,'.','markerfacecolor',pc5,'markeredgecolor',pc5)
%%%% phytane
% errorbar(phytane_age,phytane_co2,phytane_low,phytane_high,'color',[0.7 0.7 0.4],'linestyle','none')
plot(phytane_age, phytane_co2,'.','markerfacecolor',pc6,'markeredgecolor',pc6)
%%%% plot this model
plot((sens.time_myr),mean(sens.CO2ppm,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.CO2ppm,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.CO2ppm,[],2),'linewidth',0.5,'color',c_range)

%%%% TEMP
subplot(5,2,8)
hold on
box on
ylim([5 40])
xlabel('Time (Ma)')
ylabel('GAST (C)')
%%%% plot data comparison
patch(T_x,T_y,pc2,'edgecolor','none')
%%%% plot this model
plot((sens.time_myr),mean(sens.T_gast,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.T_gast,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.T_gast,[],2),'linewidth',0.5,'color',c_range)
%%%% plot this model torpical T
plot((sens.time_myr),mean(sens.SAT_tropical,2),'linewidth',1,'color',c_mean,'linestyle',':')
plot((sens.time_myr),max(sens.SAT_tropical,[],2),'linewidth',0.5,'color',c_range,'linestyle',':')
plot((sens.time_myr),min(sens.SAT_tropical,[],2),'linewidth',0.5,'color',c_range,'linestyle',':')
%%%% plot this model equatorial T
plot((sens.time_myr),mean(sens.SAT_equator,2),'linewidth',1,'color',c_mean,'linestyle',':')
plot((sens.time_myr),max(sens.SAT_equator,[],2),'linewidth',0.5,'color',c_range,'linestyle',':')
plot((sens.time_myr),min(sens.SAT_equator,[],2),'linewidth',0.5,'color',c_range,'linestyle',':')

%%%% ICE LINE
subplot(5,2,9)
hold on
box on
xlabel('Time (Ma)')
ylabel('Ice line')
%%%% plot iceline proxy
plot(paleolat_x,paleolat_y,'color' ,pc2) ;
%%%% plot this model
plot((sens.time_myr),mean(sens.iceline,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.iceline,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.iceline,[],2),'linewidth',0.5,'color',c_range)

%%%% P and N
subplot(5,2,10)
hold on
box on
xlabel('Time (Ma)')
ylabel('P (-), N (--)')
%%%% plot this model
plot((sens.time_myr),mean(sens.P,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.P,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.P,[],2),'linewidth',0.5,'color',c_range)
plot((sens.time_myr),mean(sens.N,2),'--','linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.N,[],2),'--','linewidth',0.5,'color',c_range)
plot((sens.time_myr),min(sens.N,[],2),'--','linewidth',0.5,'color',c_range)


%%%%% plotting script finished
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Cleanup workspace   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear stepnumber
clear u
clear numfields
clear trecords
clear finalrecord
clear field_names
clear n
clear veclength
clear xvec
clear yvec
clear endtime


%%%%% plotting script finished
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )
