%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SCION - Spatial Continuous Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Earth Evolution Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Coded by BJW Mills
%%%% b.mills@leeds.ac.uk
%%%%
%%%% plot model fluxes

%%%% output to screen
fprintf('running plotting script... \t')
tic
global state

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   define colorbars   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% IPCC precip colorbar modified
IPCC_pre = [ 223 194 125 ;
246 232 195 ;
245 245 245 ;
199 234 229 ;
128 205 193 ;
53 151 143 ;
1 102 94 ;
0 60 48 ] ./ 255 ;

%%%% IPCC temp colorbar
IPCC_temp = flipud( [103 0 31 ;
178 24 43 ;
214 96 77 ;
244 165 130 ;
253 219 199 ;
247 247 247 ;
209 229 240 ;
146 197 222 ;
67 147 195 ;
33 102 172 ;
5 48 97 ]./ 255 ) ;

%%%% IPCC sequential
IPCC_seq = [255 255 204 ;
161 218 180 ;
65 182 196 ;
44 127 184 ;
37 52 148] ./ 255 ;

%%%% IPCC sequential 2
IPCC_seq_2 = [ 237 248 251 ;
179 205 227 ;
140 150 198 ;
136 86 167 ;
129 15 124 ] ./ 255 ;

%%%% Proxy color chart
pc1 = [65 195 199]./255 ;
pc2 = [73 167 187]./255 ;
pc3 = [82 144 170]./255 ;
pc4 = [88 119 149]./255 ;
pc5 = [89 96 125]./255 ;
pc6 = [82 56 100]./255 ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Plot global variables   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% load geochem data
load('data\geochem_data_2020.mat')
load('data\Scotese_GAT_2021.mat')

%%%%%%% make figure
figure('Color',[1 0.98 0.95])


%%%% GLOBAL FORCINGS
subplot(4,4,1)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
ylim([0 2.5])
xlabel('Time (Ma)')
ylabel('Relative forcing')
%%%% plot this model
plot(state.time_myr,state.DEGASS,'r')
plot(state.time_myr,state.BAS_AREA,'k')
plot(state.time_myr,state.EVO,'g')
plot(state.time_myr,state.W,'b')
plot(state.time_myr,state.Bforcing,'m')
plot(state.time_myr,state.GRAN_AREA,'color',[0.8 0.8 0.8])
%%%% Legend
text(-590,2.4,'D','color','r')
text(-590,2.2,'E','color','g')
text(-590,2,'W','color','b')
text(-590,1.8,'B','color','m')
text(-590,1.6,'BA','color','k')
text(-590,1.4,'GA','color',[0.8 0.8 0.8])
%%%% Title
title('Forcings')

%%% Corg fluxes
subplot(4,4,2)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Flux (mol/yr)')
%%%% plot this model
plot(state.time_myr,state.mocb,'b')
plot(state.time_myr,state.locb,'g')
plot(state.time_myr,state.oxidw,'r')
plot(state.time_myr,state.ocdeg,'k') 
%%%% Legend
text(-590,5e12,'mocb','color','b')
text(-590,4e12,'locb','color','g')
text(-590,3e12,'oxidw','color','r')
text(-590,2e12,'ocdeg','color','k')
%%%% Title
title('C_{org} fluxes')

%%% Ccarb fluxes
subplot(4,4,3)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Flux (mol/yr)')
%%%% plot this model
plot(state.time_myr,state.silw,'r')
plot(state.time_myr,state.carbw,'c')
plot(state.time_myr,state.sfw,'b')
plot(state.time_myr,state.mccb,'k') 
%%%% Legend
text(-590,28e12,'silw','color','r')
text(-590,24e12,'carbw','color','c')
text(-590,20e12,'sfw','color','b')
text(-590,16e12,'mccb','color','k')
%%%% Title
title('C_{carb} fluxes')

%%% S fluxes
subplot(4,4,4)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
% ylim([0 5e12])
ylabel('Fluxes (mol/yr)')
%%%% plot this model
plot(state.time_myr,state.mpsb,'k')
plot(state.time_myr,state.mgsb,'c')
plot(state.time_myr,state.pyrw,'r')
plot(state.time_myr,state.pyrdeg,'m') 
plot(state.time_myr,state.gypw,'b')
plot(state.time_myr,state.gypdeg,'g') 
%%%% Legend
text(-590,1.9e12,'mpsb','color','k')
text(-590,1.7e12,'mgsb','color','c')
text(-590,1.5e12,'pyrw','color','r')
text(-590,1.2e12,'pyrdeg','color','m')
text(-590,1e12,'gypw','color','b')
text(-590,0.8e12,'gypdeg','color','g')
%%%% Title
title('S fluxes')

%%%% C SPECIES
subplot(4,4,5)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Relative size')
%%%% plot this model
plot(state.time_myr,state.G/pars.G0,'k')
plot(state.time_myr,state.C/pars.C0,'c')
plot(state.time_myr,state.VEG,'g--')
%%%% Legend
text(-590,1.5,'VEG','color','g')
text(-590,1.25,'G','color','k')
text(-590,1,'C','color','c')
%%%% Title
title('C reservoirs')

%%%% S SPECIES
subplot(4,4,6)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Relative size')
%%%% plot this model
plot(state.time_myr,state.PYR/pars.PYR0,'k')
plot(state.time_myr,state.GYP/pars.GYP0,'c')
%%%% Legend
text(-590,1,'PYR','color','k')
text(-590,0.9,'GYP','color','c')
%%%% Title
title('S reservoirs')

%%% NUTRIENTS P N
subplot(4,4,7)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
ylim([0 3])
xlabel('Time (Ma)')
ylabel('Relative size')
%%%% plot this model
plot(state.time_myr,state.P/pars.P0,'b')
plot(state.time_myr,state.N/pars.N0,'g')
%%%% Legend
text(-590,1.5,'P','color','b')
text(-590,1,'N','color','g')
%%%% Title
title('Nutrient reservoirs')

%%%% Forg and Fpy ratos
subplot(4,4,8)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('f_{org}, f_{py}')
%%%% plot this model
plot(state.time_myr,state.mocb ./ (state.mocb + state.mccb),'k')
%%%% plot fpy
plot(state.time_myr, state.mpsb ./ (state.mpsb + state.mgsb),'m')

%%%% d13C record
subplot(4,4,9)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('\delta^{13}C_{carb}')
%%%% plot data comparison
plot(d13c_x,d13c_y,'.','color',pc2)
%%%% plot this model
plot(state.time_myr,state.delta_mccb,'k')

%%%% d34S record
subplot(4,4,10)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('\delta^{34}S_{sw}')
%%%% plot data comparison
plot(d34s_x,d34s_y,'.','color',pc2)
%%%% plot this model
plot(state.time_myr,state.d34s_S,'k')

%%%% Ocean 87Sr/86Sr 
subplot(4,4,11)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
ylim([0.706 0.71])
xlabel('Time (Ma)')
ylabel('^{87}Sr/^{86}Sr seawater')
%%%% plot data comparison
plot(sr_x,sr_y,'color',pc2)
%%%% plot this model
plot(state.time_myr,state.delta_OSr,'k')

%%%% SO4
subplot(4,4,12)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
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
plot(state.time_myr,(state.S./pars.S0)*28,'k')

%%%% O2 (%) 
subplot(4,4,13)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Atmospheric O_{2} (%)')
%%%% plot data comparison
for u = 1:2:length(O2_x) - 1
   plot( [O2_x(u) O2_x(u)] , [O2_y(u) O2_y(u+1)] , 'color' , pc2  ) ;     
end
%%%% plot this model
plot(state.time_myr,state.mrO2.*100,'k')

%%%% CO2ppm
subplot(4,4,14)
set(gca, 'YScale', 'log')
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
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
plot(state.time_myr,state.RCO2.*280,'k')


%%%% TEMP
subplot(4,4,15)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
ylim([5 40])
xlabel('Time (Ma)')
ylabel('GAST (C)')
%%%% plot data comparison
% patch(T_x,T_y,pc1,'edgecolor','none')
plot(Scotese_2021_age,Scotese_2021_GAT,'color',pc1)
%%%% plot this model
plot(state.time_myr,state.tempC,'k')
plot(state.time_myr,state.SAT_equator,'r')

%%%% ICE LINE
subplot(4,4,16)
hold on
box on
xlim([pars.whenstart/1e6 pars.whenend/1e6])
xlabel('Time (Ma)')
ylabel('Ice line')
%%%% plot iceline proxy
plot(paleolat_x,paleolat_y,'color' ,pc1) ;
%%%% plot this model
plot(state.time_myr,state.iceline,'k') ;
ylim([0 90])
colormap(gca,'gray')



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
