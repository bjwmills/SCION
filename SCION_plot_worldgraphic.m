%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                               %
%              110111010                                                                        %
%           111010-1-----101                                                                    %
%        1011111---------101111                                                                 %
%      11011------------------101         SCION: Spatial Continuous Integration                 %
%     111-----------------10011011        Earth Evolution Model                                 %
%    1--10---------------1111011111                                                             %
%    1---1011011---------1010110111       Lead developer: Benjamin J. W. Mills                  %
%    1---1011000111----------010011       email: b.mills@leeds.ac.uk                            %
%    1----1111011101----------10101                                                             %
%     1----1001111------------0111        'World graphic' plotting script                       %
%      1----1101-------------1101         Plots model grid states at keyframe timepoints        %
%        1--111----------------1          CALLED IN SCRIPT - DO NOT RUN DIRECTLY                %
%           1---------------1                                                                   %
%               111011011                                                                       %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% output to screen
fprintf('running plotting script... \t')
tic
global gridstate

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
%%%%%%%%%%%%%%%%%   Plot gridstates   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for f = [1 2]

    %%%% use multiple figures to plot all time slices
    if f == 1
        choose_gridsubs = [1:14] ;
    else
        choose_gridsubs = [12:28] ;
    end

    %%%%%%% make figure
    figure('Color',[1 0.98 0.95])
    ha = tight_subplot(length(choose_gridsubs),5,0,0.01,0.01) ;

    %%%% loop over number of gridstates
    subnumber = 1 ;
    for gridsub = choose_gridsubs

        %%%% make topography
        this_TOPO = gridstate.TOPO(:,:,gridsub) ;
        this_TOPO(this_TOPO<1) = NaN ;
        cmap = ttcmap([-5000 5000],'cmap','gmtrelief');


        %%%% plot simplified topography
        axes(ha( 1 + 5*( length(choose_gridsubs) - subnumber)  )) ;
        hold on
        m_proj('robinson','longitude',[-180 172.5],'latitude',[-86.6 86.6])
        h = m_pcolor(INTERPSTACK.lon - 180,INTERPSTACK.lat,this_TOPO) ;
        m_grid('box','on','xticklabels','','yticklabels','','linestyle','-','gridcolor',[0.7 0.7 0.7],'ticklength',0)
        colormap(gca,cmap)
        clim([-5000 5000])
        axis off
        text(4,0,num2str( round(gridstate.time_myr(gridsub) ) ))


        %%%% Nan out the ocean on Tair
        thisfield = gridstate.Tair(:,:,gridsub);
        thisfield(gridstate.land(:,:,gridsub) == 0 ) = NaN ;

        %%%% plot Tair
        axes(ha( 2 + 5*( length(choose_gridsubs) - subnumber)  )) ;
        hold on
        m_proj('robinson','longitude',[-180 172.5],'latitude',[-86.6 86.6])
        h = m_pcolor(INTERPSTACK.lon - 180,INTERPSTACK.lat,thisfield) ;
        m_grid('box','on','xticklabels','','yticklabels','','linestyle','-','gridcolor',[0.7 0.7 0.7],'ticklength',0)
        clim([-40 40])
        axis off
        colormap(gca, IPCC_temp )
        if gridsub == 1
            text(4,0,'Air Temp (C)')
        end


        %%%% Nan out the ocean on Runoff
        thisfield = gridstate.Q(:,:,gridsub);
        thisfield(gridstate.land(:,:,gridsub) == 0 ) = NaN ;
 
        %%%% plot Q
        axes(ha( 3 + 5*( length(choose_gridsubs) - subnumber)  )) ;
        hold on
        m_proj('robinson','longitude',[-180 172.5],'latitude',[-86.6 86.6])
        h = m_pcolor(INTERPSTACK.lon - 180,INTERPSTACK.lat,log10(thisfield)) ;
        m_grid('box','on','xticklabels','','yticklabels','','linestyle','-','gridcolor',[0.7 0.7 0.7],'ticklength',0)
        clim([0 4])
        axis off
        colormap(gca, IPCC_pre )
        if gridsub == 1
            text(4,0,'Runoff (log mm/yr)')
        end


        %%%% Nan out the ocean on epsilon
        thisfield = gridstate.EPSILON(:,:,gridsub);
        thisfield(gridstate.land(:,:,gridsub) == 0 ) = NaN ;
 
        %%%% plot epsilon
        axes(ha( 4 + 5*( length(choose_gridsubs) - subnumber)  )) ;
        hold on
        m_proj('robinson','longitude',[-180 172.5],'latitude',[-86.6 86.6])
        h = m_pcolor(INTERPSTACK.lon - 180,INTERPSTACK.lat,log10(thisfield)) ;
        m_grid('box','on','xticklabels','','yticklabels','','linestyle','-','gridcolor',[0.7 0.7 0.7],'ticklength',0)
        clim([0 4])
        axis off
        colormap(gca, IPCC_seq_2 )
        if gridsub == 1
            text(4,0,'Erosion (log t/km2/yr)')
        end


        %%%% Nan out the ocean on silw
        thisfield = gridstate.CW(:,:,gridsub);
        thisfield(gridstate.land(:,:,gridsub) == 0 ) = NaN ;

        %%%% plot silw
        axes(ha( 5 + 5*( length(choose_gridsubs) - subnumber)  )) ;
        hold on
        m_proj('robinson','longitude',[-180 172.5],'latitude',[-86.6 86.6])
        h = m_pcolor(INTERPSTACK.lon - 180,INTERPSTACK.lat,log10(thisfield)) ;
        m_grid('box','on','xticklabels','','yticklabels','','linestyle','-','gridcolor',[0.7 0.7 0.7],'ticklength',0)
        clim([0 2])
        axis off
        colormap(gca, IPCC_seq )
        if gridsub == 1
            text(4,0,'Silw (log t/km2/yr)')
        end


        subnumber = subnumber + 1 ;
    end

end

%%%%% plotting script finished
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )
