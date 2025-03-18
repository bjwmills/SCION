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
%     1----1001111------------0111        'coastlines' utility plotting script                  %
%      1----1101-------------1101         Plots coastlines at keyframe timepoints               %
%        1--111----------------1          Call SCION_plot_coastlines(topo,lon,lat)              %
%           1---------------1                                                                   %
%               111011011                                                                       %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function coastlines = SCION_plot_coastlines(lon,lat,topo)

    %%%% alter topogrpahy to binary matrix for land/sea and trace with line segments
    coastline = topo ;
    coastline(topo>0) = 1 ;
    coastline(isnan(topo)==1) = 0 ;
    lonfix = 172.5;
    latfix = 4.5 ;

    %%%% draw line segments
    for i = 2:size(coastline, 1)
        for j = 2:size(coastline, 2)
            if coastline(i, j) == 1
                % Check edges for transition to 0
                if i > 1 && coastline(i-1, j) == 0 % Top edge
                    m_line([lon(j-1)-lonfix, lon(j)-lonfix], [lat(i-1)-latfix, lat(i-1)-latfix], 'Color', 'k', 'LineWidth', 0.5);
                end
                if i < size(coastline, 1) && coastline(i+1, j) == 0 % Bottom edge
                    m_line([lon(j-1)-lonfix, lon(j)-lonfix], [lat(i)-latfix, lat(i)-latfix], 'Color', 'k', 'LineWidth', 0.5);
                end
                if j > 1 && coastline(i, j-1) == 0 % Left edge
                    m_line([lon(j-1)-lonfix, lon(j-1)-lonfix], [lat(i-1)-latfix, lat(i)-latfix], 'Color', 'k', 'LineWidth', 0.5);
                end
                if j < size(coastline, 2) && coastline(i, j+1) == 0 % Right edge
                    m_line([lon(j)-lonfix, lon(j)-lonfix], [lat(i-1)-latfix, lat(i)-latfix], 'Color', 'k', 'LineWidth', 0.5);
                end
            end
        end
    end

end