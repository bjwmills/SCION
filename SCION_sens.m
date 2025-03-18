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
%     1----1001111------------0111        Model sensitivity initialiser                         %
%      1----1101-------------1101         Run this script to start a sensitivity analysis       %
%        1--111----------------1                                                                %
%           1---------------1                                                                   %
%               111011011                                                                       %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% number of runs
sensruns = 100 ;

global tuning

%%%%%%% multiple runs
parfor N = 1:sensruns
       
    %%%%%%% run model
    fprintf('Run # %d', N )
    fprintf(' of %d \n', sensruns )
    run(N) = SCION_initialise(N) ;

end

%%%%%% define standard time grid for outputs using first model run's grid.
tgrid =  run(1).state.time ;

%%%%%% option for standard spaced grid
% tgrid = ( run(1).state.time(1) : 1e6 : run(1).state.time(end) ) ;

%%%%%% sens analysis states mapped to tgrid
for N = sensruns:-1:1 %% iterate in reverse to avoid index shifting
    field_names = fieldnames(run(N).state);
    for numfields = 1:length(field_names)
        field_name = field_names{numfields};
        sens.(field_name)(:, N) = interp1(run(N).state.time, run(N).state.(field_name), tgrid);
    end

    %%%% Check if run completed and remove if not
    if any(isnan(sens.ANOX(:, N)))
        for numfields = 1:length(field_names)
            field_name = field_names{numfields};
            sens.(field_name)(:, N) = []; 
        end
    end
end


%%%%%% plotting
SCION_plot_sens

%%%%%% write output file
save('SCION_results.mat','sens','-mat')





