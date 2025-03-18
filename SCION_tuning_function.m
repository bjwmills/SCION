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
%     1----1001111------------0111        Model tuning function                                 %
%      1----1101-------------1101         Runs model and calculates cost function for tuning    %
%        1--111----------------1          CALLED IN SCRIPT - DO NOT RUN DIRECTLY                %
%           1---------------1                                                                   %
%               111011011                                                                       %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function costfunc = SCION_tuning_function(params)
    global state
    global tuning
    
    %%% set up starting reservoir modifiers
    tuning.Gtune = params(1) ;
    tuning.Ctune = params(2) ;
    tuning.PYRtune = params(3) ;
    tuning.GYPtune = params(4) ;
    tuning.Otune = params(5) ;
    tuning.Stune = params(6) ;
    tuning.Atune = params(7) ;
        
    %%% run model
    SCION_initialise(1) ;
    
    %%%% error and const function
    costfunc = ( state.Orel(end) - 1 )^2 + ( state.Arel(end) - 1 )^2 + ( state.Srel(end) - 1 )^2 + ( state.PYRrel(end) - 1 )^2 ...
        + ( state.GYPrel(end) - 1 )^2 + ( state.Crel(end) - 1 )^2 + ( state.Grel(end) - 1 )^2 ;
    
    %%%% print update
    fprintf('Reservoirs: G    C    PYR  GYP  O    S    A %d \n') ; 
    fprintf('\n')
    fprintf('Parameters: ')
    fprintf('%d    ',params)
    fprintf('\n')
    fprintf('Final vals: ')
    fprintf('%1.1f  ',[state.Grel(end) state.Crel(end) state.PYRrel(end) state.GYPrel(end) state.Orel(end) state.Srel(end) state.Arel(end)])
    fprintf('\n')
    fprintf('chisquared: %d \n', costfunc )

end



