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
%     1----1001111------------0111        Tuning script                                         %
%      1----1101-------------1101         Run this script to perform tuning for SCION model     %
%        1--111----------------1          runs model in parallel pattern search to determine    %
%           1---------------1             best choice for intiial reservoir values              %
%               111011011                                                                       %
%                                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% initial guess for [ G C PYR GYP O S A ] 
inputs = [0.45 1 1.1 1 0.1 0.05 3] ;

%%%% single core
% outputs = fminsearch(@SCION_tuning_function,inputs) ;

%%%% single core
% outputs = patternsearch(@SCION_tuning_function,inputs,'','','','',0.2.*inputs,2.*inputs) ;

%%%% parellel
options = optimoptions('patternsearch','UseParallel',true,'MeshTolerance',5e-7) ;
[X1,Fval,Exitflag,Output] = patternsearch(@SCION_tuning_function,inputs,'','','','',0.05.*inputs,3.*inputs,'',options) ;

%%%% display output
fprintf('The number of iterations is: %d\n', Output.iterations);
fprintf('The number of function evaluations is: %d\n', Output.funccount);
fprintf('The best function value found is: %g\n', Fval);
