classdef RunParams < LFADS.RunParams
   properties
       align = 'MoveOnsetOnline';
       preKeep = 400;
       postKeep = 500; 
       nTrialsKeep = 0;
       
       % for alignement matrix
       pcsKeep = 8;
       pcTrialAvg = true;
       
       
   end
   
   methods
       function p = RunParams()
           % set default values here
           p.spikeBinMs = 2;
           p.batchSize = 200; % trials
           p.regularizerIncreaseSteps = 900;
           p.learningRateDecayFactor = 0.98;
           p.keepProb = 0.95; % dropout of units in the network
       end
       
       function suffix = generateSuffix(p)
           % only include a parameter if it differs from the default
           if p.nTrialsKeep == 0
               nTrialsStr = 'All';
           else
               nTrialsStr = num2str(p.nTrialsKeep);
           end
           if p.pcTrialAvg
               pcStr = sprintf('trialAvgPCs%d', p.pcsKeep);
           else
               pcStr = sprintf('singleTrialPCs%', p.pcsKeep);
           end
           if p.regularizerIncreaseSteps == 900
               regIncStr = '';
           else
               regIncStr = sprintf('_regIncSteps%d', p.regularizerIncreaseSteps);
           end
           if p.learningRateDecayFactor == 0.98
               learnRateStr = '';
           else
               learnRateStr = sprintf('_learnRateFac%g', p.learningRateDecayFactor);
           end
           if p.keepProb == 0.95
               keepProbStr = '';
           else
               keepProbStr = sprintf('_keepProb%g', p.learningRateDecayFactor);
           end
           
           switch p.align
                case 'GoCue'
                    if p.preKeep ~= 200 || ppostKeep ~= 700 
                        alignStr = sprintf('GoCue_pre%d_post%d', p.preKeep, p.postKeep);
                    else
                        alignStr = 'GoCue';
                    end
                case 'MoveOnsetOnline'
                    if p.preKeep ~= 400 || p.postKeep ~= 500 
                        alignStr = sprintf('MoveOnsetOnline_pre%d_post%d', p.preKeep, p.postKeep);
                    else
                        alignStr = 'MoveOnsetOnline';
                    end
                otherwise
                    error('Unknown align %s', r.params.align);
            end
           
           suffix = sprintf('bin%02d_batch%03d_align%s_nTrialsKeep%s_%s%s%s%s', ...
               p.spikeBinMs, p.batchSize, alignStr, nTrialsStr, pcStr, regIncStr,learnRateStr, keepProbStr);
       end
   end
end