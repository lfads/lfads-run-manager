classdef ModelTrainedParams < matlab.mixin.CustomDisplay
    
    properties(SetAccess=protected)
        datasetNames
        
        % Read-in from spikes to factors
        x_to_infac_W % readin alignment weights, mapping from counts to input factors; nDatasets x 1 cell of nNeurons x c_factors_dim
        x_to_infac_b % readin alignment biases to input factors; nDatasets x 1 cell of 1 x c_factors_dim
        
        % Initial condition encoder (forward)
        ic_enc_fwd_t0 % forward IC encoder prior on t0; 1 x c_ic_enc_dim
        ic_enc_fwd_gru_xh_to_gates_ru_W % forward IC encoder GRU, mapping input+hiddens to gates r and u, weights; (c_ic_enc_dim + factors_dim) x (2 * c_ic_enc_dim)
        ic_enc_fwd_gru_xh_to_gates_ru_b % forward IC encoder GRU bmapping input+hiddens to gates r and u, biases; 1 x (2*c_ic_enc_dim)
        ic_enc_fwd_gru_xrh_to_c_W % forward IC encoder GRU mapping input, r, and hidden to candidates (weights); (c_ic_enc_dim + c_factors_dim) x c_ic_enc_dim
        ic_enc_fwd_gru_xrh_to_c_b % forward IC encoder GRU mapping input, r, and hidden to candidates (bias); 1 x c_ic_enc_dim
        
        % Initial condition encoder (reverse)
        ic_enc_rev_t0 % reverse IC encoder prior on t0; 1 x c_ic_enc_dim
        ic_enc_rev_gru_xh_to_gates_ru_W % reverse IC encoder GRU, mapping input+hidden to gates r and u, weights; (c_factors_dim + c_ic_enc_dim) x (2*c_ic_enc_dim)
        ic_enc_rev_gru_xh_to_gates_ru_b % reverse IC encoder GRU bmapping input+hidden to gates r and u, biases; 1 x (2*c_ic_enc_dim)
        ic_enc_rev_gru_xrh_to_c_W % reverse IC encoder GRU mapping input+r+hidden to candidates (weights); (c_ic_enc_dim + c_factors_dim) x c_ic_enc_dim
        ic_enc_rev_gru_xrh_to_c_b % reverse IC encoder GRU mapping input+r+hidden to candidates (bias); 1 x c_ic_enc_dim
        
        % Initial condition g0
        prior_g0_mean % Mean parameter in prior on initial condition g0; 1 x c_ic_dim
        prior_g0_logvar % Logvar parameter in prior on initial condition g0; 1 x c_ic_dim
        ic_enc_to_posterior_g0_mean_W % Weights for mean parameter in posterior of the initial condition g0; (2 * c_ic_enc_dim) x c_ic_dim
        ic_enc_to_posterior_g0_mean_b  % Bias for mean parameter in posterior of the initial condition g0; 1 x c_ic_dim
        ic_enc_to_posterior_g0_logvar_W  % Weights for logvar parameter in posterior of the initial condition g0; (2 * c_ic_enc_dim) x c_ic_dim
        ic_enc_to_posterior_g0_logvar_b  % Bias for logvar parameter in posterior of the initial condition g0;  1 x c_ic_dim
        g0_to_gen_ic_W % mapping from g0 to generator initial condition, weights; c_ic_dim x c_gen_dim
        g0_to_gen_ic_b % mapping from g0 to generator initial condition, bias; 1 x c_gen_dim
        
        % Controller encoder (forward)
        ci_enc_fwd_t0 % forward controller prior on t0; 1 x c_ci_enc_dim
        ci_enc_fwd_gru_xh_to_ru_W % forward controller encoder GRU, mapping input+hidden to gates r and u, weights; (ci_enc_dim + factors_dim) x (2*c_ci_enc_dim)
        ci_enc_fwd_gru_xh_to_ru_b % forward controller encoder GRU, mapping input+hidden to gates r and u, bias; 1 x (2*c_ci_enc_dim)
        ci_enc_fwd_gru_xrh_to_c_W % forward controller encoder GRU mapping input, r, and hidden to candidates (weights); (c_ci_enc_dim + c_factors_dim) x c_ci_enc_dim)
        ci_enc_fwd_gru_xrh_to_c_b % forward controller encoder GRU mapping input, r, and hidden to candidates (bias); 1 x c_ci_enc_dim
        
        % Controller encoder (reverse)
        ci_enc_rev_t0 % reverse controller prior on t0; 1 x c_ci_enc_dim
        ci_enc_rev_gru_xh_to_ru_W % reverse controller encoder GRU, mapping input+hidden to gates r and u, weights; (ci_enc_dim + factors_dim) x (2*c_ci_enc_dim)
        ci_enc_rev_gru_xh_to_ru_b % reverse controller encoder GRU, mapping input+hidden to gates r and u, bias; 1 x (2*c_ci_enc_dim)
        ci_enc_rev_gru_xrh_to_c_W % reverse controller encoder GRU mapping input, r, and hidden to candidates (weights); (c_ci_enc_dim + c_factors_dim) x c_ci_enc_dim)
        ci_enc_rev_gru_xrh_to_c_b % reverse controller encoder GRU mapping input, r, and hidden to candidates (bias); 1 x c_ci_enc_dim
        
        % Controlller RNN
        con_gengru_x_to_ru_W % controller GenGRU, mapping input to gates r+u, weights; (c_ci_enc_dim * 2 + c_factors_dim) x (2*c_con_dim)
        con_gengru_h_to_ru_W % controller GenGRU, mapping hidden to gates r+u, weights; c_con_dim x (2*c_con_dim)
        con_gengru_h_to_ru_b % controller GenGRU, mapping hidden to gates r+u, weights; 1 x (2*c_con_dim)
        con_gengru_x_to_c_W % controller GenGRU, mapping input to candidates, weights; (c_ci_enc_dim * 2 + c_factors_dim) x c_con_dim
        con_gengru_rh_to_c_W % controller GenGRU, mapping r+hidden to candidates, weights; 1 x c_con_dim
        con_gengru_rh_to_c_b % controller GenGRU, mapping r+hidden to candidates, bias; 1 x c_con_dim
        
        % Controller output co
        prior_ar1_logevars % autoregressive prior on controller outputs; 1 x c_co_dim
        prior_ar1_logatau % autoregressive time constant prior on controller outputs; 1 x c_co_dim
        con_co % prior on controller output; 1 x c_con_dim
        con_to_posterior_co_mean_W % mapping from controller to mean parameter of co, weights; c_con_dim x c_co_dim
        con_to_posterior_co_mean_b % mapping from controller to mean parameter of co, biases; 1 x c_co_dim
        con_to_posterior_co_logvar_W % mapping from controller to logvar parameter of co, weights; c_con_dim x c_co_dim
        con_to_posterior_co_logvar_b % mapping from controller to logvar parameter of co, biases; 1 x c_co_dim
        
        % Generator RNN
        gen_gengru_x_to_ru_W % generator GRU, mapping from input to gates r+u, weights; c_co_dim x (2*c_gen_dim)
        gen_gengru_h_to_ru_W % generator GRU, mapping from input to gates r+u, weights; c_gen_dim x (2*c_gen_dim)
        gen_gengru_h_to_ru_b % generator GRU, mapping from input to gates r+u, biases; 1 x (2*c_gen_dim)
        gen_gengru_x_to_c_W % generator GRU, mapping from input to candidates, weights; c_co_dim x c_gen_dim
        gen_gengru_rh_to_c_W % generator GRU, mapping from r+hidden to candidates, weights; c_gen_dim x c_gen_dim
        gen_gengru_rh_to_c_b % generator GRU, mapping from r+hidden to candidates, biases; 1 x c_gen_dim
        
        % Generator output
        gen_to_factors_W % mapping from generator to factors, weights; c_gen_dim x c_factors_dim
        factors_to_logrates_W % readout alignment weights; nDatasets x 1 cell of c_factors_dim x nNeurons
        factors_to_logrates_b % readout alignment biases; nDatasets x 1 cell of 1 x nNeurons
    end
    
    properties(Dependent)
        nDatasets 
    end
   
    methods
        function mtp = ModelTrainedParams(fname, datasetNames,varargin)
            p = inputParser();
            p.addParameter('verbose', false, @islogical);
            p.parse(varargin{:});
            verbose = p.Results.verbose;
            
            assert(exist(fname, 'file') ~= 0, 'File does not exist');
            
            info = h5info(fname);
            keys = {info.Datasets.Name}';
            readList = {};
            
            nD = numel(datasetNames);
            mtp.datasetNames = datasetNames;
            
            function val = read(key)
                if ismember(key, keys)
                    val = h5read(fname, ['/' key]);
                else
                    if verbose
                        fprintf('Key %s not found\n', key);
                    end
                    val = [];
                end
                
                if verbose && ismember(key, readList)
                    fprintf('Key %s read twice\n', key);
                end
                readList = union(readList, key);
            end
            
            function valCell = readByDataset(key)
                valCell = cell(nD, 1);
                for iD = 1:nD
                    valCell{iD} = read(sprintf(key, datasetNames{iD}));
                    if isempty(valCell{iD})
                        valCell = {};
                        return;
                    end
                end
            end
            
            % feeling grateful for regular expressions =)
            mtp.x_to_infac_W = readByDataset('LFADS_x_2_infac_%s.h5_W:0');
            mtp.x_to_infac_b = readByDataset('LFADS_x_2_infac_%s.h5_b:0');
            mtp.factors_to_logrates_W = readByDataset('LFADS_glm_fac_2_logrates_%s.h5_W:0');
            mtp.factors_to_logrates_b = readByDataset('LFADS_glm_fac_2_logrates_%s.h5_b:0');
                
            mtp.ic_enc_fwd_t0 = read('LFADS_ic_enc_fwd_ic_enc_t0_fwd:0');
            mtp.ic_enc_fwd_gru_xh_to_gates_ru_W = read('LFADS_ic_enc_fwd_GRU_Gates_xh_2_ru_W:0');
            mtp.ic_enc_fwd_gru_xh_to_gates_ru_b = read('LFADS_ic_enc_fwd_GRU_Gates_xh_2_ru_b:0');
            mtp.ic_enc_fwd_gru_xrh_to_c_W = read('LFADS_ic_enc_fwd_GRU_Candidate_xrh_2_c_W:0');
            mtp.ic_enc_fwd_gru_xrh_to_c_b = read('LFADS_ic_enc_fwd_GRU_Candidate_xrh_2_c_b:0');
            
            mtp.ic_enc_rev_t0 = read('LFADS_ic_enc_rev_ic_enc_t0_rev:0');
            mtp.ic_enc_rev_gru_xh_to_gates_ru_W = read('LFADS_ic_enc_rev_GRU_Gates_xh_2_ru_W:0');
            mtp.ic_enc_rev_gru_xh_to_gates_ru_b = read('LFADS_ic_enc_rev_GRU_Gates_xh_2_ru_b:0');
            mtp.ic_enc_rev_gru_xrh_to_c_W = read('LFADS_ic_enc_rev_GRU_Candidate_xrh_2_c_W:0');
            mtp.ic_enc_rev_gru_xrh_to_c_b = read('LFADS_ic_enc_rev_GRU_Candidate_xrh_2_c_b:0');
            
            mtp.prior_g0_mean = read('LFADS_z_prior_g0_mean:0');
            mtp.prior_g0_logvar = read('LFADS_z_prior_g0_logvar:0');
            mtp.ic_enc_to_posterior_g0_mean_W = read('LFADS_z_ic_enc_2_post_g0_mean_W:0');
            mtp.ic_enc_to_posterior_g0_mean_b = read('LFADS_z_ic_enc_2_post_g0_mean_b:0');
            mtp.ic_enc_to_posterior_g0_logvar_W = read('LFADS_z_ic_enc_2_post_g0_logvar_W:0');
            mtp.ic_enc_to_posterior_g0_logvar_b = read('LFADS_z_ic_enc_2_post_g0_logvar_b:0');
            mtp.g0_to_gen_ic_W = read('LFADS_gen_g0_2_gen_ic_W:0');
            mtp.g0_to_gen_ic_b = read('LFADS_gen_g0_2_gen_ic_b:0');
            
            mtp.ci_enc_fwd_t0 = read('LFADS_ci_enc_fwd_ci_enc_t0_fwd:0');
            mtp.ci_enc_fwd_gru_xh_to_ru_W =  read('LFADS_ci_enc_fwd_GRU_Gates_xh_2_ru_W:0');
            mtp.ci_enc_fwd_gru_xh_to_ru_b = read('LFADS_ci_enc_fwd_GRU_Gates_xh_2_ru_b:0');
            mtp.ci_enc_fwd_gru_xrh_to_c_W = read('LFADS_ci_enc_fwd_GRU_Candidate_xrh_2_c_W:0');
            mtp.ci_enc_fwd_gru_xrh_to_c_b = read('LFADS_ci_enc_fwd_GRU_Candidate_xrh_2_c_b:0');
            
            mtp.ci_enc_rev_t0 = read('LFADS_ci_enc_rev_ci_enc_t0_rev:0');
            mtp.ci_enc_rev_gru_xh_to_ru_W = read('LFADS_ci_enc_rev_GRU_Gates_xh_2_ru_W:0');
            mtp.ci_enc_rev_gru_xh_to_ru_b = read('LFADS_ci_enc_rev_GRU_Gates_xh_2_ru_b:0');
            mtp.ci_enc_rev_gru_xrh_to_c_W = read('LFADS_ci_enc_rev_GRU_Candidate_xrh_2_c_W:0');
            mtp.ci_enc_rev_gru_xrh_to_c_b = read('LFADS_ci_enc_rev_GRU_Candidate_xrh_2_c_b:0');
            
            mtp.prior_ar1_logevars = read('LFADS_z_u_prior_ar1_logevars:0');
            mtp.prior_ar1_logatau = read('LFADS_z_u_prior_ar1_logatau:0');
            mtp.con_gengru_x_to_ru_W = read('LFADS_con_GenGRU_Gates_x_2_ru_W:0');
            mtp.con_gengru_h_to_ru_W = read('LFADS_con_GenGRU_Gates_h_2_ru_W:0');
            mtp.con_gengru_h_to_ru_b = read('LFADS_con_GenGRU_Gates_h_2_ru_b:0');
            mtp.con_gengru_x_to_c_W = read('LFADS_con_GenGRU_Candidate_x_2_c_W:0');
            mtp.con_gengru_rh_to_c_W = read('LFADS_con_GenGRU_Candidate_rh_2_c_W:0');
            mtp.con_gengru_rh_to_c_b = read('LFADS_con_GenGRU_Candidate_rh_2_c_b:0');
            
            mtp.con_co = read('LFADS_con_c0:0');
            mtp.con_to_posterior_co_mean_W = read('LFADS_con_con_to_post_co_mean_W:0');
            mtp.con_to_posterior_co_mean_b = read('LFADS_con_con_to_post_co_mean_b:0');
            mtp.con_to_posterior_co_logvar_W = read('LFADS_con_con_to_post_co_logvar_W:0');
            mtp.con_to_posterior_co_logvar_b = read('LFADS_con_con_to_post_co_logvar_b:0');
            
            mtp.gen_gengru_x_to_ru_W = read('LFADS_gen_GenGRU_Gates_x_2_ru_W:0');
            mtp.gen_gengru_h_to_ru_W = read('LFADS_gen_GenGRU_Gates_h_2_ru_W:0');
            mtp.gen_gengru_h_to_ru_b = read('LFADS_gen_GenGRU_Gates_h_2_ru_b:0');
            mtp.gen_gengru_x_to_c_W = read('LFADS_gen_GenGRU_Candidate_x_2_c_W:0');
            mtp.gen_gengru_rh_to_c_W = read('LFADS_gen_GenGRU_Candidate_rh_2_c_W:0');
            mtp.gen_gengru_rh_to_c_b = read('LFADS_gen_GenGRU_Candidate_rh_2_c_b:0');

            mtp.gen_to_factors_W = read('LFADS_gen_gen_2_fac_W:0');
            
            if verbose
                unread = setdiff(keys, readList);
                unread = unread(~contains(unread, 'Adam'));
                if ~isempty(unread)
                    fprintf('Unread keys:\n');
                    disp(unread);
                end
            end
        end
    end
    
    methods(Access = protected)
      function groups = getPropertyGroups(obj)
         if isscalar(obj)
            groupTitle = 'Read-in from spikes to input factors';
            propList = {'x_to_infac_W', 'x_to_infac_b'};
            groups(1) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

            groupTitle = 'Initial condition encoder (forward)';
            propList = {'ic_enc_fwd_t0', 'ic_enc_fwd_gru_xh_to_gates_ru_W', 'ic_enc_fwd_gru_xh_to_gates_ru_b', 'ic_enc_fwd_gru_xrh_to_c_W', 'ic_enc_fwd_gru_xrh_to_c_b'};
            groups(2) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

            groupTitle = 'Initial condition encoder (reverse)';
            propList = {'ic_enc_rev_t0', 'ic_enc_rev_gru_xh_to_gates_ru_W', 'ic_enc_rev_gru_xh_to_gates_ru_b', 'ic_enc_rev_gru_xrh_to_c_W', 'ic_enc_rev_gru_xrh_to_c_b'};
            groups(3) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

            groupTitle = 'Initial condition g0';
            propList = {'prior_g0_mean', 'prior_g0_logvar', 'ic_enc_to_posterior_g0_mean_W', 'ic_enc_to_posterior_g0_mean_b', 'ic_enc_to_posterior_g0_logvar_W', 'ic_enc_to_posterior_g0_logvar_b', 'g0_to_gen_ic_W', 'g0_to_gen_ic_b'};
            groups(4) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

            groupTitle = 'Controller encoder (forward)';
            propList = {'ci_enc_fwd_t0', 'ci_enc_fwd_gru_xh_to_ru_W', 'ci_enc_fwd_gru_xh_to_ru_b', 'ci_enc_fwd_gru_xrh_to_c_W', 'ci_enc_fwd_gru_xrh_to_c_b'};
            groups(5) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

            groupTitle = 'Controller encoder (reverse)';
            propList = {'ci_enc_rev_t0', 'ci_enc_rev_gru_xh_to_ru_W', 'ci_enc_rev_gru_xh_to_ru_b', 'ci_enc_rev_gru_xrh_to_c_W', 'ci_enc_rev_gru_xrh_to_c_b'};
            groups(6) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

            groupTitle = 'Controlller RNN';
            propList = {'con_gengru_x_to_ru_W', 'con_gengru_h_to_ru_W', 'con_gengru_h_to_ru_b', 'con_gengru_x_to_c_W', 'con_gengru_rh_to_c_W', 'con_gengru_rh_to_c_b'};
            groups(7) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

            groupTitle = 'Controller output co';
            propList = {'prior_ar1_logevars', 'prior_ar1_logatau', 'con_co', 'con_to_posterior_co_mean_W', 'con_to_posterior_co_mean_b', 'con_to_posterior_co_logvar_W', 'con_to_posterior_co_logvar_b'};
            groups(8) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

            groupTitle = 'Generator RNN';
            propList = {'gen_gengru_x_to_ru_W', 'gen_gengru_h_to_ru_W', 'gen_gengru_h_to_ru_b', 'gen_gengru_x_to_c_W', 'gen_gengru_rh_to_c_W', 'gen_gengru_rh_to_c_b'};
            groups(9) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

            groupTitle = 'Generator output';
            propList = {'gen_to_factors_W', 'factors_to_logrates_W', 'factors_to_logrates_b'};
            groups(10) = matlab.mixin.util.PropertyGroup(propList,groupTitle);

         else
            % Nonscalar case: call superclass method
            groups = getPropertyGroups@matlab.mixin.CustomDisplay(obj);
         end
      end
   end
end