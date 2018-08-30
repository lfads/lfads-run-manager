classdef Run < LFADS.Run
    properties(SetAccess=protected)
        twoStageAlignmentTool
    end

    methods
        function v = generateExtraLFADSInputsByDataset(r, seqData, regenerate, maskDatasetsGenerate) %#ok<INUSD>
            % override this method to include extra fields in the LFADS input h5 files (e.g. for modified versions of LFADS)
            if r.params.c_do_subpop_readin || r.params.c_do_subpop_readout
                % generate two stage alignment matrices
                if regenerate || isempty(r.twoStageAlignmentTool) || isempty(r.twoStageAlignmentTool.readinMatrices_cxs)
                    r.setupTwoStageAlignment(seqData);
                end
                tool = r.twoStageAlignmentTool;
                v = struct();
                v.readin_matrix_cxs = tool.readinMatrices_cxs;
                v.readin_matrix_cxs_mask = tool.readinMatrices_cxs_mask;
                v.readin_subfactor_subpopulation_ids = tool.readin_subfactor_subpopulation_ids;
                v.bias_c = tool.biases_c;
                v.out_bias_c = tool.out_biases_c;
                v.readin_matrix_sxf = tool.readinMatrices_sxf;

                v.readout_matrix_sxc = tool.readoutMatrices_sxc;
                v.readout_matrix_sxc_mask = tool.readoutMatrices_sxc_mask;
            else
                v = struct([]);
            end
        end

        function generateExtraLFADSInputFiles(r, seqData, regenerate, maskGenerate)
            % override this moethods to generate any additional LFADS input files
            % this will be called at the end of makeLFADSInput()
            if ~any(maskGenerate) && ~regenerate
                return;
            end
            if r.params.c_do_subpop_readin || r.params.c_do_subpop_readout
                % generate two stage alignment matrices
                if isempty(r.twoStageAlignmentTool) || isempty(r.twoStageAlignmentTool.readinMatrices_cxs)
                    r.setupTwoStageAlignment(seqData);
                end
                tool = r.twoStageAlignmentTool;
                shared_readoutMatrix_fxs = tool.shared_readoutMatrix_fxs;

                outfile = fullfile(GetFullPath(r.pathLFADSInput), 'shared_data.h5');
                if exist(outfile, 'file')
                    delete(outfile);
                end
                %% permute to deal with matlab v python (column-major v row-major)
                ndim = numel(size(shared_readoutMatrix_fxs));
                shared_readoutMatrix_fxs = permute(shared_readoutMatrix_fxs,[ndim:-1:1]);
                h5create(outfile, '/readout_matrix_fxs', size(shared_readoutMatrix_fxs), 'Datatype','double');
                h5write(outfile, '/readout_matrix_fxs', shared_readoutMatrix_fxs);
            end
        end
        
        function tool = setupTwoStageAlignment(r, seqData)
            if nargin < 2
                seqData = r.loadSequenceData();
            end
            r.twoStageAlignmentTool = TwoStageLFADS.TwoStageMultisessionAlignmentTool(r, seqData);
            r.twoStageAlignmentTool.computeAlignmentMatrices();
            tool = r.twoStageAlignmentTool;
        end

    end

    methods % needed for LFADS prep
        function r = Run(varargin)
           r@LFADS.Run(varargin{:});
        end
    end
end
