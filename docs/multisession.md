# Multi-dataset Stitching Models

If you specify multiple datasets to be included in an LFADS run by selecting multiple datasets in a `RunSpec`, the resulting model will stitch together the multiple datasets. The concept is to generate the spiking data in all of the included datasets using the same encoder and generator RNNs, but to interface to the separate neural datasets through _read-in_ and _readout_ alignment matrices.

Below is a schematic of the readout side. Here, the generator RNN and readout from generator units to factors is the same for all datasets. Therefore, one intends that the factor trajectories would be similar for similar trials / conditions across the datasets. Going from factors to rates, however, the recorded neurons are, in general, not the same across datasets, and the cardinality may differ. Thus, dataset-specific readout matrices are used to combine the factors to produce each of the recorded neurons' rates on each dataset.

![Multi-session stitching schematic](images/stitching_schematic.png)

A similar set of dataset specific read-in matrices are used to connect the spiking data to the encoder RNN in order to produce initial conditions and inferred inputs for each trial.

## Generating alignment matrices

These read-in and readout alignment matrices are learned from the data along with the other parameters. However, it's useful to seed the alignment matrices with an initial guess that suggests the correspondence between the datasets. If you have multiple datasets in a `RunSpec`, and the hyperparameter `useAlignmentMatrix` is set to `true` in the `RunParams`, then `lfads-run-manager` will automatically generate read-in alignment matrices from your data using a principal components regression algorithm that proceeds as:

* Generate condition-averaged firing rates for each neuron for each condition for each dataset
* Concatenate all of neurons from all datasets together to build a matrix which is (`nTime * nConditions`) `nNeuronsTotal`
* Perform PCA on this matrix and keep the projections of the data into the top `nFactors` components. These represent the global shared structure of the data across all datasets.
* For each dataset individually, regress these projection scores onto the condition-averaged rates from that dataset alone. The regression coefficients thus transform from that dataset's neurons to the global shared structure, and consequently, we take this matrix of regression coefficients as the readout matrix.

These matrices will be computed for you automatically by `run.prepareForLFADS()` and exported in the LFADS input folder. LFADS will generate an initial guess for the readout alignment matrix, which transforms from common factors back to dataset-specific rates, using the pseudo-inverse of the read-in alignment matrix computed by `lfads-run-manager`.

!!! tip "Alignment biases"
    In addition to this alignment read-in matrix, there is also an alignment bias vector which will be added to each neuron's counts befor projecting through the matrix. Consequently, `lfads-run-manager` seeds this bias with the negative mean of the rates of each neuron.


## Verifying the alignment matrices

To visualize how well these initial alignment matrices are working, we can compare the common global PCs from all datasets against the projection of each dataset through the read-in matrices. That is, we can plot the regression target (global PCs) against the best possible reconstruction from each dataset.

Under the hood, the alignment matrix calculations are performed by an instance of `LFADS.MutlisessionAlignmentTool`. To plot the reconstruction quality, you can call `tool.plotAlignmentReconstruction(numberOrIndicesOfFactorsToPlot, numberOrIndicesOfConditionsToPlot)`, like so:

```matlab
run.doMultisessionAlignment();
tool = run.multisessionAlignmentTool;

nFactorsPlot = 3;
conditionsToPlot = [1 20 40];
tool.plotAlignmentReconstruction(nFactorsPlot, conditionsToPlot);
```

![Alignment matrix reconstruction quality](images/alignment_reconstruction.png)

In this example, the single-dataset predictions look quite similar to the global target, especially in the first 2 principal components which capture most of the variance.

The actual alignment matrices can be accessed using:
```matlab
tool.alignmentMatrices % nDatasets x 1 cell array of read-in matrices
```
