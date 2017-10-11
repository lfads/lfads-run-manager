# Key Concepts

The run manager defines a small set of key classes that encapsulate specific concepts within LFADS and enable you to organize datasets and LFADS models effectively within Matlab.

## LFADS.Dataset

A `Dataset` instance represents a collection of trials with associated neural spiking channels. One or more datasets will be used by LFADS to train and evaluate the model. An individual dataset would include simultaneously recorded neural signals collected during one day or one experimental session.

## LFADS.DatasetCollection

A `DatasetCollection` is a set or array of one or more datasets.

## LFADS.RunSpec

A `RunSpec`, short for run specification, defines which of the datasets within a dataset collection will be included as input to a specific LFADS model. Typically, only one dataset is specified. If multiple datasets are specified, the resulting LFADS model will be a _stitched_ model which uses alignment matrices. Stitched models share a common generator to generate spiking data collected in different experimental sessions. Refer to the LFADS paper for more information on stitched models.

## LFADS.RunParams

A `RunParams` encapsulates the hyperparameters of an LFADS run. Most of these hyperparameters are fed directly to the Python+Tensorflow LFADS code and are defined in the LFADS paper. Examples are the size of the generator RNN and the dropout probability during training. Another key parameter is the bin width used to convert spike times into time-varying spike rates.

When adapting the run manager to work with your datasets, you are encouraged to include your own hyperparameters that can be used to specify the way data is extracted and processed from your datasets. For example, you might wish to define a `timeWindowPre` and `timeWindowPost` that specify the window of time from each trial in which neural spiking data is extracted. Or you might wish to define hyperparameters that affect which trials are included, e.g. `keepSuccessTrialsOnly` or `includePerturbationTrials`.

The advantage of including your dataset-specific hyperparameters in your `RunParams` subclass is that the values of these fields will then affect the hash value that is used to uniquely define individual LFADS runs on disk, enabling you to easily compare across sweeps of these hyperparameter settings just as you would with the built-in LFADS hyperparameters.

## LFADS.Run

A `Run` encapsulates an actual LFADS model that will be trained using Python+Tensorflow. An `Run` is defined by the combination of an `RunSpec` instance (which specifies the datasets included) and an `RunParams` instance (which specifies the hyperparameters). Each `Run` will be associated with a run of the Python+Tensorflow code that defines and trains the LFADS model.

## LFADS.RunCollection

A `RunCollection` is a set of one or more `LFADS.Run`s. This collection is organized as a two-dimensional matrix of runs.

The first dimension of this matrix is specified by an array of `LFADS.RunSpec` instances. This enables different datasets or sets of datasets to be used within each model. For example, if you had 10 datasets, you could run LFADS on each dataset individually by having 10 `LFADS.RunSpec`s, each specifying an individual dataset to be included.

The second dimension of this matrix is specified by an array of `LFADS.RunParams` instances. This enables you to vary hyperparameter settings across the runs.

Each cell of this matrix, defined by a particular `RunSpec` and `RunParams` combination defines a specific `Run` which can then be generated and trained using Tensorflow.
