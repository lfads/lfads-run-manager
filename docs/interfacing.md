# Using LFADS Run Manager with your datasets

## Copying the `LorenzExperiment` working example code
Below we describe how to use the run manager code with datasets from a specific experiment. The recommended way to begin this process is to copy the folder `+LorenzExperiment` inside the` lfads-run-manager` repository to some other folder on your Matlab path, and then to rename it to something related to the experiment. Below, we'll stick with the name `LorenzExperiment`.

Each of the classes you have just created are defined to inherit from the corresponding `LFADS.ClassName` inside the `lfads-run-manager` repo. Consequently, only a small amount of code is present in each file; the rest of the properties and methods for each class are define inside the `+LFADS` folder in the repo.

## Editing the core classes

Here we walk through each of the classes that you've just copied. Most of the classes can be left as is to get started, but you may find it helpful to add utility methods and addtional metadata in certain locations. However, the only required edits are:

* `loadData` in `Dataset.m` - specify how to load a dataset from disk. The default implementation assumes that the data live in a `.mat` file that can be loaded using `load`.
* `generateCountsForDataset` in `Run.m`  - preprocess data and perform spike binning

### Editing `DatasetCollection.m` _(Optional)_

Edit the file `+LorenzExperiment/DatasetCollection.m`. Recall that a dataset collection refers to a set of multiple individual datasets. Note the definition of the constructor:

```matlab
function ds = DatasetCollection(path)
    ds = ds@LFADS.DatasetCollection(path);
end
```

You may edit this to fit your needs, but the default approach is to create a new dataset collection by specifying a path on disk where the data live. For example, you could run:

```matlab
dc = LorenzExperiment.DatasetCollection('/path/to/experimentData');
```

This path will then be used as the parent folder by all of the datasets that are added to this collection.

Below is a function included as an example of how to filter or down-select datasets within a collection. A typical approach might be to add all of the datasets that were collected, and then filter by those having a sufficiently high trial count (or satisfying some other set of criteria). You can use the utility function `filterDatasets` to specify the indices or mask over datasets to keep.

```matlab
function filterHavingMinimumTrials(dc, minTrials)
    % example of a function that will filter down datasets based on
    % their metadata.
    nTrials = cat(1, dc.datasets.nTrials);

    % filterDatasets is provided by DatasetCollection
    dc.filterDatasets(nTrials >= minTrials);
end
```

**No edits are necessary to `DatasetCollection.m` to get up and running**, but feel free to add any additional methods or properties as needed for your application.

### Editing `Dataset.m` _(Required)_

Edit the file `+LorenzExperiment/Dataset.m`. Recall that a dataset encapsulates a collection of trials with simultaneously recorded neural data from an individual experimental session. Here, we will make a few light edits to specify metadata about each dataset.

First, look at the constructor.

```matlab
function ds = Dataset(collection, relPath)
    ds = ds@LFADS.Dataset(collection, relPath);
end
```

In order to encapsulate a particular dataset on disk, you will create a new `LorenzExperiment.Dataset` instance in Matlab. The first argument `collection` is the DatasetCollection to add this dataset to, which will provide the parent path. The second argument `relPath` specifies the path to this dataset relative to the collection. For example, if the dataset were stored in `/path/to/experimentalData/dataset001.mat`, you might run:

```matlab
ds1 = LorenzExperiment.Dataset(dc, 'dataset001.mat');
```

You may need to specify how to load the actual data into Matlab in order to facilitate preprocessing. The default simply calls Matlab's `load` method and assumes that `ds.path` points to a `.mat` file. `ds.path` will be equal to the dataset collection path joined to `relPath`. If your data is stored differently, you will need to replace the implementation of `loadData`:

```matlab
function data = loadData(ds)
    % load this dataset's data file from .path
    in = load(ds.path);
    data = in.data;
end
```

You can then specify how to determine certain metadata about the dataset, simply for display and organizational purposes. These metadata will then be assigned into specific properties of the `Dataset` class. The simplest approach is to simply load the data and copy or compute the values from the data. However, if loading the data is expensive, you might store the metadata in a separate file to save time. This implementation is up to you, and you can simply specify empty values for metadata fields you do not care about. Note that the logical property `infoLoaded` can be used to determine if the metadata has already been loaded, since this method will be called several times to ensure the metadata is loaded when needed.

```matlab
function loadInfo(ds)
    % Load this Dataset's metadata if not already loaded

    if ds.infoLoaded, return; end

    % modify this to extract the metadata loaded from the data file
    data = ds.loadData();
    ds.subject = data.subject;
    ds.saveTags = data.saveTags;
    ds.datenum  = data.datenum;
    ds.nChannels = data.nChannels;
    ds.nTrials = numel(data.trials);

    ds.infoLoaded = true;
end
```

The metadata fields you might assign are as follows:

* `subject` - dataset subject or participant name
* `datenum` - a Matlab datenum identifying the collection time of the dataset
* `nChannels` - the number of unique spiking channels recorded in this dataset
* `nTrials` - the number of trials included in this dataset
* `saveTags` - a vector of numbers indicating within-day blocks of trials included

!!! note "Metadata are optional"
    We note that none of these fields is used for subsequent processing, and are defined only for convenience and consistency. Feel free to ignore these, and to add additional fields as properties directly to your `Dataset` class.

### Editing `RunParams.m` _(Optional)_

Edit the file `+LorenzExperiment/RunParams.m`. Recall that `RunParams` encapsulates all of the hyperparameters used by LFADS but can also be used to specify any experiment specific hyperparameters you wish to add.

You can add these additional properties anywhere in the file:
```matlab
classdef RunParams < LFADS.RunParams
   properties
       % Add additional custom parameters here. The default you assign to
       % them will be used when computing the hash value. Any params whose value
       % differs from the default will be included in the hash value, to allow new
       % parameters to be added without invalidating old hashes. So choose
       % the default once and don't change it. If you decide to use another
       % value later by default, override it in the constructor instead.
   end
```

For example, you might add:
```matlab
   properties
       timeWindowPre = 500; % ms before to include
       timeWindowPost = 500; % ms after to include
   end
```
!!! warning "Pick default values carefully"
    **The default values you assign next to each property should be chosen carefully and never changed once added.** The reason for this is that when generating the hash of the hyperparameters (which specifies where LFADS-related files live on disk), each property is compared against this default value. The current value of a particular property is only included in the hashing process if it differs from this default value. This design ensures that it is always safe to add new hyperparameters; previously performed LFADS runs will still have the same hash value and will be assigned the default hyperparameter. However, if you change the default value here, all of the hash values for all previously performed runs will change, which will require directories to be manually renamed on disk and symbolic links to be corrected. If you wish to change the default value that a property takes for new runs, you can change its value in the `RunParams` constructor without affecting the hash. However, you will then want to manually assign this property to its _old value_ in any drive scripts you used to setup previous LFADS runs, in order to correctly specify the hyperparameters used and the corresponding hash values.

**No changes are required to `RunParams.m` to get up and running.**

### Editing `RunCollection.m` _(Optional)_

Edit the file `+LorenzExperiment/RunCollection.m`. Recall that a `RunCollection` specifies a set of individual LFADS runs defined by an array of `RunSpec`s crossed with an array of `RunParams`.

```matlab
classdef RunCollection < LFADS.RunCollection
    % no need to modify anything here, but feel free to add useful methods
    % and properties as useful

    methods
        function rc = RunCollection(varargin)
            rc@LFADS.RunCollection(varargin{:});
        end
    end
end
```

**No changes are required to `RunCollection.m` to get up and running**, but you can add any utility methods to facilitate analysis for your specific application.

### Editing `Run.m` _(Required)_

Edit the file `+LorenzExperiment/Run.m`. Recall that a `Run` represents a specific LFADS model training run. The main function you will need to provide a definition for is `generateCountsForDataset`. This is where you will actually need to process your datasets and return a structure array containing binned spike counts. The function signature looks like this:

```matlab
function out = generateCountsForDataset(r, dataset, mode, varargin)
```

Here, `r` refers to the `LorenzExperiment.Run` instance. It may be particularly helpful to refer to the `RunParams` instance assigned to this run via `r.params`, especially if you have defined any additional hyperparameters that affect the way in which neural data should be extracted, e.g. which trials and what time window are included.

#### Inputs to `generateCountsForDataset`:

**`r`**:
: The `Run` instance. The current `RunParams` instance can be accessed through `r.params`.

**`dataset`**:
: `LorenzExperiment.Dataset` instance that is to be processed. If this is a single-dataset run, this will be the dataset used. If this is a multi-dataset stitched run, `generateCountsForDataset` will be called once for each dataset, one at a time. You might use `dataset.loadData()` to load the actual data, as you defined above.

**`mode`**:
: String that indicates the intended purpose of the output data. You may ignore this and simply return the same sequence struct regardless of the mode, or you may process the data differently according to the context. Currently two modes are defined:

    * `export` - indicates that the sequence data will be exported as the actual input to the Python+Tensorflow LFADS run
    * `alignment` - for multi-dataset stitched runs, indicates that the output data will be used only to construct the alignment matrices that translate between the spiking channels across different datasets. For example, you might wish to include a subset of trials or a different time window for fitting the alignment matrices, but include all trials for the actual LFADS run.

    If you decide you do wish to handle the `alignment` case differently, you will need to override the `usesDifferentDataForAlignment` method in your `Run` class to return `true`, by adding:

    ```matlab
    function tf = usesDifferentDataForAlignment(r)
        tf = true;
    end
    ```

**`varargin`**:
: Currently not being used, but this enables additional arguments to be passed as named-parameter value pairs (e.g. `:::matlab 'paramName', paramValue, ...`) in the future without breaking existing implementations.

#### Outputs to `generateCountsForDataset`:

**`out`**:
: A scalar struct which holds the following fields:

    **`:::matlab .counts`** (Required):
    : A tensor of binned spike counts (not rates) with size `nTrials` x `nChannels` x `nTime`. These should be total counts, not normalized rates, as they will be added together during re-binning.

    **`:::matlab .timeVecMs`** (Optional):
    : A vector of timepoints with length `nTime` in milliseconds associated with each time bin in `counts`. You can start this wherever you like, but timeVecMs(2) - timeVecMs(1) will be treated as the _raw_ spike bin width used when the data are later rebinned to match `r.params.spikeBinMs`. Default is `1:size(counts, 3)`.

    **`:::matlab .conditionId`** (Optional):
    : Vector with length `nTrials` identifying the condition to which each trial belongs. This can either be a cell array of strings or a numeric vector. Default is `[]`.

    **`:::matlab .truth`** (Optional):
    : For synthetic datasets, provides the ground-truth counts for each trial. Same size as `.counts`. Default is `[]`.

    **`:::matlab externalInputs`** (Optional):
    : Specifies the observed, external inputs which will be passed either to the generator directly or to the encoder. Default is `[]`.

    !!! note "A note on bin widths"
        There are two different bin widths in `lfads-run-manager`. First is this `binWidthMs` within `seq`, which is the spike binning that you will do to the data inside `generateCountsForDataset`. **We recommend binning here at 1 ms or the smallest bin width you might wish to use.** Second is the field `spikeBinMs` inside the `RunParams` class. The expectation is that you will bin using a small bin width inside `generateCountsForDataset`, and then **the run manager code will automatically re-bin the data at the larger bin width set by `r.params.spikeBinMs`** for you. However, you are responsible for ensuring that the larger spike bin width is an integer multiple of the smaller bin width, otherwise an error will be generated.
