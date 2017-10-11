# File organization used by LFADS Run Manager

You won't likely need to dive into the raw files produced by lfads-run-manager and LFADS, but we include a brief overview here to help understand what's stored where. You can safely skip to [running your model](running.md).

## Run collection organization

After running `MyExperiment.drive_script` and calling `rc.prepareForLFADS()`, you'll see the following directory tree on your hard drive:

```bash
$ tree -L 4 ~/lorenz_example/
.
├── datasets
│   ├── dataset001.mat
│   ├── dataset002.mat
│   ├── dataset003.mat
│   └── dataset004.mat
└── runs
    └── exampleRun
        ├── data_IR3OQV
        │   ├── inputInfo_dataset001.mat
        │   ├── inputInfo_dataset002.mat
        │   ├── inputInfo_dataset003.mat
        │   ├── lfads_dataset001.h5
        │   ├── lfads_dataset002.h5
        │   └── lfads_dataset003.h5
        ├── launch_tensorboard.sh
        ├── param_pqQbzB
        │   ├── all
        │   ├── single_dataset001
        │   ├── single_dataset002
        │   └── single_dataset003
        ├── run_lfadsqueue.py
        └── summary.txt

9 directories, 13 files
```

**`datasets`**:
:    Inside live the raw datasets, unrelated to lfads-run-manager. These were created by you when you called `LFADS.Utils.generateDemoDatasets(...)`.

**`runs`**:
:   The root runs folder you specified in the constructor to `RunCollection`.
    ```matlab
    runRoot = '~/lorenz_example/runs';
    rc = MyExperiment.RunCollection(runRoot, 'exampleRun', dc);
    ```

    **`exampleRun`**:
    :   The location of this specific `RunCollection`, based on the name you passed to constructor of `RunCollection`
        ```matlab
        rc = MyExperiment.RunCollection(runRoot, 'exampleRun', dc);
        ```

        **`data_IR3OQV`**:
        :   The location of the exported datasets for `RunParams` whose data hash is `IR3OQV`. The data hash includes properties of `RunParams` that affect the exported data, as described [here](drive.md#runparams-data-and-param-hashes).

            **`inputInfo_datasetName.mat`**
            :   Contains data collected when generating the LFADS input, including the raw spike counts, condition ids, time vector, and trial indices assigned to the training and validation sets.

            **`lfads_datasetName.h5`**:
            :   The spike counts data directly read by LFADS

        **`param_pqQbzB`**:
            :   Location of the individual runs generated with the `RunParams` instance whose param hash is `pqQbzB`. The subfolders correspond to the run names passed to `RunSpec`, and their contents will be discussed below.
                ```matlab
                rc.addRunSpec(MyExperiment.RunSpec('all', dc, 1:dc.nDatasets));
                ```

        **`launch_tensorboard.sh`**:
        :   Shell script which will launch TensorBoard, optionally on a specific port
            ```bash
            ./launch_tensorboard.sh 50000
            ```

        **`run_lfadsqueue.py`**:
        :   Python script which will launch the LFADS Run Queue on all runs within the `exampleRun` `RunCollection`:
            ```bash
            python run_lfadsqueue.py
            ```

## Individual run folders

Within an run folder, we find:

```bash
$ tree ~/lorenz_example/runs/exampleRuns/param_pqQbzB/all
.
├── lfads.done
├── lfads.out
├── lfadsInput
│   ├── inputInfo_dataset001.mat -> ../../../data_IR3OQV/inputInfo_dataset001.mat
│   ├── inputInfo_dataset002.mat -> ../../../data_IR3OQV/inputInfo_dataset002.mat
│   ├── inputInfo_dataset003.mat -> ../../../data_IR3OQV/inputInfo_dataset003.mat
│   ├── lfads_dataset001.h5 -> ../../../data_IR3OQV/lfads_dataset001.h5
│   ├── lfads_dataset002.h5 -> ../../../data_IR3OQV/lfads_dataset002.h5
│   └── lfads_dataset003.h5 -> ../../../data_IR3OQV/lfads_dataset003.h5
├── lfadsOutput
│   ├── checkpoint
│   ├── checkpoint_lve
│   ├── fitlog.csv
│   ├── hyperparameters-0.txt
│   ├── hyperparameters-38740.txt
│   ├── lfads_log
│   │   └── events.out.tfevents.1507613307.photon
│   ├── lfads_vae.ckpt-37206.data-00000-of-00001
│   ├── lfads_vae.ckpt-37206.index
│   ├── lfads_vae.ckpt-37206.meta
│   ├── ...
│   ├── model_params
│   ├── model_runs_dataset001.h5_train_posterior_sample_and_average
│   ├── model_runs_dataset001.h5_valid_posterior_sample_and_average
│   ├── model_runs_dataset002.h5_train_posterior_sample_and_average
│   ├── model_runs_dataset002.h5_valid_posterior_sample_and_average
│   ├── model_runs_dataset003.h5_train_posterior_sample_and_average
│   └── model_runs_dataset003.h5_valid_posterior_sample_and_average
└── lfads_train.sh
```

**`lfadsInput`**:
:    Contains relative symbolic links to the contents of `data_IR3OQV`, enabling multiple runs to share data without duplication. The `.h5` files will be read in by LFADS.

**`lfadsOutput`**:
:   The directory to which LFADS will write generated output. Some of the key files within are:

    **`checkpoint_lve`**:
    :   Contains the saved checkpoint with the lowest validation error.

    **`fitlog.csv`**:
    :   Contains information about the various costs through training

    **`hyperparameters-0.txt`**:
    :   Records the hyperparameters used by LFADS

    **`model_runs_datasetName.h5_train_posterior_sample_and_average`**:
    :   Contains the posterior mean samples and averages for the training trials

    **`model_runs_datasetName.h5_valid_posterior_sample_and_average`**:
    :   Contains the posterior mean samples and averages for the validation trials

**`lfads_train.sh`**:
:   Shell script which will launch LFADS to train the model. This script may potentially chain performing posterior mean sampling and writing the model parameters, depending on how it was generated by lfads-run-manager.

**`lfads.done`**:
:   Empty text file indicating to the LFADS Run Queue that this model has already completed successfully

**`lfads.out`**:
:   Logged output of LFADS generated by the LFADS Run Queue

## Clearing LFADS Output

If you wanted to re-train a model from scratch, you can call `run.deleteLFADSOutput()` from Matlab, or you could manually delete the `lfadsOutput` folder, `lfads.done`, and `lfads.out`.
