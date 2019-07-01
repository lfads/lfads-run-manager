# Installation

These instructions will walk you through the basic setup process to get you up and running with LFADS.

!!! warning "Use Python 2.7"
    While TensorFlow fully supports Python 3, the LFADS code itself does not yet. We expect to fix the few incompatibilities soon, but for now, use Python 2.7.

## Install TensorFlow

You'll need to install TensorFlow to run LFADS. Please note that LFADS will run much faster on a GPU, and **Tensorflow on GPU is not supported on MacOS**.

!!! tip "Recommended installation with Anaconda"
    The simplest way to get up and running is to install [Anaconda](https://www.anaconda.com/download/) and then to install Tensorflow in its own conda environment (here named `tensorflow`) using:

    ```bash
    conda create --name tensorflow python=2.7 tensorflow-gpu h5py matplotlib
    ```

    This approach was suggested by Harveen Singh [here](https://towardsdatascience.com/tensorflow-gpu-installation-made-easy-use-conda-instead-of-pip-52e5249374bc) and should take care of installing the compatible version of the NVIDIA dependencies, including the CUDA toolkit and cuDNN.

Tensorflow updates quickly, and the recommended approach is to install the latest version of Tensorflow (which will be used by the `tensorflow-gpu` conda package). If you encounter compatibility issues with the LFADS code, please let us know by [filing an issue on Github](https://github.com/lfads/lfads-run-manager/issues).

### Alternative manual installation

To install Tensorflow manually along with the dependencies, follow the [documentation for installing Tensorflow](https://www.tensorflow.org/install/) and be sure to install the version for GPUs if you wish to take advantage of the LFADS run queue. You may wish to install everything in a Python `virtualenv` or inside a `conda` environment, both of which are supported by `lfads-run-manager`.

### Test Tensorflow install
Test your installation by trying to `import tensorflow` in Python. If you're using a conda environment, be sure to activate it first:

```bash
source activate tensorflow
```

```python
# Python
import tensorflow as tf
hello = tf.constant('Hello, TensorFlow!')
sess = tf.Session()
print(sess.run(hello))
```

Which should output:
```
Hello, TensorFlow!
```

If you're getting errors, check this [helpful list of common error messages](https://www.tensorflow.org/install/errors).

## Install LFADS

You'll then need to clone the [Tensorflow models repo containing LFADS](https://github.com/lfads/models/tree/master/research/lfads) somewhere convenient on your system.

```bash
git clone https://github.com/lfads/models.git
```

Then add this LFADS folder both to your `PYTHONPATH` and system `PATH`. The LFADS Run Manager src folder should also be added to your `PYTHONPATH`. Add the following to your `.bashrc`:
```bash
export PYTHONPATH=$PYTHONPATH:/path/to/models/research/lfads/:/path/to/lfads-run-manager/src
export PATH=$PATH:/path/to/models/research/lfads/
```

Ensure that typing `which run_lfads.py` at your terminal prompt shows the path to `run_lfads.py`.

LFADS depends on the Python libraries `h5py` and `matplotlib` being installed as well:

```bash
conda install h5py matplotlib
```

## Install tmux

LFADS Run Manager uses tmux to run LFADS within to enable queuing many runs across the available GPUs and to facilitate online monitoring. Fortunately, installing tmux is pretty straightforward on most distributions.

* Ubuntu: `sudo apt-get install tmux`
* Mac: `brew install tmux` using [Homebrew](https://brew.sh).

There are many nice guides to using `tmux`:

* [A Quick and Easy Guide to tmux - Ham Vocke](http://www.hamvocke.com/blog/a-quick-and-easy-guide-to-tmux/)
* [A tmux Crash Course - Josh Clayton](https://robots.thoughtbot.com/a-tmux-crash-course)


## Install LFADS Run Manager

Finally, clone the lfads-run-manager repository somewhere convenient on your system.

```bash
git clone https://github.com/djoshea/lfads-run-manager.git
```

You'll need to have Matlab installed. Then you can add the root folder of the lfads-run-manager to your Matlab path, either using `pathtool` or by running:

```matlab
addpath('/path/to/lfads-run-manager/src')
```

No need to add the subfolders to the path recursively.

## Having issues?   

If you're having issues, please let us know by [filing an issue on Github](https://github.com/lfads/lfads-run-manager/issues). Unless you have a scientific question or question about the paper, Github issues works better than emailing us directly as other people can use the thread as a resource in the future, as well as creating a to-do list for us to ensure that everything get fixed. Thanks!
