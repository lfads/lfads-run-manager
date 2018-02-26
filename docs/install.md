# Installation

## Install TensorFlow

You'll need to install TensorFlow to run LFADS. Follow the [documentation for installing Tensorflow](https://www.tensorflow.org/install/) and be sure to install the version for GPUs if you wish to take advantage of the LFADS run queue. You may wish to install everything in a Python `virtualenv` or inside a `conda` environment, both of which are supported by `lfads-run-manager`. Be sure to test that you can `import tensorflow` in Python correctly:

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

## Install LFADS

You'll then need to clone the [Tensorflow models repo containing LFADS](https://github.com/tensorflow/models/tree/master/research/lfads) somewhere convenient on your system.

```bash
git clone https://github.com/tensorflow/models.git
```

Then add this LFADS folder both to your `PYTHONPATH` and system `PATH`. Add the following to your `.bashrc`:
```bash
export PYTHONPATH=$PYTHONPATH:/path/to/models/research/lfads/
export PATH=$PATH:/path/to/models/research/lfads/
```

Ensure that typing `which run_lfads.py` at your terminal prompt shows the path to `run_lfads.py`.

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
