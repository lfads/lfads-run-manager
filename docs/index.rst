.. LFADS Run Manager documentation master file, created by
   sphinx-quickstart on Sun Dec 18 12:57:00 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LFADS Run Manager documentation
=============================================

Each of these classes can be subclassed by you to meet the needs of your dataset or scientific goals. Two of these classes, :ref:`LFADS_Run` and LFADS_Dataset are abstract classes, which means you *must* subclass them and fill in the missing abstract methods in order to use them. Besides these mandatory methods which you must write yourself, you may find it helpful to add additional functionality to the classes themselves, e.g. to facilitate screening of datasets or post-hoc analyses of LFADS runs.

LFADS Run Manager Classes
-------------------------

.. toctree::
   :maxdepth: 2

   LFADS_DatasetCollection
   LFADS_Dataset
   LFADS_RunSpec
   LFADS_Run
   LFADS_RunCollection
   LFADS_RunParams
   LFADS_PosteriorMeans

Index and Search
-----------------

* :ref:`genindex`
* :ref:`search`
