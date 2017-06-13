Parameters files
================

parameters.yaml
---------------
Defines paths for both local and remote binaries and databases. A template is provided in the ``examples`` directory.


step.yaml
---------
Defines the steps that the pipeline will execute. A template is provided in the ``/examples`` directory.

Step names correspond to a python module that will launch the step. Step names are split based on the '_' character so you can launch multiple instance. For example you might want to launch blastx and blastn, so step names could be 'Blast_N' and 'Blast_X'. What is after the underscore do not matters, it is just used to differanciate the two steps.

Special words in bracket are used as substitution string.
- (file), (file1) and (file2)
- (SampleID)
- (library)

.. literalinclude:: ../../examples/step.yaml


.. _map.txt:

map.txt
-------
The map file describe the experiment. It is a tabulated file with the first line containing headers starting with '#'. It must contain at least two column: SampleID and file.
A template is provided in the ``examples`` directory.
This is a minimum map.txt file:

.. literalinclude:: ../../examples/map.txt

You can add categories for each sample so they can be used when coloring sequences in trees from the Rps2tree module.
