Options controlling the overall behaviour
=========================================


verbose
-------

If set to 'YES' or 'True', will run in verbose mode. 'NO' or 'False' will print minimal information. 

debug
-----

If set to 'YES' or 'True', will print additional information to help locating problems

no_opt
------

If set to 'YES' or 'True', will stop execution before doing the inversion. This option is usually run together with the save option

save
----

If save is a string, then all model information is saved as a pickle dump of model, with the observation normal matrix (N) and and observation normal vector (Nd) separately saved as numpy npy files. This option is used when the same observations and same green functions will be used for different regularizations.

print_result
------------

If set to 'False', will only write a pickle dump of model after the inversion. This option is useful when runs are performed on a remote machine and results analysis on a local machine. All results information can then be generated using pyek_print_results.py. If set to 'True', will generate all the results information.

plot
----

If set to 'False', no plot is generated.

tar
---

performs a compress tarball of the result directory. 



  


 