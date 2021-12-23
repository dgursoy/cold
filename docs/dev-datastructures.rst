***************
Data structures
***************

Measurement data
================

.. code-block:: python

    import cold 
    path = 'myfolder/'
    data = cold.read(path)
    print (data.dtype, data.shape)


Coded-apertures
================

.. code-block:: python

    import cold 
    path = 'myfolder/mymask.npy'
    mask = cold.mask(path)
    print (mask.dtype, mask.shape)