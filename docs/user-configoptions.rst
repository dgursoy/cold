****************************
Runtime cofiguration options
****************************

Data readout options
====================

.. code-block:: yaml

    file:
        path: myfolder'
        range: [0, 100]
        threshold: 10 
        frame: [0, 2048, 0, 2048] 
        ext: 'h5'
        chunks: 8
        type: 'stacked'
        h5:
            key: '/entry1/data/data'


Geometry options
================

Aperture
--------

.. code-block:: yaml

    geo:
        mask: 
            material: 'Au'
            path: 'codes/code-debruijn-2-8-000.npy'
            pad: 300
            bitsizes: [15, 7.5] # [mu]
            resolution: 0.5
            thickness: 7.5 # [mu]
            smoothness: 0 # [mu]
            widening: 2.5 # [mu]
            dist: 1.16 # [mu]
            tiltx: 5 # 
            cenx: 290 # 
            tilty: 18 # 
            ceny: 100 #
            step: 1 # [mu]
            calibrate: 
                dist: [1.15, 1.25, 0.01] # [mm]
                tiltx: [0, 8, 1]  # -6
                cenx: [0, 0, 20] # -180
                tilty: [0, 23, 1] # 23
                ceny: [-100, 100, 20] 

Detector
--------

.. code-block:: yaml

    geo:
        detector: 
            shape: [2048, 2048] # [pixels]
            size: [409.6, 409.6] # [mm]
            rot: [-1.20139958, -1.21416739, -1.21878591] # [radian]
            pos: [28.871, 2.786, 513.140] # [mm]


Source
------

.. code-block:: yaml

    geo:
        source: 
            grid: [-1.365, -1.156, 0.001] # [mm]


Algorithm options
=================

Decoding
--------

.. code-block:: yaml

    algo:
        pos: 
            method: 'lsqr' 
            init: 0
        sig: 
            method: 'nnls' 
            init:
                maxsize: 120 # [mu]
                avgsize: 20 # [mu]
                atol: 4


Peak searching
--------------

Indexing
--------