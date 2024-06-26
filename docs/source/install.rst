Installation
=============

Download
-----------

To obtain the code you can clone the git repository:

.. code-block:: bash

    git clone https://github.com/henriquemiranda/yambopy.git

Or download the `zip` file from:

.. code-block:: bash

    wget https://github.com/henriquemiranda/yambo-py/archive/master.zip

Install
--------

`yambopy` uses distutils for the instalation. To install it run:

.. code-block:: bash

    sudo python setup.py install

If you do not have root permisisons (when you want to install in your cluster for example):

.. code-block:: bash

    python setup.py install --user

Another option is to use developer installation:

.. code-block:: bash

    sudo python setup.py develop
