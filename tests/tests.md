Tests
=====

There are two types of tests:
    - Integration tests
    - Unit tests

Integration tests
==========================
The integration tests are in the folder tests in the root of the git repository.

To run the basic tests:

    :::bash
    python test_si.py -i

This will extract a reference database, generate input files and compare them to a reference.

To run the advanced tests:

    :::bash
    python test_si.py -f

This will generate the databases using Quantum espresso and p2y, generate the input files, run the calculations and plot the result.
For this test you need to have Yambo and Quantum Espresso installed. 

To clean the folder:

    :::bash
    python test_si.py -c

Unit tests
============================
The unit tests test specific routines and can always be found in a folder
called tests in the same folder as the python source code.

To run all the tests use pytest.

Packages
============================
To obtain the coverage when running pytest we use 'pytest-cov':

    :::python
    pip install pytest-cov
