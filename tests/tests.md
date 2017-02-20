Tests
=====

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
