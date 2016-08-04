# Author:
# Tests for the yambopy library
# 
#
import unittest
import sys
import os
import argparse
import subprocess
import filecmp
from yamboparser import YamboFile, YamboFolder

folder = os.path.dirname(os.path.realpath(__file__))+'/testdata/t1_errors_warnings'


fold = YamboFolder(folder)
fold.yambofiles


fl = YamboFile('o-GW_run.10.720.qp',folder)
print 
print fl.type
print fl.data.keys()

fl = YamboFile('l-GW_run.8.480_em1d_ppa_HF_and_locXC_gw0_rim_cut_CPU_1',folder)
print 
print "type", fl.type
print "data", fl.data 
print "kpts", fl.kpoints
print "warn", fl.warnings
print "erro", fl.errors

fl = YamboFile('r-GW_run.8.480_em1d_ppa_HF_and_locXC_gw0_rim_cut',folder)
print 
print fl.type
print fl.kpoints
print fl.data




folder = os.path.dirname(os.path.realpath(__file__))+'/testdata/t2_parse_qps'

fl = YamboFile('o-yambo.qp',folder)
print 
print fl.type
print fl.data.keys()

fl = YamboFile('l-yambo_em1d_HF_and_locXC_gw0',folder)
print 
print "type", fl.type
print "data", fl.data 
print "kpts", fl.kpoints
print "warn", fl.warnings
print "erro", fl.errors

fl = YamboFile('r-yambo_em1d_HF_and_locXC_gw0',folder)
print 
print "type", fl.type
print "kpts", fl.kpoints
print "data", fl.data
