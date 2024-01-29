import os
from netCDF4 import *
from glob import glob
import argparse

"""
Script to update serial numbers of yambo ndb.* databases in order to import them to new calculations.

Inputs:
 1. --new_serial='path/to/folder/with/new/dbs' [e.g., the new SAVE]
 2. --old_serial='path/to/folder/with/old/dbs' [e.g., an old ndb.em1s]

This script will prompt the user to go through with updating the dbs.
"""

def get_serials(fldr):
    """
    Identify the databases involved in the change and read their serial numbers
    """
    ndbs_tmp = glob('%s/ndb.*'%fldr)
    ns_tmp   = glob('%s/ns.*'%fldr)
    ndbs = []
    for ndb in ndbs_tmp:
        if "fragment" not in ndb: ndbs.append(ndb)
    for ns in ns_tmp:
        if "fragment" not in ns: ndbs.append(ns)

    ndbs_target = []
    ndbs_values = []
    for ndb in ndbs:
        dbs = Dataset(ndb)
        if 'SERIAL_NUMBER' in dbs.variables:
            ndbs_target.append(ndb)
            serial = int(dbs['SERIAL_NUMBER'][0])
            ndbs_values.append(str(serial))
        dbs.close()
    return ndbs_target, ndbs_values

def prompt_user(new,old,check=False):
    """
    Print info for the user and ask if they want to go through with the change
    """
    new_fn,new_sn = get_serials(new)
    old_fn,old_sn = get_serials(old)
    SERIAL_NUMBER = new_sn[0]
    print("Serial numbers found:")
    print("=====================")
    print(">> New databases:")
    if '%s/ndb.gops'%new in new_fn: is_gops_there = True
    else: is_gops_there = False
    are_serials_equal = True
    for sn in new_sn:
        if sn != new_sn[0]: are_serials_equal = False
    if not are_serials_equal:
        if not is_gops_there: raise ValueError('Conflicting serial numbers in newer folder!')
        else: print('Conflicting serial numbers in newere folder. The value of ndb.gops will be used.')
    for i in range(len(new_fn)): 
        print('%s: %s'%(new_fn[i],new_sn[i]))
        if is_gops_there and new_fn[i]=='%s/ndb.gops'%new: SERIAL_NUMBER = new_sn[i]
    print("--------------")
    print(">> Old databases:")
    for i in range(len(old_fn)): print('%s: %s'%(old_fn[i],old_sn[i]))
    print("--------------")
    print("=====================")
    if not check:
        usr_inp = input("Serial number %s will be put in the old databases. Proceed ['y','n']? "%SERIAL_NUMBER)
        if usr_inp != 'y':
            print("Serial numbers not updated.")
            exit(0)
        else:
            update_serial_numbers(SERIAL_NUMBER,old_fn)
            print("Serial numbers updated. Check:")
            prompt_user(new,old,check=True)

def update_serial_numbers(SN,db_names):
    """
    Edit the SERIAL_NUMBER variables in the old dbs
    """
    for db_name in db_names:
        dbs = Dataset(db_name,'r+')
        dbs['SERIAL_NUMBER'][0] = float(SN)
        dbs.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Updated serial numbers in yambo databases')
    parser.add_argument('-new','--new_serial', type=str, default="./SAVE", help='<Optional> Path to folder with the newer databases (Default is ./SAVE)')
    parser.add_argument('-old','--old_serial', type=str,help='<Required> Path to folder with the older databases',required=True)
    args = parser.parse_args()

    new_dir = args.new_serial
    old_dir = args.old_serial

    prompt_user(new_dir,old_dir)
