# Copyright (c) 2016, Henrique Miranda
# All rights reserved.
#
# This file is part of the yambopy project
#
from yambopy import *
from bse_cutoff import *



def plot(read_folder,layer_separations):
    #create plot
    ax = plt.gca()

    for layer_separation in layer_separations:
        
        print "layer separation:", layer_separation

        folder = "%s/%d/%d"%(read_folder,layer_separation,layer_separation)
        ys = YamboStaticScreeningDB(save=folder)
        print ys 

        #plot static screening 
        ys.plot(ax,marker='o',label=layer_separation)

    #final plot
    plt.legend()
    plt.show()


if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='Convergence test of the colomb cutoff')
    parser.add_argument('folder', help='Folder where to read em1s from')
    args = parser.parse_args()
    
    plot(args.folder,layer_separations)
