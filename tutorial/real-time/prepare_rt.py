from yambopy import *
import sys
import argparse

if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='RT setup')
    parser.add_argument('-f', '--field_direction',nargs='+',type=float,help='<Required> Set field direction',required=True)    
    parser.add_argument('-p', '--prefix',type=str,help='<Required> QE prefix',required=True)

    args = parser.parse_args()

    if len(sys.argv)<=2:
        parser.print_help()
        sys.exit(1)

    CreateYamboSave(args.prefix,field_dir=args.field_direction)
