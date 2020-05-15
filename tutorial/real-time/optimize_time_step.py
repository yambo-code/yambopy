from yambopy import *
import sys
import argparse

if __name__ == "__main__":
    #parse options
    parser = argparse.ArgumentParser(description='RT Time step optimization')
    parser.add_argument('-F', '--input_file',type=str,help='<Required> RT input file',required=True)
    parser.add_argument('-D', '--directory', type=str,help='RT conv directory',required=False)

    args = parser.parse_args()

    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    SAVE_path = 'database/FixSymm/SAVE'    
    #YamboRTStep_Optimize(input_path=args.input_file,SAVE_path=SAVE_path)
    YamboRTStep_Optimize(input_path=args.input_file,SAVE_path=SAVE_path,RUN_path=args.directory)
    #YamboRTStep_Optimize(input_path=args.input_file,SAVE_path=SAVE_path,RUN_path='RT_test_QSSIN_NAN',TSteps_min_max=[100,2000],NSimulations=8)
