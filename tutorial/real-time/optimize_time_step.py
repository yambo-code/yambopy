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

    yambo_rt = '/Users/fulvio.paleari/software/yambo-devel/bin/yambo_rt' 
    SAVE_path = './database/FixSymm/SAVE'    
    #YamboRTStep_Optimize(input_path=args.input_file,SAVE_path=SAVE_path)
    #YamboRTStep_Optimize(input_path=args.input_file,SAVE_path=SAVE_path,RUN_path=args.directory,TStep_MAX=80,TStep_increase=20,NSimulations=4,ref_time=60,tol_pol=5e-3)
    #YamboRTStep_Optimize(input_path=args.input_file,SAVE_path=SAVE_path,RUN_path=args.directory,TStep_MAX=12,TStep_increase=2,NSimulations=6,ref_time=60,tol_pol=5e-3)
    YamboRTStep_Optimize(input_path=args.input_file,SAVE_path=SAVE_path,RUN_path=args.directory,TStep_MAX=12,TStep_increase=2,NSimulations=6,yambo_rt=yambo_rt)
    #YamboRTStep_Optimize(input_path=args.input_file,SAVE_path=SAVE_path,RUN_path=args.directory,TStep_MAX=10,TStep_increase=2.5,NSimulations=4,yambo_rt=yambo_rt)
