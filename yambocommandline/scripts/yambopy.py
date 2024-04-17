#!/usr/bin/env python
#TODO: delete/move any other scripts in this directory
from yambopy import *
from yambocommandline.commands import *
import argparse
import sys

class Cmd():
    """
    Define some generic functions for a command
    """
    def info(self):
        """
        display the available commands
        """
        print('yambopy')
        print('Available commands are:\n')
        for cmd,c in list(self._commands.items()):
            print("%15s -> %s"%(cmd, c.__doc__.split('\n')[1]))
    
    def run(self,cmds,args):
        """
        generic run command
        cmds is a dictionary that maps the command to the function to run
        """
        cmd = args[0]
        if cmd in list(cmds.keys()):
            cmds[cmd](args[1:])

class PlotExcitons(Cmd):
    """
    Plot excitons calculation
        
        possible arguments are:
        
        Arguments:
        filename -> json file containing the absorption spectra. Default: 'absorptionspectra.json' 
        -s       -> Size of the materis in the plot
    """
    def __init__(self,args):
        import matplotlib

        #check for args
        if len(args) < 1:
            print((self.__doc__))
            exit(0)

        parser = argparse.ArgumentParser(description='Study convergence on BS calculations using ypp calls.')
        pa = parser.add_argument
        pa('filename',     help='json file containing the absorption spectra. Default: \'absorptionspectra.json\'' )
        pa('-s','--size',  help='Size of the markers in the plot', default=20, type=int)
        args = parser.parse_args(args)

        if os.path.isfile(args.filename):
            #create plot
            recipes.plot_excitons(args.filename,size=args.size) 
        else:
            print('file %s is invalid'%filename)

class PlotEm1sCmd(Cmd):
    """
    Plot em1s calculation

        possible arguments are:
        
        Arguments:
        folders        -> Folders containing the ndb.ems1 files
        -w, --write    -> Write data file in a text file
        -p, --plot     -> Save a file with the plot
        -v, --verbose  -> Print which files are not folder
        --fontsize     -> Choose the font size of the plot
    """
    def __init__(self,args):
        import matplotlib
        import matplotlib.pyplot as plt        
        from glob import glob

        #check for args
        if len(args) < 1:
            print((self.__doc__))
            exit(0)

        parser = argparse.ArgumentParser(description='Plot em1s calculation.')
        pa = parser.add_argument
        pa('folders',        nargs='+', help='json file containing the absorption spectra. Default: \'absorptionspectra.json\'' )
        pa('-w','--write',   help='Write data file in a text file', default='em1s.dat', type=str)
        pa('-p','--plot',    help='Save a file with the plot', default='em1s.pdf', type=str)
        pa('-v','--verbose', help='Print which files are not folder', action='store_true')
        pa('--fontsize',     help='Choose the font size of the plot', default=10, type=int)
        args = parser.parse_args(args)
        folders = args.folders
        #print(glob('bse_cutoff/*/*'))
        #print(folders)
        #exit() 
        #create plot
        fig = plt.figure(figsize=(6,5))
        matplotlib.rcParams.update({'font.size': args.fontsize})
        ax = plt.gca()

        epsilons = []
        for folder in folders: 
            if os.path.isdir(folder) and os.path.isfile("%s/ndb.em1s"%folder):
                if not os.path.isfile("%s/ns.db1"%folder):
                    if not os.path.isfile("%s/../SAVE/ns.db1"%folder):
                        if args.verbose: print("SAVE folder not found")
                    else:
                        os.system('cp %s/../SAVE/ns.db1 %s/'%(folder,folder))
                print(folder)
                ys = YamboStaticScreeningDB(save=folder,em1s=folder)
                #plot epsilon^-1_{00} = [(1+vX)]_{00}
                ys.plot_epsm1(ax,marker='o',markersize=2,label=folder)
                #get epsilon^{-1}_{00} = [1+vX]_{00}
                x,vX = ys._getvxq()
                epsilons.append([folder,x,1+vX])
            else:
                if args.verbose:
                    print("path %s is not a folder"%folder)

        if args.write:
            #write a text file with the data
            f=open(args.write,'w')
            f.write('# Real and imaginary part of \\epsilon^{-1}_{00}(\\omega=0,q) = [1+vX]_{00} as a funciton of |q|')
            for folder,x,y in epsilons:
                f.write('#%s\n'%folder)
                for xi,yi in zip(x,y):
                    f.write("%12.8lf %12.8lf %12.8lf\n"%(xi,yi.real,yi.imag))
                f.write('\n\n')
            f.close()

        plt.legend(frameon=False)
        plt.tight_layout()

        #final plot
        if args.write: plt.savefig(args.plot) 
        plt.show()

class AnalyseGWCmd(Cmd):
    """
    Study the convergence of GW calculations by looking at the change in band-gap value.

    The script reads from <folder> all results from <variable> calculations and display them.

    Use the band and k-point options according to the size of your k-grid
    and the location of the band extrema.

        Mandatory arguments are:
            folder   -> Folder containing SAVE and convergence runs.
            var      -> Variable tested (e.g. FFTGvecs)

        Optional variables are:
            -bc, --bandc   (int)  -> Lowest conduction band number
            -kc, --kpointc (int)  -> k-point index for conduction band
            -bv, --bandv   (int)  -> Highest valence band number
            -kv, --kpointv (int)  -> k-point index for valence band
            -np, --nopack  (flag) -> Do not call 'pack_files_in_folder'
            -nt, --notext  (flag) -> Do not print a text file
            -nd, --nodraw  (flag) -> Do not draw (plot) the result
    """

    def __init__(self,args):

        #check for args
        if len(args) <= 1:
            print((self.__doc__))
            exit(0)

        #all the other arguments are passed to the analysegw function

        parser = argparse.ArgumentParser(description='Study GW convergence with regards to the band-gap value.')
        pa = parser.add_argument
        pa('folder'            , help='Folder containing SAVE and convergence runs.')
        pa('variable'          , help='Variable tested (e.g. FFTGvecs)' )
        pa('-bc','--bandc'     , help='Lowest conduction band number'    , default=53, type=int)
        pa('-kc','--kpointc'   , help='K-point index for conduction band', default=19, type=int)
        pa('-bv','--bandv'     , help='Highest valence band number'      , default=52, type=int)
        pa('-kv','--kpointv'   , help='K-point index for valence band'   , default=1 , type=int)
        pa('-np','--no-pack'   , help='Skip the packing of output files'        , dest='pack', action='store_false')
        pa('-nt','--no-text'   , help='Skip the writing of the analysis result' , dest='text', action='store_false')
        pa('-nd','--no-draw'   , help='Skip the plotting of the analysis result', dest='draw', action='store_false')
        pa('-v', '--verbose'   , action='store_false')
        parser.set_defaults(pack=True,text=True,draw=True)
        args = parser.parse_args(args)

        folder = args.folder ; var     = args.variable
        bandc  = args.bandc  ; kpointc = args.kpointc
        bandv  = args.bandv  ; kpointv = args.kpointv
        pack   = args.pack   ; text    = args.text
        draw   = args.draw   ; verbose = args.verbose
        
        #call analyse_gw from recipes.py
        recipes.analyse_gw(folder,var,bandc,kpointc,bandv,kpointv,pack,text,draw,verbose)

    def info(self):
        """
        display help to use this command
        """
        print(self.__doc__)

class AnalyseBSECmd(Cmd):
    """
    Using ypp, you can study the convergence of BSE calculations in 2 ways:
      Create a .png of all absorption spectra relevant to the variable you study
      Look at the eigenvalues of the first n "bright" excitons (given a threshold intensity)

    The script reads from <folder> all results from <variable> calculations for processing.
    The resulting pictures and data files are saved in the ./analyse_bse/ folder.

    Mandatory arguments are:
        folder   -> Folder containing SAVE and convergence runs.
        var      -> Variable tested (e.g. FFTGvecs)

    Optional arguments are:
        -ne,--numbexc  (int)   -> Number of excitons to read beyond threshold (default=2)
        -ie,--intexc   (float) -> Minimum intensity for excitons to be considered bright (default=0.05)
        -de,--degenexc (float) -> Energy threshold under which different peaks are merged (eV) (default=0.01)
        -me,--maxexc   (float) -> Energy threshold after which excitons are not read anymore (eV) (default=8.0)
        -np,--nopack   (flag)  -> Skips packing o- files into .json files
        -nt,--notext   (flag)  -> Skips writing the .dat file
        -nd,--nodraw   (flag)  -> Skips drawing (plotting) the abs spectra
    """
    def __init__(self,args):
        #check for args
        if len(args) < 2:
            print((self.__doc__))
            exit(0)

        parser = argparse.ArgumentParser(description='Study convergence on BS calculations using ypp calls.')
        pa = parser.add_argument
        pa('folder',           help='Folder containing SAVE and convergence runs.' )
        pa('variable',         help='Variable tested (e.g. FFTGvecs)' )
        pa('-ne','--numbexc',  help='Number of excitons to read beyond threshold', default=2,type=int)
        pa('-ie','--intexc',   help='Minimum intensity for excitons to be considered bright', default=0.05,type=float)
        pa('-de','--degenexc', help='Energy threshold under which different peaks are merged (eV)', default=0.01,type=float)
        pa('-me','--maxexc',   help='Energy threshold after which excitons are not read anymore (eV)', default=8.0,type=float)
        pa('-nt','--notext',   help='Skips writing the .dat file', action='store_false')
        pa('-nd','--nodraw',   help='Skips drawing (plotting) the abs spectra', action='store_false')
        args = parser.parse_args(args)

        folder    = args.folder
        var       = args.variable
        exc_n     = args.numbexc
        exc_int   = args.intexc
        exc_degen = args.degenexc
        exc_max_E = args.maxexc
        text      = args.notext
        draw      = args.nodraw

        #all the other arguments are passed to the analyse bse function
        recipes.analyse_bse( folder, var, exc_n, exc_int, exc_degen, exc_max_E, text=text, draw=draw )

class TestCmd(Cmd):
    """
    Run yambopy tests
    
        possible arguments are:

        basic -> fast test where input/output is compared with reference files
        full  -> requires yambo and quantum espresso to be installed
    """
    
    def __init__(self,args):
        #check for args
        if len(args) < 1:
            print((self.__doc__))
            exit(0)

        cmds = {'basic':self.basic,
                'full':self.full}
        self.run(cmds,args)

    def basic(self,*args):
        os.system('py.test --cov-config=.coveragerc --cov')
        
    def full(self,*args):
        print(args)
    

class MergeQPCmd(Cmd):
    """
    Merge QP databases
    
        possible arguments are:

           <QP files>    -> list of QP files produced by yambo
        -o <output file> -> output file where to save the merged db
    """

    def __init__(self,args):
        """ 
        possible arguments are:
        """ 
        #check for args
        if len(args) <= 1:
            print((self.__doc__))
            exit(0)
        
        #all the other arguments are passed to the merge_qp function
        parser = argparse.ArgumentParser(description='Join different NetCDF quasi-particle databases')
        parser.add_argument('files', nargs='+', type=argparse.FileType('r'))
        parser.add_argument('-o','--output',                       help='Output filename', default='ndb_out.QP')
        args = parser.parse_args(args)

        output  = args.output
        files   = args.files
        
        #call merge_qp from recipes.py
        recipes.merge_qp(output,files)

    def info(self):
        """
        display help to use this command 
        """
        print(self.__doc__)


class AddQPCmd(Cmd):
    """
    Add corrections from QP databases.
    This function reads the QP correction from Yambo databases and add them together. Pass a list of files (separated by a space) after each flag to get the total correction you need.
    The k-grid must be the same, the number of bands can be different. There will always be at least the LDA value. The Z factor is taken as one.

    Optional arguments are:

        -a,  --add       -> Add the real-part of a QP correction
        -s,  --subtract -> Substract the real-part of a QP correction
        -ai, --addimg    -> Add the imaginary-part of a QP correction
        -o,  --output    -> Output filename
        -v,  --verbose   -> Increased verbosity

    """

    def __init__(self,args):
        #check for args
        if len(args) <= 1:
            print((self.__doc__))
            exit(0)

        parser = argparse.ArgumentParser(description='Join different NetCDF quasi-particle databases.')
        parser.add_argument('-a', '--add',       nargs='+', type=argparse.FileType('r'), help="Add the real part correction E-Eo to the final db",default=[])
        parser.add_argument('-s', '--subtract', nargs='+', type=argparse.FileType('r'), help="Subtract the real part correction E-Eo part to the final db", default=[])
        parser.add_argument('-ai','--addimg',    nargs='+', type=argparse.FileType('r'), help="Add the imaginary part Im(E) to the final db",default=[])
        parser.add_argument('-o', '--output',  default='ndb_out.QP', help='Output filename')
        parser.add_argument('-v', '--verbose', action="store_true",  help='Verbose mode')
        args = parser.parse_args(args)


        output  = args.output
        add     = args.add
        subtract = args.subtract
        addimg  = args.addimg
        verbose = args.verbose


        #call add_qp from recipes.py
        recipes.add_qp(output,add,subtract,addimg,verbose)

    def info(self):
        """
        display help
        """
        print(self.__doc__)

class GkkpCmd(Cmd):
    """
    Produce a SAVE folder including elph_gkkp databases
    
    Arguments are:
        -nscf, --nscf_dir  -> <Optional> Path to nscf save folder
        -elph, --elph_dir  -> Path to elph_dir folder
        -y, --yambo_dir    -> <Optional> Path to yambo executables
        -e, --expand       -> <Optional> Expand gkkp databases
    """
    def __init__(self,args):
        #check for args
        if len(args) < 2:
            print((self.__doc__))
            exit(0)
            
        parser = argparse.ArgumentParser(description='Generate SAVE folder including gkkp databases')
        parser.add_argument('-nscf','--nscf_dir', type=str, default="", help='<Optional> Path to nscf save folder')
        parser.add_argument('-elph','--elph_dir', type=str,help='<Required> Path to elph_dir folder',required=True)
        parser.add_argument('-y','--yambo_dir', default="", type=str,help='<Optional> Path to yambo executables')
        parser.add_argument('-e','--expand', action="store_true", help="Expand GKKP")
        args = parser.parse_args(args)
        
        nscf_dir  = args.nscf_dir
        elph_dir  = args.elph_dir
        yambo_dir = args.yambo_dir
        expand    = args.expand
        
        database = './'
        scheduler = Scheduler.factory

        #call gkkp
        gkkp.generate_gkkp(database,nscf_dir,elph_dir,yambo_dir,expand,scheduler)
        
class SaveCmd(Cmd):
    """
    Produce a SAVE folder
    
    Arguments are:
        -nscf, --nscf_dir  -> Path to nscf save folder
        -y, --yambo_dir    -> <Optional> Path to yambo executables
    """
    def __init__(self,args):
        #check for args
        if len(args) < 1:
            print((self.__doc__))
            exit(0)
            
        parser = argparse.ArgumentParser(description='Generate SAVE folder including gkkp databases')
        parser.add_argument('-nscf','--nscf_dir', type=str,help='<Required> Path to nscf save folder', required=True)
        parser.add_argument('-y','--yambo_dir', default="", type=str,help='<Optional> Path to yambo executables')
        args = parser.parse_args(args)
        
        nscf_dir  = args.nscf_dir
        yambo_dir = args.yambo_dir
        
        database = './'
        scheduler = Scheduler.factory

        #call generate_save
        generate_save.generate_save(database,nscf_dir,yambo_dir,scheduler)

class UpdtSrlNmbrCmd(Cmd):
    """
    Script to update serial numbers of yambo ndb.* databases in order to import them to new calculations.

    Inputs:
    1. -new, --new_serial='path/to/folder/with/new/dbs' [e.g., the new SAVE]
    2. -old, --old_serial='path/to/folder/with/old/dbs' [e.g., an old ndb.em1s]

    This script will prompt the user to go through with updating the dbs.
    """
    def __init__(self,args):
        #check for args
        if len(args) < 2:
            print((self.__doc__))
            exit(0)

        parser = argparse.ArgumentParser(description='Updated serial numbers in yambo databases')
        parser.add_argument('-new','--new_serial', type=str, default="./SAVE", help='<Optional> Path to folder with the newer databases (Default is ./SAVE)')
        parser.add_argument('-old','--old_serial', type=str,help='<Required> Path to folder with the older databases',required=True)
        args = parser.parse_args(args)

        new_dir = args.new_serial
        old_dir = args.old_serial

        #check dbs
        update_serial.prompt_user(new_dir,old_dir)

class PlotBndStrCmd(Cmd):
    """
    Script to produce band structure data and visualization from QE.
    
    It reads data-file-schema.xml found in the save folders of newer pw versions.
    
    It requires a yaml input file (see below)

    Arguments are:
        -i, --band_input   -> Path to yaml input file
        -o, --output_name  -> Name of output file <optional>
        -e, --erange       -> Energy window for plot <optional>
        -s, --show         -> Show plot at runtime <optional>

    - YAML INPUT FILE 
    This is a yaml file in the following format, to be copy-pasted and edited:
    
    :: yaml
    
        ---
         save_dir: "PREFIX.save"
         KPTs: [[0.,0.,0.],[0.,0.5,0.],[0.5,0.5,0.],[0.5,0.,0.]]
         KPTs_labels: [G, Y, S, X]
         shift_Delta_c_v: [0.,1.,1.]
    
    :: 
    
    The input parameters are:
    - save_dir: path to folder containing 'data-file-schema.xml' (in principle the QE save folder)
    - KPTs: band circuit in reduced coordinates
    - KPTs_labels: labels for the band circuit points
    - shift_Delta_c_v: k-dependent scissor shift as a list of three values (gap shift, cond. stretch, val. stretch)
    """
    def __init__(self,args):

        #check for args
        if len(args) < 2:
            print((self.__doc__))
            exit(0)

        parser = argparse.ArgumentParser(description='Generate band plot')
        parser.add_argument('-i','--band_input', type=str,  help='<Required> Path to input file',      required=True)
        parser.add_argument('-o','--output_name', type=str,  help='<Optional> Name of output file',    required=False, default=None)
        parser.add_argument('-e','--erange',      type=float,help='<Optional> Energy window for plot', required=False, default=4.)
        parser.add_argument('-s','--show',        type=bool, help='<Optional> Show plot at runtime',   required=False, default=False)
        args = parser.parse_args(args)

        inp = args.band_input

        # Read user input file
        input_params = generate_bands.read_input(inp)

        # Read data from .xml file
        output_data = generate_bands.get_data_from_xml(input_params)

        # Plot setup
        data_to_plot = generate_bands.setup_BZ_points(input_params,output_data)

        # Produce the plot
        generate_bands.launch_plot(data_to_plot,plt_type='bands',out_name=args.output_name,erange=args.erange,show=args.show)

class GwSubspace(Cmd):
    """
    Script to calculate off-diago corrections of yambo ndb.QP databases in order to plot band structure.

    Inputs:
    1. -d,--fld_diag='path/to/folder/with/diago/dbs' [e.g., the ./diag]
    2. -o,--fld_offdiag='path/to/folder/with/offdiagoold/dbs' [e.g., the ./offdiago]

    This script will prompt the user to go through with updating the dbs.
    """
    def __init__(self,args):
        #check for args
        if len(args) < 2:
            print((self.__doc__))
            exit(0)

        parser = argparse.ArgumentParser(description='Create a new diago dbs after a GW subspace calculation numbers in yambo databases')
        parser.add_argument('-d','--fld_diag', type=str, default="./diago", help='<Required> Path to folder with the ndb.QP diago-database   (Default is ./diago)',required=True)
        parser.add_argument('-o','--fld_offdiag', type=str,default ="./offdiago",help='<Required> Path to folder with the ndb.QP offdiagodiago-database   (Default is ./diago)',required=True)
        args = parser.parse_args(args)

        fld_diag = args.fld_diag
        fld_offdiag = args.fld_offdiag
        gw_subspace.create_newdb(fld_diag,fld_offdiag)

class GetPHqInputCmd(Cmd):        
    """
    Script to update the explicit list of q-points in a ph input file (ldisp=.false., qplot=.true.).

    - Reads the output of the scf calculation and the ph input.

    Usage:
    :: -pw,--pwout='path/to/pw/output/file
    :: -ph,--phin='path/to/ph/input/file'
    """
    def __init__(self,args):

        #check for args
        if len(args) < 4:
            print((self.__doc__))
            exit(0)

        parser = argparse.ArgumentParser(description='Append explicit qpoints to ph input')
        parser.add_argument('-pw','--pwout', type=str, help='Path to pw (scf) output file',required=True)
        parser.add_argument('-ph','--phin', type=str,help='Path to ph (dvscf or elph) input file',required=True)
        args = parser.parse_args(args)

        pwout = args.pwout
        phin  = args.phin

        get_phq_input.get_phq_input(pwout,phin)

class ConvertRLtoRyCmd(Cmd):
    """
    Script to convert RL number in Ry energy units using ndb.gops.

    Inputs:
    1. -gops,--ndb_gops='path/to/folder/with/ndb.gops' [i.e., SAVE]
    2. -v,--value='value to convert with units' [e.g., 11 RL or 5 Ry]

    The script will read ndb.gops and find the nearest completed G-shell, then give the
    converted value in Ry (RL) to the one supplied in input.
    """
    def __init__(self,args):

        #check for args
        if len(args) < 4:
            print((self.__doc__))
            exit(0)

        parser = argparse.ArgumentParser(description='Convert RL number in energy units (Ry)')
        parser.add_argument('-gops','--ndb_gops', type=str, default="./SAVE", help='<Optional> Path to folder with ndb.gops (default: ./SAVE)')
        parser.add_argument('-v','--value', type=str,help="<Required> Value to be converted along with units, e.g.: '11 RL' or '5 Ry'",nargs=2,required=True)
        args = parser.parse_args(args)

        ndb_gops = args.ndb_gops
        value    = args.value

        convert_RL_to_Ry.convert(value,ndb_gops)

class ConvertLELPHCtoYAMBO(Cmd):
	"""
	Calculate gauge-invariant electron-phonon matrix elements with LetzElPhC and convert them into Yambo format

	:: Usage:

	>> yambopy l2y -ph phinp -b b1 b2 -par nq nk [--lelphc lelphc] [--debug]

	:: Input parameters:
		-ph,--ph_inp_path     : path to ph.x input file, e.g. dvscf/ph.in
		-b,--bands            : initial and final band indices (counting from 1)
		-par,--pools [OPT]    : MPI pools for q and k (needs mpirun)
		-lelphc,--lelphc [OPT]: path to lelphc executable (code will prompt)
		-D,--debug [OPT]      : won't remove LetzElPhC input and outputs

	:: Prerequisites:

	* ph.x phonon calculation must be complete, e.g. the phinp folder should contain:
		- ph.x input file
		- pw.x (scf) save directory
		- [prefix].dyn files
		- _ph* directories
	* Yambo SAVE directory must be present. We run in the directory where the SAVE is.
	* LetzElPhC must be installed
	* mpirun must be linked for parallel runs
	"""
	def __init__(self,args):

		#check for args
		if len(args) < 5:
			print((self.__doc__))
			exit(0)

		parser = argparse.ArgumentParser(description='Generate electron-phonon coupling databases via LetzElPhC')
		parser.add_argument('-ph','--ph_inp_path', type=str, help='<Required> Path to ph.x (dvscf) input file',required=True)
		parser.add_argument('-b','--bands',nargs=2,type=str,help="<Required> First and last band (counting from 1), e.g. 'b_i b_f'",required=True)
		parser.add_argument('-par','--pools',nargs=2,type=str, default=[1,1], help="<Optional> MPI tasks as 'nqpools nkpools' (default serial)")
		parser.add_argument('-lelphc','--lelphc',type=str,default='lelphc',help="<Optional> Path to lelphc executable (default assumed in Path, otherwise prompted)")
		parser.add_argument('-D','--debug', action="store_true", help="Debug mode")

		args = parser.parse_args(args)

		phinp  = args.ph_inp_path
		bands  = args.bands
		pools  = args.pools
		lelphc = args.lelphc
		debug  = args.debug

		# Check inputs
		lelphc,ph_path,inp_ph,inp_lelphc,inp_name = \
		lelph_interface.checks(phinp,lelphc,bands,pools)

		# run preprocessing
		lelph_interface.run_preprocessing(lelphc,ph_path,inp_ph)

		# run el-ph calculation and rotation
		lelph_interface.run_elph(lelphc,inp_lelphc,inp_name)

		# load database and convert to yambo format
		lelph_interface.letzelph_to_yambo()

		# clean
		lelph_interface.clean_lelphc(debug,inp_name,ph_path)

class YambopyCmd(Cmd):
    """
    class to implement commands for yambopy.
    each new command to be added should be implemented as a class inheriting from this one
    """
    _commands = {'plotem1s':     PlotEm1sCmd,
                 'analysebse':   AnalyseBSECmd,
                 'analysegw':    AnalyseGWCmd,
                 'plotexcitons': PlotExcitons,
                 'addqp':        AddQPCmd,
                 'mergeqp':      MergeQPCmd,
                 'save':         SaveCmd,
                 'gkkp':         GkkpCmd,
                 'bands':        PlotBndStrCmd,
                 'serial':       UpdtSrlNmbrCmd,
                 'gwsubspace':   GwSubspace,
                 'phinp':        GetPHqInputCmd,
                 'convert':      ConvertRLtoRyCmd,
				 'l2y':          ConvertLELPHCtoYAMBO,
                 'test':         TestCmd}

    def __init__(self,*args):
        """
        parse the command from the command line and initialize the class responsible
        for handling such command
        """
        
        #check for args
        if len(args) <= 1:
            self.info()
            exit(0)
 
        #start call graph     
        if args[1] in self._commands:
            cmdclass = self._commands[args[1]]
            self.cmd = cmdclass(args[2:]) 
        else:
            self.info()  
            print() 
            print("Command %s is not known to yambopy"%args[1])
 
#parse options
#def run_script(): return YambopyCmd(*sys.argv)
ycmd = YambopyCmd(*sys.argv)
exit()
