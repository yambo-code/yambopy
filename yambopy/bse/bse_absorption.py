# Copyright (c) 2018, Henrique Pereira Coutada Miranda
#
# All rights reserved.
#
# This file is part of the yambopy project
#
"""
This file contains a class to analyse the BSE absorption spectra

The features to be provided include: 
    1. Read the BSE diago file and plot chi including excitonic effects
    2. Obtain the excitonic states energy, degeneracies and optical intensity
    3. Create flows to run ypp and obtain the excitonic wavefunctions
    4. Produce a json file that can be used in the excitonwebsite

Long term:
    1. For the same exciton store multiple WFs changing the hole position

To initialize this structure we require instances of:
    - YamboExcitonDB (which already requires YamboLatticeDB)

The json format produced by this class for the excitonwebsite is:

structure stuff
    -lattice
    -atoms
    -atypes
bse stuff
    -eps - Absorption spectra
    -eel - Electron energy loss spectra
    -nx,ny,nz - the dimensions of the datagrid
    -datagrid
    -excitons - list of excitons to show

each exciton entry contains:
    -datagrid
    -intensity
    -index
    -energy

"""
from __future__ import print_function, division
import os
from yambopy.tools.string import marquee
from yambopy import *

class YamboBSEAbsorptionSpectra():
    """
    Create a file with information about the excitons from Yambo files
    """
    def __init__(self,excitondb):
        """
        Parameters:
            excitondb - an instance of the excitonDB class
        """
        self.excitondb = excitondb

    @classmethod
    def from_folder(cls,folder):
        raise NotImplementedError('from_folder not implemented yet')
        #try to find a excitonDB class
        #initialize this class
        return cls

    def get_excitons(self,min_intensity=0.1,max_energy=4,eps=1e-4):
        """ 
        Obtain the excitons using ypp
        Parameters:
            min_intensity - Only plot excitons with intensity larger than this value (default: 0.1)
            max_energy    - Only plot excitons with energy below this value (default: 4 eV)
            Degen_Step    - Only plot excitons whose energy is different by more that this value (default: 0.0)
        """
        excitons = self.excitondb.get_nondegenerate(eps=eps)

        #filter with energy
        excitons = excitons[excitons[:,0]<max_energy]

        #filter with intensity
        excitons = excitons[excitons[:,1]>min_intensity]

        #filter with degen
        if Degen_Step:

            #create a list with differences in energy
            new_excitons = []
            prev_exc = 0
            for exc in self.excitons:
                e,i,index = exc
                #if the energy of this exciton is too diferent then we add it to the list
                if abs(e-prev_exc)<Degen_Step:
                    new_excitons[-1][1] += i
                    continue
                new_excitons.append([e,i,index])
                intensity = 0
                prev_exc = e
            self.excitons = np.array(new_excitons)

        #create dictionary with excitons 
        excitons = self.data["excitons"]
        for e,intensity,i in self.excitons:
            exciton = {"energy": e,
                       "intensity": intensity,
                       "index": i}
            excitons.append(exciton)
        return self.excitons

    def get_wavefunctions(self, FFTGvecs=30,
                          Cells=[1,1,1], Hole=[0,0,0],
                          Direction="123", Format="x",
                          Degen_Step=0.0100,
                          MinWeight=1e-8,
                          repx=list(range(-1,2)), repy=list(range(-1,2)), repz=list(range(-1,2)),
                          wf=False):
        """
        Collect all the wavefuncitons with an intensity larger than self.threshold
        Parameters:
          FFTGvecs   - Number of FFTGvecs. Related to how accurate the representation is
          Cells      - Number of cells to plot in real space
          Hole       - Define the hole position in cartesian coordinates
          Direction  - Choose the directions to plot along
          Format     - Choose the format to plot in. Can be: x for xcrysden or g for gnuplot (default: 'x' for xcrysden)
          Degen_Step - Threshold to merge degenerate states. If the difference in energy between the states is smaller than
                       this value their wavefunctions will be plotted together
          repx       - Number or repetitions along the x direction
          repy       - Number or repetitions along the x direction
          repz       - Number or repetitions along the x direction
          wf         - Get the wavefuncitons in real space or not (default: False)
        """
        if self.excitons is None:
            raise ValueError( "Excitons not present. Run YamboBSEAbsorptionSpectra.get_excitons() first" )
        self.data["excitons"] = []

        #create a ypp file using YamboIn for reading the wavefunction
        yppwf = YamboIn('ypp -e w -V all',filename='ypp.in',folder=self.path)
        yppwf['Format'] = Format
        yppwf['Direction'] = Direction
        yppwf['FFTGvecs'] = [FFTGvecs,'Ry']
        yppwf['Degen_Step'] = [Degen_Step,'eV']
        yppwf['Hole'] = [Hole,'']
        yppwf['Cells'] = [Cells,'']

        #create a ypp file using YamboIn for reading the excitonic weights
        yppew = YamboIn('ypp -e a',filename='ypp.in',folder=self.path)
        yppew['MinWeight'] = MinWeight
        yppew['Degen_Step'] = Degen_Step

        keywords = ["lattice", "atoms", "atypes", "nx", "ny", "nz"]
        for exciton in self.excitons:
            #get info
            e,intensity,i = exciton

            if wf:
                ##############################################################
                # Excitonic Wavefunction
                ##############################################################
                #create ypp input for the wavefunction file and run
                yppwf["States"] = "%d - %d"%(i,i)
                yppwf.write("%s/yppwf_%d.in"%(self.path,i))

                filename = "o-%s.exc_%dd_%d%s"%(self.job_string,len(Direction),i,{"g":"","x":".xsf"}[Format] )
                print(filename)
                if not os.path.isfile(filename):
                    os.system("cd %s; ypp -F yppwf_%d.in -J %s"%(self.path,i,self.job_string))

                #read the excitonic wavefunction
                if Format == 'x':
                    ewf = YamboExcitonWaveFunctionXSF()
                else:
                    ewf = YamboExcitonWaveFunctionGnuplot()
                ewf.read_file("%s/%s"%(self.path,filename))
                data = ewf.get_data()
                for word in keywords:
                    if word in data:
                        self.data[word] = data[word]

                #calculate center of mass of atoms
                lat = np.array(data["lattice"])
                center_atom = np.zeros([3])
                for atype,x,y,z in data["atoms"]:
                    center_atom += np.array([x,y,z])
                center_atom /= len(data["atoms"])
                center_atom_red = car_red([center_atom],lat)[0]

                #shift wavefunctions grid to center of mass
                nx = data['nx'] 
                ny = data['ny'] 
                nz = data['nz']

                #make center_atom_red commensurate with fft
                center_atom_red = center_atom_red * np.array([nx,ny,nz])
                center_atom_red_int = [int(x) for x in center_atom_red]
                displacement = np.array([nx,ny,nz])/2-center_atom_red_int
                dx,dy,dz = displacement

                # shift grid
                # http://www.xcrysden.org/doc/XSF.html 
                dg  = np.array(data["datagrid"]).reshape([nz,ny,nx])
                dg  = np.roll(dg,dx,axis=2)
                dg  = np.roll(dg,dy,axis=1)
                dg  = np.roll(dg,dz,axis=0)
                data["datagrid"] = dg.flatten()

                #shift atoms
                atoms = []
                dx,dy,dz = red_car([displacement/np.array([nx,ny,nz],dtype=float)],lat)[0]
                for atype,x,y,z in data["atoms"]:
                    atoms.append([atype,x+dx,y+dy,z+dz])
                self.data["atoms"] = atoms

            ##############################################################
            # Excitonic Amplitudes
            ##############################################################
            #create ypp input for the amplitudes file and run
            yppew["States"] = "%d - %d"%(i,i)
            yppew.write("%s/yppew_%d.in"%(self.path,i))

            filename = "%s/o-%s.exc_weights_at_%d"%(self.path,self.job_string,i)
            if not os.path.isfile(filename):
                os.system("cd %s; ypp -F yppew_%d.in -J %s"%(self.path,i,self.job_string))

            #read the excitonic weigths
            ew = YamboExcitonWeight(filename,save=self.save,path=self.path)
            qpts, weights = ew.calc_kpts_weights(repx=repx,repy=repy,repz=repz)

            ############
            # Save data
            ############
            exciton = {"energy": e,
                       "intensity": intensity,
                       "weights": weights,
                       "qpts": qpts,
                       "index": i}
            if wf:
                exciton["hole"] = Hole
                exciton["datagrid"] = np.array(data["datagrid"])

            self.data["excitons"].append(exciton)

    def plot(self):
        pass

    def write_json(self,filename="absorptionspectra",verbose=0):
        """ 
        Write a jsonfile with the absorption spectra 
        and the wavefunctions of certain excitons
        """
        JsonDumper(self.data,"%s.json"%filename)

    def __str__(self):
        lines = []; app = lines.append
        app(marquee(self.__class__.__name__))
        app(str(self.excitondb.get_string(mark='-')))
        return "\n".join(lines)
