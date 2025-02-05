import numpy as np

class GreensFunctions():
    def __init__(self,w, E, eta, focc=1):
        self.w = w
        self.E = E
        self.eta = eta
        self.focc=focc
        self.GR = self.G_retarded()
        self.GA = self.G_advanced()
        self.A = self._spectralfunction()
        #to-do fix occupations
        self.G_less = self.G_lesser(focc,self.A)
        self.G_great = self.G_greater(focc,self.A)

    def G_retarded(self):
        GR = 1/(self.w-self.E+1j*self.eta)
        return GR
    
    def G_advanced(self):
        GA = 1/(self.w-self.E-1j*self.eta)
        return GA

    def G_lesser(self, focc, A):
        G_less = 1j*focc * A
        return G_less

    def G_greater(self, focc, A):
        G_great = -1j*(1-focc) * A
        return G_great
    
    def _spectralfunction(self):
        A = -2*np.imag(self.GR)
        return A