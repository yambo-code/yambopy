from qepy import *
#test supercell class

qe = PwIn('scf/bn.scf')
R=[4,3,2]
sup = supercell(qe,R)
print sup.qe

calculation = ''.join(sup.qe.control['calculation'].split('\''))
prefix = ''.join(sup.qe.control['prefix'].split('\''))
sup.qe.write('%s.%s'%(prefix,calculation))
print 'supercell input file written.'
