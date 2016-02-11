from yambopy import *

pack_files_in_folder('gw_par')

#plot the results using yambm analyser
ya = YamboAnalyser()
print(ya)
print('plot all qpoints')
ya.plot_gw('qp')
print('plot along a path')
path = [[   0,   0,   0],
        [ 0.5,   0,   0],
        [1./3,1./3,   0],
        [   0,   0,   0]]

ya.plot_gw_path('qp',path)

print('done!')
