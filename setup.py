from setuptools import setup

packages_yambopy = ['yambopy',
                    'yambopy.io',
                    'yambopy.dbs',
                    'yambopy.bse',
                    'yambopy.rt',
                    'yambopy.double_grid',
                    'yambopy.data',
                    'yambopy.plot',
                    'yambopy.em1s',
                    'yambopy.tools',
                    'yambopy.common',
                    'yambopy.gkkp',
                    'yambopy.flow',
                    'qepy',
                    'qepy.upf_interface',
                    'qepy.data.pseudos',
                    'schedulerpy',
                    'yamboparser',
                    'command_line']

install_requires = [
"numpy",
"scipy",
"netCDF4",
"matplotlib",
"pyyaml",
"lxml",
]

if __name__ == '__main__':
    setup(name='yambopy',
          version='0.2.0',
          description='Pre-Postprocessing and automatic workflows for Yambo (and Quantum Espresso).',
          author='Fulvio Paleari, Alejandro Molina-Sánchez, Riccardo Reho, Davide Romanin, Alexandre Morlet and Henrique Miranda',
          author_email='fulvio.paleari90@gmail.com',
          requires=['numpy','scipy','matplotlib','netCDF4','pyyaml','lxml'],
          scripts=['scripts/yambopy'],
          packages=packages_yambopy,
          install_requires=install_requires,
          )
