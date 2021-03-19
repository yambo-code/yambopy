from setuptools import setup

packages_yambopy = ['yambopy',
                    'yambopy.io',
                    'yambopy.dbs',
                    'yambopy.bse',
                    'yambopy.rt',
                    'yambopy.double_grid',
                    'yambopy.data',
                    'yambopy.plot',
                    'yambopy.tools',
                    'yambopy.common',
                    'yambopy.gkkp',
                    'qepy',
                    'materials',
                    'schedulerpy',
                    'yamboparser',
                    'command_line']

install_requires = [
"numpy",
"scipy",
"netCDF4",
"matplotlib",
]

if __name__ == '__main__':
    setup(name='yambopy',
          version='1',
          description='Pre-Postprocessing and automatic workflows for Yambo (and Quantum Espresso).',
          author='Henrique Miranda, Alejandro Molina-Sánchez, Fulvio Paleari, Alexandre Morlet',
          author_email='fulvio.paleari90@gmail.com',
          requires=['numpy','scipy','matplotlib','netCDF4'],
          scripts=['scripts/yambopy'],
          packages=packages_yambopy,
          install_requires=install_requires,
          )
