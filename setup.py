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
                    'qepy',
                    'materials',
                    'schedulerpy',
                    'yamboparser']

install_requires = [
"numpy",
"scipy",
"netCDF4",
"matplotlib",
]

if __name__ == '__main__':
    setup(name='yambopy',
          version='0.1',
          description='Automatic workflows for Yambo.',
          author='Henrique Miranda',
          author_email='miranda.henrique@gmail.com',
          requires=['numpy','matplotlib','netCDF4'],
          scripts=['scripts/yambopy'],
          packages=packages_yambopy,
          install_requires=install_requires,
          )
