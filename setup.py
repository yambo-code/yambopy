from distutils.core import setup

packages_yambopy = ['yambopy','qepy']

if __name__ == '__main__':
    setup(name='yambopy',
          version='0.1',
          description='Automatic workflows for Yambo.',
          author='Henrique Miranda',
          author_email='miranda.henrique@gmail.com',
          packages=packages_yambopy,
          )
