#/usr/bin/env sh

set -e # terminate if error is encountered
set -x # print to stdout

echo "Yambo: $YAMBO_VERSION"
echo "Quantum Espresso: $PW_VERSION"

ls /usr/lib

# Yambo 4.4 tar
wget https://www.dropbox.com/s/mi493jb67u1lcog/yambo_4.4.tar.gz?dl=0 -O yambo.tar
sudo tar -xf yambo.tar -C /bin/

# Espresso 6.3 tar
wget https://www.dropbox.com/s/bxim3jgmjgd5v3h/qe_6.4.1.tar.gz?dl=0 -O espresso.tar
sudo tar -xf espresso.tar -C /bin/
