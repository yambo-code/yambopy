#/usr/bin/env sh

set -e # terminate if error is encountered
set -x # print to stdout

echo "Yambo: $YAMBO_VERSION"
echo "Quantum Espresso: $PW_VERSION"

wget https://www.dropbox.com/s/7g2wxesc0sxu19h/yambo?dl=0 -O yambo # yambo 4.1
chmod +x yambo
sudo mv yambo /bin/
which yambo
yambo

