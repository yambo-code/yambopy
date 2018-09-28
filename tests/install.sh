#/usr/bin/env sh

set -e # terminate if error is encountered
set -x # print to stdout

echo "Yambo: $YAMBO_VERSION"
echo "Quantum Espresso: $PW_VERSION"

ls /usr/lib

# Yambo 4.4 tar
wget https://www.dropbox.com/s/7w8f0ll2zq6rwzy/yambo_4.4.tar?dl=0 -O yambo.tar
sudo tar -xf yambo.tar -C /bin/

# Espresso 6.3 tar
wget https://www.dropbox.com/s/lv82dwp434sud6d/espresso_6.3.tar?dl=0 -O espresso.tar
sudo tar -xf espresso.tar -C /bin/
