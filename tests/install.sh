#/usr/bin/env sh

set -e # terminate if error is encountered
set -x # print to stdout

echo "Yambo: $YAMBO_VERSION"
echo "Quantum Espresso: $PW_VERSION"

wget https://www.dropbox.com/s/hna9v4yw6u7ea5e/yambo_4.1.tar?dl=0 -O yambo.tar # Yambo 4.1 tar
sudo tar -xf yambo.tar -C /bin/

wget https://www.dropbox.com/s/o91hx046chfk26d/espresso_5.4.tar?dl=0 -O espresso.tar # Espresso 5.4 tar
sudo tar -xf espresso.tar -C /bin/
