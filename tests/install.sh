/usr/bin/env sh

set -e # terminate if error is encountered
set -x # print to stdout

# Parse args
while echo $1 | grep -q ^-; do
    eval $( echo $1 | sed 's/^-//' )=$2
    shift
    shift
done

echo "YAMBO_VERSION = $y"
echo "PW_VERSION    = $pw"

wget https://www.dropbox.com/s/7g2wxesc0sxu19h/yambo?dl=0 -O yambo # yambo 4.1
chmod +x yambo
mv yambo /bin/
which yambo
yambo

