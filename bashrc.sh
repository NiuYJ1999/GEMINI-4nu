# set environment
export GTOP=`pwd`
cd $GTOP
mkdir ROOT
export GINPUT="${GTOP}/gemini/"
export GROOT="${GTOP}/ROOT/"

# compile
cd $GINPUT
cmake .
make
cd $GTOP

