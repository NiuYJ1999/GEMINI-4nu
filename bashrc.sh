# set environment
export GTOP=`pwd`
cd $GTOP
export GINPUT="${GTOP}/gemini/"
export GROOT="${GTOP}/ROOT/"

# check if the executable file already exists
if [ ! -f "${GTOP}/bin/deexG" ]; then
    mkdir bin
    # compile if the executable file does not exist
    cd $GINPUT
    cmake .
    make
    cd $GTOP

    # move the executable file to the bin directory
    mv ${GINPUT}/deexGen ${GTOP}/bin/deexG
else
    echo "Executable file deexGen already exists. Skipping cmake and make operations."
fi

# create an alias for deexGen
# alias deexG="${GTOP}/bin/deexG"

export PATH=${GTOP}/bin:$PATH

