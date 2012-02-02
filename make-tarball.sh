#!/bin/sh
function usage {
    echo "Usage: $0 [-s] <exercise>"
    echo "  Build tarball for specified exercise, defaults to template code"
    echo "  -s  Tar up solution code rather than template code"
    exit 1
}

if [[ $# != 1 && $# != 2 ]]; then
    usage
fi

SRC_DIR=exercises
HANDOUTS_DIR=handouts

if [[ $# == 2 ]]; then
    if [[ $1 != "-s" ]]; then
        echo "Unrecognised argument '$1'"
        usage
    else
        CODEDIR=solution
        DOSOLUTION=1
        EXERCISE=$2
    fi
else
    CODEDIR=template
    DOSOLUTION=0
    EXERCISE=$1
fi

# Find PDF handout for specified exercise
HANDOUT=$HANDOUTS_DIR/$EXERCISE.pdf

if [[ $DOSOLUTION == 1 ]]; then
    TARBALL=$EXERCISE-solution.tar
else
    TARBALL=$EXERCISE.tar
    if [ ! -f $HANDOUT ]; then
        echo "Unable to find handout PDF for $EXERCISE, have you written it?"
        exit 1
    fi
fi

# Guard against overwriting existing tarball
if [ -f $TARBALL ]; then
    echo -n "Tarball for $EXERCISE already exists.  Overwrite? [yN] "
    while read f; do
        canon=${f,,}
        if [ "x$f" == "x" ]; then
            canon="n"
        fi
        if [ $canon == "n" ]; then
            echo "Aborting tarball creation"
            exit 1
        elif [ $canon == "y" ]; then
            break
        else
            echo "Please answer y or n"
        fi
    done
fi

TMPDIR=`mktemp -d`

# Build directory structure
TARDIR=$TMPDIR/mandelbrot
mkdir -p $TARDIR

if [[ $DOSOLUTION == 0 ]]; then
    cp $HANDOUT $TARDIR
fi

mkdir $TARDIR/C
mkdir $TARDIR/F

# Clean up executables and output files (we don't want them)
(cd $SRC_DIR/C/$EXERCISE/$CODEDIR; make clean 2>/dev/null; rm -f output.ppm)>/dev/null
(cd $SRC_DIR/F/$EXERCISE/$CODEDIR; make clean 2>/dev/null; rm -f output.ppm)>/dev/null

cp -r $SRC_DIR/C/utils $TARDIR/C/

cp -r $SRC_DIR/F/utils $TARDIR/F/

mkdir -p $TARDIR/C/$EXERCISE
mkdir -p $TARDIR/F/$EXERCISE

cp $SRC_DIR/config.mak $TARDIR/
cp -r $SRC_DIR/C/$EXERCISE/$CODEDIR $TARDIR/C/$EXERCISE
cp -r $SRC_DIR/F/$EXERCISE/$CODEDIR $TARDIR/F/$EXERCISE


# Build tarball
tar -cf $TARBALL -C $TMPDIR mandelbrot/

# Remove temp directory
rm -rf $TMPDIR

