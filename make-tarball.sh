#!/bin/sh

SRC_DIR=exercises
HANDOUTS_DIR=handouts

EXERCISE=$1

if [[ $EXERCISE == "" ]]; then
    echo "You need to tell me which exercise to tar up (exercise1-5)"
    exit 1
fi

# Find PDF handout for specified exercise
HANDOUT=$HANDOUTS_DIR/$EXERCISE.pdf

TARBALL=$EXERCISE.tar
if [ ! -f $HANDOUT ]; then
    echo "Unable to find handout PDF for $EXERCISE, have you written it?"
    exit 1
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

cp $HANDOUT $TARDIR

mkdir $TARDIR/C
mkdir $TARDIR/F

# Clean up executables and output files (we don't want them)
(cd $SRC_DIR/C/$EXERCISE/template; make clean 2>/dev/null; rm -f output.ppm)>/dev/null
(cd $SRC_DIR/F/$EXERCISE/template; make clean 2>/dev/null; rm -f output.ppm)>/dev/null

cp -r $SRC_DIR/C/utils $TARDIR/C/

cp -r $SRC_DIR/F/utils $TARDIR/F/

mkdir -p $TARDIR/C/$EXERCISE
mkdir -p $TARDIR/F/$EXERCISE

cp $SRC_DIR/config.mak $TARDIR/
cp -r $SRC_DIR/C/$EXERCISE/template $TARDIR/C/$EXERCISE
cp -r $SRC_DIR/F/$EXERCISE/template $TARDIR/F/$EXERCISE


# Build tarball
tar -cf $TARBALL -C $TMPDIR mandelbrot/

# Remove temp directory
rm -rf $TMPDIR

