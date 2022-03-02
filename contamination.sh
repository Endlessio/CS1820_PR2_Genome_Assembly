# This example shell script calls a Python script and passes in all arguments
# To locate your interpreter or compiler, run "whereis [interpreter]"
# Edit as necessary for different interpreters and/or file locations
# Note: "$@" is functionally-equivalent to "$1 $2 $3 ..."
# Come to hours or post on Piazza if you have questions about scripting!
python3 src/contamination.py $1 $2 $3
