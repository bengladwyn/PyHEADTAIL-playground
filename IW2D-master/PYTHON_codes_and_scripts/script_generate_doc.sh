#!/bin/bash

# simple script to generate the doc from the commented functions
# in .py files
#
# to launch it:
# ./script_generate_doc.sh name_doc_file
# where 'name_doc_file' is the name of the file in which to put the doc
# (NOTE: this file should not be there or be empty)

out=$1
touch $out

for i in `ls *.py`; do

  j=${i/.py/}
  echo $j
  pydoc $j >> $out
  
done

