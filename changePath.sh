#!/bin/bash

PWD_ESCAPED=$(sed 's/\//\\\//g' <<< $PWD) 
for annot in `ls output`; do sed -e "s/TAMIR_INSTALL_ROOT/$PWD_ESCAPED/g" output/$annot/${annot}_seg_correct.annotation.template > output/$annot/${annot}_seg_correct.annotation; done
