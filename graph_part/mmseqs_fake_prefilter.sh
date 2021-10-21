#!/bin/bash

# copied from https://mmseqs.com/latest/userguide.pdf

fake_pref() {
QDB="$1"
TDB="$2"
RES="$3"
LDB="$4"
# create link to data file which contains a list of all targets that should be aligned
ln -s "${LDB}.index" "${RES}" # added .0 here, seems to fix it. default makes multiple files.
# create new index repeatedly pointing to same entry
INDEX_SIZE="$(echo $(wc -c < "${TDB}.index"))"
awk -v size=$INDEX_SIZE '{ print $1"\t0\t"size; }' "${QDB}.index" > "${RES}.index"
# create dbtype (7)
awk 'BEGIN { printf("%c%c%c%c",7,0,0,0); exit; }' > "${RES}.dbtype"
}

fake_pref "$@"