#!/usr/bin/env bash

#
# Run this script to download relevant data files,
# which are not committed to the repository.
#

server="https://starfit.org"
data_dir="src/starfit/data"
data_dirs="db stars ref"

for d in $data_dirs; do

  dir="${data_dir}/${d}"

  echo "==> Downloading data to: $dir/"
  cd $dir
  while read line; do
    file=$(echo "$line" | cut -d ' ' -f 3)
    echo $file
    curl -fLO# "${server}/data/${d}/${file}"
  done < .hashlist

  echo "==> Inspecting checksums..."
  shasum -a 256 -c .hashlist
  echo; cd - > /dev/null

done
