#!/bin/bash

if [ "$#" -ne 1 ]
then
  echo "Usage: Specify version of the benchmark [debug,release,intel]"
  exit 1
fi

version=$1
directory="../cmake-build-${version}/benchmarks"

for file in data/*.nwk;
    do
    filename=$(basename "$file")
    extension="${filename##*.}"
    filename="${filename%.*}"
    for i in {1..10};
    do
        #echo "${directory}/tshlib_benchmarks ${file}"
        ${directory}/tshlib_benchmarks ${file}
    done
    #echo "mv tshlib_statistics.csv ${filename}_${version}.csv"
    mv tshlib_statistics.csv ${filename}_${version}.csv
done

./generateBenchmarkPlots.R
