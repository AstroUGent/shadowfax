#! /bin/bash

if [ $# -lt 1 ]
then
echo "Usage: ./format_script <clang-format command>"
exit
fi

files=( src/*.cpp src/io/*.cpp \
        src/io/*.hpp src/utilities/*.cpp \
        src/utilities/*.hpp \
        src/*.hpp test/*.hpp test/*.cpp )

for f in "${files[@]}"
do $1 -style=file -i $f
done
