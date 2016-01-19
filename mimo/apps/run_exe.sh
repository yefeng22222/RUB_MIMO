#!/bin/bash
INDEX=10
while [ 100 -gt $INDEX ]; do
  './main.exe' >> output-$INDEX.txt
  let INDEX=INDEX+1
done
