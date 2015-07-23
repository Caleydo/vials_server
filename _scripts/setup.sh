#!/usr/bin/env bash

#search for the right parent directory
while [ ! -f "Vagrantfile" ]
do
  cd ..
done

mkdir -p _data/
cd _data

if [ -d "vials" ]
then
  echo "vials directory already there"
else
  echo "downloading vials data"
  #TODO
  #baseurl="https://googledrive.com/host/0B7lah7E3BqlAfmNnQ3ptNUhtbG1fWklkemVGc0xnZkNyZ21lUi15aFlIb3NSZ2FWOTR3NHM/"
  #wget -O ccle.h5.gz "${baseurl}/ccle.h5.gz"
  #gunzip ccle.h5.gz
  #rm -f ccle.h5.gz
fi
