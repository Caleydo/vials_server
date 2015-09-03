#!/usr/bin/env bash

#search for the right parent directory such that we have a common start directory
while [[ ! -f "run.sh" ]] && [[ ! -f "Vagrantfile" ]]
do
  cd ..
done


function setup {
  echo "setup"
}

function update {
  echo "update"
}

function uninstall {
  echo "uninstall"
}

#command switch
case "$1" in
update)
  update
  ;;
uninstall)
  uninstall
  ;;
*)
  setup
  ;;
esac