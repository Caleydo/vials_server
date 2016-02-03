#!/usr/bin/env bash

#search for the right parent directory such that we have a common start directory
while [[ ! -f "run.sh" ]] && [[ ! -f "Vagrantfile" ]]
do
  cd ..
done


DATA_PROJECT_DIR="bodymap.vials_project"
DATA_PROJECT_FILE="${DATA_PROJECT_DIR}.zip"
#https://www.dropbox.com/s/c99p8bk65zdrsf6/bodymap.vials_project.tar.gz?dl=1
DATA_PROJECT_URL="https://www.dropbox.com/s/xrbs250tjafjvpd/${DATA_PROJECT_FILE}?dl=1"
DATA_REF_DIR="reference_genomes"
DATA_REF_FILE="${DATA_REF_DIR}.zip"
DATA_REF_URL="https://www.dropbox.com/s/zoqnihdrhony4bh/${DATA_REF_DIR}?dl=1"

mkdir -p _data/vials_projects
cd _data/vials_projects

function update_project_file {
  echo "downloading ${DATA_PROJECT_FILE} file"
  wget --timestamping -O "${DATA_PROJECT_FILE}" "${DATA_PROJECT_URL}"
  rm -rf "${DATA_PROJECT_DIR}"
  unzip "${DATA_PROJECT_FILE}"
  rm -f "${DATA_PROJECT_FILE}"
}
function update_ref_file {
  echo "downloading ${DATA_REF_FILE} file"
  wget --timestamping -O "${DATA_REF_FILE}" "${DATA_REF_URL}"
  rm -rf "${DATA_REF_DIR}"
  unzip "${DATA_REF_FILE}"
  rm -f "${DATA_REF_FILE}"
}

function setup {
  if [ -d "${DATA_PROJECT_FILE}" ]
  then
    echo "${DATA_PROJECT_FILE} already there"
  else
    update_project_file
  fi
  if [ -d "${DATA_REF_DIR}" ]
  then
    echo "${DATA_REF_DIR} already there"
  else
    update_ref_file
  fi
}

function update {
  #update_file
  echo "nothing to update"
}

function uninstall {
  rm -rf "${DATA_PROJECT_DIR}"
  rm -rf "${DATA_REF_DIR}"
}


#command switch
case "$1" in
update)
  update
  ;;
updatedata)
  update_project_file
  update_ref_file
  ;;
updateproject)
  update_project_file
  ;;
updateref)
  update_ref_file
  ;;
uninstall)
  uninstall
  ;;
*)
  setup
  ;;
esac
