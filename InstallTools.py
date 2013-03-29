#!/usr/bin/env bash


echo "Installing numerical optimization tools:"
sudo bash << EOF
apt-get install coinor-csdp
apt-get install coinor-libcbc-dev
apt-get install coinor-libcgl-dev
apt-get install coinor-libclp-dev 
apt-get install coinor-coinutils-dev 
apt-get install coinor-libdylp-dev
sudo apt-get install coinor-libflopc++-dev
apt-get install coinor-libipopt-dev
apt-get install coinor-libosi-dev
apt-get install coinor-libsymphony-dev 
apt-get install coinor-libvol-dev
apt-get install python-pip             
easy_install -U pulp
EOF
