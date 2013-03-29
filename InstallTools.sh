#!/usr/bin/env bash


echo "Installing numerical optimization tools:"
sudo bash <<EOF
apt-get --force-yes -y install coinor-csdp
apt-get --force-yes -y install coinor-libcbc-dev
apt-get --force-yes -y install coinor-libcgl-dev
apt-get --force-yes -y install coinor-libclp-dev 
apt-get --force-yes -y install coinor-coinutils-dev 
apt-get --force-yes -y install coinor-libdylp-dev
apt-get --force-yes -y install coinor-libflopc++-dev
apt-get --force-yes -y install coinor-libipopt-dev
apt-get --force-yes -y install coinor-libosi-dev
apt-get --force-yes -y install coinor-libsymphony-dev 
apt-get --force-yes -y install coinor-libvol-dev
apt-get --force-yes -y install python-pip             
easy_install -U pulp 
EOF

