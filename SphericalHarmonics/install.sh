#!/bin/bash
cd Definitions
python setup.py install
cd ..
cd NSSHT
python setup.py install
cd ..
cd SHT
./install.sh
python setup.py install
cd ..