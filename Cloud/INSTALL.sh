#!/bin/bash
# INSTALL.sh
# Install the dependencies required by MQLib to run on an Ubuntu EC2 instance

sudo apt-get update --yes > progress_A_$1.txt 2>&1
sudo apt-get install git build-essential python-pip python-paramiko --yes > progress_B_$1.txt 2>&1
# Assume .boto already copied
sudo pip install boto > progress_C_$1.txt 2>&1
rm -rf MQLib
echo "Progress report" > progress_D_$1.txt 2>&1
git clone $2 > progress_E_$1.txt 2>&1
cd MQLib
make > progress_F_$1.txt 2>&1
