#!/bin/bash

echo "Configuring bashrc to allow Python scripts to use the Impedance library in this directory"
echo "(NOTE: WORKS ONLY IF YOU USE BASH SHELL!)"

echo $'\n# paths for Impedance_lib_Python' >> ~/.bashrc
echo "PYTHONPATH="`pwd`":\$PYTHONPATH" >> ~/.bashrc
echo "export PYTHONPATH" >> ~/.bashrc
echo "YOK_PATH="`pwd` >> ~/.bashrc
echo "export YOK_PATH" >> ~/.bashrc

echo $'\nDone\n'

echo $'DON\'T FORGET TO DO:\nsource ~/.bashrc'

