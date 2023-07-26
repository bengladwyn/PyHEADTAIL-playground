#!/bin/bash

echo "Configuring bashrc to allow Python scripts to use the Python tools in this directory"
echo "(NOTE: WORKS ONLY IF YOU USE BASH SHELL!)"

echo $'\n# Python path for General_Python_tools' >> ~/.bashrc
echo "PYTHONPATH="`pwd`":\$PYTHONPATH" >> ~/.bashrc
echo "export PYTHONPATH" >> ~/.bashrc

echo $'\nDone\n'

echo $'DON\'T FORGET TO DO:\nsource ~/.bashrc'

