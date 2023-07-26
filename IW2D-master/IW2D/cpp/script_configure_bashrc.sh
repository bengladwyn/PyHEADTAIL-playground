#!/bin/bash

echo "Configuring bashrc to allow Python scripts to use IW2D library and exectuables\n and add this folder to the PATH"
echo "(NOTE: WORKS ONLY IF YOU USE BASH SHELL!)"
echo $'\n# paths for ImpedanceWake2D' >> ~/.bashrc
echo "LD_LIBRARY_PATH="`pwd`":\$LD_LIBRARY_PATH" >> ~/.bashrc
echo "export LD_LIBRARY_PATH" >> ~/.bashrc
echo "PATH="`pwd`":\$PATH" >> ~/.bashrc
echo "export PATH" >> ~/.bashrc
echo "IW2D_PATH="`pwd` >> ~/.bashrc
echo "export IW2D_PATH" >> ~/.bashrc

echo $'\nDone\n'

echo $'DON\'T FORGET TO DO:\nsource ~/.bashrc'
