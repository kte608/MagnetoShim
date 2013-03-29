This is how the program is supposed to work on ubuntu Linux:

1) You run the script 'InstallTools.sh' which supposedly installs all the pre-requisites for the program to run.
2) You run the program 'InkShim.py' with appropriate arguments. Running it without arguments or with --help should give you the information that you need.
3) The resulting .cPickle file can be fed as an argument to 'plotProgram.py' which shows the ink pattern and produces a .ps file.
4) Optionall use the program ps2pdf to make a pdf file. (apt-get install ps2pdf)
