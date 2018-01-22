# HDXSimulation
This project calculates hydrogen-deuterium exchange, optimizes parameters for a model and compares the results with the experimental data.


What is HDX:
HDX is a chemical reaction where a covalently bonded hydrogen atom is exchanged with a deuterium atom. In proteins, HDX refers to a process in which amide hydrogen atoms in the backbone are replaced by deuterium atoms. HDX-MS has been widely used in the literature as an excellent technique for understanding protein dynamics. HDX can be used to compare solvent exposure of distinct protein regions.
Flexible protein regions will undergo more frequent interruptions of secondary structure, allowing more frequent solvent exposure of backbone amide hydrogens and faster incorporation of solvent deuterium. To investigate and study HDX modeling, I used microsecond time-scale MD simulations for three isoforms of RGS proteins using CHARMM and AMBER force-fields. One of models that I used for HDX modeling for quantify HDX rates for RGS protein isoforms is based on Pesson model which the correspodning pythone codes are presented. 

Files:
There are two main files in the directroy "2015_persson_optimization.py" and "2015_persson.py" which the first one computes the parameters of Persson model using optimization algorithem and the second python code plots the comparison of experimental HDX with modeling (the resulting figures are saved in "figs" and the other outputs are saved in "outs"). The TCL code Sec. F.2 of my doctoral dissertation was used to analyze the trajectory of MD simulations and extract particular data from the trajectories (e.g. hydrogen bonds, distance to the first polar atoms, SASA values and etc.). The data were saved in "models" directory to be used by python code. Each python code reads corresponding data file to carry out further calculations and analysis. The python code for one model (out of 9 models in the dissertation) are shown in which command “np.loadtext” reads the corresponding data file.


