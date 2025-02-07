# Free metal calculator

![License](https://img.shields.io/github/license/joel-rosenbaum/free_metal)
![GitHub last commit](https://img.shields.io/github/last-commit/joel-rosenbaum/free_metal)
![GitHub issues](https://img.shields.io/github/issues/joel-rosenbaum/free_metal)
![GitHub stars](https://img.shields.io/github/stars/joel-rosenbaum/free_metal?style=social)

## Overview

The **free metal calculator** is an interactive web application built with Dash that is designed to visualize metal ion buffering under standard conditions. The computational approach was borrowed entirely from [WEBMAXC](https://somapp.ucdmc.ucdavis.edu/pharmacology/bers/maxchelator/webmaxc/webmaxcS.htm) and so the results should be very similar. This appplication extends on the capabilities of that tool by providing data visualizatioin and extensibility through configurable JSON files. 

Updated to include adjustments for metal hydrolysis.

## Notes

Please note that the limitations of **WEBMAXC** still apply:

> **NOTE:** It is possible with the **FIND FREE METAL** calculation to produce answers that are not valid. This is more likely to occur with multiple metals. If you find one of the metals shows almost no effect of the chelator presence, be suspicious. Try different concentrations of the second metal. If there is a point where the result suddenly changes the first metal, then you are having problems with a convergence error. Try the Windows version. Here double precision math is used and it is also able to warn specifically about convergence problems.

## Attribution

This project uses methods from **WEBMAXC** originally developed by Chris Patton. Any citations should use the references found [here](https://somapp.ucdmc.ucdavis.edu/pharmacology/bers/maxchelator/references.htm).

## 📂 File Structure

### `app.py`
The main entry point of the Dash application. It initializes the app, defines the layout, and sets up callbacks to handle user interactions.

### `calculations.py`
Contains the core computational logic for modeling metal-chelator equilibria.<br>
-`Calculator` class implements functions to compute equilibrium concentrations and stability constants.<br>
-`ConstantUtils` class provides utility functions related to chemical constants and conversions.

### `processing.py`
Handles data preprocessing tasks.<br>
-`Preprocessing` class includes methods to structure inputs.

### `plotting.py`
Plotting functions based on user inputs.<br>
-`Plotting` class generates Plotly graphs based on calculated data.

### `webmaxc_chelator.json` and `webmaxc_protonation.json`
Provides binding and acid dissociation constants extracted from WEBMAXC.

### `perrin_metals.json`
Provides pKa values for metal hydrolysis (from Perrin, 1969)
