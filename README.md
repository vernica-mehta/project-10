# ASTR4022 Project: Computing the emergent spectrum of an externally-irradiated atmosphere

**Project Number:** 10

**Group Members:** Griffin, Nathan, Vernica

## Repo Organisation

We have approached this project by applying a combination of analytical and numerical models. The analytical code, done by Nathan, is in the `non-grey_model` folder and focuses on retrieving a grid of values describing the temperatures of layers within the exoplanet atmosphere as a function of optical depths / opacities, exploring both semi-grey and non-grey models. The main code file here is `non_grey.ipynb`. The numerical code, done by Vernica, is in the `grey_model` folder and applies a grey atmosphere model to develop a simplified exoplanet spectrum with irradiation from a few initial conditions. The main code file here is `grey_model.py`. There is also a related subfolder `absorption_lines` which contains files related to some attempts to incorporate absorption features into the grey model spectra. The `EoS` subfolder within the grey model folder, contains code written mostly by Griffin to calculate equations of state for a molecular atmospheric composition and opacity calculations which follow on from that. The data files used within the grey model calculations are also in this subfolder.

Individual logbooks can be found within the `Logbooks` folder, under named subfolders. The file `Logbook.md` contains some notes on early group discussions from the first few weeks of the project, before the bulk of the research was divided up. General notes on some literature sources written by Griffin can be found in the `Literature notes` folder with related figures in `Image`.
