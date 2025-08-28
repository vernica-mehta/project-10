# Project 10 Group Logbook

## 26-08-2025

- Git repo has been set up, still ironing out some issues so that everyone can see everything
- Have discussed the problem and approaches, sketched out the plan and ideas as below

![scribble ](image.png)

## 28 - 08 - 2025 - minutes from class time 3-5pm

- Zotero project created so we can easily group all sources together
- Discussing possible methods for tackling project
- Textbook has grey atmospher approximation for externally irradiated atmospheres
- planet structure considered to be gas giant, Vernica provided a paper that discusses processes of externally irradiated atmosphere of gas giant
    - Model from paper takes non-LTE effects into account, might be too complex of a starting point
- Looking at provided code, specifically structure.py file. We plan too look over the code to see how it works
- We aim to neglect convection due to the problem being too complex to add, possibly might come back to it if we have time.

For pseudo-grey model, we need to know the temp, g and rho of planet. Albedo of the planet remains fixed and will come from literature. We aim to find how the temperature of the planet is affected by the external radiation field.

We discuss 3 main tasks to split between the three of this:

- Fixing the temperature code so that it works
- Finding literature values for temperature, g and density of exoplanets
- Changing code so that it becomes an exoplanet instead of a star