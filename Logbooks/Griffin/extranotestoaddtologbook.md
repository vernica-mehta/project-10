

Existing chemical models for transiting exoplanets tend to fall
into two groups: 
- thermochemical-equilibrium models
- photochemical models.

Thermochemical equilibrium is a reasonable starting point for exoplanet composition predictions, but the strong ultraviolet flux incident on close-in transiting planets like HD 189733b and HD 209458b can drive the observable upper regions of the atmosphere out of equilibrium. Moreover,



The Venus model is based on Caltech/JPL-KINETICS code which solves the continuity equation for a particular species i of a 1-D model, 
\\
$$\dfrac{\partial n_i}{\partial t} + \dfrac{\partial \Phi_i}{\partial z} = P_i -L_i$$ 
where $n_i, \Phi_i,P_i, \textrm{and} \ L_i$ are the concentrations, vertical diffusive flux, chemical production, and chemical loss terms, for species i \cite{allenVerticalTransportPhotochemistry1981} \cite{millsObservationsPhotochemicalModeling1998}.
\\
\\
The vertical diffusive flux, $\Phi_i$ is given by: 
\\
$$\Phi_i = -D_i(\dfrac{\mathrm{d}n_i}{\mathrm{d}z}+\frac{n_i}{H_i}+\frac{n_i(1+\alpha)}{T}\dfrac{\mathrm{dT}}{\mathrm{d}z})-K(\dfrac{\mathrm{d}n_i}{dz}+\frac{n_i}{H_i}+\frac{n_i}{T}\dfrac{\mathrm{dT}}{\mathrm{d}z})$$ Where $D_i$ is the molecular diffusion coefficient of species i through the background atmosphere, $H_i$ is the scale height of species i, T is the temperature, $\alpha$ is the thermal diffusion factor, K is the eddy diffusion coefficient, and H is the scale height of the background atmosphere \cite{allenVerticalTransportPhotochemistry1981} \cite{millsObservationsPhotochemicalModeling1998}. Background atmosphere is the temperature and pressure profiles and remain fixed throughout the simulations run.