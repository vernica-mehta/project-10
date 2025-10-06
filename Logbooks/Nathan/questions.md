# Progress Report for Nathan

### Initial semi-grey atmosphere model

Initial work was done on a semi-grey atmosphere model from the textbook, first was understanding the underlying equations and theory before trying to code as there were some questions with the final model given (mainly the dilution factor W). The model takes a two stream approach and uses the moment equations to get some first order DEs to solve for the mean intensity as a function of optical depth and thus a temperature profile.

Additional reading into the model and we found that the dilution factor W was found to be $W = (R/D)^2 \times f$, where f is an angular factor $f = \mu^2$

The two-stream approximation was used for the pseudo-grey model but as the model was semi-grey, it used only mean opacities which is a large limitation, so while this was coded up, we also needed to look for something more complex. Vernica and Griff did the coding and I found a analytical solution for a non-grey atmosphere using the picket fence model.

### Non-grey atmosphere

The paper below details the process for an analytical solution for a non-grey atmosphere using the picket fence model and compares it to more simple grey models that use discrete ordinance methods and to numerical solutions.

https://www.aanda.org/articles/aa/full_html/2014/02/aa22342-13/aa22342-13.html

My understanding of the picket fence model is that it introduces two more opacities, one for the lines and one for the continuum to act as the thermal opacities for the two-stream approximation. This is to simulate the effect of line blanketing without needing to model all of the different lines numerically which would be a project in itself. The paper introduces the picket fence model as numerical solutions have shown that there is cooling effects in the low optical depth regions, from 0 to 1, though the high optical depths are modelled well by the semi-grey model. They introduce $\kappa_1$ as the line opacity and $\kappa_2$ as the continuum opacity, with $\kappa_1 > \kappa_2$. They convert most of the opacities into ratios of the rosseland mean opacity which allows for some simplification of the equations used.

Another opacity they use is the opacity for stellar radiation, $\kappa_V$ which is the main part of confusion for me as I am trying to understand how this variable is determined so we can use this model. I'll discuss my issues in the questions section further down.

Similar to the semi-grey model, they use the moment equations and some approximations to get a set of 3 coupled equations, specifically two second order ODEs for the moment equations and one first order ODE for the radiative equilibrium equation. They decople the equations by introducig new variables, for brevity I will not go into the full derivation here.

The paper sets the boundary conditions as the diffuesion approximation for the high optical depths, i.e.:

$$\lim_{\tau \to \infty} J(\tau) = B$$

For the low optical depths, they set the boundary condition as:
$$J(0) = 2H(0)$$
and the system is in LTE. They solve the ODEs by using the constant coefficients method, this gives a general solution for the mean intensity as a function of optical depth and is converted to a temperature profile:
$$ T(\tau) = \frac{3}{4}T_{int}^4\left(\tau + A + Be^{-\tau/\tau_{lim}}\right) + \frac{3}{4}T_{irr}^4\mu\left(C + De^{-\tau/\tau_{lim}} + Ee^{-\gamma_v \tau}\right)$$

Where A, B, C, D, E are coefficients defined in the paper, they are complex due to the number of terms but are made up entirely of the opacity ratios and the limit of optical depth defined by Chandraskhar for the picket fence model:

$$\gamma_1 = \kappa_1/\kappa_R$$
$$\gamma_2 = \kappa_2/\kappa_R$$
$$\gamma_V = \kappa_V/\kappa_R$$
$$\gamma_P = \kappa_P/\kappa_R$$
$$\kappa_P = \frac{\int_0^\infty \kappa_\nu B_\nu d\nu}{\int_0^\infty B_\nu d\nu}$$
$$\tau_{lim} = \frac{1}{\gamma_1\gamma_2}\sqrt{\frac{\gamma_P}{3}}$$



### Questions

1. $\kappa_V$ is described in the table of quantities of the paper on page 16 as the opacity in the visible.How is $\kappa_V$ determined? We modified the opac.py code for the calculation of the opacity $\kappa_\nu$ for the semi-grey atmosphere but I am not sure if this is the same as $\kappa_V$. Vernica mentioned that it might have some issues with the model as how opac.py calculates the opacity. She mentioned that it uses a kappa matrix to calculate optical depths which then produces the kappa_nu array, so there is 
2. In your response on Tuesday, you mentioned trying a "Toy atmosphere" where we add lines to the semi-grey approximation we used. Could you elaborate on this? Would this be similar to the picket fence model but without needing to modify the equations too much such that it becomes a second order ODE like the paper mentioned?
3. For the semi-grey model, for $T_{eff}$ we used the effective surface temperature of jupiter (144K), but the paper makes $T_{eff}$ as the sum of internal and irradiation temperatures, i.e. $T_{eff}^4 = T_{int}^4 + \mu T_{irr}^4$. For the semi-grey model we made the irradiation temperature the effective temperature of the star, approximately 6500K for solar radiation. The definition in the paper seems to imply that irradiation temperature would need to be lower to provide a reasonable effective temperature, because if we used our choice of irradiation temperature for the paper, the jupiter-like planet would have a temperature of ~6500K if the internal temperature is 0K and $\mu = 1$. This obviously is way too high, so are we misunderstanding the definition of irradiation temperature or is the paper defining it differently?