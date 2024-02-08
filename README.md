# Phytoplankton population dynamics in a vertical freshwater column
Phytoplankton is a significant primary producer and at the base of the aquatic food chain, and understanding its dynamics can help form a basis for predictions about entire marine ecosystems. 

A 1D reaction-advection-diffusion equation reduces the complex list of factors which affect some population to it’s patch size, production rate, death rate, and rates of spreading (advection and diffusion).


<img src="ProjectDiagram.jpg" width="300" height="250" /> 

Phytoplankton take in dissolved nutrients (white dots) at a rate $\rho$, converting them to particulate nutrients (black dots) bound in the phytoplankton. The bound nutrients diffuse up and down at a rate d and sink at a rate $v$ with the phytoplankton cloud. When bound nutrients reach $z_\max$ they become sedimented as $R_s$. Sedimented nutrients are released back into the water column as dissolved nutrients at a rate $r$. The dissolved nutrients also diffuse up and down at a rate d but do not sink.