# _aztekas_-code 

_aztekas_ is a program that solves hyperbolic partial differential equations in conservative form using High Resolution Shock-Capturing (HRSC) schemes. This version of _aztekas_ allows to solve the non-relativistic and relativistic hydordyamic equations of motion (Euler equations for a perfect fluid. The relativistic part allows to solve these equations on a background fixed metric (Schwarzschild, Minkowski, Kerr-Schild, etc.).

## Introduction

_aztekas_ uses a conservative finite-volume approach to obtain the dicrete form of a hyperbolic partial differential system of equations (PDE).

<img src="http://www.sciweavers.org/tex2img.php?eq=%20%5Cfrac%7B%5Cpartial%20Q%7D%7B%5Cpartial%20t%7D%20%20%2B%20%5Cfrac%7B%5Cpartial%20F%5Ei%7D%7B%5Cpartial%20x%5Ei%7D%20%3D%20S&bc=White&fc=Black&im=jpg&fs=12&ff=txfonts&edit=0" align="center" border="0" alt=" \frac{\partial Q}{\partial t}  + \frac{\partial F^i}{\partial x^i} = S" width="104" height="39" />

where <img src="http://www.sciweavers.org/tex2img.php?eq=Q&bc=White&fc=Black&im=jpg&fs=12&ff=txfonts&edit=0" align="center" border="0" alt="Q" width="19" height="17" /> is the vector of conservative variables, <img src="http://www.sciweavers.org/tex2img.php?eq=F%5Ei&bc=White&fc=Black&im=jpg&fs=12&ff=txfonts&edit=0" align="center" border="0" alt="F^i" width="19" height="17" /> are the flux vectors of the conservative charges and <img src="http://www.sciweavers.org/tex2img.php?eq=S&bc=White&fc=Black&im=jpg&fs=12&ff=txfonts&edit=0" align="center" border="0" alt="S" width="17" height="14" /> is the source term vector.
