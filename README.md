# _aztekas_-code 

_aztekas_ is a program that solves hyperbolic partial differential equations in conservative form using High Resolution Shock-Capturing (HRSC) schemes. This version of _aztekas_ allows to solve the non-relativistic and relativistic hydordyamic equations of motion (Euler equations for a perfect fluid. The relativistic part allows to solve these equations on a background fixed metric (Schwarzschild, Minkowski, Kerr-Schild, etc.).

## Introduction

_aztekas_ uses a conservative finite-volume approach to obtain the dicrete form of a hyperbolic partial differential system of equations (PDE).

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{\partial&space;Q}{\partial&space;t}&space;&plus;&space;\frac{\partial&space;F^i}{\partial&space;x^i}&space;=&space;S" target="_blank"><img src="https://latex.codecogs.com/pdf.latex?\frac{\partial&space;Q}{\partial&space;t}&space;&plus;&space;\frac{\partial&space;F^i}{\partial&space;x^i}&space;=&space;S" title="\frac{\partial Q}{\partial t} + \frac{\partial F^i}{\partial x^i} = S" /></a>
