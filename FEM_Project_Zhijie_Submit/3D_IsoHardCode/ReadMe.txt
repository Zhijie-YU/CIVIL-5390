two_element.m is the example file for a 3D problem which can be run directly with Matlab to get the results.
NLFEA.m is the main program to this problem where N-R iteration are implemented.
PLAST3D.m is the function to compute and assemble global tangent stiffness matrix and residual force vector.
SHAPEL.m is the function to compute variables related to shape functions.
combHard.m is the function to implement the isotropic hardening model to get new stress and history variables given strain increment. Elastic assumption and return mapping algorithm are included here.
combHardTan.m is the function similar to comHard and is used to get new tangent operator.