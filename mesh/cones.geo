// Gmsh project created on Wed May 11 13:18:09 2022
SetFactory("OpenCASCADE");
//+
Sphere(1) = {0, -4.95368-0.5, 0, 0.3, -Pi/2, Pi/2, Pi};
//+
//+
Rotate {{0, 0, 1}, {0, -4.95368-0.5, 0}, Pi} {
  Volume{1}; Surface{1}; Curve{5}; Surface{2}; Curve{4}; Curve{2}; Surface{3}; 
}
//+
Cone(2) = {0, 0, 0, 0, 2, 0, 0.1, 0.6, 2*Pi};
//+
Cone(3) = {0, -4.95368-0.5, 0, 0, 4.95368, 0, 0.3, 0.1, 2*Pi};
//+
Cylinder(4) = {0, -0.5, 0, 0, 0.5, 0, 0.1, 2*Pi};
//+
Sphere(5) = {0, 2, 0, 0.6, -Pi/2, Pi/2, Pi};
//+//+
Coherence;
//+
Coherence;
