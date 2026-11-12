(* Mathematica visualization for G+Smo Planar Domain *)
ClearAll[gismoDomain, plotGismoDomain];

gismoDomain = <|
  "NumCorners" -> 4,
  "NumHoles" -> 1,
  "Corners" -> {{6.73556e-17, 1.1}, {-1.1, 1.34711e-16}, {-2.02067e-16, -1.1}, {1.1, -2.69422e-16}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{1.1, 0}, {1.38399, 0.433867}, {0.433867, 0.816014}, {6.73556e-17, 1.1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{6.73556e-17, 1.1}, {-0.433867, 1.38399}, {-0.816014, 0.433867}, {-1.1, 1.34711e-16}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1.1, 1.34711e-16}, {-1.38399, -0.433867}, {-0.433867, -0.816014}, {-2.02067e-16, -1.1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-2.02067e-16, -1.1}, {0.433867, -1.38399}, {0.816014, -0.433867}, {1.1, -2.69422e-16}}|>
  }
|>;

bernsteinPoly[p_, j_, t_] := Binomial[p, j] * t^j * (1 - t)^(p - j);
evalBSpline[curve_, t_] := Module[{p = curve["Degree"], cps = curve["ControlPoints"]},
  Sum[cps[[j + 1]] * bernsteinPoly[p, j, t], {j, 0, p}]
];

Graphics[{
  {Thick, Blue, Table[ParametricPlot[evalBSpline[c, t], {t, 0, 1}][[1]], {c, gismoDomain["OuterLoop"]}]},
  {Red, PointSize[0.025], Point[gismoDomain["Corners"]]}
}, Axes -> True, GridLines -> Automatic, AspectRatio -> Automatic]
