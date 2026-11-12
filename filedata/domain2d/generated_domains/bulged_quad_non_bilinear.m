(* Mathematica visualization for G+Smo Planar Domain *)
ClearAll[gismoDomain, plotGismoDomain];

gismoDomain = <|
  "NumCorners" -> 4,
  "NumHoles" -> 0,
  "Corners" -> {{6.12323e-17, 1}, {-1, 1.22465e-16}, {-1.83697e-16, -1}, {1, 0}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{1, 0}, {0.509167, 0.175833}, {0.175833, 0.509167}, {6.12323e-17, 1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{6.12323e-17, 1}, {-0.175833, 0.509167}, {-0.509167, 0.175833}, {-1, 1.22465e-16}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1, 1.22465e-16}, {-0.509167, -0.175833}, {-0.175833, -0.509167}, {-1.83697e-16, -1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1.83697e-16, -1}, {0.175833, -0.509167}, {0.509167, -0.175833}, {1, 0}}|>
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
