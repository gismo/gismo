(* Mathematica visualization for G+Smo Planar Domain *)
ClearAll[gismoDomain, plotGismoDomain];

gismoDomain = <|
  "NumCorners" -> 2,
  "NumHoles" -> 3,
  "Corners" -> {{1.2, 0}, {-1.2, 0}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{-1.2, 0}, {-0.48, 0.36}, {0.48, 0.36}, {1.2, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1.2, 0}, {0.6, 1.04}, {-0.6, 1.04}, {-1.2, 0}}|>
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
