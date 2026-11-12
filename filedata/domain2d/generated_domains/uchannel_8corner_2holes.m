(* Mathematica visualization for G+Smo Planar Domain *)
ClearAll[gismoDomain, plotGismoDomain];

gismoDomain = <|
  "NumCorners" -> 8,
  "NumHoles" -> 2,
  "Corners" -> {{-1.2, -1}, {1.2, -1}, {1.2, 1}, {0.5, 1}, {0.5, 0}, {-0.5, 0}, {-0.5, 1}, {-1.2, 1}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{-1.2, 1}, {-1.2, 0.333333}, {-1.2, -0.333333}, {-1.2, -1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1.2, -1}, {-0.4, -1}, {0.4, -1}, {1.2, -1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1.2, -1}, {1.2, -0.333333}, {1.2, 0.333333}, {1.2, 1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1.2, 1}, {0.966667, 1}, {0.733333, 1}, {0.5, 1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.5, 1}, {0.5, 0.666667}, {0.5, 0.333333}, {0.5, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.5, 0}, {0.166667, 0}, {-0.166667, 0}, {-0.5, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.5, 0}, {-0.5, 0.333333}, {-0.5, 0.666667}, {-0.5, 1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.5, 1}, {-0.733333, 1}, {-0.966667, 1}, {-1.2, 1}}|>
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
