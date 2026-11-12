(* Mathematica visualization for G+Smo Planar Domain *)
ClearAll[gismoDomain, plotGismoDomain];

gismoDomain = <|
  "NumCorners" -> 6,
  "NumHoles" -> 0,
  "Corners" -> {{1, -1}, {1, 0}, {0, 0}, {0, 1}, {-1, 1}, {-1, -1}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{-1, -1}, {-0.333333, -0.825}, {0.333333, -0.825}, {1, -1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, -1}, {0.9125, -0.666667}, {0.9125, -0.333333}, {1, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, 0}, {0.666667, -0.0875}, {0.333333, -0.0875}, {0, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0, 0}, {-0.0875, 0.333333}, {-0.0875, 0.666667}, {0, 1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0, 1}, {-0.333333, 0.9125}, {-0.666667, 0.9125}, {-1, 1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1, 1}, {-0.825, 0.333333}, {-0.825, -0.333333}, {-1, -1}}|>
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
