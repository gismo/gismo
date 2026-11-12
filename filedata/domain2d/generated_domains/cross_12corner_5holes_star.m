(* Mathematica visualization for G+Smo Planar Domain *)
ClearAll[gismoDomain, plotGismoDomain];

gismoDomain = <|
  "NumCorners" -> 12,
  "NumHoles" -> 5,
  "Corners" -> {{-0.4, 0.4}, {-1, 0.4}, {-1, -0.4}, {-0.4, -0.4}, {-0.4, -1}, {0.4, -1}, {0.4, -0.4}, {1, -0.4}, {1, 0.4}, {0.4, 0.4}, {0.4, 1}, {-0.4, 1}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{-0.4, 1}, {-0.4, 0.8}, {-0.4, 0.6}, {-0.4, 0.4}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.4, 0.4}, {-0.6, 0.4}, {-0.8, 0.4}, {-1, 0.4}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1, 0.4}, {-1, 0.133333}, {-1, -0.133333}, {-1, -0.4}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1, -0.4}, {-0.8, -0.4}, {-0.6, -0.4}, {-0.4, -0.4}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.4, -0.4}, {-0.4, -0.6}, {-0.4, -0.8}, {-0.4, -1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.4, -1}, {-0.133333, -1}, {0.133333, -1}, {0.4, -1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.4, -1}, {0.4, -0.8}, {0.4, -0.6}, {0.4, -0.4}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.4, -0.4}, {0.6, -0.4}, {0.8, -0.4}, {1, -0.4}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, -0.4}, {1, -0.133333}, {1, 0.133333}, {1, 0.4}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, 0.4}, {0.8, 0.4}, {0.6, 0.4}, {0.4, 0.4}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.4, 0.4}, {0.4, 0.6}, {0.4, 0.8}, {0.4, 1}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.4, 1}, {0.133333, 1}, {-0.133333, 1}, {-0.4, 1}}|>
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
