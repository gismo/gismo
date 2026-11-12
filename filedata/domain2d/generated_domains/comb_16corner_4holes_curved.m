(* Mathematica visualization for G+Smo Planar Domain *)
ClearAll[gismoDomain, plotGismoDomain];

gismoDomain = <|
  "NumCorners" -> 16,
  "NumHoles" -> 4,
  "Corners" -> {{-1.2, -0.8}, {1.2, -0.8}, {1.2, 0.8}, {0.857143, 0.8}, {0.857143, 0}, {0.514286, 0}, {0.514286, 0.8}, {0.171429, 0.8}, {0.171429, 0}, {-0.171429, 0}, {-0.171429, 0.8}, {-0.514286, 0.8}, {-0.514286, 0}, {-0.857143, 0}, {-0.857143, 0.8}, {-1.2, 0.8}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{-1.2, 0.8}, {-1.2, 0.266667}, {-1.2, -0.266667}, {-1.2, -0.8}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1.2, -0.8}, {-0.4, -0.8}, {0.4, -0.8}, {1.2, -0.8}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1.2, -0.8}, {1.2, -0.266667}, {1.2, 0.266667}, {1.2, 0.8}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1.2, 0.8}, {1.08571, 0.8}, {0.971429, 0.8}, {0.857143, 0.8}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.857143, 0.8}, {0.857143, 0.533333}, {0.857143, 0.266667}, {0.857143, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.857143, 0}, {0.742857, 0}, {0.628571, 0}, {0.514286, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.514286, 0}, {0.514286, 0.266667}, {0.514286, 0.533333}, {0.514286, 0.8}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.514286, 0.8}, {0.4, 0.8}, {0.285714, 0.8}, {0.171429, 0.8}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.171429, 0.8}, {0.171429, 0.533333}, {0.171429, 0.266667}, {0.171429, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.171429, 0}, {0.0571429, 0}, {-0.0571429, 0}, {-0.171429, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.171429, 0}, {-0.171429, 0.266667}, {-0.171429, 0.533333}, {-0.171429, 0.8}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.171429, 0.8}, {-0.285714, 0.8}, {-0.4, 0.8}, {-0.514286, 0.8}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.514286, 0.8}, {-0.514286, 0.533333}, {-0.514286, 0.266667}, {-0.514286, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.514286, 0}, {-0.628571, 0}, {-0.742857, 0}, {-0.857143, 0}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.857143, 0}, {-0.857143, 0.266667}, {-0.857143, 0.533333}, {-0.857143, 0.8}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.857143, 0.8}, {-0.971429, 0.8}, {-1.08571, 0.8}, {-1.2, 0.8}}|>
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
