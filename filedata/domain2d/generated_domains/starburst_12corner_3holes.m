(* Mathematica visualization for G+Smo Planar Domain *)
ClearAll[gismoDomain, plotGismoDomain];

gismoDomain = <|
  "NumCorners" -> 12,
  "NumHoles" -> 3,
  "Corners" -> {{0.33775, 0.195}, {0.58, 1.00459}, {0, 0.39}, {-0.58, 1.00459}, {-0.33775, 0.195}, {-1.16, 1.38778e-16}, {-0.33775, -0.195}, {-0.58, -1.00459}, {-1.11022e-16, -0.39}, {0.58, -1.00459}, {0.33775, -0.195}, {1.16, 0}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{1.16, 0}, {0.885917, 0.065}, {0.611833, 0.13}, {0.33775, 0.195}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.33775, 0.195}, {0.4185, 0.464863}, {0.49925, 0.734726}, {0.58, 1.00459}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.58, 1.00459}, {0.386667, 0.799726}, {0.193333, 0.594863}, {0, 0.39}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{2.38806e-17, 0.39}, {-0.193333, 0.594863}, {-0.386667, 0.799726}, {-0.58, 1.00459}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.58, 1.00459}, {-0.49925, 0.734726}, {-0.4185, 0.464863}, {-0.33775, 0.195}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.33775, 0.195}, {-0.611833, 0.13}, {-0.885917, 0.065}, {-1.16, 1.38778e-16}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1.16, 1.42059e-16}, {-0.885917, -0.065}, {-0.611833, -0.13}, {-0.33775, -0.195}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.33775, -0.195}, {-0.4185, -0.464863}, {-0.49925, -0.734726}, {-0.58, -1.00459}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.58, -1.00459}, {-0.386667, -0.799726}, {-0.193333, -0.594863}, {-1.11022e-16, -0.39}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-7.16418e-17, -0.39}, {0.193333, -0.594863}, {0.386667, -0.799726}, {0.58, -1.00459}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.58, -1.00459}, {0.49925, -0.734726}, {0.4185, -0.464863}, {0.33775, -0.195}}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.33775, -0.195}, {0.611833, -0.13}, {0.885917, -0.065}, {1.16, 0}}|>
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
