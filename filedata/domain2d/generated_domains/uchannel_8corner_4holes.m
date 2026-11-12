(* ========================================================================= *)
(* G+Smo Planar Domain Definition & Mathematica Visualization Script        *)
(* Domain: corners=8, convex=False, star=False, holes=4 *)
(* ========================================================================= *)

ClearAll[gismoDomain, plotGismoDomain];

(* 1. Domain Specification *)
gismoDomain = <|
  "NumCorners" -> 8,
  "IsConvex" -> false,
  "IsStarShaped" -> false,
  "NumHoles" -> 4,
  "Corners" -> {{-1, -1}, {1, -1}, {1, 1}, {0.6, 1}, {0.6, -0.6}, {-0.6, -0.6}, {-0.6, 1}, {-1, 1}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{-1, -1}, {-0.33333333, -1}, {0.33333333, -1}, {1, -1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, -1}, {1, -0.33333333}, {1, 0.33333333}, {1, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, 1}, {0.86666667, 1}, {0.73333333, 1}, {0.6, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.6, 1}, {0.6, 0.46666667}, {0.6, -0.066666667}, {0.6, -0.6}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.6, -0.6}, {0.2, -0.6}, {-0.2, -0.6}, {-0.6, -0.6}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.6, -0.6}, {-0.6, -0.066666667}, {-0.6, 0.46666667}, {-0.6, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.6, 1}, {-0.73333333, 1}, {-0.86666667, 1}, {-1, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1, 1}, {-1, 0.33333333}, {-1, -0.33333333}, {-1, -1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
  },
  "Holes" -> {
    {
      <|"Degree" -> 3, "ControlPoints" -> {{-0.828, 0.53}, {-0.80933333, 0.53}, {-0.79066667, 0.53}, {-0.772, 0.53}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.772, 0.53}, {-0.772, 0.51266667}, {-0.772, 0.49533333}, {-0.772, 0.478}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.772, 0.478}, {-0.75466667, 0.478}, {-0.73733333, 0.478}, {-0.72, 0.478}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.72, 0.478}, {-0.72, 0.45933333}, {-0.72, 0.44066667}, {-0.72, 0.422}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.72, 0.422}, {-0.73733333, 0.422}, {-0.75466667, 0.422}, {-0.772, 0.422}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.772, 0.422}, {-0.772, 0.40466667}, {-0.772, 0.38733333}, {-0.772, 0.37}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.772, 0.37}, {-0.79066667, 0.37}, {-0.80933333, 0.37}, {-0.828, 0.37}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.828, 0.37}, {-0.828, 0.38733333}, {-0.828, 0.40466667}, {-0.828, 0.422}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.828, 0.422}, {-0.84533333, 0.422}, {-0.86266667, 0.422}, {-0.88, 0.422}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.88, 0.422}, {-0.88, 0.44066667}, {-0.88, 0.45933333}, {-0.88, 0.478}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.88, 0.478}, {-0.86266667, 0.478}, {-0.84533333, 0.478}, {-0.828, 0.478}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.828, 0.478}, {-0.828, 0.49533333}, {-0.828, 0.51266667}, {-0.828, 0.53}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
    },
    {
      <|"Degree" -> 3, "ControlPoints" -> {{-0.72, -0.05}, {-0.73587977, -0.057837137}, {-0.75175955, -0.065674273}, {-0.76763932, -0.07351141}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.76763932, -0.07351141}, {-0.77018576, -0.09103578}, {-0.7727322, -0.10856015}, {-0.77527864, -0.12608452}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.77527864, -0.12608452}, {-0.78763932, -0.11340377}, {-0.8, -0.10072301}, {-0.81236068, -0.088042261}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.81236068, -0.088042261}, {-0.82981424, -0.09103578}, {-0.8472678, -0.0940293}, {-0.86472136, -0.09702282}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.86472136, -0.09702282}, {-0.85648091, -0.081348547}, {-0.84824045, -0.065674273}, {-0.84, -0.05}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.84, -0.05}, {-0.84824045, -0.034325727}, {-0.85648091, -0.018651453}, {-0.86472136, -0.0029771798}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.86472136, -0.0029771798}, {-0.8472678, -0.0059706997}, {-0.82981424, -0.0089642195}, {-0.81236068, -0.011957739}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.81236068, -0.011957739}, {-0.8, 0.0007230142}, {-0.78763932, 0.013403768}, {-0.77527864, 0.026084521}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.77527864, 0.026084521}, {-0.7727322, 0.0085601509}, {-0.77018576, -0.0089642195}, {-0.76763932, -0.02648859}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.76763932, -0.02648859}, {-0.75175955, -0.034325727}, {-0.73587977, -0.042162863}, {-0.72, -0.05}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
    },
    {
      <|"Degree" -> 3, "ControlPoints" -> {{0.8, 0.53}, {0.82309401, 0.49}, {0.84618802, 0.45}, {0.86928203, 0.41}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.86928203, 0.41}, {0.82309401, 0.41}, {0.77690599, 0.41}, {0.73071797, 0.41}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.73071797, 0.41}, {0.75381198, 0.45}, {0.77690599, 0.49}, {0.8, 0.53}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
    },
    {
      <|"Degree" -> 3, "ControlPoints" -> {{0.8, 0.03}, {0.82, 0.0033333333}, {0.84, -0.023333333}, {0.86, -0.05}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.86, -0.05}, {0.84, -0.076666667}, {0.82, -0.10333333}, {0.8, -0.13}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.8, -0.13}, {0.78, -0.10333333}, {0.76, -0.076666667}, {0.74, -0.05}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.74, -0.05}, {0.76, -0.023333333}, {0.78, 0.0033333333}, {0.8, 0.03}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
    }
  }
|>;

(* 2. Visualization Function *)
bernsteinPoly[p_, j_, t_] := Binomial[p, j] * t^j * (1 - t)^(p - j);
evalBSpline[curve_, t_] := Module[{p = curve["Degree"], cps = curve["ControlPoints"]},
  Sum[cps[[j + 1]] * bernsteinPoly[p, j, t], {j, 0, p}]
];

plotGismoDomain[domain_] := Module[{
  outerCurves, holeCurves, cornerPts, cnetLines, curvePlots, plotRange, bbox
},
  outerCurves = domain["OuterLoop"];
  holeCurves = Flatten[domain["Holes"]];
  cornerPts = domain["Corners"];

  (* Control polygons *)
  cnetLines = Join[
    Line /@ (outerCurves[[All, "ControlPoints"]]),
    Line /@ (holeCurves[[All, "ControlPoints"]])
  ];

  (* Parametric boundary curves *)
  curvePlots = Join[
    Table[ParametricPlot[evalBSpline[c, t], {t, 0, 1}, PlotStyle -> Directive[Thick, Blue]], {c, outerCurves}],
    Table[ParametricPlot[evalBSpline[c, t], {t, 0, 1}, PlotStyle -> Directive[Thick, Red]], {c, holeCurves}]
  ];

  Show[
    Graphics[{
      {Directive[LightGray, Opacity[0.3]], EdgeForm[None]},
      {Directive[Dashed, Gray, Opacity[0.6]], cnetLines},
      {Directive[PointSize[0.015], Darker[Green]], Point /@ (outerCurves[[All, "ControlPoints"]][[All, 1]])},
      {Directive[PointSize[0.025], Red], Point[cornerPts]},
      MapIndexed[Text[Style[Row[{"C", #2[[1]] - 1}], 12, Bold, Darker[Red]], #1 + {0.05, 0.05}] &, cornerPts]
    }],
    curvePlots,
    Axes -> True, GridLines -> Automatic, AspectRatio -> Automatic,
    PlotLabel -> Style[StringTemplate["Planar Domain: `` Corners, Convex: ``, Holes: ``"][
      domain["NumCorners"], domain["IsConvex"], domain["NumHoles"]], 14, Bold]
  ]
];

(* 3. Render Plot *)
domainVisualization = plotGismoDomain[gismoDomain];
Print["Domain loaded successfully. Displaying plot:"];
domainVisualization
