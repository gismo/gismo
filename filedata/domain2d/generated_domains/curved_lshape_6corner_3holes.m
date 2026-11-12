(* ========================================================================= *)
(* G+Smo Planar Domain Definition & Mathematica Visualization Script        *)
(* Domain: corners=6, convex=False, star=False, holes=3 *)
(* ========================================================================= *)

ClearAll[gismoDomain, plotGismoDomain];

(* 1. Domain Specification *)
gismoDomain = <|
  "NumCorners" -> 6,
  "IsConvex" -> false,
  "IsStarShaped" -> false,
  "NumHoles" -> 3,
  "Corners" -> {{-1, -1}, {1, -1}, {1, 0}, {0, 0}, {0, 1}, {-1, 1}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{-1, -1}, {-0.33333333, -1.06}, {0.33333333, -1.06}, {1, -1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, -1}, {1.03, -0.66666667}, {1.03, -0.33333333}, {1, 0}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, 0}, {0.66666667, 0.03}, {0.33333333, 0.03}, {0, 0}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0, 0}, {0.03, 0.33333333}, {0.03, 0.66666667}, {0, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0, 1}, {-0.33333333, 1.03}, {-0.66666667, 1.03}, {-1, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1, 1}, {-1.06, 0.33333333}, {-1.06, -0.33333333}, {-1, -1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
  },
  "Holes" -> {
    {
      <|"Degree" -> 3, "ControlPoints" -> {{0.508, -0.38}, {0.536, -0.38}, {0.564, -0.38}, {0.592, -0.38}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.592, -0.38}, {0.592, -0.406}, {0.592, -0.432}, {0.592, -0.458}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.592, -0.458}, {0.618, -0.458}, {0.644, -0.458}, {0.67, -0.458}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.67, -0.458}, {0.67, -0.486}, {0.67, -0.514}, {0.67, -0.542}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.67, -0.542}, {0.644, -0.542}, {0.618, -0.542}, {0.592, -0.542}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.592, -0.542}, {0.592, -0.568}, {0.592, -0.594}, {0.592, -0.62}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.592, -0.62}, {0.564, -0.62}, {0.536, -0.62}, {0.508, -0.62}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.508, -0.62}, {0.508, -0.594}, {0.508, -0.568}, {0.508, -0.542}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.508, -0.542}, {0.482, -0.542}, {0.456, -0.542}, {0.43, -0.542}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.43, -0.542}, {0.43, -0.514}, {0.43, -0.486}, {0.43, -0.458}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.43, -0.458}, {0.456, -0.458}, {0.482, -0.458}, {0.508, -0.458}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.508, -0.458}, {0.508, -0.432}, {0.508, -0.406}, {0.508, -0.38}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
    },
    {
      <|"Degree" -> 3, "ControlPoints" -> {{-0.38, 0.55}, {-0.40381966, 0.53824429}, {-0.42763932, 0.52648859}, {-0.45145898, 0.51473288}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.45145898, 0.51473288}, {-0.45527864, 0.48844633}, {-0.4590983, 0.46215977}, {-0.46291796, 0.43587322}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.46291796, 0.43587322}, {-0.48145898, 0.45489435}, {-0.5, 0.47391548}, {-0.51854102, 0.49293661}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.51854102, 0.49293661}, {-0.54472136, 0.48844633}, {-0.5709017, 0.48395605}, {-0.59708204, 0.47946577}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.59708204, 0.47946577}, {-0.58472136, 0.50297718}, {-0.57236068, 0.52648859}, {-0.56, 0.55}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.56, 0.55}, {-0.57236068, 0.57351141}, {-0.58472136, 0.59702282}, {-0.59708204, 0.62053423}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.59708204, 0.62053423}, {-0.5709017, 0.61604395}, {-0.54472136, 0.61155367}, {-0.51854102, 0.60706339}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.51854102, 0.60706339}, {-0.5, 0.62608452}, {-0.48145898, 0.64510565}, {-0.46291796, 0.66412678}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.46291796, 0.66412678}, {-0.4590983, 0.63784023}, {-0.45527864, 0.61155367}, {-0.45145898, 0.58526712}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.45145898, 0.58526712}, {-0.42763932, 0.57351141}, {-0.40381966, 0.56175571}, {-0.38, 0.55}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
    },
    {
      <|"Degree" -> 3, "ControlPoints" -> {{-0.5, -0.38}, {-0.46535898, -0.44}, {-0.43071797, -0.5}, {-0.39607695, -0.56}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.39607695, -0.56}, {-0.46535898, -0.56}, {-0.53464102, -0.56}, {-0.60392305, -0.56}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.60392305, -0.56}, {-0.56928203, -0.5}, {-0.53464102, -0.44}, {-0.5, -0.38}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
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
