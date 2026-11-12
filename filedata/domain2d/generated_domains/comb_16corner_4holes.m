(* ========================================================================= *)
(* G+Smo Planar Domain Definition & Mathematica Visualization Script        *)
(* Domain: corners=16, convex=False, star=False, holes=4 *)
(* ========================================================================= *)

ClearAll[gismoDomain, plotGismoDomain];

(* 1. Domain Specification *)
gismoDomain = <|
  "NumCorners" -> 16,
  "IsConvex" -> false,
  "IsStarShaped" -> false,
  "NumHoles" -> 4,
  "Corners" -> {{-1, -1}, {1, -1}, {1, 1}, {0.6, 1}, {0.6, -0.4}, {0.45, -0.4}, {0.45, 1}, {0.05, 1}, {0.05, -0.4}, {-0.05, -0.4}, {-0.05, 1}, {-0.45, 1}, {-0.45, -0.4}, {-0.6, -0.4}, {-0.6, 1}, {-1, 1}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{-1, -1}, {-0.33333333, -1}, {0.33333333, -1}, {1, -1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, -1}, {1, -0.33333333}, {1, 0.33333333}, {1, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1, 1}, {0.86666667, 1}, {0.73333333, 1}, {0.6, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.6, 1}, {0.6, 0.53333333}, {0.6, 0.066666667}, {0.6, -0.4}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.6, -0.4}, {0.55, -0.4}, {0.5, -0.4}, {0.45, -0.4}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.45, -0.4}, {0.45, 0.066666667}, {0.45, 0.53333333}, {0.45, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.45, 1}, {0.31666667, 1}, {0.18333333, 1}, {0.05, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.05, 1}, {0.05, 0.53333333}, {0.05, 0.066666667}, {0.05, -0.4}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.05, -0.4}, {0.016666667, -0.4}, {-0.016666667, -0.4}, {-0.05, -0.4}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.05, -0.4}, {-0.05, 0.066666667}, {-0.05, 0.53333333}, {-0.05, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.05, 1}, {-0.18333333, 1}, {-0.31666667, 1}, {-0.45, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.45, 1}, {-0.45, 0.53333333}, {-0.45, 0.066666667}, {-0.45, -0.4}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.45, -0.4}, {-0.5, -0.4}, {-0.55, -0.4}, {-0.6, -0.4}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.6, -0.4}, {-0.6, 0.066666667}, {-0.6, 0.53333333}, {-0.6, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.6, 1}, {-0.73333333, 1}, {-0.86666667, 1}, {-1, 1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1, 1}, {-1, 0.33333333}, {-1, -0.33333333}, {-1, -1}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
  },
  "Holes" -> {
    {
      <|"Degree" -> 3, "ControlPoints" -> {{-0.828, 0.38}, {-0.80933333, 0.38}, {-0.79066667, 0.38}, {-0.772, 0.38}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.772, 0.38}, {-0.772, 0.36266667}, {-0.772, 0.34533333}, {-0.772, 0.328}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.772, 0.328}, {-0.75466667, 0.328}, {-0.73733333, 0.328}, {-0.72, 0.328}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.72, 0.328}, {-0.72, 0.30933333}, {-0.72, 0.29066667}, {-0.72, 0.272}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.72, 0.272}, {-0.73733333, 0.272}, {-0.75466667, 0.272}, {-0.772, 0.272}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.772, 0.272}, {-0.772, 0.25466667}, {-0.772, 0.23733333}, {-0.772, 0.22}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.772, 0.22}, {-0.79066667, 0.22}, {-0.80933333, 0.22}, {-0.828, 0.22}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.828, 0.22}, {-0.828, 0.23733333}, {-0.828, 0.25466667}, {-0.828, 0.272}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.828, 0.272}, {-0.84533333, 0.272}, {-0.86266667, 0.272}, {-0.88, 0.272}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.88, 0.272}, {-0.88, 0.29066667}, {-0.88, 0.30933333}, {-0.88, 0.328}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.88, 0.328}, {-0.86266667, 0.328}, {-0.84533333, 0.328}, {-0.828, 0.328}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.828, 0.328}, {-0.828, 0.34533333}, {-0.828, 0.36266667}, {-0.828, 0.38}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
    },
    {
      <|"Degree" -> 3, "ControlPoints" -> {{-0.17, 0.3}, {-0.18587977, 0.29216286}, {-0.20175955, 0.28432573}, {-0.21763932, 0.27648859}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.21763932, 0.27648859}, {-0.22018576, 0.25896422}, {-0.2227322, 0.24143985}, {-0.22527864, 0.22391548}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.22527864, 0.22391548}, {-0.23763932, 0.23659623}, {-0.25, 0.24927699}, {-0.26236068, 0.26195774}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.26236068, 0.26195774}, {-0.27981424, 0.25896422}, {-0.2972678, 0.2559707}, {-0.31472136, 0.25297718}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.31472136, 0.25297718}, {-0.30648091, 0.26865145}, {-0.29824045, 0.28432573}, {-0.29, 0.3}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.29, 0.3}, {-0.29824045, 0.31567427}, {-0.30648091, 0.33134855}, {-0.31472136, 0.34702282}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.31472136, 0.34702282}, {-0.2972678, 0.3440293}, {-0.27981424, 0.34103578}, {-0.26236068, 0.33804226}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.26236068, 0.33804226}, {-0.25, 0.35072301}, {-0.23763932, 0.36340377}, {-0.22527864, 0.37608452}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.22527864, 0.37608452}, {-0.2227322, 0.35856015}, {-0.22018576, 0.34103578}, {-0.21763932, 0.32351141}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.21763932, 0.32351141}, {-0.20175955, 0.31567427}, {-0.18587977, 0.30783714}, {-0.17, 0.3}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
    },
    {
      <|"Degree" -> 3, "ControlPoints" -> {{0.25, 0.38}, {0.27309401, 0.34}, {0.29618802, 0.3}, {0.31928203, 0.26}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.31928203, 0.26}, {0.27309401, 0.26}, {0.22690599, 0.26}, {0.18071797, 0.26}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.18071797, 0.26}, {0.20381198, 0.3}, {0.22690599, 0.34}, {0.25, 0.38}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
    },
    {
      <|"Degree" -> 3, "ControlPoints" -> {{0.8, 0.38}, {0.82, 0.35333333}, {0.84, 0.32666667}, {0.86, 0.3}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.86, 0.3}, {0.84, 0.27333333}, {0.82, 0.24666667}, {0.8, 0.22}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.8, 0.22}, {0.78, 0.24666667}, {0.76, 0.27333333}, {0.74, 0.3}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.74, 0.3}, {0.76, 0.32666667}, {0.78, 0.35333333}, {0.8, 0.38}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
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
