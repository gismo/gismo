(* ========================================================================= *)
(* G+Smo Planar Domain Definition & Mathematica Visualization Script        *)
(* Domain: corners=5, convex=False, star=True, holes=1 *)
(* ========================================================================= *)

ClearAll[gismoDomain, plotGismoDomain];

(* 1. Domain Specification *)
gismoDomain = <|
  "NumCorners" -> 5,
  "IsConvex" -> false,
  "IsStarShaped" -> true,
  "NumHoles" -> 1,
  "Corners" -> {{1.14, 0}, {0.12669697, 0.38993317}, {-0.92227937, 0.67007519}, {-0.33169697, -0.24099195}, {0.35227937, -1.0842044}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{1.14, 0}, {0.81782965, 0.17050985}, {0.48006197, 0.30048757}, {0.12669697, 0.38993317}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.12669697, 0.38993317}, {-0.23136607, 0.45184455}, {-0.58102485, 0.54522523}, {-0.92227937, 0.67007519}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.92227937, 0.67007519}, {-0.76186126, 0.34276284}, {-0.56500046, 0.039073797}, {-0.33169697, -0.24099195}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-0.33169697, -0.24099195}, {-0.07840848, -0.50154349}, {0.14958363, -0.78261431}, {0.35227937, -1.0842044}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{0.35227937, -1.0842044}, {0.65822109, -0.75431178}, {0.92079464, -0.3929103}, {1.14, 0}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
  },
  "Holes" -> {
    {
      <|"Degree" -> 3, "ControlPoints" -> {{0.22, 0}, {0.17633062, -0.021552126}, {0.13266125, -0.043104252}, {0.088991869, -0.064656378}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.088991869, -0.064656378}, {0.081989159, -0.1128484}, {0.074986449, -0.16104041}, {0.067983739, -0.20923243}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.067983739, -0.20923243}, {0.033991869, -0.17436036}, {1.3877788e-17, -0.13948829}, {-0.033991869, -0.10461622}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.033991869, -0.10461622}, {-0.081989159, -0.1128484}, {-0.12998645, -0.12108058}, {-0.17798374, -0.12931276}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.17798374, -0.12931276}, {-0.15532249, -0.086208504}, {-0.13266125, -0.043104252}, {-0.11, -1.3471115e-17}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.11, -1.3471115e-17}, {-0.13266125, 0.043104252}, {-0.15532249, 0.086208504}, {-0.17798374, 0.12931276}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.17798374, 0.12931276}, {-0.12998645, 0.12108058}, {-0.081989159, 0.1128484}, {-0.033991869, 0.10461622}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.033991869, 0.10461622}, {-2.0816682e-17, 0.13948829}, {0.033991869, 0.17436036}, {0.067983739, 0.20923243}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.067983739, 0.20923243}, {0.074986449, 0.16104041}, {0.081989159, 0.1128484}, {0.088991869, 0.064656378}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0.088991869, 0.064656378}, {0.13266125, 0.043104252}, {0.17633062, 0.021552126}, {0.22, 0}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
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
