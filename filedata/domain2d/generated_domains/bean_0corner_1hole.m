(* ========================================================================= *)
(* G+Smo Planar Domain Definition & Mathematica Visualization Script        *)
(* Domain: corners=0, convex=False, star=True, holes=1 *)
(* ========================================================================= *)

ClearAll[gismoDomain, plotGismoDomain];

(* 1. Domain Specification *)
gismoDomain = <|
  "NumCorners" -> 0,
  "IsConvex" -> false,
  "IsStarShaped" -> true,
  "NumHoles" -> 1,
  "Corners" -> {},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{1.1, 0}, {0.96757243, 0.38334296}, {0.4055723, 0.52}, {3.1840817e-17, 0.52}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{3.1840817e-17, 0.52}, {-0.4055723, 0.52}, {-0.96757243, 0.38334296}, {-1.1, 1.3471115e-16}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-1.1, 1.3471115e-16}, {-1.2836916, -0.53173892}, {-0.56257345, -1.28}, {-2.3513219e-16, -1.28}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{-2.3513219e-16, -1.28}, {0.56257345, -1.28}, {1.2836916, -0.53173892}, {1.1, -2.694223e-16}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
  },
  "Holes" -> {
    {
      <|"Degree" -> 3, "ControlPoints" -> {{0.22, 0}, {0.22, -0.12150264}, {0.12150264, -0.22}, {0, -0.22}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0, -0.22}, {-0.12150264, -0.22}, {-0.22, -0.12150264}, {-0.22, 0}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{-0.22, 0}, {-0.22, 0.12150264}, {-0.12150264, 0.22}, {0, 0.22}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
      <|"Degree" -> 3, "ControlPoints" -> {{0, 0.22}, {0.12150264, 0.22}, {0.22, 0.12150264}, {0.22, 0}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
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
