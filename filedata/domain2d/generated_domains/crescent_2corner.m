(* ========================================================================= *)
(* G+Smo Planar Domain Definition & Mathematica Visualization Script        *)
(* Domain: corners=2, convex=False, star=True, holes=0 *)
(* ========================================================================= *)

ClearAll[gismoDomain, plotGismoDomain];

(* 1. Domain Specification *)
gismoDomain = <|
  "NumCorners" -> 2,
  "IsConvex" -> false,
  "IsStarShaped" -> true,
  "NumHoles" -> 0,
  "Corners" -> {{-1.2, 0}, {1.2, 0}},
  "OuterLoop" -> {
    <|"Degree" -> 3, "ControlPoints" -> {{-1.2, 0}, {-0.48, 0.36}, {0.48, 0.36}, {1.2, 0}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>,
    <|"Degree" -> 3, "ControlPoints" -> {{1.2, 0}, {0.6, 1.04}, {-0.6, 1.04}, {-1.2, 0}}, "Knots" -> {0, 0, 0, 0, 1, 1, 1, 1}|>
  },
  "Holes" -> {
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
