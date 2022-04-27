(* ::Package:: *)

(* ::Title:: *)
(*Visualization of Complex Functions on the Riemann Sphere*)


(* ::Section:: *)
(*Info*)


(*:Title: Visualization of Complex Functions on the Riemann Sphere *)

(*:Context: complexVisualize` *) (*new*)

(*:Author:
  Antonio Hernandez-Garduno (UAM-I),
  Angeles Sandoval-Romero (UNAM),
  December 2013
*)

(*:Summary: Produces a domain colored Riemann sphere which gives 
information about the complex function defined on it.  For Mobius
transformations, a 'fast rendering' version allows to dynamically 
manipulate the parameters. *)

(*:Mathematica Version: 8.0 *)

(*:Package Version: 2.0 *)

(*:Name: Complex Visualize *)

(*:Keywords: domain coloring, complex functions, Riemann sphere, 
Mobius transformation *)


(* ::Section::Closed:: *)
(*Startup*)


BeginPackage["complexVisualize`"];


(* ::Section::Closed:: *)
(*Usage*)


(* ::Subsection::Closed:: *)
(*commands defined*)


complexVisualize::usage = "complexVisualize[f[z], z] provides a visualization of a complex function f[z] on the Riemann sphere (compactified complex plane) through domain coloring.  It admits the options:
	'colorScheme'  (default is \"azimuth\")
	'targetMesh'  (default is {15,15})
	'targetMeshColors'  (default is {Gray, Black})
	'referenceMeshColors'  (default is {Red, Blue, Green})
	'referenceMeshThickness'  (default is Thickness[0.0045])
    (Further details available within each option's help.)";


mobiusVisualize::usage = "mobiusVisualize[f[z], z] provides a visualization of the M\[ODoubleDot]bius transformation f[z] on the Riemann sphere.  Its option fastRender (True/False) controls whether only the preimages of the reference curves (the real and imaginary axes, and the unit circle) are drawn.";


azimuthLatitude::usage = "azimuthLatitude[a,r] is a setting for option colorScheme. Parameters must satisfy 0<a<1, 0<r<1. With 'a' near 1, 'r' controlls the concentration of black/white around zeroes/poles. Setting azimuthLatitude[] or \"azimuthLatitude\" is equivalent to azimuthLatitude[0.9,0.3].";


cylC::usage =
"cylC[{\[Theta],z}] gives the complex number represented by the point on the Riemann sphere with cylindrical coordinates {\[Theta],z}."; 


drawStereographic::usage = "drawStereographic[z] draws the stereographic projection associated to a complex number z.";


(* ::Subsection::Closed:: *)
(*options defined*)


(* ::Subsubsection::Closed:: *)
(*options for 'complexVisualize'*)


colorScheme::usage = "colorScheme is an option for complexVisualize that determines how the target Riemann sphere is identified with the hue-saturation-brightness (HSB) space.

      With the default setting \"azimuth\", only the hue is accounted for and is identified with the azimuthal angle \[Theta].

      With colorScheme -> \"azimuthLatitude\" the colatitude \[CurlyPhi] is taken into account in the saturation and brightness, so that zero (infinity) on the target sphere is totally dark (bright).  Further control is achieved with colorScheme -> azimuthLatitude[a,r].

      A general color scheme is specified as colorScheme->cs, with cs[\[Theta], \[CurlyPhi]] a Hue function of the azimuthal and colatitudinal angles in the target Riemann sphere.

      colorScheme -> \"azimuth\" is equivalent to
	colorScheme -> Function[{\[Theta], \[CurlyPhi]}, Hue[\[Theta]/(2\[Pi])]]
      colorScheme -> \"azimuthLatitude\" is equivalent to
	colorScheme -> azimuthLatitude[0.9,0.3]
      colorScheme -> \"metallic\" is equivalent to
	colorScheme -> Function[{\[Theta], \[CurlyPhi]}, Hue[\[Theta]/(2\[Pi]), (\[CurlyPhi]/\[Pi])^(1/2), ((\[Pi]-\[CurlyPhi])/\[Pi])^(1/2)]]";


targetMesh::usage = "targetMesh is an option for complexVisualize that specifies the number of mesh lines on the target Riemann sphere.  targetMesh -> {n1,n2} corresponds to n1 azimuthal and n2 latitudinal mesh lines (default is {15,15}).  Setting 'None' is equivalent to {0,0}.";


targetMeshColors::usage = "targetMeshColors is an option for complexVisualize that specifies the coloring of the pullback of the mesh on the target Riemann sphere.  Default is targetMeshColors -> {Gray, Black}.";


referenceMeshColors::usage = "referenceMeshColors is an option for complexVisualize that specifies the coloring of the pullback of the reference curves (i.e. the geodesics representing the unit circle, real axis, and imaginary axis) on the target Riemann sphere.  Default is referenceMeshColors -> {Red, Blue, Green}.";


referenceMeshThickness::usage = "referenceMeshThickness is an option for complexVisualize that specifies the thickness of the pullback of the reference curves (i.e. the geodesics representing the unit circle, real axis, and imaginary axis) on the target Riemann sphere.  Default is referenceMeshThickness -> Thickness[0.0045].";


(* ::Subsubsection::Closed:: *)
(*options for 'mobiusVisualize'*)


fastRender::usage = "fastRender is an option for mobiusVisualize.  With the setting True, only the circles corresponding to the preimages of the reference curves (the real and imaginary axes, and the unit circle) are drawn. With the setting False, complexVisualize is invoked. Default setting is True.";


lightingFR::usage = 
"lightingFR is an option for mobiusVisualize that specifies the simulated lighting to use incoloring the sphere when fastRender->True.";


referenceMeshThicknessFR::usage = "referenceMeshThicknessFR is an option for mobiusVisualize that specifies the thickness of the reference circles (i.e. the geodesics representing the unit circle, and the real and imaginary axes) when fastRender->True.";


(* ::Subsubsection::Closed:: *)
(*options for 'drawStereographic'*)


showLabels::usage = "showLabels (True/False) is an option for drawStereographic that specifies whether to draw labels for graphic elements.";


(* ::Section::Closed:: *)
(*Private*)


Begin["`Private`"];


(* ::Subsection::Closed:: *)
(*Options*)


(* ::Subsubsection::Closed:: *)
(*prolog*)


(* standardPlotOptions1 and standardPlotOptions2 are lists of 
standard ParametricPlot3D option names *)


standardPlotOptions1={AlignmentPoint,AspectRatio,AxesEdge,AxesOrigin,AxesStyle,Background,BoxStyle,ControllerLinking,ControllerMethod,ControllerPath,DisplayFunction,Epilog,ImageMargins,ImagePadding,Lighting,MaxRecursion,Method,PlotRegion,PreserveImageOptions,Prolog,RegionFunction,RotationAction,SphericalRegion,Ticks,TicksStyle,TouchscreenAutoZoom,ViewAngle,ViewCenter,ViewMatrix,ViewPoint,ViewRange,ViewVector,ViewVertical,WorkingPrecision};


standardPlotOptions2={AlignmentPoint,AspectRatio,AxesEdge,AxesOrigin,AxesStyle,Background,BoxStyle,ControllerLinking,ControllerMethod,ControllerPath,DisplayFunction,Epilog,ImageMargins,ImagePadding,Lighting,MaxRecursion,Method,PlotRegion,PreserveImageOptions,Prolog,RegionFunction,RotationAction,SphericalRegion,Ticks,TicksStyle,TouchscreenAutoZoom,WorkingPrecision};


standardPlotOptions3={AlignmentPoint,AspectRatio,AxesEdge,AxesOrigin,AxesStyle,Background,BoxStyle,DisplayFunction,Epilog,ImageMargins,ImagePadding,ImageSize,Lighting,PlotRegion,PreserveImageOptions,Prolog,RotationAction,SphericalRegion,Ticks,TicksStyle,TouchscreenAutoZoom};


(* ::Subsubsection::Closed:: *)
(*default options*)


Options[complexVisualize] = Join[
{
colorScheme -> "azimuth",
(*targetMesh \[Rule] {15,15},*)
targetMesh->{{-((7 \[Pi])/8),-((3 \[Pi])/4),-((5 \[Pi])/8),-((3 \[Pi])/8),-(\[Pi]/4),-(\[Pi]/8),\[Pi]/8,\[Pi]/4,(3 \[Pi])/8,(5 \[Pi])/8,(3 \[Pi])/4,(7 \[Pi])/8},{\[Pi]/16,\[Pi]/8,(3 \[Pi])/16,\[Pi]/4,(5 \[Pi])/16,(3 \[Pi])/8,(7 \[Pi])/16,(9 \[Pi])/16,(5 \[Pi])/8,(11 \[Pi])/16,(3 \[Pi])/4,(13 \[Pi])/16,(7 \[Pi])/8,(15 \[Pi])/16}},
targetMeshColors -> {Gray,Black},
referenceMeshColors -> {Red,Blue,Green},
referenceMeshThickness -> Thickness[0.0045],
PlotPoints -> 70, ImageSize->500, Lighting->"Neutral", ViewPoint->{2.05,-1.8,2}, RotationAction->"Clip", Axes->False, 
Boxed -> False, AxesLabel -> {"x","y","z"}
},
FilterRules[Options[ParametricPlot3D],standardPlotOptions1]
];


Options[mobiusVisualize] = Join[
{
fastRender->True,
lightingFR->"Neutral",
colorScheme -> "azimuth",
(*targetMesh \[Rule]{15,15},*)
targetMesh->{{-((7 \[Pi])/8),-((3 \[Pi])/4),-((5 \[Pi])/8),-((3 \[Pi])/8),-(\[Pi]/4),-(\[Pi]/8),\[Pi]/8,\[Pi]/4,(3 \[Pi])/8,(5 \[Pi])/8,(3 \[Pi])/4,(7 \[Pi])/8},{\[Pi]/16,\[Pi]/8,(3 \[Pi])/16,\[Pi]/4,(5 \[Pi])/16,(3 \[Pi])/8,(7 \[Pi])/16,(9 \[Pi])/16,(5 \[Pi])/8,(11 \[Pi])/16,(3 \[Pi])/4,(13 \[Pi])/16,(7 \[Pi])/8,(15 \[Pi])/16}},
targetMeshColors -> {Gray,Black},
referenceMeshColors -> {Red,Blue,Green},
referenceMeshThickness -> Thickness[0.0045],
referenceMeshThicknessFR->AbsoluteThickness[2],
PlotPoints -> 70, ImageSize->340,Lighting->"Neutral", RotationAction->"Clip",
Axes->False,Boxed->False,AxesLabel -> {"x","y","z"}
},
FilterRules[Options[ParametricPlot3D],standardPlotOptions2]
];


Options[drawStereographic] = Join[
{
showLabels -> True,
Boxed -> False,
ViewVector -> {9,2,3},
ImageSize -> 410
},
FilterRules[Options[ParametricPlot3D],standardPlotOptions3]
];


(* ::Subsection::Closed:: *)
(*Coordinate transformations*)


(* ::Subsubsection:: *)
(*Stereographic projection*)


(* stereo[{x,y,z}] gives the complex number \[Zeta] that is the
stereographic projection of the point on the sphere {x,y,z} *)


stereo[{x_,y_,z_}]= {x ,y}/(1 - z);


(* stereoInv[\[Zeta]] gives the point on the sphere that
stereographically projects to the complex number \[Zeta] *)


stereoInv[z_] = {2u, 2v, \[Rho]sq-1}/(\[Rho]sq+1)/.{\[Rho]sq->u^2+v^2}/.{u->Re[z], v->Im[z]};


stereoInv[ComplexInfinity]={0,0,1};


(* ::Subsubsection:: *)
(*Spherical and projective coordinates*)


(*
Let {u,v} be the stereographic projection of a point {\[Theta],\[CurlyPhi]} 
on the unit sphere.  Then
	sphericalToProjective[{\[Theta],\[CurlyPhi]}] = {u,v}
	projectiveToSpherical[{u,v}] = {\[Theta],\[CurlyPhi]}
*)


sphericalToProjective[{\[Theta]_, \[CurlyPhi]_}] = {Cos[\[Theta]] Cot[\[CurlyPhi]/2], Cot[\[CurlyPhi]/2] Sin[\[Theta]]};


projectiveToSpherical[{u_, v_}] = {ArcTan[u,v], 2ArcTan[1/Sqrt[u^2+v^2]]};


(* ::Subsubsection:: *)
(*Cylindrical to complex coordinates*)


(* cylC transforms from cilindrical to complex *)


cylinderToSphere[{\[Theta]_,z_}]:={Sqrt[1-z^2]Cos[\[Theta]],Sqrt[1-z^2]Sin[\[Theta]],z}


cylC[{\[Theta]_,z_}]:=#[[1]]+I #[[2]]&[stereo[cylinderToSphere[{\[Theta],z}]]]


(* ::Subsection::Closed:: *)
(*Color schemes*)


(* coloring[f, cs] gives the composition of color scheme cs[\[Theta],\[CurlyPhi]] 
with function f in spherical coordinates.  Common color schemes are
'azimuthColorScheme', 'azimuthLatitudeColorScheme' and 'metallic' *)


coloring[f_, cs_] := Composition[cs@@#&, projectiveToSpherical, f, sphericalToProjective]


(* 'azimuth' color scheme *)


azimuthColorScheme[\[Theta]_, \[CurlyPhi]_] := Hue[\[Theta]/(2\[Pi])]


coloring[f_, "azimuth"] := coloring[f, azimuthColorScheme]


(* 'azimuthLatitude' color scheme *)


azimuthLatitude[a_,r_][\[Theta]_,\[CurlyPhi]_] = Hue[\[Theta]/(2\[Pi]),Exp[Log[a](r \[Pi]/\[CurlyPhi])^2],Exp[Log[a](r \[Pi]/(\[Pi]-\[CurlyPhi]))^2]];


azimuthLatitude[] = azimuthLatitude[0.9,0.3];


coloring[f_, "azimuthLatitude"] := coloring[f, azimuthLatitude[]]


(* 'metallic' color scheme *)


metallicColorScheme[\[Theta]_, \[CurlyPhi]_] := Hue[\[Theta]/(2\[Pi]), (\[CurlyPhi]/\[Pi])^(1/2), ((\[Pi]-\[CurlyPhi])/\[Pi])^(1/2)]


coloring[f_, "metallic"] := coloring[f, metallicColorScheme]


(* ::Subsection::Closed:: *)
(*The command 'complexVisualize'*)


(* ::Subsubsection:: *)
(*vector form of a complex function*)


(* complexToVectorFunction[f[z],z] gives the vector form of the 
complex function f[z] *)


complexToVectorFunction[f_, z_][{x_,y_}] := ComplexExpand[{Re[#],Im[#]}&[f]/.{z->x+I y}]


(* ::Subsubsection:: *)
(*polar and Cartesian projections*)


(* Cartesian and polar coordinate- projections *)


polarR[{x_,y_}] := Sqrt[x^2+y^2];
cartesianX[{x_,y_}] := x;
cartesianY[{x_,y_}] := y;


(* ::Subsubsection:: *)
(*sphere meshing*)


(* sphereFun[f] gives the spherical-coordinates representation of 
vector function f *)


sphereFun[f_] := Composition[projectiveToSpherical, f, sphericalToProjective]


(* ::Subsubsection:: *)
(*reference curves*)


(* By 'reference curves' we mean the Real Axis, the Imaginary Axis, 
and the Unit Circle *)


(* {MeshFunction \[Rule] referenceR[f], Mesh \[Rule] {{1}}} draws the 
preimage of the UNIT CIRCLE by function f *)


(* {MeshFunction \[Rule] referenceX[f], Mesh \[Rule] {{0}}} draws the 
preimage of the REAL AXIS by function f *)


(* {MeshFunction \[Rule] referenceY[f], Mesh \[Rule] {{0}}} draws the 
preimage of the IMAGINARY AXIS by function f *)


referenceR[f_] := Composition[polarR, f, sphericalToProjective]


referenceX[f_] := Composition[cartesianX, f, sphericalToProjective]


referenceY[f_] := Composition[cartesianY, f, sphericalToProjective];


(* ::Subsubsection:: *)
(*Cartesian axes*)


(* theCartesianAxes draws reference Cartesian axes. *)


theCartesianAxes = MapThread[Graphics3D[{#1, AbsoluteThickness[2], Line[{{0,0,0},#2}]}]&, {{Blue,Green,Black}, {{1.5,0,0},{0,1.5,0},{0,0,1.5}}}];


(* ::Subsubsection:: *)
(*complexVisualize*)


complexVisualize[f_, \[Zeta]_, opts:OptionsPattern[]] := Block[
 { fn, vecFun, sphereMeshing, referenceMeshing, referenceMeshStyle, theSphere, \[Epsilon]=0.01},
 fn = Which[Re[f]===0, f+0.0001, Im[f]===0, f+0.0001I, True, f]; (* trick to deal with purely real or imaginary functions *)
 vecFun[{x_, y_}] = complexToVectorFunction[fn, \[Zeta]][{x,y}]; (* vecFun is the vector form of f[\[Zeta]] *)
 (* 'sphereMeshing' is the pullback of the standard spherical-coordinates-meshing determined by function f[\[Zeta]] *)
 sphereMeshing = {Composition[cartesianX,sphereFun[vecFun]][{#4,#5}]&, Composition[cartesianY,sphereFun[vecFun]][{#4,#5}]&};
 (* 'referenceMeshing' is the pullback of the reference curves determined by function f[\[Zeta]] *)
 referenceMeshing = {referenceR[vecFun][{#4,#5}]&, referenceY[vecFun][{#4,#5}]&, referenceX[vecFun][{#4,#5}]&};
 referenceMeshStyle = Map[{OptionValue[referenceMeshThickness],#}&, OptionValue[referenceMeshColors]];
 (* 'theSphere' is the domain-colored and meshed Riemann Sphere *)
 theSphere = ParametricPlot3D[
{Cos[\[Theta]]Sin[\[CurlyPhi]],Sin[\[Theta]]Sin[\[CurlyPhi]],Cos[\[CurlyPhi]]},
{\[Theta],0-\[Epsilon],2\[Pi]-\[Epsilon]},{\[CurlyPhi],\[Epsilon],\[Pi]-\[Epsilon]},
Evaluate[FilterRules[{opts}, Options[ParametricPlot3D]]],
Evaluate[FilterRules[Options[complexVisualize], Options[ParametricPlot3D]]], 
ColorFunction -> Function[{x,y,z,\[Theta],\[CurlyPhi]}, coloring[vecFun,OptionValue[colorScheme]][{\[Theta],\[CurlyPhi]}]],
ColorFunctionScaling -> False,
MeshFunctions -> Join[sphereMeshing, referenceMeshing],
MeshStyle -> Join[OptionValue[targetMeshColors], referenceMeshStyle],
Mesh -> Join[OptionValue[targetMesh]/.{None->{0,0}}, {{1},{0},{0}}]
];
Show[theSphere, theCartesianAxes, PlotRange->All]
]


(* ::Subsection::Closed:: *)
(*The command 'mobiusVisualize'*)


(* ::Subsubsection:: *)
(*preimages of reference curves for a Mobius Transformation*)


(* Assuming f[z]=(a z + b)/(c z + d) *)


(* coefficientsMobius[f[z],z] gives {a,b,c,d} *)


coefficientsMobius[expr_, z_] := Map[Reverse[Table[Coefficient[#,z,i],{i,0,1}]]&,{Numerator[#], Denominator[#]}&[Together[expr]]]//Flatten


(* Notation: {z0,z1,zInf,z\[ImaginaryI],zm1} are the preimages under f[z] of 
{0,1,Infinity,\[ImaginaryI],-1} *)


(* fivePreImages[{a,b,c,d}] gives {z0,z1,zInf,z\[ImaginaryI],zm1} *)


fivePreImages[{a_,b_,c_,d_}]:=Block[{fpi},
Off[Power::infy];
fpi={-(b/a),(d-b)/(a-c),-(d/c),(I d-b)/(a-I c),-((b+d)/(a+c))} ;
On[Power::infy];
fpi
]


referenceTriads[{z0_,z1_,zInf_,zI_,zm1_}] := {{z1,zI,zm1},{z0,z1,zInf},{z0,zI,zInf}}


normalize[vec:{x_,y_,z_}] = vec/Sqrt[vec . vec];


(* ::Subsubsection:: *)
(*circle through three points*)


(* circleThroughThreePoints draws a circle passing 
through three given points on the unit sphere *)


circleThroughThreePoints[{pt1_,pt2_,pt3_},opts___] := Block[
{center, v1, v2, v3, r},
v3 = Cross[pt1-pt2, pt3-pt2];
center = pt2 . v3/v3 . v3 v3;
r = Sqrt[1-center . center];
v1 = r normalize[pt2-center];
v2 = r normalize[v3\[Cross]v1];
ParametricPlot3D[
center+Cos[\[Theta]]v1+Sin[\[Theta]]v2, 
{\[Theta],0,2\[Pi]},
Evaluate[FilterRules[{opts},Options[ParametricPlot3D]]]
]
]


(* referenceCircles[f,\[Zeta]] draws the three reference geodesics *)


referenceCircles[f_,\[Zeta]_,opts:OptionsPattern[mobiusVisualize]] :=
Block[{zpts,g},
zpts = stereoInv/@fivePreImages[coefficientsMobius[f,\[Zeta]]];
g[1] = Graphics3D[Sphere[]];
g[2] = MapThread[
circleThroughThreePoints[#1, PlotStyle->{OptionValue[referenceMeshThicknessFR], #2}]&,
{referenceTriads[zpts], OptionValue[referenceMeshColors]}
];
Show[
g[1], g[2], theCartesianAxes,
PlotRange->All,Lighting->OptionValue[lightingFR],Evaluate[FilterRules[Join[{opts},Options[mobiusVisualize]],Options[Graphics3D]]]
]
]


(* ::Subsubsection:: *)
(*mobiusVisualize*)


(* mobiusVisualize[f,\[Zeta]] switches to referenceCircles[f,\[Zeta]] 
or complexVisualize[f,\[Zeta]] depending on whether option 
fastRender is set to True or False *)


mobiusVisualize[f_,\[Zeta]_,opts:OptionsPattern[]] := If[OptionValue[fastRender],
referenceCircles[f, \[Zeta],opts],
complexVisualize[f, \[Zeta],
Evaluate[FilterRules[
Join[{opts},Options[mobiusVisualize]],
Options[complexVisualize]
]]
], Message[mobiusVisualize::mderr,fR];
]


mobiusVisualize::mderr = "Mode specification `1` is ambiguous.";


(* ::Subsection::Closed:: *)
(*The command 'drawStereographic'*)


drawStereographic[zComplex_, opts:OptionsPattern[]] := Block[
{orig={0,0,0},north={0,0,1},arcRadius=0.3,zPoint,\[Zeta]Point,zDir,plX,plY,radius,arc,ax,ay,b,\[Theta],\[CurlyPhi],axesEP,axes,legs,gr,tt},
(* Graphic elements *)
zPoint = {Re[#],Im[#],0}&[zComplex];
{ax,ay,b} = \[Zeta]Point=stereoInv[zComplex];
zDir = zPoint/Max[Sqrt[zPoint . zPoint],0.01];
\[CurlyPhi] = ArcCos[b]; \[Theta] = If[ay==0,0.001,ArcTan[ax,ay]];
axesEP = 1.2{{1,0,0},{0,1,0},{0,0,1}};
axes = Line[{orig,#}]&/@axesEP;
legs = {Line[{\[Zeta]Point,{0,0,b}}],Line[{\[Zeta]Point,{ax,ay,0}}],Line[{orig,zPoint}]};
radius = Line[{orig,\[Zeta]Point}];
arc[1] = ParametricPlot3D[arcRadius{Cos[t],Sin[t],0},{t,0,\[Theta]}];
arc[2] = ParametricPlot3D[arcRadius(Cos[t]north+Sin[t]zDir),{t,0,\[CurlyPhi]}];
(* u-v plane dimensions *)
plX[1] = Min[-1.7,Re[zComplex]-0.5];plX[2]=Max[1.7,Re[zComplex]+0.5];
plY[1] = Min[-1.7,Im[zComplex]-0.5];plY[2]=Max[1.7,Im[zComplex]+0.5];
(* Labels *)
tt["uvplane"] = Text[Style["u-v plane","Label",9],0.8{plX[1],plY[2],0}];
tt["axes"] = MapThread[Text[Style[#1,Larger,Italic],1.1#2]&,{{"u","v"},Most[axesEP]}];
tt["rho"] = Text[Style["\[Rho]",Larger],0.72zPoint-{0,0,0.1}];
tt["a"] = Text[Style["a",Larger,Italic],{0.45ax,0.45ay,b+0.1}];
tt["b"] = Text[Style["b",Larger,Italic],{0,0,0.5b},{2.5,0}];
tt["\[Theta]"] = Text[Style["\[Theta]",Larger],1.3arcRadius{Cos[0.5\[Theta]],Sin[0.5\[Theta]],-0.15}];
tt["\[CurlyPhi]"] = Text[Style["\[CurlyPhi]",Larger],1.4 arcRadius{Cos[\[Theta]]Sin[0.6\[CurlyPhi]],Sin[\[Theta]]Sin[0.6\[CurlyPhi]],Cos[0.6\[CurlyPhi]]}];
(* Drawing elements *)
gr["axes"] = Graphics3D[axes];
gr["plane"] = Plot3D[0,{x,plX[1],plX[2]},{y,plY[1],plY[2]},Mesh->None,PlotStyle->Opacity[0.7]];
gr["sphere"] = Graphics3D[{Opacity[0.5],Sphere[]}];
gr["points"] = Graphics3D[{AbsolutePointSize[4]}~Join~(Point/@{north,zPoint,\[Zeta]Point})];
gr["projection"] = Graphics3D[{AbsoluteThickness[1],Line[{north,zPoint}]}];
gr["auxLines"] = Graphics3D[{legs,radius}];
gr["arcs"] = Show[{arc[1],arc[2]}];
gr["labels"] = If[OptionValue[showLabels]==True,
Graphics3D[tt/@{"uvplane","axes","rho","a","b","\[Theta]","\[CurlyPhi]"}], 
Graphics3D[]];
Show[
gr/@{
"axes","plane","sphere","points","projection","auxLines","arcs","labels"
},
Evaluate[FilterRules[{opts}, Options[ParametricPlot3D]]],
Evaluate[FilterRules[Options[drawStereographic], Options[ParametricPlot3D]]]
]
]


(* ::Section::Closed:: *)
(*End*)


End[];
EndPackage[];
