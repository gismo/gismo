% Input for Isogemetric Analysis Volume Decomposition 
% Object: the unit cube

% vetices
nvert = 8;
vert = [
  0 0 1
  1 0 1
  1 1 1
  0 1 1
  0 0 0
  1 0 0
  1 1 0
  0 1 0
  ];

% faces
nface = 6;

% -------------------------------------------------------------------------
% face 1
face{1}.pCorner = [1 2 3 4]; % corners of the face
face{1}.pDeg1 = 3; % patch spline degree
face{1}.pDeg2 = 3; % patch spline degree
face{1}.pKnotv1 = [0 0 0 0 1 1 1 1];
face{1}.pKnotv2 = [0 0 0 0 1 1 1 1];
face{1}.pCp = [
    0 0 0
    1/3 0 0
    2/3 0 0
    1 0 0
    0 1/3 0
    1/3 1/3 0
    2/3 1/3 0
    1 1/3 0    
    0 2/3 0
    1/3 2/3 0
    2/3 2/3 0
    1 2/3 0       
    0 1 0
    1/3 1 0
    2/3 1 0
    1 1 0           
];
face{1}.nLoop = 1; % number of loops
% - the trimming curves of the first loop
face{1}.pLoop{1}.deg = 2; 
face{1}.pLoop{1}.knotv = [0 0 0 1 1 1];
face{1}.pLoop{1}.cp = [
    0 0   
    .5 0
    1 0
    ];
face{1}.pLoop{2}.deg = 2; 
face{1}.pLoop{2}.knotv = [0 0 0 1 1 1];
face{1}.pLoop{2}.cp = [
    1 0   
    1 .5
    1 1
    ];
% - - note that the order of the cps for the third trimming curve is
% opposite to that of the patch control points
face{1}.pLoop{3}.deg = 2; 
face{1}.pLoop{3}.knotv = [0 0 0 1 1 1];
face{1}.pLoop{3}.cp = [
    1 1   
    .5 1
    0 1
    ];
face{1}.pLoop{4}.deg = 2; 
face{1}.pLoop{4}.knotv = [0 0 0 1 1 1];
face{1}.pLoop{4}.cp = [
    0 1   
    0 .5
    0 0
    ];

% -------------------------------------------------------------------------
% face 2
face{2}.pCorner = [1 2 6 5]; 
face{2}.pDeg1 = 3; 
face{2}.pDeg2 = 3; 
face{2}.pKnotv1 = [0 0 0 0 1 1 1 1];
face{2}.pKnotv2 = [0 0 0 0 1 1 1 1];
face{2}.pCp = [
    0 0 0
    1/3 0 0
    2/3 0 0
    1 0 0
    0 0 1/3
    1/3 0 1/3
    2/3 0 1/3
    1 0 1/3 
    0 0 2/3
    1/3 0 2/3
    2/3 0 2/3
    1 0 2/3    
    0 0 1
    1/3 0 1
    2/3 0 1
    1 0 1               
];
face{2}.nLoop = 1; 
face{2}.pLoop{1}.deg = 2; 
face{2}.pLoop{1}.knotv = [0 0 0 1 1 1];
face{2}.pLoop{1}.cp = [
    0 0   
    .5 0
    1 0
    ];
face{2}.pLoop{2}.deg = 2; 
face{2}.pLoop{2}.knotv = [0 0 0 1 1 1];
face{2}.pLoop{2}.cp = [
    1 0   
    1 .5
    1 1
    ];
face{2}.pLoop{3}.deg = 2; 
face{2}.pLoop{3}.knotv = [0 0 0 1 1 1];
face{2}.pLoop{3}.cp = [
    1 1   
    .5 1
    0 1
    ];
face{2}.pLoop{4}.deg = 2; 
face{2}.pLoop{4}.knotv = [0 0 0 1 1 1];
face{2}.pLoop{4}.cp = [
    0 1   
    0 .5
    0 0
    ];

% -------------------------------------------------------------------------
% face 3
face{3}.pCorner = [2 3 7 6]; 
face{3}.pDeg1 = 3; 
face{3}.pDeg2 = 3; 
face{3}.pKnotv1 = [0 0 0 0 1 1 1 1];
face{3}.pKnotv2 = [0 0 0 0 1 1 1 1];
face{3}.pCp = [
    1 0 0
    1 1/3 0
    1 2/3 0
    1 1 0
    1 0 1/3
    1 1/3 1/3
    1 2/3 1/3
    1 1 1/3    
    1 0 2/3
    1 1/3 2/3
    1 2/3 2/3
    1 1 2/3     
    1 0 1
    1 1/3 1
    1 2/3 1
    1 1 1      
];
face{3}.nLoop = 1; 
face{3}.pLoop{1}.deg = 2; 
face{3}.pLoop{1}.knotv = [0 0 0 1 1 1];
face{3}.pLoop{1}.cp = [
    0 0   
    .5 0
    1 0
    ];
face{3}.pLoop{2}.deg = 2; 
face{3}.pLoop{2}.knotv = [0 0 0 1 1 1];
face{3}.pLoop{2}.cp = [
    1 0   
    1 .5
    1 1
    ];
face{3}.pLoop{3}.deg = 2; 
face{3}.pLoop{3}.knotv = [0 0 0 1 1 1];
face{3}.pLoop{3}.cp = [
    1 1   
    .5 1
    0 1
    ];
face{3}.pLoop{4}.deg = 2; 
face{3}.pLoop{4}.knotv = [0 0 0 1 1 1];
face{3}.pLoop{4}.cp = [
    0 1   
    0 .5
    0 0
    ];

% -------------------------------------------------------------------------
% face 4
face{4}.pCorner = [4 3 7 8]; 
face{4}.pDeg1 = 3; 
face{4}.pDeg2 = 3; 
face{4}.pKnotv1 = [0 0 0 0 1 1 1 1];
face{4}.pKnotv2 = [0 0 0 0 1 1 1 1];
face{4}.pCp = [
        0 1 0
        1/3 1 0
        2/3 1 0
        1 1 0
        0 1 1/3
        1/3 1 1/3
        2/3 1 1/3
        1 1 1/3       
        0 1 2/3
        1/3 1 2/3
        2/3 1 2/3
        1 1 2/3      
        0 1 1
        1/3 1 1
        2/3 1 1
        1 1 1          
];
face{4}.nLoop = 1; 
face{4}.pLoop{1}.deg = 2; 
face{4}.pLoop{1}.knotv = [0 0 0 1 1 1];
face{4}.pLoop{1}.cp = [
    0 0   
    .5 0
    1 0
    ];
face{4}.pLoop{2}.deg = 2; 
face{4}.pLoop{2}.knotv = [0 0 0 1 1 1];
face{4}.pLoop{2}.cp = [
    1 0   
    1 .5
    1 1
    ];
face{4}.pLoop{3}.deg = 2; 
face{4}.pLoop{3}.knotv = [0 0 0 1 1 1];
face{4}.pLoop{3}.cp = [
    1 1   
    .5 1
    0 1
    ];
face{4}.pLoop{4}.deg = 2; 
face{4}.pLoop{4}.knotv = [0 0 0 1 1 1];
face{4}.pLoop{4}.cp = [
    0 1   
    0 .5
    0 0
    ];

% -------------------------------------------------------------------------
% face 5
face{5}.pCorner = [1 4 8 5]; 
face{5}.pDeg1 = 3; 
face{5}.pDeg2 = 3; 
face{5}.pKnotv1 = [0 0 0 0 1 1 1 1];
face{5}.pKnotv2 = [0 0 0 0 1 1 1 1];
face{5}.pCp = [
        0 0 0 
        0 1/3 0
        0 2/3 0
        0 1 0
        0 0 1/3 
        0 1/3 1/3
        0 2/3 1/3
        0 1 1/3       
        0 0 2/3 
        0 1/3 2/3
        0 2/3 2/3
        0 1 2/3  
        0 0 1 
        0 1/3 1
        0 2/3 1
        0 1 1         
];
face{5}.nLoop = 1; 
face{5}.pLoop{1}.deg = 2; 
face{5}.pLoop{1}.knotv = [0 0 0 1 1 1];
face{5}.pLoop{1}.cp = [
    0 0   
    .5 0
    1 0
    ];
face{5}.pLoop{2}.deg = 2; 
face{5}.pLoop{2}.knotv = [0 0 0 1 1 1];
face{5}.pLoop{2}.cp = [
    1 0   
    1 .5
    1 1
    ];
face{5}.pLoop{3}.deg = 2; 
face{5}.pLoop{3}.knotv = [0 0 0 1 1 1];
face{5}.pLoop{3}.cp = [
    1 1   
    .5 1
    0 1
    ];
face{5}.pLoop{4}.deg = 2; 
face{5}.pLoop{4}.knotv = [0 0 0 1 1 1];
face{5}.pLoop{4}.cp = [
    0 1   
    0 .5
    0 0
    ];

% -------------------------------------------------------------------------
% face 6
face{6}.pCorner = [5 6 7 8]; 
face{6}.pDeg1 = 3; 
face{6}.pDeg2 = 3; 
face{6}.pKnotv1 = [0 0 0 0 1 1 1 1];
face{6}.pKnotv2 = [0 0 0 0 1 1 1 1];
face{6}.pCp = [
        0 0 1
        1/3 0 1
        2/3 0 1
        1 0 1
        0 1/3 1
        1/3 1/3 1
        2/3 1/3 1
        1 1/3 1    
        0 2/3 1
        1/3 2/3 1
        2/3 2/3 1
        1 2/3 1   
        0 1 1
        1/3 1 1
        2/3 1 1
        1 1 1           
];
face{6}.nLoop = 1; 
face{6}.pLoop{1}.deg = 2; 
face{6}.pLoop{1}.knotv = [0 0 0 1 1 1];
face{6}.pLoop{1}.cp = [
    0 0   
    .5 0
    1 0
    ];
face{6}.pLoop{2}.deg = 2; 
face{6}.pLoop{2}.knotv = [0 0 0 1 1 1];
face{6}.pLoop{2}.cp = [
    1 0   
    1 .5
    1 1
    ];
face{6}.pLoop{3}.deg = 2; 
face{6}.pLoop{3}.knotv = [0 0 0 1 1 1];
face{6}.pLoop{3}.cp = [
    1 1   
    .5 1
    0 1
    ];
face{6}.pLoop{4}.deg = 2; 
face{6}.pLoop{4}.knotv = [0 0 0 1 1 1];
face{6}.pLoop{4}.cp = [
    0 1   
    0 .5
    0 0
    ];
