eps = 1.e-6;

function linspace( x0, x1, N ) = [ for( i=[0:N] ) x0 + (i/N) * (x1-x0) ];

module tubular( intervals, profile, path, dpath_t, 
                scale=false, rotation=false, prin_axis=[0,0,1] ) {
    I = len( path ) - 1;
    J = len( profile );

    // Set defaults
    _rot = [ for(i=[0:I]) (rotation==false) ?  [ [0,1,0], [0,0,1], [1,0,0] ] :
             [ [ 0., cos( rotation[i] ), -sin( rotation[i] ) ],
               [ 0., sin( rotation[i] ),  cos( rotation[i] ) ],
               [ 1., 0.,                  0.                 ] ]] ;    
    _s =   (scale==false) ? [ for(i=[0:I]) 1. ] : scale;

    // Set orientation along path
    X =      [ for( v = dpath_t ) cross( v, prin_axis ) ];
    Y =      [ for( i = [0:I] ) cross( X[i], dpath_t[i] ) ];
    orient = [ for( i = [0:I]) [ dpath_t[i]/norm(dpath_t[i]),
                                 X[i]/norm(X[i]), 
                                 Y[i]/norm(Y[i]) ] ];

    // Compute points on surface
    all_pts = [ for( i = [0:I] ) 
                for( v = profile * _s[i] * _rot[i] * orient[i] ) path[i] + v ];

    // Compute facets for main body
    body_faces = [ for( i = [0:I-1] ) 
                   for( j = [0:J-2] )  
                   [ i*J+j, i*J+(j+1), (i+1)*J+(j+1), (i+1)*J+j ] ];
                   
    // Add end caps
    end_faces =  [ [ for( j = [J-1:-1:0] ) j ], 
                   [ for( j = [0:J-1] ) j + I*J ] ];

    polyhedron( points=all_pts, 
                faces=concat( body_faces, end_faces ) );    
}
