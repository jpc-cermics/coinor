// This is the builder.sce 
// must be run from this directory 

// [1] generate Path.incl 
if exists('%nsp') then 
  ilib_path_incl()
end 

// [3] the part devoted to shared lib generation 
ilib_name  = 'libcoinor' 		// interface library name 

// objects files (but do not give mexfiles here)
files = ['nspcoinor-IN.o';'nspcoinmp.o';'nspcoinread.o';'nspclp.o'];

// other libs needed for linking (must be shared library names)
// libs  = ['lib/liblpsolve55']; 				
libs =[];

// table of (scilab_name,interface-name or mexfile-name, type) 
table =[]; // 'lpsolve', 'int_lpsolve',	'cnsp'];

// we need  gtk+-2.0 because of imatrices 
ldflags = "`pkg-config gtk+-2.0 coinmp --libs`";
cflags = "`pkg-config gtk+-2.0 coinmp --cflags`"

//function ilib_build(ilib_name,table,files,libs,...
//		    makename='Makelib',ldflags="",cflags="",fflags="",verbose=%t)
// do not modify below 
// ----------------------------------------------

ilib_build(ilib_name,table,files,libs,ldflags = ldflags,cflags = cflags );


