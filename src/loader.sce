// this loader is not generated 
// we have added the global=%t option 
// because some other dynamic libraries 
// may want to use symbols defined internally here.

libcoinor_path=file('join',['.','libcoinor'+%shext]);
addinter(libcoinor_path,'libcoinor',global=%t);
clear libcoinor_path;
