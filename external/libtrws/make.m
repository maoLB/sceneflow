cmdLine = ['mex ', 'trwsMex.cpp', ' -output trwsMex -largeArrayDims CXXFLAGS="\$CXXFLAGS -std=c++0x -fpermissive"'];
eval(cmdLine);
disp('done!');