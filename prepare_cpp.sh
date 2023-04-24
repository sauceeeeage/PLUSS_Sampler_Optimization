# get prepared
# brew install binutils


cd c_lib
g++ -c pluss.cpp
/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ar rvs libpluss.a pluss.o 
cd ..

rustc src/main.rs --extern pluss=libpluss.a -Lc_lib