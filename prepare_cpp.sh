#!/usr/bin/env bash

cd c_lib
/opt/homebrew/Cellar/gcc/12.2.0/bin/g++-12 -c pluss.cpp
/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ar rvs libpluss.a pluss.o -v
cd ..

#RUSTFLAGS='-L /opt/homebrew/include'
#rustc src/.rs --extern pluss=libpluss.a -Lc_lib

cargo build --bin RUST_PLUSS

export RUST_BACKTRACE=1
export RUST_LOG=info

./target/debug/RUST_PLUSS