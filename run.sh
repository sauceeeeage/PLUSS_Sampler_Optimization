#!/usr/bin/env bash

METHOD=$1

cd c_lib/test
make clean
make "$METHOD"

cd ../..

cargo build --release --bin RUST_PLUSS
./target/release/RUST_PLUSS "$METHOD" >> output.txt

