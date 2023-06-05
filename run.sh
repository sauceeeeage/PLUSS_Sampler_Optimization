#!/usr/bin/env bash

METHOD=$1

cd c_lib/test
make clean
make "$METHOD"

cd ../..

env RUSTFLAGS="-C target-cpu=native" cargo build --release --bin RUST_PLUSS
./target/release/RUST_PLUSS "$METHOD" >> output.txt

