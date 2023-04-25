#!/usr/bin/env bash

# need to set RUSTFLAGS to include the path to the homebrew gsl library
#RUSTFLAGS='-L /opt/homebrew/include'

cargo build --release --bin RUST_PLUSS

export RUST_BACKTRACE=1
export RUST_LOG=info

./target/release/RUST_PLUSS >> output.txt