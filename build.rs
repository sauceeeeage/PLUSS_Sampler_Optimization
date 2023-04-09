// use std::env;
// use std::path::PathBuf;
// use cbindgen::Config;
//
// fn main() {
//     let crate_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
//
//     let package_name = env::var("CARGO_PKG_NAME").unwrap();
//     let output_file = target_dir()
//         .join(format!("{}.hpp", package_name))
//         .display()
//         .to_string();
//
//     let config = Config {
//         namespace: Some(String::from("ffi")),
//         ..Default::default()
//     };
//
//     cbindgen::generate_with_config(&crate_dir, config)
//         .unwrap()
//         .write_to_file(&output_file);
// }
//
// /// Find the location of the `target/` directory. Note that this may be
// /// overridden by `cmake`, so we also need to check the `CARGO_TARGET_DIR`
// /// variable.
// fn target_dir() -> PathBuf {
//     if let Ok(target) = env::var("CARGO_TARGET_DIR") {
//         PathBuf::from(target)
//     } else {
//         PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap()).join("target")
//     }
// }

extern crate bindgen;
extern crate cc;

fn main() {
    cc::Build::new().file("c_lib/playground.c").compile("playground");

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .generate()
        .expect("Couldn't generate bindings!");

    let out_path = std::path::PathBuf::from(std::env::var("OUT_DIR").unwrap());

    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}