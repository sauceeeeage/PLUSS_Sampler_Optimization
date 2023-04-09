include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

#[test]
fn test() {
    unsafe {
        assert_eq!(something(), 0);
    }
}
