#[link(name = "pluss")]
extern "C" {
    fn pluss_timer_start();
    fn pluss_timer_stop();
    fn pluss_timer_print();
}

fn main() {
    unsafe {
        pluss_timer_start();
        let mut a = 1;
        a = a + 1;
        pluss_timer_stop();
        pluss_timer_print();
    }
}
