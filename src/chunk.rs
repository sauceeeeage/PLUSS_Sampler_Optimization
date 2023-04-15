pub(crate) struct Chunk {
    lb: i32,
    ub: i32,
}

impl Chunk {
    pub fn new(lb: i32, ub: i32) -> Chunk {
        Chunk {
            lb,
            ub,
        }
    }
    pub fn first(&self) -> i32 {
        self.lb
    }
    pub fn second(&self) -> i32 {
        self.ub
    }
}