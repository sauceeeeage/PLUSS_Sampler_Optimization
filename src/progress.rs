use crate::chunk::Chunk;
pub (crate)struct Progress {
    pub(crate) refs: &'static str, // cannot name it ref because it's a keyword
    pub(crate) chunk: Chunk,
    pub(crate) iteration: Vec<i32>,
}

impl Progress {
    pub fn new() -> Self {
        Progress {
            refs: "",
            chunk: Chunk::new(0, 0),
            iteration: Vec::new(),
        }
    }
    pub fn new_with_ref(refs: &'static str, iteration: Vec<i32>, chunk: Chunk) -> Self {
        Progress {
            refs,
            chunk,
            iteration,
        }
    }
    pub fn get_iteration(&self) -> String {
        let mut ret = String::from("(");
        for i in 0..self.iteration.len() {
            ret += &self.iteration[i].to_string();
            if i != self.iteration.len() - 1 {
                ret += ",";
            }
        }
        ret += ")";
        ret
    }
    pub fn get_reference(&self) -> &'static str {
        self.refs.clone()
    }
    pub fn increment(&mut self, refs: &'static str, iteration: Vec<i32>) {
        self.refs = refs;
        self.iteration = iteration;
    }
    pub fn increment_with_ref(&mut self, refs: &'static str) {
        self.refs = refs;
    }
    pub fn is_in_bound(&self) -> bool {
        assert!(self.iteration[0] >= self.chunk.first());
        // println!("{} <= {}", self.iteration[0], self.chunk.second());
        self.iteration[0] <= self.chunk.second()
    }
    pub fn to_string(&self) -> String {
        let mut ss = String::new();
        ss += "(";
        ss += &self.refs;
        ss += ",";
        ss += &self.get_iteration();
        ss += ")";
        ss
    }
}