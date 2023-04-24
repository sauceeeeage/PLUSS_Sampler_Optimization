// class Progress {
// public:
// string ref;
// Chunk chunk;
// vector<int> iteration;
// Progress() { }
// Progress(string ref, vector<int> iteration, Chunk c) {
// this->ref = ref;
// this->iteration = iteration;
// this->chunk = c;
// }
// string getIteration() {
// string ret = "(";
// for (unsigned i = 0; i < this->iteration.size(); i++) {
// ret += to_string(this->iteration[i]);
// if (i != this->iteration.size() - 1)
// ret += ",";
// }
// ret += ")";
// return ret;
// }
// string getReference() {
// return this->ref;
// }
// void increment(string ref, vector<int> iteration) {
// this->ref = ref;
// this->iteration = iteration;
// }
// void increment(string ref) {
// this->ref = ref;
// }
// bool isInBound() {
// assert(this->iteration[0] >= chunk.first);
// return this->iteration[0] <= chunk.second;
// }
// string toString() {
// stringstream ss;
// ss << "(" << this->ref << ",";
// ss << getIteration();
// ss << ")";
// return ss.str();
// }
// };
use crate::chunk::Chunk;
pub (crate)struct Progress {
    pub(crate) refs: String, // cannot name it ref because it's a keyword
    pub(crate) chunk: Chunk,
    pub(crate) iteration: Vec<i32>,
}

impl Progress {
    pub fn new() -> Self {
        Progress {
            refs: String::new(),
            chunk: Chunk::new(0, 0),
            iteration: Vec::new(),
        }
    }
    pub fn new_with_ref(refs: String, iteration: Vec<i32>, chunk: Chunk) -> Self {
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
    pub fn get_reference(&self) -> String {
        self.refs.clone()
    }
    pub fn increment(&mut self, refs: String, iteration: Vec<i32>) {
        self.refs = refs;
        self.iteration = iteration;
    }
    pub fn increment_with_ref(&mut self, refs: String) {
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