pub(crate) struct cParams {
    iInsertMax: i64,
    iUpdateMax: i64,
    iDeleteMax: i64,
    iInstanceMax: i64,
    tFirstInstance: bool,
    tCreateTables: bool,
    tContinue: bool,
}

impl cParams {
    pub fn new() -> cParams {
        cParams {
            iInsertMax: -1,
            iUpdateMax: -1,
            iDeleteMax: -1,
            iInstanceMax: -1,
            tFirstInstance: false,
            tCreateTables: false,
            tContinue: false,
        }
    }
}

fn main() {
    let p = cParams::new();
}