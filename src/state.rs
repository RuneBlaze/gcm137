use ahash::AHashMap;

pub struct AlnState {
    pub names: Vec<String>,
    pub names2id: AHashMap<String, usize>, // TODO: there is no need to keep 2 copies of the same string
    pub s: Vec<Vec<(u32, u32)>>,
    pub column_counts: Vec<usize>,
}

impl AlnState {
    pub fn request_name(&mut self, name: &str) -> usize {
        let id = self
            .names2id
            .entry(name.to_string())
            .or_insert(self.names.len());
        if *id == self.names.len() {
            self.names.push(name.to_string());
        };
        *id
    }

    pub fn retrieve(&self, name: &str) -> usize {
        *self.names2id.get(name).unwrap()
    }

    pub fn new() -> Self {
        Self {
            names: vec![],
            names2id: AHashMap::default(),
            s: vec![],
            column_counts: vec![],
        }
    }

    pub fn ncols(&self) -> usize {
        self.column_counts.len()
    }
}
