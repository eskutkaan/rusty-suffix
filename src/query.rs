use crate::fasta::Sequence;

#[derive(Debug, Clone)]
pub struct QueryBatch {
    pub queries: Vec<Sequence>,
    pub batch_size: usize,
}

impl QueryBatch {
    pub fn new(queries: Vec<Sequence>, batch_size: usize) -> Self {
        QueryBatch { queries, batch_size }
    }

    pub fn from_sequences(sequences: Vec<Sequence>) -> Self {
        let batch_size = sequences.len();
        QueryBatch {
            queries: sequences,
            batch_size,
        }
    }

    pub fn batches(&self) -> Vec<QueryBatch> {
        self.queries
            .chunks(self.batch_size)
            .map(|chunk| QueryBatch {
                queries: chunk.to_vec(),
                batch_size: chunk.len(),
            })
            .collect()
    }

    pub fn len(&self) -> usize {
        self.queries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.queries.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &Sequence> {
        self.queries.iter()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_query_batch_creation() {
        let queries = vec![
            Sequence {
                id: "seq1".to_string(),
                sequence: b"ACGT".to_vec(),
            },
            Sequence {
                id: "seq2".to_string(),
                sequence: b"TGCA".to_vec(),
            },
        ];

        let batch = QueryBatch::new(queries.clone(), 1);
        assert_eq!(batch.len(), 2);
        assert_eq!(batch.batch_size, 1);

        let batches = batch.batches();
        assert_eq!(batches.len(), 2);
    }
}
