use thiserror::Error;

#[derive(Debug, Error)]
pub enum CoreError {
    #[error("Invalid space group number: {0}. Must be between 1 and 230.")]
    InvalidSpaceGroupNumber(u16),
}
