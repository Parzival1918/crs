use moyo::data::Setting;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum CoreError {
    #[error("Invalid space group number: {0}. Must be between 1 and 230.")]
    InvalidSpaceGroupNumber(u16),
    #[error("Invalid space group settings: {0:?}, primitive: {1}.")]
    InvalidSpaceGroupSettings(Setting, bool),
    #[error("Invalid symmetry operation string: {0}")]
    InvalidSymOpString(String),
    #[error("Space group not found for the given symmetry operations: {0:?}")]
    SpaceGroupNotFound(Vec<String>),
}
