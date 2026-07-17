use thiserror::Error;
use moyo::data::Setting;

#[derive(Debug, Error)]
pub enum CoreError {
    #[error("Invalid space group number: {0}. Must be between 1 and 230.")]
    InvalidSpaceGroupNumber(u16),
    #[error("Invalid space group settings: {0:?}, primitive: {1}.")]
    InvalidSpaceGroupSettings(Setting, bool),
}
