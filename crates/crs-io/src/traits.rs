use std::fs::File;
use std::io::{BufRead, BufReader, Result as IoResult, Write};
use std::path::Path;

/// A trait for parsing atomic data (like `Molecule` or `Crystal`) from a data source.
pub trait Parser<T> {
    type E: From<std::io::Error>;

    /// Parse a single item from a reader. Returns `Ok(None)` if end-of-file is reached cleanly.
    fn parse_from_reader<R: BufRead>(reader: &mut R) -> Result<Option<T>, Self::E>;

    /// Return an iterator over all items in the reader.
    fn parse_many_from_reader<R: BufRead>(
        mut reader: R,
    ) -> impl Iterator<Item = Result<T, Self::E>> {
        std::iter::from_fn(move || match Self::parse_from_reader(&mut reader) {
            Ok(Some(item)) => Some(Ok(item)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        })
    }

    /// Parse a single item from a file. Provides a default implementation.
    fn parse_from_file<P: AsRef<Path>>(path: P) -> Result<T, Self::E> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        Self::parse_from_reader(&mut reader).and_then(|res| {
            res.ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "Empty file").into()
            })
        })
    }

    /// Parse a single item from a string. Provides a default implementation.
    fn parse_from_string(s: &str) -> Result<T, Self::E> {
        let mut reader = s.as_bytes();
        Self::parse_from_reader(&mut reader).and_then(|res| {
            res.ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "Empty string").into()
            })
        })
    }

    /// Parse all items from a file into a Vec. Provides a default implementation.
    fn parse_many_from_file<P: AsRef<Path>>(
        path: P,
    ) -> Result<impl Iterator<Item = Result<T, Self::E>>, Self::E> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        Ok(Self::parse_many_from_reader(reader))
    }

    /// Parse all items from a string into a Vec. Provides a default implementation.
    fn parse_many_from_string<'a>(s: &'a str) -> impl Iterator<Item = Result<T, Self::E>> + 'a
    where
        Self: 'a,
        T: 'a,
    {
        let reader = s.as_bytes();
        Self::parse_many_from_reader(reader)
    }
}

/// A trait for writing atomic data (like `Molecule` or `Crystal`) to a data sink.
pub trait Writer {
    /// Write a single item to a writer.
    fn write_to_writer<W: Write>(&self, writer: &mut W) -> IoResult<()>;

    /// Write a single item to a file. Provides a default implementation.
    fn write_to_file<P: AsRef<Path>>(&self, path: P) -> IoResult<()> {
        let mut file = File::create(path)?;
        self.write_to_writer(&mut file)
    }

    /// Write a single item to a string. Provides a default implementation.
    fn write_to_string(&self) -> IoResult<String> {
        let mut buf = Vec::new();
        self.write_to_writer(&mut buf)?;
        String::from_utf8(buf).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }
}
