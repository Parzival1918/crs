use std::fs::File;
use std::io::{BufRead, BufReader, Result as IoResult, Write};
use std::path::Path;

/// A trait for parsing atomic data (like `Molecule` or `Crystal`) from a data source.
pub trait Parser<T> {
    type E: From<std::io::Error>;

    /// Parse a single item from a reader. Returns `Ok(None)` if end-of-file is reached cleanly.
    fn parse_from_reader<R: BufRead>(&self, reader: &mut R) -> Result<Option<T>, Self::E>;

    /// Return an iterator over all items in the reader.
    fn parse_many_from_reader<'a, R: BufRead + 'a>(
        &'a self,
        mut reader: R,
    ) -> impl Iterator<Item = Result<T, Self::E>> + 'a {
        std::iter::from_fn(move || match self.parse_from_reader(&mut reader) {
            Ok(Some(item)) => Some(Ok(item)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        })
    }

    /// Parse a single item from a file. Provides a default implementation.
    fn parse_from_file<P: AsRef<Path>>(&self, path: P) -> Result<T, Self::E> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        self.parse_from_reader(&mut reader).and_then(|res| {
            res.ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "Empty file").into()
            })
        })
    }

    /// Parse a single item from a string. Provides a default implementation.
    fn parse_from_string(&self, s: &str) -> Result<T, Self::E> {
        let mut reader = s.as_bytes();
        self.parse_from_reader(&mut reader).and_then(|res| {
            res.ok_or_else(|| {
                std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "Empty string").into()
            })
        })
    }

    /// Parse all items from a file into a Vec. Provides a default implementation.
    fn parse_many_from_file<'a, P: AsRef<Path>>(
        &'a self,
        path: P,
    ) -> Result<impl Iterator<Item = Result<T, Self::E>> + 'a, Self::E> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        Ok(self.parse_many_from_reader(reader))
    }

    /// Parse all items from a string into a Vec. Provides a default implementation.
    fn parse_many_from_string<'a>(
        &'a self,
        s: &'a str,
    ) -> impl Iterator<Item = Result<T, Self::E>> + 'a {
        let reader = s.as_bytes();
        self.parse_many_from_reader(reader)
    }
}

/// A trait for writing atomic data (like `Molecule` or `Crystal`) to a data sink.
pub trait Writer<T> {
    /// Write a single item to a writer.
    fn write_to_writer<W: Write>(&self, data: &T, writer: &mut W) -> IoResult<()>;

    /// Write multiple items to a writer.
    fn write_many_to_writer<'a, W: Write, I: IntoIterator<Item = &'a T>>(
        &self,
        data: I,
        mut writer: W,
    ) -> IoResult<()>
    where
        T: 'a,
    {
        for item in data {
            self.write_to_writer(item, &mut writer)?;
        }
        Ok(())
    }

    /// Write a single item to a file. Provides a default implementation.
    fn write_to_file<P: AsRef<Path>>(&self, data: &T, path: P) -> IoResult<()> {
        let mut file = File::create(path)?;
        self.write_to_writer(data, &mut file)
    }

    /// Write multiple items to a file.
    fn write_many_to_file<'a, P: AsRef<Path>, I: IntoIterator<Item = &'a T>>(
        &self,
        data: I,
        path: P,
    ) -> IoResult<()>
    where
        T: 'a,
    {
        let mut file = File::create(path)?;
        for item in data {
            self.write_to_writer(item, &mut file)?;
        }
        Ok(())
    }

    /// Write a single item to a string. Provides a default implementation.
    fn write_to_string(&self, data: &T) -> IoResult<String> {
        let mut buf = Vec::new();
        self.write_to_writer(data, &mut buf)?;
        String::from_utf8(buf).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }

    /// Write multiple items to a string.
    fn write_many_to_string<'a, I: IntoIterator<Item = &'a T>>(&self, data: I) -> IoResult<String>
    where
        T: 'a,
    {
        let mut buf = Vec::new();
        for item in data {
            self.write_to_writer(item, &mut buf)?;
        }
        String::from_utf8(buf).map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }
}
