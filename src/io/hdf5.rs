//! A minimal HDF5 writer, with no dependency on libhdf5.
//!
//! This exists for one file: CellRanger's `raw_feature_bc_matrix.h5`. That file
//! uses a small, fixed corner of the format — one group tree, nine 1-D
//! datasets, three scalar types, a handful of root attributes — and writing
//! that corner directly is a few hundred lines, against a C library that would
//! need either a system `libhdf5` or `cmake` on every build (see
//! `DIVERGENCE.md` §3.5 and the discussion in the HDF5 dependency issue).
//!
//! # What it writes
//!
//! Version 0 superblock, version 1 object headers, old-style groups (symbol
//! table + version 1 B-tree + local heap), and contiguous, uncompressed
//! datasets. Every one of those is the oldest and most widely-read variant of
//! its structure, which is deliberate: this writer's only job is to produce
//! files other people's readers open.
//!
//! # What it does not write
//!
//! No chunking, no compression, no variable-length types, no references, no
//! links other than the hard links in a group's symbol table, no deletion, no
//! rewriting. Every dataset is written once, in full, from a slice already in
//! memory. Groups hold at most `2 * GROUP_LEAF_K` = 8 entries, since a second
//! B-tree node would need node splitting and nothing here comes close.
//!
//! The layout constants and structure layouts follow the *HDF5 File Format
//! Specification Version 3.0* (superblock §II.A, B-trees §II.A.1, local heaps
//! §II.D, object headers §IV.A, messages §IV.A.2).

use std::io::Write;
use std::path::Path;

use crate::error::Error;

/// `\x89HDF\r\n\x1a\n` — the format signature, chosen (like PNG's) so that a
/// transfer that mangles high bits or line endings is detectable.
const SIGNATURE: [u8; 8] = [0x89, b'H', b'D', b'F', b'\r', b'\n', 0x1a, b'\n'];

/// Sizes of offsets and lengths in the file. 8 bytes each, so a file may exceed
/// 4 GB; a raw 10x matrix will not, but a genome-scale one might.
const SIZEOF_ADDR: usize = 8;

/// Group node K values recorded in the superblock. A leaf symbol-table node
/// holds up to `2 * GROUP_LEAF_K` entries and a B-tree node up to
/// `2 * GROUP_INTERNAL_K` children; we only ever produce one of each.
const GROUP_LEAF_K: u16 = 4;
const GROUP_INTERNAL_K: u16 = 16;

/// The "undefined address" HDF5 uses where a pointer is absent: all bits set.
const UNDEF_ADDR: u64 = u64::MAX;

/// "No free block" in a local heap's free list. Not the undefined address: the
/// library uses the literal 1, which cannot be a real offset because free
/// blocks are 8-byte aligned (`H5HL_FREE_NULL` in `H5HLpkg.h`). Writing the
/// undefined address here makes `H5Ovisit` fail with "bad heap free list"
/// while the file still opens, which is a memorable way to learn the
/// difference.
const HEAP_FREE_NULL: u64 = 1;

/// Superblock v0 is 56 bytes, followed by the 40-byte root symbol table entry.
const SUPERBLOCK_LEN: usize = 96;

/// Object header message types used here.
const MSG_DATASPACE: u16 = 0x0001;
const MSG_DATATYPE: u16 = 0x0003;
const MSG_FILL_VALUE: u16 = 0x0005;
const MSG_DATA_LAYOUT: u16 = 0x0008;
const MSG_ATTRIBUTE: u16 = 0x000C;
const MSG_SYMBOL_TABLE: u16 = 0x0011;

/// One dataset's contents. Every variant is a 1-D array; HDF5 scalars appear
/// only as attributes.
///
/// `FixedString` is the padded, fixed-width byte string CellRanger uses for
/// barcodes and feature fields (`|S45`, `|S256` in numpy terms). Values longer
/// than `width` are an error rather than a silent truncation, because a
/// truncated barcode is a wrong answer that looks like a right one.
pub enum Data<'a> {
    I32(&'a [i32]),
    I64(&'a [i64]),
    FixedString { width: usize, values: Vec<&'a str> },
}

impl Data<'_> {
    fn len(&self) -> usize {
        match self {
            Data::I32(v) => v.len(),
            Data::I64(v) => v.len(),
            Data::FixedString { values, .. } => values.len(),
        }
    }

    fn element_size(&self) -> usize {
        match self {
            Data::I32(_) => 4,
            Data::I64(_) => 8,
            Data::FixedString { width, .. } => *width,
        }
    }

    /// The datatype message body for this data (see `datatype_message`).
    fn datatype(&self) -> Vec<u8> {
        match self {
            Data::I32(_) => fixed_point_datatype(4),
            Data::I64(_) => fixed_point_datatype(8),
            Data::FixedString { width, .. } => string_datatype(*width),
        }
    }

    /// The raw bytes of the array, little-endian, in element order.
    fn bytes(&self) -> Result<Vec<u8>, String> {
        let mut out = Vec::with_capacity(self.len() * self.element_size());
        match self {
            Data::I32(v) => v
                .iter()
                .for_each(|x| out.extend_from_slice(&x.to_le_bytes())),
            Data::I64(v) => v
                .iter()
                .for_each(|x| out.extend_from_slice(&x.to_le_bytes())),
            Data::FixedString { width, values } => {
                for s in values {
                    if s.len() > *width {
                        return Err(format!(
                            "value {s:?} is {} bytes, longer than the {width}-byte field",
                            s.len()
                        ));
                    }
                    out.extend_from_slice(s.as_bytes());
                    out.resize(out.len() + (width - s.len()), 0);
                }
            }
        }
        Ok(out)
    }
}

/// A scalar or 1-D attribute attached to a group.
pub enum AttrValue {
    /// Fixed-width string, written at exactly its own length.
    Str(String),
    /// 1-D array of fixed-width strings, each padded to `width`.
    StrArray {
        width: usize,
        values: Vec<String>,
    },
    I64(i64),
    I64Array(Vec<i64>),
}

pub struct Attr {
    pub name: String,
    pub value: AttrValue,
}

pub struct DatasetSpec<'a> {
    pub name: &'a str,
    pub data: Data<'a>,
}

/// A group: named sub-groups, named datasets, and attributes. Children are
/// sorted by name at write time, as a symbol table requires.
pub struct GroupSpec<'a> {
    pub name: &'a str,
    pub groups: Vec<GroupSpec<'a>>,
    pub datasets: Vec<DatasetSpec<'a>>,
    pub attrs: Vec<Attr>,
}

impl<'a> GroupSpec<'a> {
    #[must_use]
    pub fn new(name: &'a str) -> Self {
        Self {
            name,
            groups: Vec::new(),
            datasets: Vec::new(),
            attrs: Vec::new(),
        }
    }

    #[must_use]
    pub fn dataset(mut self, name: &'a str, data: Data<'a>) -> Self {
        self.datasets.push(DatasetSpec { name, data });
        self
    }

    #[must_use]
    pub fn group(mut self, g: GroupSpec<'a>) -> Self {
        self.groups.push(g);
        self
    }

    #[must_use]
    pub fn attr(mut self, name: &str, value: AttrValue) -> Self {
        self.attrs.push(Attr {
            name: name.to_string(),
            value,
        });
        self
    }

    fn child_count(&self) -> usize {
        self.groups.len() + self.datasets.len()
    }
}

// ---------------------------------------------------------------------------
// Byte-level helpers
// ---------------------------------------------------------------------------

/// The file image under construction. Blocks are appended and never moved, so
/// an address handed out stays valid; everything is 8-byte aligned, which every
/// structure here wants anyway.
struct Image {
    buf: Vec<u8>,
}

impl Image {
    fn new() -> Self {
        // The superblock is written last, once the root group's address and the
        // end-of-file address are known, so leave room for it.
        Self {
            buf: vec![0u8; SUPERBLOCK_LEN],
        }
    }

    fn align(&mut self) {
        while !self.buf.len().is_multiple_of(8) {
            self.buf.push(0);
        }
    }

    /// Append a block and return its address.
    fn alloc(&mut self, bytes: &[u8]) -> u64 {
        self.align();
        let addr = self.buf.len() as u64;
        self.buf.extend_from_slice(bytes);
        addr
    }

    /// Reserve `len` zeroed bytes and return the address, for a block that is
    /// filled in after its own contents are known.
    fn reserve(&mut self, len: usize) -> u64 {
        self.align();
        let addr = self.buf.len() as u64;
        self.buf.resize(self.buf.len() + len, 0);
        addr
    }

    fn patch(&mut self, addr: u64, bytes: &[u8]) {
        let start = addr as usize;
        self.buf[start..start + bytes.len()].copy_from_slice(bytes);
    }
}

fn pad_to_8(v: &mut Vec<u8>) {
    while !v.len().is_multiple_of(8) {
        v.push(0);
    }
}

/// The first byte of a datatype message: version in the high nibble, class in
/// the low one.
fn datatype_tag(class: u8) -> u8 {
    (1 << 4) | class
}

/// Datatype message body for a signed little-endian integer of `size` bytes.
///
/// Class 0 (fixed-point). The class bit field's bit 3 marks it signed; byte
/// order (bit 0) and padding (bits 1-2) are zero, i.e. little-endian with zero
/// padding. The properties are the bit offset and precision.
fn fixed_point_datatype(size: u32) -> Vec<u8> {
    const CLASS_FIXED_POINT: u8 = 0;
    let mut m = Vec::with_capacity(12);
    m.push(datatype_tag(CLASS_FIXED_POINT));
    m.extend_from_slice(&[0x08, 0x00, 0x00]); // signed
    m.extend_from_slice(&size.to_le_bytes());
    m.extend_from_slice(&0u16.to_le_bytes()); // bit offset
    m.extend_from_slice(&((size * 8) as u16).to_le_bytes()); // precision
    m
}

/// Datatype message body for a fixed-width, null-padded ASCII string.
///
/// Class 3 (string), padding type 1 (null pad) in bits 0-3 and character set 0
/// (ASCII) in bits 4-7 — the combination numpy reads back as `|S<width>`.
fn string_datatype(width: usize) -> Vec<u8> {
    let mut m = Vec::with_capacity(8);
    const CLASS_STRING: u8 = 3;
    m.push(datatype_tag(CLASS_STRING));
    m.extend_from_slice(&[0x01, 0x00, 0x00]); // null pad, ASCII
    m.extend_from_slice(&(width as u32).to_le_bytes());
    m
}

/// Dataspace message body: version 1, `dims.len()` dimensions, no maximum
/// dimensions and no permutation index. An empty `dims` is a scalar.
fn dataspace_message(dims: &[u64]) -> Vec<u8> {
    let mut m = Vec::with_capacity(8 + dims.len() * 8);
    m.push(1); // version
    m.push(dims.len() as u8);
    m.push(0); // flags: no max dims
    m.extend_from_slice(&[0, 0, 0, 0, 0]); // reserved
    for d in dims {
        m.extend_from_slice(&d.to_le_bytes());
    }
    m
}

/// Fill value message, version 2, with no fill value defined: readers use the
/// type's default. Nothing here is ever read before it is written.
fn fill_value_message() -> Vec<u8> {
    vec![2, 2, 0, 0] // version, space allocation time (early), write time, undefined
}

/// Data layout message, version 3, contiguous storage at `addr` for `size`
/// bytes.
fn data_layout_message(addr: u64, size: u64) -> Vec<u8> {
    let mut m = Vec::with_capacity(18);
    m.push(3); // version
    m.push(1); // class 1: contiguous
    m.extend_from_slice(&addr.to_le_bytes());
    m.extend_from_slice(&size.to_le_bytes());
    m
}

/// One object header message, prefixed and padded as the header expects.
fn header_message(kind: u16, body: &[u8]) -> Vec<u8> {
    let mut padded = body.to_vec();
    pad_to_8(&mut padded);
    let mut m = Vec::with_capacity(8 + padded.len());
    m.extend_from_slice(&kind.to_le_bytes());
    m.extend_from_slice(&(padded.len() as u16).to_le_bytes());
    m.push(0); // flags
    m.extend_from_slice(&[0, 0, 0]); // reserved
    m.extend_from_slice(&padded);
    m
}

/// A version 1 object header wrapping already-encoded messages.
fn object_header(messages: &[Vec<u8>]) -> Vec<u8> {
    let body_len: usize = messages.iter().map(Vec::len).sum();
    let mut h = Vec::with_capacity(16 + body_len);
    h.push(1); // version
    h.push(0); // reserved
    h.extend_from_slice(&(messages.len() as u16).to_le_bytes());
    h.extend_from_slice(&1u32.to_le_bytes()); // reference count
    h.extend_from_slice(&(body_len as u32).to_le_bytes());
    h.extend_from_slice(&[0, 0, 0, 0]); // pad the prefix to 16 bytes
    for m in messages {
        h.extend_from_slice(m);
    }
    h
}

/// An attribute message, version 1: name, then the datatype and dataspace
/// messages inline, then the value. Each of the three is padded to 8 bytes.
fn attribute_message(name: &str, datatype: &[u8], dataspace: &[u8], data: &[u8]) -> Vec<u8> {
    let mut name_b = name.as_bytes().to_vec();
    name_b.push(0);
    let name_size = name_b.len();
    pad_to_8(&mut name_b);
    let mut dt = datatype.to_vec();
    let dt_size = dt.len();
    pad_to_8(&mut dt);
    let mut ds = dataspace.to_vec();
    let ds_size = ds.len();
    pad_to_8(&mut ds);

    let mut m = Vec::new();
    m.push(1); // version
    m.push(0); // reserved
    m.extend_from_slice(&(name_size as u16).to_le_bytes());
    m.extend_from_slice(&(dt_size as u16).to_le_bytes());
    m.extend_from_slice(&(ds_size as u16).to_le_bytes());
    m.extend_from_slice(&name_b);
    m.extend_from_slice(&dt);
    m.extend_from_slice(&ds);
    m.extend_from_slice(data);
    m
}

/// One 40-byte symbol table entry.
///
/// `scratch` is the 16-byte scratch pad: for a group (cache type 1) it holds
/// the B-tree and local heap addresses, letting a reader find a subgroup's
/// contents without reading its object header first. Datasets use cache type 0
/// and leave it zero.
fn symbol_entry(name_offset: u64, header_addr: u64, cache_type: u32, scratch: [u8; 16]) -> Vec<u8> {
    let mut e = Vec::with_capacity(40);
    e.extend_from_slice(&name_offset.to_le_bytes());
    e.extend_from_slice(&header_addr.to_le_bytes());
    e.extend_from_slice(&cache_type.to_le_bytes());
    e.extend_from_slice(&[0, 0, 0, 0]); // reserved
    e.extend_from_slice(&scratch);
    e
}

// ---------------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------------

/// What writing one group produced: where its object header, B-tree and local
/// heap live. The parent needs all three — the header for the symbol table
/// entry, the other two for that entry's scratch pad.
struct WrittenGroup {
    header: u64,
    btree: u64,
    heap: u64,
}

/// Write `root`'s tree into `img` and return its addresses.
///
/// Children are written first so that their addresses are known when the
/// parent's symbol table is built, which is why this recurses before it
/// allocates anything of its own.
fn write_group(img: &mut Image, group: &GroupSpec<'_>) -> Result<WrittenGroup, String> {
    if group.child_count() > 2 * GROUP_LEAF_K as usize {
        return Err(format!(
            "group {:?} has {} children; this writer emits a single symbol table node, \
             which holds at most {}",
            group.name,
            group.child_count(),
            2 * GROUP_LEAF_K
        ));
    }

    // (name, object header address, cache type, scratch pad)
    let mut children: Vec<(String, u64, u32, [u8; 16])> = Vec::new();

    for sub in &group.groups {
        let w = write_group(img, sub)?;
        let mut scratch = [0u8; 16];
        scratch[..8].copy_from_slice(&w.btree.to_le_bytes());
        scratch[8..].copy_from_slice(&w.heap.to_le_bytes());
        children.push((sub.name.to_string(), w.header, 1, scratch));
    }

    for ds in &group.datasets {
        let bytes = ds.data.bytes()?;
        let data_addr = if bytes.is_empty() {
            // A zero-length dataset still needs an address; point it at the
            // current end of file rather than at the undefined address, which
            // some readers treat as "not allocated" and then refuse to read.
            img.reserve(0)
        } else {
            img.alloc(&bytes)
        };
        let messages = vec![
            header_message(MSG_DATASPACE, &dataspace_message(&[ds.data.len() as u64])),
            header_message(MSG_DATATYPE, &ds.data.datatype()),
            header_message(MSG_FILL_VALUE, &fill_value_message()),
            header_message(
                MSG_DATA_LAYOUT,
                &data_layout_message(data_addr, bytes.len() as u64),
            ),
        ];
        let hdr = img.alloc(&object_header(&messages));
        children.push((ds.name.to_string(), hdr, 0, [0u8; 16]));
    }

    // A symbol table is searched by binary search on the name, so entries must
    // be in lexicographic order.
    children.sort_by(|a, b| a.0.cmp(&b.0));

    // Local heap data segment: an empty string at offset 0 (the B-tree's first
    // key points at it), then each child's name, null-terminated.
    let mut heap_data: Vec<u8> = vec![0u8; 8];
    let mut name_offsets: Vec<u64> = Vec::with_capacity(children.len());
    for (name, _, _, _) in &children {
        let off = heap_data.len() as u64;
        heap_data.extend_from_slice(name.as_bytes());
        heap_data.push(0);
        pad_to_8(&mut heap_data);
        name_offsets.push(off);
    }
    let heap_data_addr = img.alloc(&heap_data);

    let mut heap = Vec::with_capacity(32);
    heap.extend_from_slice(b"HEAP");
    heap.push(0); // version
    heap.extend_from_slice(&[0, 0, 0]); // reserved
    heap.extend_from_slice(&(heap_data.len() as u64).to_le_bytes());
    heap.extend_from_slice(&HEAP_FREE_NULL.to_le_bytes()); // no free block
    heap.extend_from_slice(&heap_data_addr.to_le_bytes());
    let heap_addr = img.alloc(&heap);

    // Symbol table node. Allocated at full capacity (2K entries) whatever the
    // occupancy, as HDF5 does.
    let capacity = 2 * GROUP_LEAF_K as usize;
    let mut snod = Vec::with_capacity(8 + capacity * 40);
    snod.extend_from_slice(b"SNOD");
    snod.push(1); // version
    snod.push(0); // reserved
    snod.extend_from_slice(&(children.len() as u16).to_le_bytes());
    for (i, (_, hdr, cache, scratch)) in children.iter().enumerate() {
        snod.extend_from_slice(&symbol_entry(name_offsets[i], *hdr, *cache, *scratch));
    }
    snod.resize(8 + capacity * 40, 0);
    let snod_addr = img.alloc(&snod);

    // A single-node B-tree over that one leaf. The keys bracket the node: the
    // first is the empty string at heap offset 0, the last is the greatest
    // name present.
    let last_key = name_offsets.last().copied().unwrap_or(0);
    let mut btree = Vec::new();
    btree.extend_from_slice(b"TREE");
    btree.push(0); // node type 0: group node
    btree.push(0); // level 0: leaf
    btree.extend_from_slice(&1u16.to_le_bytes()); // entries used
    btree.extend_from_slice(&UNDEF_ADDR.to_le_bytes()); // left sibling
    btree.extend_from_slice(&UNDEF_ADDR.to_le_bytes()); // right sibling
    btree.extend_from_slice(&0u64.to_le_bytes()); // key 0
    btree.extend_from_slice(&snod_addr.to_le_bytes()); // child 0
    btree.extend_from_slice(&last_key.to_le_bytes()); // key 1
    // Allocated at full capacity, like the symbol table node.
    let btree_capacity = 24 + (2 * GROUP_INTERNAL_K as usize) * (SIZEOF_ADDR * 2) + SIZEOF_ADDR;
    btree.resize(btree_capacity, 0);
    let btree_addr = img.alloc(&btree);

    let mut messages = vec![header_message(MSG_SYMBOL_TABLE, &{
        let mut b = Vec::with_capacity(16);
        b.extend_from_slice(&btree_addr.to_le_bytes());
        b.extend_from_slice(&heap_addr.to_le_bytes());
        b
    })];
    for a in &group.attrs {
        messages.push(header_message(MSG_ATTRIBUTE, &encode_attr(a)));
    }
    let header_addr = img.alloc(&object_header(&messages));

    Ok(WrittenGroup {
        header: header_addr,
        btree: btree_addr,
        heap: heap_addr,
    })
}

fn encode_attr(a: &Attr) -> Vec<u8> {
    let (dt, ds, data) = match &a.value {
        AttrValue::Str(s) => (string_datatype(s.len().max(1)), dataspace_message(&[]), {
            let mut v = s.as_bytes().to_vec();
            if v.is_empty() {
                v.push(0);
            }
            v
        }),
        AttrValue::StrArray { width, values } => {
            let mut data = Vec::with_capacity(values.len() * width);
            for s in values {
                let b = s.as_bytes();
                let n = b.len().min(*width);
                data.extend_from_slice(&b[..n]);
                data.resize(data.len() + (width - n), 0);
            }
            (
                string_datatype(*width),
                dataspace_message(&[values.len() as u64]),
                data,
            )
        }
        AttrValue::I64(v) => (
            fixed_point_datatype(8),
            dataspace_message(&[]),
            v.to_le_bytes().to_vec(),
        ),
        AttrValue::I64Array(vs) => {
            let mut data = Vec::with_capacity(vs.len() * 8);
            for v in vs {
                data.extend_from_slice(&v.to_le_bytes());
            }
            (
                fixed_point_datatype(8),
                dataspace_message(&[vs.len() as u64]),
                data,
            )
        }
    };
    attribute_message(&a.name, &dt, &ds, &data)
}

/// Serialise `root` as a complete HDF5 file image.
pub fn build(root: &GroupSpec<'_>) -> Result<Vec<u8>, String> {
    let mut img = Image::new();
    let w = write_group(&mut img, root)?;
    img.align();
    let eof = img.buf.len() as u64;

    let mut sb = Vec::with_capacity(SUPERBLOCK_LEN);
    sb.extend_from_slice(&SIGNATURE);
    sb.push(0); // superblock version
    sb.push(0); // free space storage version
    sb.push(0); // root group symbol table entry version
    sb.push(0); // reserved
    sb.push(0); // shared header message format version
    sb.push(SIZEOF_ADDR as u8);
    sb.push(SIZEOF_ADDR as u8); // size of lengths
    sb.push(0); // reserved
    sb.extend_from_slice(&GROUP_LEAF_K.to_le_bytes());
    sb.extend_from_slice(&GROUP_INTERNAL_K.to_le_bytes());
    sb.extend_from_slice(&0u32.to_le_bytes()); // file consistency flags
    sb.extend_from_slice(&0u64.to_le_bytes()); // base address
    sb.extend_from_slice(&UNDEF_ADDR.to_le_bytes()); // free space info
    sb.extend_from_slice(&eof.to_le_bytes());
    sb.extend_from_slice(&UNDEF_ADDR.to_le_bytes()); // driver information block

    let mut scratch = [0u8; 16];
    scratch[..8].copy_from_slice(&w.btree.to_le_bytes());
    scratch[8..].copy_from_slice(&w.heap.to_le_bytes());
    sb.extend_from_slice(&symbol_entry(0, w.header, 1, scratch));

    debug_assert_eq!(sb.len(), SUPERBLOCK_LEN);
    img.patch(0, &sb);
    Ok(img.buf)
}

/// Serialise `root` and write it to `path`.
pub fn write_file(path: &Path, root: &GroupSpec<'_>) -> Result<(), Error> {
    let image = build(root).map_err(|e| Error::io(std::io::Error::other(e), path))?;
    let file = std::fs::File::create(path).map_err(|e| Error::io(e, path))?;
    let mut w = std::io::BufWriter::new(file);
    w.write_all(&image).map_err(|e| Error::io(e, path))?;
    w.flush().map_err(|e| Error::io(e, path))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The file image starts with the format signature and declares 8-byte
    /// offsets, an end-of-file address equal to its own length, and a root
    /// group whose object header is inside the file.
    #[test]
    fn superblock_is_well_formed() {
        let root = GroupSpec::new("/").dataset("x", Data::I32(&[1, 2, 3]));
        let img = build(&root).unwrap();

        assert_eq!(&img[..8], &SIGNATURE);
        assert_eq!(img[8], 0, "superblock version");
        assert_eq!(img[13], 8, "size of offsets");
        assert_eq!(img[14], 8, "size of lengths");

        let eof = u64::from_le_bytes(img[40..48].try_into().unwrap());
        assert_eq!(eof as usize, img.len(), "end-of-file address");

        let root_hdr = u64::from_le_bytes(img[64..72].try_into().unwrap());
        assert!(root_hdr >= SUPERBLOCK_LEN as u64 && (root_hdr as usize) < img.len());
        assert_eq!(img[root_hdr as usize], 1, "object header version");
    }

    /// Every structure this writer emits is 8-byte aligned, which the format
    /// requires of object headers and which readers assume of the rest.
    #[test]
    fn every_allocation_is_eight_byte_aligned() {
        let mut img = Image::new();
        let a = img.alloc(&[1, 2, 3]);
        let b = img.alloc(&[4]);
        let c = img.reserve(5);
        let d = img.alloc(&[]);
        for addr in [a, b, c, d] {
            assert_eq!(addr % 8, 0, "address {addr} is not 8-byte aligned");
        }
    }

    /// A string longer than its field is rejected. Truncating a barcode would
    /// produce a file that looks valid and holds a different barcode.
    #[test]
    fn an_overlong_string_is_an_error_not_a_truncation() {
        let root = GroupSpec::new("/").dataset(
            "bc",
            Data::FixedString {
                width: 4,
                values: vec!["ACGTACGT"],
            },
        );
        let err = build(&root).unwrap_err();
        assert!(err.contains("longer than"), "unexpected error: {err}");
    }

    /// Group children are written in lexicographic order, as a symbol table's
    /// binary search requires, whatever order they were added in.
    #[test]
    fn symbol_table_entries_are_sorted_by_name() {
        let root = GroupSpec::new("/")
            .dataset("zeta", Data::I32(&[1]))
            .dataset("alpha", Data::I32(&[2]))
            .dataset("mu", Data::I32(&[3]));
        let img = build(&root).unwrap();

        // The heap data segment holds an 8-byte empty string then the names in
        // the order they were written.
        let names: Vec<&str> = ["alpha", "mu", "zeta"].into();
        let mut last = 0usize;
        for n in names {
            let pos = find(&img, n.as_bytes()).unwrap_or_else(|| panic!("{n} not in image"));
            assert!(pos > last, "{n} is out of order in the local heap");
            last = pos;
        }
    }

    /// More children than one symbol table node holds is refused rather than
    /// written as a file with entries silently dropped.
    #[test]
    fn too_many_children_for_one_node_is_refused() {
        let mut root = GroupSpec::new("/");
        for i in 0..9 {
            root.datasets.push(DatasetSpec {
                name: match i {
                    0 => "a",
                    1 => "b",
                    2 => "c",
                    3 => "d",
                    4 => "e",
                    5 => "f",
                    6 => "g",
                    7 => "h",
                    _ => "i",
                },
                data: Data::I32(&[0]),
            });
        }
        let err = build(&root).unwrap_err();
        assert!(err.contains("at most 8"), "unexpected error: {err}");
    }

    fn find(hay: &[u8], needle: &[u8]) -> Option<usize> {
        hay.windows(needle.len()).position(|w| w == needle)
    }

    /// Integers land in the file little-endian, in element order.
    #[test]
    fn integer_data_is_little_endian_in_element_order() {
        let d = Data::I32(&[1, 256]);
        assert_eq!(d.bytes().unwrap(), vec![1, 0, 0, 0, 0, 1, 0, 0]);
        let d = Data::I64(&[-1]);
        assert_eq!(d.bytes().unwrap(), vec![0xff; 8]);
    }

    /// The image is byte-deterministic: no timestamps, no allocator addresses,
    /// nothing that varies between runs. Determinism is a project-wide property
    /// and this is the file format most likely to smuggle in a violation.
    ///
    /// The checksum also locks the layout. It was produced from an image that
    /// HDF5 1.14.6 accepted — `h5dump -H`, `h5ls -r`, and an `h5repack` +
    /// `h5diff` round-trip reporting no difference (`test/h5_conformance.sh`).
    /// If a change to this module moves it, re-run that script before updating
    /// the constant, because only libhdf5 can say whether the new bytes are
    /// still a valid file.
    #[test]
    fn the_image_is_deterministic_and_matches_the_validated_layout() {
        let build_it = || {
            let root = GroupSpec::new("/")
                .group(GroupSpec::new("matrix").dataset("data", Data::I32(&[1, 2])))
                .attr("filetype", AttrValue::Str("matrix".into()));
            build(&root).unwrap()
        };
        let a = build_it();
        assert_eq!(a, build_it(), "two builds of the same tree differ");
        assert_eq!(a.len(), 2192, "image length");
        assert_eq!(fnv1a(&a), 0x67f1_1bdc_8bcd_42e4, "image checksum");
    }

    /// FNV-1a, for the golden check above. Not a security property: a
    /// dependency-free way to notice that the bytes moved.
    fn fnv1a(bytes: &[u8]) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        for b in bytes {
            h ^= *b as u64;
            h = h.wrapping_mul(0x1000_0000_01b3);
        }
        h
    }

    /// Strings are padded with NULs to the field width, which is what makes
    /// numpy read the dataset back as `|S<width>`.
    #[test]
    fn strings_are_null_padded_to_the_field_width() {
        let d = Data::FixedString {
            width: 4,
            values: vec!["AC", "ACGT"],
        };
        assert_eq!(d.bytes().unwrap(), b"AC\0\0ACGT");
    }
}
