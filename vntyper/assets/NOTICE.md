# Vendored third-party assets

Everything in this directory is a third-party file VNtyper redistributes verbatim.
Nothing here is hand-edited: the bytes are what the digests below say they are, and
`vntyper/scripts/report_assets.py` refuses to embed them if they are not.

## igv.js 3.0.2

| | |
| --- | --- |
| Upstream | <https://github.com/igvteam/igv.js> |
| Version | 3.0.2 |
| Fetched from | `https://cdn.jsdelivr.net/npm/igv@3.0.2/dist/igv.min.js` |
| Licence | MIT (text below) |
| File here | `igv-3.0.2.min.js.gz` |

| Form | Bytes | SHA-256 |
| --- | --- | --- |
| `igv.min.js` as fetched | 1,310,337 | `ab1aa79c514ee3a0d66a0ffc788b6d37803910e62cf6d114d9b2909d96b5e790` |
| `igv-3.0.2.min.js.gz` (this file) | 372,690 | `0d8b512654b2ef588009453c403c8f1329dce88eedca90ba9e60888af6b2f79f` |

The first digest is the provenance claim: it is what a reader can check against
upstream, and it is the one printed in every report's Provenance section. The second is
what the file on disk carries. `report_assets.igv_payload` verifies **both** — the file
it read, and then the bytes that come back out of it, which are exactly the source the
report later evaluates.

The gzip is reproducible: run `gzip.compress(raw, compresslevel=9, mtime=0)`, then set
byte 9 of the result to `3` (RFC 1952's Unix OS value). `mtime=0` fixes the timestamp;
canonicalising the OS byte avoids a Python-version difference where supported
interpreters emit either Unix (`3`) or unknown (`255`). Both are required for the
same source to retain the same file digest across the supported Python matrix.

### Why 3.0.2 and not the newest release

Size, measured rather than assumed. 3.0.2 costs 496,920 base64 characters in the
document. 3.8.5 costs 575,932 — 79 KB *worse* for a report whose whole point is being
small enough to archive and forward. Upgrade when there is a reason other than
recency, and re-measure when you do.

### Licence

igv.js is redistributed under the MIT licence, reproduced here in full as that licence
requires:

```text
The MIT License (MIT)

Copyright (c) 2014-2019 Broad Institute of MIT and Harvard, Regents of the University of California

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

The minified bundle additionally carries its own dependencies' notices inline — jQuery
Foundation, Gliffy Inc., Robert Buels' `two-bit.js`, and the `vanilla-picker` authors
are each named in comments inside `igv.min.js` itself. Those comments are part of the
bytes digested above and are redistributed with them unchanged.
