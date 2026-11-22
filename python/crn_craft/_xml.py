"""The one XML reader of the package.

Every XML document that CRN-craft reads (SBML, SED-ML and the COMBINE
manifest) goes through :func:`_read_xml`, and no other module of the package
parses XML. The function holds all safety checks, so that the importers and
the experiment reader cannot drift apart.

A document is read only if it is at most ``limit`` bytes of well-formed UTF-8
XML without a document type declaration. Nothing in the standard library
parser then has anything to expand or fetch: a DTD, internal or external
entities and parameter entities cannot occur, and no external resource is
ever opened. A byte order mark, a processing instruction, a comment and CDATA
are ordinary XML and read normally; any of ``<!DOCTYPE`` or ``<!ENTITY``
anywhere in the bytes, even inside a comment or CDATA section, is refused.
"""

import re
import xml.etree.ElementTree as ET

# A document that declares an encoding other than UTF-8 (the parser would
# reinterpret the bytes). Only the XML declaration at the start is looked at.
_ENCODING = re.compile(rb"""<\?xml[^>]*?encoding\s*=\s*["']([^"']*)["']""")
_UTF8 = (b"utf-8", b"utf8")
_PARSE_ERRORS = (ET.ParseError, UnicodeDecodeError, LookupError, ValueError)


class _XmlRefused(Exception):
    """A document the guard refuses.

    ``reason`` is one of ``"type"`` (not bytes), ``"size"`` (over the limit),
    ``"nul"`` (NUL bytes: UTF-16 or UTF-32), ``"dtd"`` (a document type or
    entity declaration), ``"encoding"`` (a declared encoding other than
    UTF-8), ``"utf8"`` (bytes that are not UTF-8) or ``"syntax"`` (not
    well-formed XML). ``message`` is a short human message. ``detail`` is the
    text the parser itself gives when it also objects to a document with a
    declared encoding, and ``None`` otherwise.
    """

    def __init__(self, reason, message, detail=None):
        super().__init__(message)
        self.reason = reason
        self.message = message
        self.detail = detail


def _parse(raw):
    raw.decode("utf-8")
    return ET.fromstring(raw)


def _read_xml(raw, limit):
    """Return the root element of one XML document given as ``bytes``.

    ``raw`` may be ``bytes``, ``bytearray`` or ``memoryview``; ``limit`` is the
    largest size in bytes. Raises :class:`_XmlRefused` for everything the
    package does not read. The checks run in this order: type, size, NUL
    bytes, document type or entity declaration, declared encoding, UTF-8,
    well-formedness; a document that breaks several reports the first.
    """
    if not isinstance(raw, bytes | bytearray | memoryview):
        raise _XmlRefused("type", "XML must be given as bytes")
    size = memoryview(raw).nbytes
    if size > limit:
        raise _XmlRefused("size", f"XML document exceeds {limit} bytes")
    raw = bytes(raw)
    # NUL bytes: UTF-16 or UTF-32 would hide a DTD from the byte check below.
    if b"\x00" in raw:
        raise _XmlRefused("nul", "XML must be UTF-8, not UTF-16 or UTF-32")
    # The standard parser expands entities, so no document type is read. The
    # check covers comments and CDATA as well, on purpose.
    upper = raw.upper()
    if b"<!DOCTYPE" in upper or b"<!ENTITY" in upper:
        raise _XmlRefused("dtd", "XML may not contain a document type or entities")
    declared = _ENCODING.match(raw.removeprefix(b"\xef\xbb\xbf"))
    if declared and declared[1].lower() not in _UTF8:
        # No DTD can be present here, so the parse is safe; it only tells
        # what the parser itself says about this document.
        try:
            _parse(raw)
            detail = None
        except _PARSE_ERRORS as error:
            detail = str(error)
        raise _XmlRefused("encoding", "XML must be UTF-8", detail)
    try:
        return _parse(raw)
    except UnicodeDecodeError as error:
        raise _XmlRefused("utf8", str(error)) from error
    except _PARSE_ERRORS as error:
        raise _XmlRefused("syntax", str(error)) from error
