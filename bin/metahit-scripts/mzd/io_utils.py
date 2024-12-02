import bz2
import pickle
import gzip
import json
import io
import yaml

# default buffer for incremental read/write
DEF_BUFFER = 16384


def save_object(file_name, obj):
    """
    Serialize an object to a file with gzip compression. .gz will automatically be
    added if missing.

    :param file_name: output file name
    :param obj: object to serialize
    """
    with open_output(file_name, compress='gzip') as out_h:
        pickle.dump(obj, out_h)


def load_object(file_name):
    """
    Deserialize an object from a file with automatic support for compression.

    :param file_name: input file name
    :return: deserialized object
    """
    with open_input(file_name) as in_h:
        return pickle.load(in_h)


def open_input(file_name):
    """
    Open a text file for input. The filename is used to indicate if it has been
    compressed. Recognizing gzip and bz2.

    :param file_name: the name of the input file
    :return: open file handle, possibly wrapped in a decompressor
    """
    suffix = file_name.split('.')[-1].lower()
    if suffix == 'bz2':
        return bz2.open(file_name, 'rb')
    elif suffix == 'gz':
        return gzip.open(file_name, 'rb')
    else:
        return open(file_name, 'r')


def open_output(file_name, append=False, compress=None, gzlevel=6):
    """
    Open a text stream for reading or writing. Compression can be enabled
    with either 'bzip2' or 'gzip'. Additional option for gzip compression
    level. Compressed filenames are only appended with suffix if not included.

    :param file_name: file name of output
    :param append: append to any existing file
    :param compress: gzip, bzip2
    :param gzlevel: gzip level (default 6)
    :return:
    """
    mode = 'wb' if not append else 'ab'

    if compress == 'bzip2':
        if not file_name.endswith('.bz2'):
            file_name += '.bz2'
        return bz2.open(file_name, mode, compresslevel=9)
    elif compress == 'gzip':
        if not file_name.endswith('.gz'):
            file_name += '.gz'
        return gzip.open(file_name, mode, compresslevel=gzlevel)
    else:
        return open(file_name, mode)


def multicopy_tostream(file_name, *ostreams, **kwargs):
    bufsize = DEF_BUFFER if 'bufsize' not in kwargs else kwargs['bufsize']

    with open(file_name, 'rb') as in_h:
        while True:
            buf = in_h.read(bufsize)
            if not buf:
                break
            for oi in ostreams:
                oi.write(buf)


def multicopy_tofile(file_name, *onames, **kwargs):
    bufsize = DEF_BUFFER if 'bufsize' not in kwargs else kwargs['bufsize']
    write_mode = "wb" if 'write_mode' not in kwargs else kwargs['write_mode']
    compress = None if 'compress' not in kwargs else kwargs['compress']

    try:
        with open(file_name, 'rb') as in_h:
            out_h = [open_output(oi, write_mode, compress) for oi in onames]

            while True:
                buf = in_h.read(bufsize)
                if not buf:
                    break
                for oi in out_h:
                    oi.write(buf)
    finally:
        if out_h:
            for oi in out_h:
                if oi:
                    oi.close()


def write_to_stream(stream, data, fmt='plain'):
    if fmt == 'yaml':
        yaml.dump(data, stream, default_flow_style=False)
    elif fmt == 'json':
        json.dump(data, stream, indent=1)
    elif fmt == 'plain':
        stream.write(f'{data}\n')


def read_from_stream(stream, fmt='yaml'):
    if fmt == 'yaml':
        return yaml.safe_load(stream)
    elif fmt == 'json':
        return json_load_byteified(stream)


def json_loads_byteified(json_text):
    return _byteify(
        json.loads(json_text, object_hook=_byteify),
        ignore_dicts=True
    )


def json_load_byteified(file_handle):
    return _byteify(
        json.load(file_handle, object_hook=_byteify),
        ignore_dicts=True
    )


def _byteify(data, ignore_dicts=False):
    if isinstance(data, str):
        return data
    if isinstance(data, list):
        return [_byteify(item, ignore_dicts=True) for item in data]
    if isinstance(data, dict) and not ignore_dicts:
        return {
            _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True)
            for key, value in data.items()
        }
    return data
