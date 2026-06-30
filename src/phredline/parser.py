import gzip
import sys

CHUNK_SIZE = 16


class ParseError(Exception):
    """Error for when FATQ format is malformed"""

    pass


def read_chunks(file_path, chunk_size):
    if file_path is None or file_path == "-":
        stream = sys.stdin.buffer

        while True:
            bchunk = stream.read(chunk_size)
            if bchunk == b"":
                break
            yield bchunk
    else:
        with open(file_path, "rb") as f:
            magic = f.read(2)

            f.seek(0)

            if magic == b"\x1f\x8b":
                stream = gzip.open(f, "rb")
            else:
                stream = f

            while True:
                bchunk = stream.read(chunk_size)

                if bchunk == b"":
                    break
                else:
                    yield bchunk


def read_lines(chunk_stream):
    tail = b""

    for chunk in chunk_stream:
        combined_data = tail + chunk

        chunk_split = combined_data.split(b"\n")

        for split in chunk_split[:-1]:
            yield split

        tail = chunk_split[-1]

        pass

    if len(tail) > 0:
        yield tail


def parse_records(line_stream):
    record = []

    for line in line_stream:
        record.append(line)

        if len(record) == 4:
            if not record[0].startswith(b"@") or not record[2].startswith(b"+"):
                raise ParseError("Error when parsing FASTQ. Malformed line")

            yield (record[0], record[1], record[2], record[3])

            record = []

def decode_phred(quality_string):
    quality_array = []

    for c in quality_string:
        quality_array.append(c-33)

    return quality_array

def decode_error_probability(q_array):
    error_probability_array = []

    for q in q_array:
        error_probability_array.append(10**(-q/10))

    return error_probability_array

if __name__ == "__main__":
    print("--- Chunk Reader ---")

    file_input = sys.argv[1] if len(sys.argv) > 1 else None

    chunks = read_chunks(file_input, CHUNK_SIZE)
    lines = read_lines(chunks)
    records = parse_records(lines)

    for i, rec in enumerate(records):
        header, seq, plus, qual = rec
        print(f"Record {i}:")
        print(f"  Header: {header}")
        print(f"  Seq:    {seq}")
        print(f"  Qual:   {qual}")

        phred_qualities = decode_phred(qual)

        print(f"  Phred: ", end=" ")
        for phred in phred_qualities: 
            print(f"{phred}", end="-")

        error_probabilities = decode_error_probability(phred_qualities)

        print(f"\n  Error probability: ", end=" ")
        for prob in error_probabilities: 
            print(f"{prob}", end="-")


        if i >= 0:
            print("\nStopping test early")
            break

            
