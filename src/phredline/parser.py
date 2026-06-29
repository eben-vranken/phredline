FILE_PATH = "data\ERR1410268.fastq\ERR1410268.fastq"
CHUNK_SIZE = 16

def read_chunks(file_path, chunk_size):
    with open(file_path, "rb") as stream:
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

if __name__ == "__main__":
    print("--- Chunk Reader ---")

    chunks = read_chunks(FILE_PATH, CHUNK_SIZE)
    lines = read_lines(chunks)

    for i, chunk in enumerate(lines):
        print(f"Chunk {i}: {chunk} (Length: {len(chunk)})")

        if i >= 100:
           print("Stopping test early")
           break