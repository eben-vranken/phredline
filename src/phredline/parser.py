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


if __name__ == "__main__":
     print("--- Chunk Reader ---")

     chunk_stream = read_chunks(FILE_PATH, CHUNK_SIZE)

     for i, chunk in enumerate(chunk_stream):
         print(f"Chunk {i}: {chunk} (Length: {len(chunk)})")

         if i >= 100:
            print("Stopping test early")
            break